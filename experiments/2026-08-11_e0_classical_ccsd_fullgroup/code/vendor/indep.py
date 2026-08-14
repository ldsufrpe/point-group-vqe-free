"""
Independent from-scratch reimplementation for auditing SymUCCSD claims.
No code inherited from He et al. (MindQuantum) or from any prior local script.

Conventions
-----------
spin-orbital index  = 2*p + s   (p spatial, s=0 alpha, s=1 beta)
determinant         = python int bitmask over spin orbitals,
                      |D> = prod_{k in occ, ascending} a^dag_k |vac>
Sector              = fixed N electrons, fixed Sz = 0.
"""
import itertools, math
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from pyscf import gto, scf, ao2mo, fci, symm


# ----------------------------------------------------------------- determinants
def popcount(x):
    return bin(x).count('1')


def build_det_basis(norb, nalpha, nbeta):
    dets = []
    for oa in itertools.combinations(range(norb), nalpha):
        ma = 0
        for p in oa:
            ma |= 1 << (2 * p)
        for ob in itertools.combinations(range(norb), nbeta):
            mb = 0
            for p in ob:
                mb |= 1 << (2 * p + 1)
            dets.append(ma | mb)
    dets.sort()
    return dets, {d: i for i, d in enumerate(dets)}


def _sign_below(det, k):
    """(-1)^{number of occupied spin-orbitals with index < k}"""
    return -1 if (popcount(det & ((1 << k) - 1)) & 1) else 1


def apply_annih(det, k):
    if not (det >> k) & 1:
        return None, 0
    return det & ~(1 << k), _sign_below(det, k)


def apply_creat(det, k):
    if (det >> k) & 1:
        return None, 0
    return det | (1 << k), _sign_below(det, k)


def apply_string(det, ops):
    """ops = list of ('c'|'a', spin_orbital), applied RIGHT to LEFT."""
    s = 1
    d = det
    for kind, k in reversed(ops):
        d, sg = (apply_creat(d, k) if kind == 'c' else apply_annih(d, k))
        if d is None:
            return None, 0
        s *= sg
    return d, s


def op_matrix(ops, dets, index):
    """Sparse matrix of the fermionic operator string in the determinant basis."""
    rows, cols, vals = [], [], []
    for j, d in enumerate(dets):
        nd, s = apply_string(d, ops)
        if nd is not None:
            i = index.get(nd)
            if i is not None:
                rows.append(i); cols.append(j); vals.append(float(s))
    n = len(dets)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))


# ----------------------------------------------------------------- Hamiltonian
class SpinFreeHam:
    """H |psi> without ever forming H, via E_pq = sum_sigma a^dag_{p sigma} a_{q sigma}."""

    def __init__(self, h1, eri, ecore, dets, index):
        self.norb = h1.shape[0]
        self.ecore = ecore
        self.eri = eri                                   # chemist (pq|rs)
        # h1_eff = h1 - 1/2 sum_r (p r | r q)
        self.h1eff = h1 - 0.5 * np.einsum('prrq->pq', eri)
        n = self.norb
        self.E = [[None] * n for _ in range(n)]
        for p in range(n):
            for q in range(n):
                M = op_matrix([('c', 2 * p), ('a', 2 * q)], dets, index) \
                    + op_matrix([('c', 2 * p + 1), ('a', 2 * q + 1)], dets, index)
                self.E[p][q] = M.tocsr()

    def matvec(self, v):
        n = self.norb
        nd = v.shape[0]
        chi = np.empty((n, n, nd))
        for r in range(n):
            for s in range(n):
                chi[r, s] = self.E[r][s].dot(v)
        out = self.ecore * v + np.einsum('pq,pqd->d', self.h1eff, chi, optimize=True)
        acc = np.einsum('pqrs,rsd->pqd', self.eri, chi, optimize=True)
        for p in range(n):
            for q in range(n):
                out += 0.5 * self.E[p][q].dot(acc[p, q])
        return out

    def as_linop(self, n):
        return sp.linalg.LinearOperator((n, n), matvec=self.matvec, dtype=float)


# ----------------------------------------------------------------- molecule
def make_molecule(atom, basis='sto-3g', symmetry=True):
    mol = gto.M(atom=atom, basis=basis, symmetry=symmetry, verbose=0)
    mf = scf.RHF(mol).run()
    norb = mol.nao
    nelec = mol.nelectron
    na = nb = nelec // 2
    orbsym = np.array(symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb, mf.mo_coeff))
    h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), norb)
    ecore = mol.energy_nuc()
    cisolver = fci.FCI(mol, mf.mo_coeff)
    efci = cisolver.kernel()[0]
    return dict(mol=mol, mf=mf, norb=norb, nelec=nelec, na=na, nb=nb,
                orbsym=orbsym, h1=h1, eri=eri, ecore=ecore,
                ehf=mf.e_tot, efci=efci, groupname=mol.groupname,
                topgroup=mol.topgroup, mo_energy=mf.mo_energy)


# ----------------------------------------------------------------- UCCSD pool
def uccsd_pool(no, nv):
    """
    Singlet-UCCSD parameterisation (OpenFermion convention):
      - singles  : one parameter per spatial (i->a)
      - doubles  : one parameter per unordered pair {(i->a),(j->b)} with repetition
    Returns list of parameters; each parameter carries its list of elementary
    spin-orbital excitation operators (as index tuples).
    """
    sing = [(i, a) for i in range(no) for a in range(nv)]
    params = []
    for (i, a) in sing:                                     # singles
        elems = []
        for s in (0, 1):
            elems.append((( 2 * i + s,), (2 * (no + a) + s,)))
        params.append(dict(kind='S', spatial=((i, a),), elems=elems))
    for p in range(len(sing)):                              # doubles
        for q in range(p, len(sing)):
            i, a = sing[p]; j, b = sing[q]
            elems = []
            seen = set()
            for s in (0, 1):
                for t in (0, 1):
                    oi = 2 * i + s;  oa = 2 * (no + a) + s
                    oj = 2 * j + t;  ob = 2 * (no + b) + t
                    if oi == oj or oa == ob:
                        continue
                    key = (tuple(sorted((oi, oj))), tuple(sorted((oa, ob))))
                    if key in seen:
                        continue
                    seen.add(key)
                    elems.append((key[0], key[1]))
            params.append(dict(kind='D', spatial=((i, a), (j, b)), elems=elems))
    return params


def param_irrep(prm, occ_sym, virt_sym):
    """XOR of the irrep ids of all spatial orbitals involved (PySCF D2h-subgroup convention)."""
    r = 0
    for (i, a) in prm['spatial']:
        r ^= int(occ_sym[i]) ^ int(virt_sym[a])
    return r


def kappa_matrix(elem, dets, index):
    """kappa = tau - tau^dag for one elementary excitation ((occ...),(virt...))."""
    occ, virt = elem
    ops = [('c', k) for k in reversed(virt)] + [('a', k) for k in occ]
    tau = op_matrix(ops, dets, index)
    return (tau - tau.T).tocsr()


# ----------------------------------------------------------------- ansatz
class TrotterAnsatz:
    """
    U(theta) = prod_k exp(theta_{p(k)} kappa_k), applied left-to-right on |HF>.
    kappa^3 = -kappa  =>  exp(t k) v = v + sin(t) k v + (1-cos t) k^2 v   (exact).
    Analytic gradient by the adjoint method.
    """

    def __init__(self, kappas, pidx, coef, nparam, psi0, ham):
        self.K = kappas            # list of sparse kappa_k (unit elementary generators)
        self.pidx = pidx           # parameter index of each gate
        self.coef = coef           # gate k implements exp(theta_{p(k)} * coef_k * kappa_k)
        self.np = nparam
        self.psi0 = psi0
        self.ham = ham

    def _gate(self, k, t, v):
        Kv = self.K[k].dot(v)
        K2v = self.K[k].dot(Kv)
        return v + math.sin(t) * Kv + (1.0 - math.cos(t)) * K2v

    def state(self, theta):
        v = self.psi0.copy()
        for k in range(len(self.K)):
            v = self._gate(k, theta[self.pidx[k]] * self.coef[k], v)
        return v

    def energy(self, theta):
        v = self.state(theta)
        return float(v @ self.ham.matvec(v))

    def energy_grad(self, theta):
        v = self.state(theta)                        # v = v_L
        Hv = self.ham.matvec(v)
        e = float(v @ Hv)
        w = 2.0 * Hv                                 # w_L
        g = np.zeros(self.np)
        for k in range(len(self.K) - 1, -1, -1):
            t = theta[self.pidx[k]] * self.coef[k]
            g[self.pidx[k]] += self.coef[k] * float(w @ self.K[k].dot(v))
            v = self._gate(k, -t, v)                 # v_{k-1} = G_k^{-1} v_k
            w = self._gate(k, -t, w)                 # w_{k-1} = G_k^T w_k
        return e, g
