"""
group_compression.py — compressão do ansatz pelo grupo pontual COMPLETO.

Consolida o que originalmente estava em c3v.py, c3v3.py, td.py e abelian_ctrl.py
(perdidos na limpeza do scratchpad). Self-contained e reproduz os mesmos números.

Constrói a representação do grupo pontual na base de MOs a partir da ação cartesiana de
cada elemento sobre as funções de base, induz em
    V        = excitações simples espaciais (occ (x) virt)
    Sym^2 V  = duplos do singlet-UCCSD
e projeta no subespaço totalmente simétrico com P = (1/|G|) sum_g Rep(g).

Resultados esperados (STO-3G):
    H2O  C2v (Abeliano)   UCCSD   65  -> subgrupo  26 -> grupo pleno  26   ganho 0
    C2H4 D2h (Abeliano)   UCCSD 1224  -> subgrupo 231 -> grupo pleno 231   ganho 0
    NH3  C3v              UCCSD  135  -> subgrupo  75 -> grupo pleno  30
    CH4  Td               UCCSD  230  -> subgrupo  65 -> grupo pleno  21
todos com o mesmo erro contra o FCI.

Uso:  python group_compression.py [h2o|c2h4|nh3|ch4|all]
"""
import sys, math, time, itertools, json
import numpy as np
import scipy.optimize as opt
import scipy.sparse.linalg as sla
from pyscf import gto, scf, ao2mo, symm

HERE = __import__('os').path.dirname(__import__('os').path.abspath(__file__))  # PATCHED: was a hard-coded absolute path
sys.path.insert(0, HERE)
from indep import build_det_basis, SpinFreeHam, kappa_matrix, TrotterAnsatz
from audit import gates_of

PMAP = {'px': 0, 'py': 1, 'pz': 2}


# --------------------------------------------------------------- geometrias
def geo_nh3(d=1.09, ang=107.3):
    t = math.radians(ang); r = d * math.sin(t / 2); z = -d * math.cos(t / 2)
    at = [('N', (0., 0., 0.))]
    for k in range(3):
        p = 2 * math.pi * k / 3
        at.append(('H', (r * math.cos(p), r * math.sin(p), z)))
    return at


def geo_ch4(d=1.09):
    c = d / math.sqrt(3)
    return [('C', (0., 0., 0.))] + [('H', (c * x, c * y, c * z))
                                    for x, y, z in ((1, 1, 1), (1, -1, -1),
                                                    (-1, 1, -1), (-1, -1, 1))]


def geo_h2o(d=1.02, ang=104.5):
    t = math.radians(ang) / 2
    return [('O', (0., 0., 0.)), ('H', (d * math.sin(t), 0., d * math.cos(t))),
            ('H', (-d * math.sin(t), 0., d * math.cos(t)))]


def geo_c2h4(d1=1.34, d2=1.08, a=120.0):
    r = math.radians(a); x, y = d2 * math.cos(r), d2 * math.sin(r); hc = d1 / 2
    return [('C', (-hc, 0., 0.)), ('C', (hc, 0., 0.)),
            ('H', (-hc - x, y, 0.)), ('H', (-hc - x, -y, 0.)),
            ('H', (hc + x, y, 0.)), ('H', (hc + x, -y, 0.))]


def to_str(at):
    return "; ".join(f"{s} {p[0]:.12f} {p[1]:.12f} {p[2]:.12f}" for s, p in at)


# --------------------------------------------------------------- grupos
def Rz(a):
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])


def group(name):
    """(elementos do grupo pleno, elementos do subgrupo Abeliano usado pelo PySCF)."""
    if name == 'C3v':
        MX = np.diag([1., -1., 1.])                      # sigma_v contendo H1
        G = [np.eye(3), Rz(2 * math.pi / 3), Rz(4 * math.pi / 3),
             MX, Rz(2 * math.pi / 3) @ MX, Rz(4 * math.pi / 3) @ MX]
        return G, [G[0], G[3]]                            # Cs = {E, sigma_v}
    if name == 'Td':
        # as 24 matrizes de permutação com sinal cujo produto de sinais é +1
        G = []
        for perm in itertools.permutations(range(3)):
            for sg in itertools.product([1, -1], repeat=3):
                if np.prod(sg) != 1:
                    continue
                M = np.zeros((3, 3))
                for i in range(3):
                    M[i, perm[i]] = sg[i]
                G.append(M)
        D2 = [M for M in G if np.allclose(M, np.diag(np.diag(M)))]
        return G, D2
    if name == 'C2v':
        G = [np.diag([1., 1., 1.]), np.diag([-1., -1., 1.]),
             np.diag([1., -1., 1.]), np.diag([-1., 1., 1.])]
        return G, G                                       # já é Abeliano
    if name == 'D2h':
        G = [np.diag([float(a), float(b), float(c)])
             for a in (1, -1) for b in (1, -1) for c in (1, -1)]
        return G, G
    raise ValueError(name)


# --------------------------------------------------------------- máquina
def mo_rep(mol, mo, G):
    """D(g) na base de MOs, com validação."""
    S = mol.intor('int1e_ovlp'); coords = mol.atom_coords(); nao = mol.nao
    parsed = [(int(L.split()[0]), L.split()[2]) for L in mol.ao_labels()]

    def U_of(R):
        U = np.zeros((nao, nao))
        for mu, (k, sh) in enumerate(parsed):
            p = R @ coords[k]
            cand = [q for q in range(mol.natm)
                    if np.linalg.norm(p - coords[q]) < 1e-6
                    and mol.atom_symbol(q) == mol.atom_symbol(k)]
            if not cand:
                raise RuntimeError("elemento não mapeia átomo em átomo")
            j = cand[0]
            if sh[-2:] in PMAP:
                i = PMAP[sh[-2:]]
                for nu, (k2, s2) in enumerate(parsed):
                    if k2 == j and s2[:-2] == sh[:-2] and s2[-2:] in PMAP:
                        U[nu, mu] = R[PMAP[s2[-2:]], i]
            else:
                for nu, (k2, s2) in enumerate(parsed):
                    if k2 == j and s2 == sh:
                        U[nu, mu] = 1.0
        return U

    Us = [U_of(R) for R in G]
    err = max(np.abs(U.T @ S @ U - S).max() for U in Us)
    D = [mo.T @ S @ U @ mo for U in Us]
    return D, err


def invariant_basis(Ms, Sym2proj):
    """bases ortonormais dos subespaços invariantes de V e de Sym^2 V."""
    nS = Ms[0].shape[0]
    Ps = sum(Ms) / len(Ms)
    w, V = np.linalg.eigh((Ps + Ps.T) / 2)
    Bs = V[:, np.abs(w - 1) < 1e-6]
    Pd = (sum(np.kron(M, M) for M in Ms) / len(Ms)) @ Sym2proj
    w2, V2 = np.linalg.eigh((Pd + Pd.T) / 2)
    Bd = V2[:, np.abs(w2 - 1) < 1e-6]
    return Bs, Bd


def build_ops(Bs, Bd, Ep, nS):
    """Operadores anti-hermitianos invariantes.
       Singles: sum_k c_k (E_k - E_k^T).
       Doubles: T_pq = E_ai E_bj (comutam: i ocupado, b virtual anulam os dois deltas).
       O fator 2 fora da diagonal vem de sum_rs c_rs v_r (x) v_s
       = sum_{r<s} 2 c_rs T_rs + sum_r c_rr T_rr."""
    ops = []
    for j in range(Bs.shape[1]):
        M = sum(Bs[k, j] * Ep[k] for k in range(nS) if abs(Bs[k, j]) > 1e-10)
        ops.append((M - M.T).tocsr())
    for j in range(Bd.shape[1]):
        c = Bd[:, j].reshape(nS, nS); O = 0
        for p in range(nS):
            for q in range(p, nS):
                w = (2.0 * c[p, q] if p < q else c[p, p])
                if abs(w) > 1e-10:
                    X = Ep[p] @ Ep[q]
                    O = O + w * (X - X.T)
        ops.append(O.tocsr())
    return ops


def vqe_expm(ops, psi0, ham, egs, label):
    L = len(ops)

    def eg(th):
        v = psi0.copy()
        for k in range(L):
            v = sla.expm_multiply(th[k] * ops[k], v)
        Hv = ham.matvec(v); e = float(v @ Hv); w = 2.0 * Hv; g = np.zeros(L)
        for k in range(L - 1, -1, -1):
            g[k] = float(w @ ops[k].dot(v))
            v = sla.expm_multiply(-th[k] * ops[k], v)
            w = sla.expm_multiply(-th[k] * ops[k], w)
        return e, g
    r = opt.minimize(eg, np.zeros(L), jac=True, method='L-BFGS-B',
                     options=dict(maxiter=4000, maxfun=8000, ftol=1e-16, gtol=1e-11))
    print(f"    {label:36s} P={L:5d}  E={r.fun:.10f}  err={(r.fun-egs)*1000:8.4f} mHa",
          flush=True)
    return float(r.fun)


def run(tag, atoms, gname, do_energy=True):
    t0 = time.time()
    mol = gto.M(atom=to_str(atoms), basis='sto-3g', symmetry=True, verbose=0)
    mf = scf.RHF(mol).run()
    nao, nelec = mol.nao, mol.nelectron; no = nelec // 2
    G, H = group(gname)
    D, err = mo_rep(mol, mf.mo_coeff, G)
    mix = max(abs(D[g][q, p]) for g in range(len(G)) for p in range(nao)
              for q in range(nao) if abs(mf.mo_energy[p] - mf.mo_energy[q]) > 1e-6)
    print(f"\n{'='*82}\n{tag}  grupo={gname} (|G|={len(G)}, |H|={len(H)})  "
          f"topgroup={mol.topgroup} abelian={mol.groupname}")
    print(f"  validação: |U^T S U - S|={err:.2e}   máx elemento fora de camada={mix:.2e}")
    sing = [(i, a) for i in range(no) for a in range(no, nao)]; nS = len(sing)
    idx = {R.tobytes(): k for k, R in enumerate(G)}
    Ms = {k: np.array([[D[k][j, i] * D[k][b, a] for (i, a) in sing]
                       for (j, b) in sing]) for k in range(len(G))}
    Hk = [idx[R.tobytes()] for R in H]
    I = np.eye(nS * nS); SW = np.zeros((nS * nS, nS * nS))
    for r_ in range(nS):
        for s_ in range(nS):
            SW[r_ * nS + s_, s_ * nS + r_] = 1.0
    PI = 0.5 * (I + SW)
    res = {'molecule': tag, 'group': gname, 'nparam_uccsd': nS + nS * (nS + 1) // 2}
    for nm, sub in (('subgrupo Abeliano', Hk), ('grupo pleno', list(range(len(G))))):
        Bs, Bd = invariant_basis([Ms[k] for k in sub], PI)
        print(f"  {nm:20s}: simples {Bs.shape[1]:4d} + duplos {Bd.shape[1]:5d} "
              f"= {Bs.shape[1]+Bd.shape[1]:5d} parâmetros")
        res['n_' + ('sub' if 'sub' in nm else 'full')] = int(Bs.shape[1] + Bd.shape[1])
    print(f"  UCCSD completo: {res['nparam_uccsd']} parâmetros")
    if not do_energy:
        print("  (energia omitida: espaço de determinantes grande demais)")
        return res
    dets, index = build_det_basis(nao, no, no); n = len(dets)
    psi0 = np.zeros(n); psi0[index[sum((1 << (2*p)) | (1 << (2*p+1)) for p in range(no))]] = 1.0
    h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), nao)
    ham = SpinFreeHam(h1, eri, mol.energy_nuc(), dets, index)
    egs = float(sla.eigsh(sla.LinearOperator((n, n), matvec=ham.matvec, dtype=float),
                          k=1, which='SA', tol=1e-12)[0][0])
    print(f"  determinantes={n}  E_FCI={egs:.10f}")
    npar, gates = gates_of(nao, nelec)
    K = [kappa_matrix((g['occ'], g['virt']), dets, index) for g in gates]
    ans = TrotterAnsatz(K, [g['param'] for g in gates],
                        [g['coef'] for g in gates], npar, psi0, ham)
    r = opt.minimize(lambda x: ans.energy_grad(x), np.zeros(npar), jac=True,
                     method='L-BFGS-B',
                     options=dict(maxiter=20000, maxfun=60000, ftol=1e-16, gtol=1e-11))
    print(f"    {'UCCSD (Trotter, referência)':36s} P={npar:5d}  E={r.fun:.10f}  "
          f"err={(r.fun-egs)*1000:8.4f} mHa", flush=True)
    res['E_uccsd'] = float(r.fun); res['efci'] = egs
    Ep = [ham.E[a][i].tocsr() for (i, a) in sing]
    for nm, sub, key in (('invariante sob o subgrupo [CONTROLE]', Hk, 'E_sub'),
                         ('invariante sob o grupo pleno', list(range(len(G))), 'E_full')):
        Bs, Bd = invariant_basis([Ms[k] for k in sub], PI)
        res[key] = vqe_expm(build_ops(Bs, Bd, Ep, nS), psi0, ham, egs, nm)
    print(f"  [{time.time()-t0:.0f}s]")
    return res


CASES = {'h2o':  ('H2O / STO-3G',  geo_h2o(),  'C2v', True),
         'c2h4': ('C2H4 / STO-3G', geo_c2h4(), 'D2h', False),
         'nh3':  ('NH3 / STO-3G',  geo_nh3(),  'C3v', True),
         'ch4':  ('CH4 / STO-3G',  geo_ch4(),  'Td',  True)}

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    keys = list(CASES) if which == 'all' else [which]
    out = [run(*CASES[k]) for k in keys]
    json.dump(out, open(f'{HERE}/group_compression.json', 'w'), indent=2)
    print("\nsalvo em group_compression.json")
