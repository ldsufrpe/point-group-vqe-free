"""
E1 — a compressão pelo grupo pleno sobrevive à Trotterização?

Motivação (ressalva R1 do handoff). A compressão 135->30 e 230->21 foi medida com
exponenciais EXATAS, exp(theta_m * O_m) via expm_multiply, onde O_m é um operador
invariante que é soma de dezenas de excitações elementares. Isso reduz PARÂMETROS, não
necessariamente profundidade de circuito.

Uma implementação em hardware Trotteriza de volta:
    prod_m prod_k exp(theta_m * w_mk * kappa_k)
com kappa_k excitações elementares (kappa^3 = -kappa, portão de forma fechada) e o
parâmetro theta_m AMARRADO entre todas as portas do operador m. O número de portas volta
a ser o do SymUCCSD; o que permanece reduzido é a contagem de parâmetros.

A pergunta: Trotterizar quebra a exatidão de exp(soma). Os ~0,14 mHa sobrevivem?

Decomposição: as elementares kappa_k são ortogonais no produto de Frobenius (conectam pares
distintos de determinantes), então w_mk = <kappa_k, O_m>_F / <kappa_k, kappa_k>_F. O resíduo
||O_m - sum_k w_mk kappa_k|| é reportado como verificação.
"""
import sys, math, time, json
import numpy as np
import scipy.optimize as opt
import scipy.sparse.linalg as sla
from pyscf import gto, scf, ao2mo

HERE = __import__('os').path.dirname(__import__('os').path.abspath(__file__))  # PATCHED: was a hard-coded absolute path
sys.path.insert(0, HERE)
from indep import build_det_basis, SpinFreeHam, kappa_matrix, TrotterAnsatz
from audit import gates_of
from group_compression import (CASES, group, mo_rep, invariant_basis, build_ops, to_str)


def sparse_invariant_basis(Ms, nS):
    """Base ESPARSA do subespaço invariante.

    `np.linalg.eigh` devolve uma base ortonormal arbitrária do autoespaço de autovalor 1,
    que é genericamente DENSA: cada operador invariante vira combinação de centenas de
    excitações elementares, e o circuito Trotterizado explode (18709 portas para o NH3).

    A base natural é {P e_j}: aplicar o projetor aos vetores coordenada. Como
    P = (1/|G|) sum_g Rep(g), cada P e_j tem no máximo |G| termos. Escolhem-se colunas
    independentes por QR com pivotamento.

    Retorna (Bs, Bd) nas mesmas convenções de invariant_basis().
    """
    def pivot(cols, tol=1e-9):
        A = np.array(cols).T
        if A.shape[1] == 0:
            return A
        import scipy.linalg as _sl; Q, R, piv = _sl.qr(A, pivoting=True, mode='economic')
        rank = int((np.abs(np.diag(R)) > tol * max(1.0, abs(R[0, 0]))).sum())
        return A[:, sorted(piv[:rank])]

    Ps = sum(Ms) / len(Ms)
    Bs = pivot([Ps[:, j] for j in range(nS)])
    Pd_full = sum(np.kron(M, M) for M in Ms) / len(Ms)
    cols = []
    for p in range(nS):
        for q in range(p, nS):
            v = np.zeros(nS * nS)
            v[p * nS + q] = 1.0
            v[q * nS + p] = 1.0 if p != q else 0.0
            if p == q:
                v[p * nS + p] = 1.0
            cols.append(Pd_full @ v)
    Bd = pivot(cols)
    return Bs, Bd


def run(tag, atoms, gname):
    t0 = time.time()
    mol = gto.M(atom=to_str(atoms), basis='sto-3g', symmetry=True, verbose=0)
    mf = scf.RHF(mol).run()
    nao, nelec = mol.nao, mol.nelectron; no = nelec // 2
    G, H = group(gname)
    D, _ = mo_rep(mol, mf.mo_coeff, G)
    sing = [(i, a) for i in range(no) for a in range(no, nao)]; nS = len(sing)
    Ms = [np.array([[Dg[j, i] * Dg[b, a] for (i, a) in sing] for (j, b) in sing]) for Dg in D]
    I = np.eye(nS * nS); SW = np.zeros((nS * nS, nS * nS))
    for r_ in range(nS):
        for s_ in range(nS):
            SW[r_ * nS + s_, s_ * nS + r_] = 1.0
    PI = 0.5 * (I + SW)

    dets, index = build_det_basis(nao, no, no); n = len(dets)
    psi0 = np.zeros(n)
    psi0[index[sum((1 << (2*p)) | (1 << (2*p+1)) for p in range(no))]] = 1.0
    h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), nao)
    ham = SpinFreeHam(h1, eri, mol.energy_nuc(), dets, index)
    egs = float(sla.eigsh(sla.LinearOperator((n, n), matvec=ham.matvec, dtype=float),
                          k=1, which='SA', tol=1e-12)[0][0])
    print(f"\n{'='*84}\n{tag}  grupo={gname}   determinantes={n}   E_FCI={egs:.10f}")

    # ---- base de excitações elementares (distintas) ----
    _, gates = gates_of(nao, nelec)
    elems = sorted({(g['occ'], g['virt']) for g in gates})
    Kel = [kappa_matrix(e, dets, index) for e in elems]
    nrm = np.array([float((k.multiply(k)).sum()) for k in Kel])
    print(f"  excitações elementares distintas: {len(elems)}   [{time.time()-t0:.0f}s]",
          flush=True)

    Ep = [ham.E[a][i].tocsr() for (i, a) in sing]
    res = {'molecule': tag, 'group': gname, 'efci': egs}

    for nm, sub, key in (('subgrupo Abeliano', [i for i, R in enumerate(G)
                                                if any(np.allclose(R, h) for h in H)], 'sub'),
                         ('grupo pleno', list(range(len(G))), 'full')):
        Bs, Bd = sparse_invariant_basis([Ms[k] for k in sub], nS)
        ops = build_ops(Bs, Bd, Ep, nS)
        M = len(ops)
        # ---- projeta cada operador invariante na base elementar ----
        pidx, coef, kap, resid = [], [], [], 0.0
        for m, O in enumerate(ops):
            w = np.array([float(O.multiply(Kel[k]).sum()) / nrm[k] for k in range(len(Kel))])
            rec = sum(w[k] * Kel[k] for k in np.nonzero(np.abs(w) > 1e-12)[0])
            resid = max(resid, float(abs(O - rec).max()) if (O - rec).nnz else 0.0)
            for k in np.nonzero(np.abs(w) > 1e-12)[0]:
                pidx.append(m); coef.append(float(w[k])); kap.append(Kel[k])
        print(f"  {nm:20s}: {M} parâmetros, {len(kap)} portas elementares, "
              f"resíduo máx da decomposição = {resid:.2e}")

        # ---- referência: exponenciais exatas ----
        def eg_exact(th):
            v = psi0.copy()
            for k in range(M):
                v = sla.expm_multiply(th[k] * ops[k], v)
            Hv = ham.matvec(v); e = float(v @ Hv); w_ = 2.0 * Hv; g = np.zeros(M)
            for k in range(M - 1, -1, -1):
                g[k] = float(w_ @ ops[k].dot(v))
                v = sla.expm_multiply(-th[k] * ops[k], v)
                w_ = sla.expm_multiply(-th[k] * ops[k], w_)
            return e, g
        r0 = opt.minimize(eg_exact, np.zeros(M), jac=True, method='L-BFGS-B',
                          options=dict(maxiter=4000, maxfun=8000, ftol=1e-16, gtol=1e-11))
        print(f"    {'exponenciais exatas':32s} P={M:4d}  err={(r0.fun-egs)*1000:8.4f} mHa",
              flush=True)

        # ---- Trotterizado, parâmetros amarrados, em duas ordens ----
        out = {'nparam': M, 'ngates': len(kap), 'residual': resid,
               'err_exact_mha': (r0.fun - egs) * 1000}
        for oname, order in (('ordem natural', list(range(len(kap)))),
                             ('ordem invertida', list(range(len(kap)))[::-1])):
            ans = TrotterAnsatz([kap[i] for i in order], [pidx[i] for i in order],
                                [coef[i] for i in order], M, psi0, ham)
            r1 = opt.minimize(lambda x: ans.energy_grad(x), np.zeros(M), jac=True,
                              method='L-BFGS-B',
                              options=dict(maxiter=20000, maxfun=60000,
                                           ftol=1e-16, gtol=1e-11))
            print(f"    {'Trotter amarrado, '+oname:32s} P={M:4d}  "
                  f"err={(r1.fun-egs)*1000:8.4f} mHa", flush=True)
            out['err_trotter_' + oname.split()[-1] + '_mha'] = (r1.fun - egs) * 1000
        res[key] = out
    print(f"  [{time.time()-t0:.0f}s]")
    return res


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'nh3'
    keys = ['nh3', 'ch4'] if which == 'all' else [which]
    out = []
    for k in keys:
        tag, atoms, gname, _ = CASES[k]
        out.append(run(tag, atoms, gname))
    json.dump(out, open(f'{HERE}/e1_trotter.json', 'w'), indent=2)
    print("\nsalvo em e1_trotter.json")
