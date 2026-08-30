"""
he_pipeline.py — reprodução ponta a ponta do pipeline de He et al., em código independente.

tools.get_mole_info      -> openfermionpyscf.run_pyscf -> gto.Mole(symmetry=False)
                            => Hamiltoniano e operadores de excitação vivem na base de
                               orbitais moleculares do RHF SEM simetria.
algorithms.get_symuccsd  -> um SEGUNDO gto.M(atom=geo, symmetry=True) + RHF, separado
                            => os rótulos de irrep vêm de OUTRO conjunto de orbitais,
                               e são aplicados aos operadores POR ÍNDICE.

Em camada degenerada os dois SCF diferem por uma rotação arbitrária interna, e o filtro
passa a remover correlação genuína.

Alvos publicados (JCTC 2026, 22, 6008, Tabela 2):
    NH3  UCCSD 1.09e-4 Ha   SymUCCSD 2.78e-2 Ha
    CH4  UCCSD 2.17e-4 Ha   SymUCCSD 4.08e-2 Ha

NOTA: o erro do SymUCCSD NÃO é reprodutível entre execuções. O LAPACK resolve o subespaço
degenerado num ponto diferente a cada processo, e o erro varia (medidos 37,5 e 28,6 mHa para
o CH4 em duas execuções). Isso é o achado, não um defeito do teste.
"""
import sys, math, json
import numpy as np
import scipy.optimize as opt
import scipy.sparse.linalg as sla
from pyscf import gto, scf, ao2mo, symm

HERE = __import__('os').path.dirname(__import__('os').path.abspath(__file__))  # PATCHED: was a hard-coded absolute path
sys.path.insert(0, HERE)
from indep import build_det_basis, SpinFreeHam, kappa_matrix, TrotterAnsatz
from audit import gates_of, irrep_of


def he_geom_nh3(dist=1.09, angle=107.3):
    th = math.radians(angle); r = dist * math.sin(th / 2); z = -dist * math.cos(th / 2)
    return [['N', [0., 0., 0.]],
            ['H', [r, 0., z]],
            ['H', [r * math.cos(2 * math.pi / 3), r * math.sin(2 * math.pi / 3), z]],
            ['H', [r * math.cos(4 * math.pi / 3), r * math.sin(4 * math.pi / 3), z]]]


def he_geom_ch4(dist=1.09):
    v = [np.array(x, dtype=float) for x in ([1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1])]
    return [['C', [0., 0., 0.]]] + [['H', list(dist * x / np.linalg.norm(x))] for x in v]


def tostr(g):
    return "; ".join(f"{a} {x[0]:.10f} {x[1]:.10f} {x[2]:.10f}" for a, x in g)


def run(label, geom):
    gs = tostr(geom)
    # --- ramo A: o que get_mole_info() realmente constrói (symmetry = False) ---
    molA = gto.M(atom=gs, basis='sto-3g', symmetry=False, verbose=0)
    mfA = scf.RHF(molA).run()
    # --- ramo B: o que get_symuccsd() constrói (symmetry = True) ---
    molB = gto.M(atom=gs, basis='sto-3g', symmetry=True, verbose=0)
    mfB = scf.RHF(molB).run()
    orbsym = np.array(symm.label_orb_symm(molB, molB.irrep_id, molB.symm_orb, mfB.mo_coeff))
    norb, nelec = molA.nao, molA.nelectron; na = nelec // 2

    print(f"\n{'='*84}\n{label}   topgroup={molB.topgroup}  abeliano={molB.groupname}")
    print(f"  RHF(sym=False)={mfA.e_tot:.10f}   RHF(sym=True)={mfB.e_tot:.10f}   "
          f"dif={abs(mfA.e_tot-mfB.e_tot):.2e}")
    print(f"  energias orbitais coincidem: {np.allclose(mfA.mo_energy, mfB.mo_energy, atol=1e-8)}")
    print(f"  rótulos tomados do ramo B: {[int(x) for x in orbsym]}")

    S = molA.intor('int1e_ovlp')
    O = mfA.mo_coeff.T @ S @ mfB.mo_coeff
    print(f"  |<phi^A_p|phi^B_p>| por orbital: {np.round(np.abs(np.diag(O)), 4)}")
    print(f"     (1 = orbitais idênticos; <1 = os dois SCF discordam nessa camada)")

    dets, index = build_det_basis(norb, na, na); n = len(dets)
    psi0 = np.zeros(n)
    psi0[index[sum((1 << (2 * p)) | (1 << (2 * p + 1)) for p in range(na))]] = 1.0
    h1 = mfA.mo_coeff.T @ mfA.get_hcore() @ mfA.mo_coeff
    eri = ao2mo.restore(1, ao2mo.kernel(molA, mfA.mo_coeff), norb)
    ham = SpinFreeHam(h1, eri, molA.energy_nuc(), dets, index)
    egs = float(sla.eigsh(sla.LinearOperator((n, n), matvec=ham.matvec, dtype=float),
                          k=1, which='SA', tol=1e-11)[0][0])
    ehf = float(psi0 @ ham.matvec(psi0))
    npar, gates = gates_of(norb, nelec)
    K = [kappa_matrix((g['occ'], g['virt']), dets, index) for g in gates]
    ans = TrotterAnsatz(K, [g['param'] for g in gates],
                        [g['coef'] for g in gates], npar, psi0, ham)
    sym_ids = sorted({g['param'] for g in gates if irrep_of(g, orbsym) == 0})
    print(f"  determinantes={n}  E_HF={ehf:.10f}  E_FCI={egs:.10f}  "
          f"E_corr={(ehf-egs)*1000:.3f} mHa")
    print(f"  UCCSD {npar} parâmetros ; filtrado {len(sym_ids)}  <-- Tabela II de He et al.")

    def go(active, lbl):
        act = np.array(sorted(active))

        def f(x):
            th = np.zeros(npar); th[act] = x
            e, g = ans.energy_grad(th)
            return e, g[act]
        r = opt.minimize(f, np.zeros(len(act)), jac=True, method='BFGS',
                         options=dict(maxiter=20000, gtol=1e-8))
        print(f"    {lbl:36s} P={len(act):4d}  E={r.fun:.10f}  "
              f"err={(r.fun-egs):.4e} Ha  ({(r.fun-egs)*1000:8.3f} mHa)", flush=True)
        return float(r.fun)

    EF = go(range(npar), "UCCSD (todos os operadores)")
    ES = go(sym_ids, "SymUCCSD (pipeline de He et al.)")
    return dict(label=label, ehf=ehf, efci=egs, npar=npar, nsym=len(sym_ids),
                err_uccsd=EF - egs, err_sym=ES - egs,
                overlap_diag=[float(x) for x in np.abs(np.diag(O))])


if __name__ == '__main__':
    out = [run("NH3/STO-3G  (geometria de He et al.: d=1,09 A, ang=107,3)", he_geom_nh3()),
           run("CH4/STO-3G  (geometria de He et al.: d=1,09 A)", he_geom_ch4())]
    json.dump(out, open(f'{HERE}/he_pipeline.json', 'w'), indent=2)
    print("\nsalvo em he_pipeline.json")
