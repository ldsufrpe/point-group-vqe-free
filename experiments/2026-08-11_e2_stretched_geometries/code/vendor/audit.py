"""
Independent audit: does the point-group (Abelian-subgroup) filter of SymUCCSD
actually cost energy for NH3/STO-3G (C3v) and CH4/STO-3G (Td)?

Two ansatz conventions are compared under IDENTICAL Hamiltonian, reference,
optimiser and convergence criteria:

  convention 'of'  : OpenFermion singlet-UCCSD parameterisation.
                     Reproduces He et al.'s counts exactly
                     (NH3 135/75 params, 315/163 operators;
                      CH4 230/65 params, 560/146 operators).

  convention 'ds'  : the enumeration used in da Silva & Santos's public code
                     (nh3_symuccsd_confinement.py, build_pool):
                     singles spin-complemented, same-spin doubles
                     combinations(occ,2) x combinations(vir,2) spin-complemented,
                     opposite-spin doubles cwr(occ,2) x cwr(vir,2), ONE spin case.
"""
import sys, math, time, json, itertools
import numpy as np
import scipy.optimize as opt
import scipy.sparse.linalg as sla
from itertools import combinations, combinations_with_replacement as cwr
from openfermion.circuits import uccsd_singlet_generator, uccsd_singlet_paramsize

HERE = __import__('os').path.dirname(__import__('os').path.abspath(__file__))  # PATCHED: was a hard-coded absolute path
sys.path.insert(0, HERE)
from indep import (build_det_basis, apply_string, SpinFreeHam, make_molecule,
                   kappa_matrix, TrotterAnsatz)


def nh3_geom(r=1.0124, ang_deg=106.7):
    a = math.radians(ang_deg)
    ct = math.sqrt((math.cos(a) + 0.5) / 1.5); st = math.sqrt(1 - ct * ct)
    at = [('N', (0.0, 0.0, 0.0))]
    for k in range(3):
        p = 2 * math.pi * k / 3
        at.append(('H', (r * st * math.cos(p), r * st * math.sin(p), r * ct)))
    return "; ".join(f"{s} {x:.10f} {y:.10f} {z:.10f}" for s, (x, y, z) in at)


NH3_DS = """N 0.0000 0.0000 0.1173; H 0.0000 0.9377 -0.2738;
            H 0.8126 -0.4689 -0.2738; H -0.8126 -0.4689 -0.2738"""


def ch4_geom(r=1.087):
    c = r / math.sqrt(3.0)
    return (f"C 0 0 0; H {c} {c} {c}; H {c} {-c} {-c}; "
            f"H {-c} {c} {-c}; H {-c} {-c} {c}")


# ---------------- convention 'of' : OpenFermion singlet UCCSD ----------------
def gates_of(norb, nelec):
    nq, nocc_so = 2 * norb, nelec
    npar = uccsd_singlet_paramsize(nq, nelec)
    gates = []
    for k in range(npar):
        amps = np.zeros(npar); amps[k] = 1.0
        T = uccsd_singlet_generator(amps, nq, nelec)
        acc = {}
        for term, coef in T.terms.items():
            if abs(coef) < 1e-12:
                continue
            cre = [i for i, d in term if d == 1]
            ann = [i for i, d in term if d == 0]
            if all(i >= nocc_so for i in cre) and all(i < nocc_so for i in ann):
                occ, virt, sgn = tuple(sorted(ann)), tuple(sorted(cre)), 1.0
                probe = sum(1 << i for i in occ)
                can = [('c', v) for v in reversed(virt)] + [('a', o) for o in occ]
            else:
                occ, virt, sgn = tuple(sorted(cre)), tuple(sorted(ann)), -1.0
                probe = sum(1 << i for i in virt)
                can = [('c', o) for o in reversed(occ)] + [('a', v) for v in virt]
            raw = [('c' if d == 1 else 'a', i) for i, d in term]
            _, sr = apply_string(probe, raw)
            _, sc = apply_string(probe, can)
            assert sr and sc
            acc[(occ, virt)] = acc.get((occ, virt), 0.0) + sgn * coef * sr * sc
        for (occ, virt), c in acc.items():
            if abs(c) > 1e-12:
                gates.append(dict(param=k, occ=occ, virt=virt, coef=float(c)))
    return npar, gates


# ------------- convention 'ds' : da Silva & Santos build_pool ----------------
def gates_ds(norb, nelec):
    """spin-orbital index here: alpha p -> 2p, beta p -> 2p+1 (our convention).
       Their code uses p / p+N; the operator content is identical."""
    nocc = nelec // 2
    occ, vir = list(range(nocc)), list(range(nocc, norb))
    A = lambda p: 2 * p
    B = lambda p: 2 * p + 1
    gates, k = [], 0
    for i in occ:                                     # singles (spin-complemented)
        for a in vir:
            gates.append(dict(param=k, occ=(A(i),), virt=(A(a),), coef=1.0,
                              spatial=((i, a),)))
            gates.append(dict(param=k, occ=(B(i),), virt=(B(a),), coef=1.0,
                              spatial=((i, a),)))
            k += 1
    for i, j in combinations(occ, 2):                 # same-spin doubles
        for a, b in combinations(vir, 2):
            gates.append(dict(param=k, occ=tuple(sorted((A(i), A(j)))),
                              virt=tuple(sorted((A(a), A(b)))), coef=1.0,
                              spatial=((i, a), (j, b))))
            gates.append(dict(param=k, occ=tuple(sorted((B(i), B(j)))),
                              virt=tuple(sorted((B(a), B(b)))), coef=1.0,
                              spatial=((i, a), (j, b))))
            k += 1
    for i, j in cwr(occ, 2):                          # opposite-spin doubles
        for a, b in cwr(vir, 2):
            gates.append(dict(param=k, occ=tuple(sorted((A(i), B(j)))),
                              virt=tuple(sorted((A(a), B(b)))), coef=1.0,
                              spatial=((i, a), (j, b))))
            k += 1
    return k, gates


def irrep_of(g, orbsym):
    r = 0
    for s in g['occ'] + g['virt']:
        r ^= int(orbsym[s // 2])
    return r


def prep(name, geom, conv, symmetry=True):
    M = make_molecule(geom, symmetry=symmetry)
    norb, na = M['norb'], M['na']
    dets, index = build_det_basis(norb, na, na)
    ham = SpinFreeHam(M['h1'], M['eri'], M['ecore'], dets, index)
    hf = sum((1 << (2 * p)) | (1 << (2 * p + 1)) for p in range(na))
    psi0 = np.zeros(len(dets)); psi0[index[hf]] = 1.0
    npar, gates = (gates_of if conv == 'of' else gates_ds)(norb, M['nelec'])
    K = [kappa_matrix((g['occ'], g['virt']), dets, index) for g in gates]
    ans = TrotterAnsatz(K, [g['param'] for g in gates],
                        [g['coef'] for g in gates], npar, psi0, ham)
    symf = [irrep_of(g, M['orbsym']) == 0 for g in gates]
    sym_ids = sorted({g['param'] for g, f in zip(gates, symf) if f})
    for g, f in zip(gates, symf):
        assert f == (g['param'] in sym_ids), "a parameter mixes irreps"
    return M, dets, index, ham, psi0, npar, gates, symf, sym_ids, ans


def optimize(ans, npar, active, x0, label, egs, maxiter=20000):
    act = np.array(sorted(active))

    def f(x):
        th = np.zeros(npar); th[act] = x
        e, g = ans.energy_grad(th)
        return e, g[act]
    r = opt.minimize(f, x0[act], jac=True, method='L-BFGS-B',
                     options=dict(maxiter=maxiter, maxfun=maxiter * 3,
                                  ftol=1e-16, gtol=1e-12))
    th = np.zeros(npar); th[act] = r.x
    print(f"    {label:38s} P={len(act):4d} E={r.fun:.10f} "
          f"err={ (r.fun-egs)*1000:9.4f} mHa  |g|={np.max(np.abs(r.jac)):.1e} "
          f"nit={r.nit}", flush=True)
    return th, r.fun


def main(mol):
    t0 = time.time()
    if mol == 'nh3':
        name, geom = "NH3/STO-3G", nh3_geom()
    elif mol == 'nh3ds':
        name, geom = "NH3/STO-3G (da Silva&Santos geometry)", NH3_DS
    else:
        name, geom = "CH4/STO-3G", ch4_geom()

    out = {}
    for conv in ('of', 'ds'):
        M, dets, index, ham, psi0, npar, gates, symf, sym_ids, ans = prep(name, geom, conv)
        n = len(dets)
        if conv == 'of':
            lo = sla.LinearOperator((n, n), matvec=ham.matvec, dtype=float)
            egs = sla.eigsh(lo, k=1, which='SA', tol=1e-12)[0][0]
            ehf = float(psi0 @ ham.matvec(psi0))
            print(f"\n{'='*86}\n{name}   topgroup={M['topgroup']}  "
                  f"Abelian subgroup used={M['groupname']}")
            print(f"  dets={n}  E(HF)={ehf:.10f} [pyscf {M['ehf']:.10f}]  "
                  f"E(FCI)={egs:.10f} [pyscf {M['efci']:.10f}]")
            print(f"  E_corr = {(ehf-egs)*1000:.3f} mHa    orbsym={list(M['orbsym'])}")
            out['ehf'], out['efci'] = ehf, egs
            EGS = egs
        egs = EGS
        allp = set(range(npar))
        print(f"\n  --- convention '{conv}' : {npar} parameters, {len(gates)} "
              f"elementary excitation operators ---")
        print(f"      symmetric subset: {len(sym_ids)} parameters, "
              f"{sum(symf)} elementary operators")

        rng = np.random.default_rng(0)
        th = rng.normal(0, 0.05, npar)
        e, g = ans.energy_grad(th)
        ii = rng.choice(npar, 4, replace=False)
        num = [ (ans.energy(th + np.eye(npar)[i]*1e-5) -
                 ans.energy(th - np.eye(npar)[i]*1e-5))/2e-5 for i in ii]
        print(f"      gradient check: max|analytic-FD| = "
              f"{np.max(np.abs(g[ii]-np.array(num))):.2e}")

        z = np.zeros(npar)
        thF, EF = optimize(ans, npar, allp, z, f"UCCSD    [{conv}] init 0", egs)
        thS, ES = optimize(ans, npar, sym_ids, z, f"SymUCCSD [{conv}] init 0", egs)
        bestF, bestS = EF, ES
        for seed in (1, 2):
            r2 = np.random.default_rng(seed); x0 = r2.normal(0, 0.15, npar)
            _, e1 = optimize(ans, npar, allp, x0, f"UCCSD    [{conv}] rand{seed}", egs)
            _, e2 = optimize(ans, npar, sym_ids, x0, f"SymUCCSD [{conv}] rand{seed}", egs)
            bestF, bestS = min(bestF, e1), min(bestS, e2)

        asym = sorted(allp - set(sym_ids))
        thP = thF.copy(); thP[asym] = 0.0
        EP = ans.energy(thP)
        print(f"      full optimum: max|theta| on the {len(asym)} REMOVED params "
              f"= {np.max(np.abs(thF[asym])):.3e}")
        print(f"      E(full optimum, removed params zeroed) = {EP:.10f}  "
              f"err={(EP-egs)*1000:.4f} mHa")
        _, gF = ans.energy_grad(thS)
        print(f"      at the SymUCCSD optimum: max|dE/dtheta| over REMOVED "
              f"directions = {np.max(np.abs(gF[asym])):.2e}")

        out[conv] = dict(npar=npar, npool=len(gates), nsym=len(sym_ids),
                         npool_sym=int(sum(symf)),
                         E_uccsd=bestF, err_uccsd_mha=(bestF - egs) * 1000,
                         E_sym=bestS, err_sym_mha=(bestS - egs) * 1000,
                         E_proj=EP, err_proj_mha=(EP - egs) * 1000,
                         max_removed_theta=float(np.max(np.abs(thF[asym]))))
    json.dump(out, open(f'{HERE}/audit_{mol}.json', 'w'), indent=2)
    print(f"\n[{time.time()-t0:.0f}s] saved audit_{mol}.json")


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'nh3')
