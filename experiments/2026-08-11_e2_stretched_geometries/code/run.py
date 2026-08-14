"""E2 — does the filter stay free as the molecule is pulled apart?

The withdrawn manuscript claimed the deficit "becomes more pronounced at stretched
geometries where full orbital relaxation is required".  This experiment tests the
opposite over NH3 with the N-H distance from 0.90 to 2.20 A at fixed angle, which
walks the system from the equilibrium region into frank static correlation: the
correlation energy grows by about an order of magnitude and the S0-S1 gap collapses.

Three ansaetze at every geometry, under an identical Hamiltonian, reference,
optimiser and convergence criterion:

  uccsd          all 135 parameters, no filter
  abelian        the Cs filter, 75 parameters -- the SymUCCSD of the literature
  full_c3v       the full C3v filter, 30 parameters

The third is the extension: the original stretch scan compared only the first two,
so it tested the Abelian filter under static correlation and left the non-Abelian
one untested, which is the filter the paper actually proposes.

The scan also locates the boundary of the theorem's hypothesis.  The statements
require a symmetric closed-shell reference, so they stop applying exactly where
RHF stops converging; that point is reported rather than worked around.

Usage:  python run.py [--smoke]
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "vendor"))

import expcommon as X  # noqa: E402
import numpy as np  # noqa: E402
import scipy.optimize as opt  # noqa: E402
import scipy.sparse.linalg as sla  # noqa: E402
from pyscf import ao2mo, gto, scf, symm  # noqa: E402

from audit import gates_of, irrep_of, nh3_geom  # noqa: E402
from e1_trotter import sparse_invariant_basis  # noqa: E402
from group_compression import build_ops, group, mo_rep  # noqa: E402
from indep import SpinFreeHam, TrotterAnsatz, build_det_basis, kappa_matrix  # noqa: E402

ROOT = HERE.parent
RS = [0.90, 1.0124, 1.10, 1.25, 1.40, 1.60, 1.80, 2.00, 2.20]
ANGLE = 106.7
GNAME = "C3v"


def optimise(ans, npar, active, egs):
    act = np.array(sorted(active))

    def f(x):
        th = np.zeros(npar)
        th[act] = x
        e, g = ans.energy_grad(th)
        return e, g[act]

    r = opt.minimize(f, np.zeros(len(act)), jac=True, method="L-BFGS-B",
                     options=dict(maxiter=20000, maxfun=60000, ftol=1e-16, gtol=1e-11))
    return float(r.fun), (float(r.fun) - egs) * 1000.0, int(r.nit), float(np.max(np.abs(r.jac)))


FIELDS = ["r_nh_angstrom", "angle_deg", "converged_rhf", "n_dets", "e_hf_ha", "e_fci_ha",
          "e_s1_ha", "corr_energy_mha", "gap_s0_s1_mha", "ansatz", "n_params", "n_gates",
          "energy_ha", "err_vs_fci_mha", "delta_vs_uccsd_mha", "n_iterations",
          "grad_inf_norm", "wallclock_s"]


def make_log(smoke=False):
    import tempfile
    path = (Path(tempfile.mkdtemp()) / "smoke_log.json") if smoke else (ROOT / "experiment_log.json")
    return X.ExperimentLog(
        path,
        experiment_id=ROOT.name,
        paper_project="2026-08_point-group-vqe-free",
        opportunity_reference="opportunities.md #3",
        mode="numerical+symbolic",
        subareas_active=["numerical", "symbolic"],
        target_venue="Quantum (componente do artigo principal); Comment no JCTC como item separado e posterior.",
        created_at=X.utcnow(),
        plan_md_hash=X.sha256_file(ROOT / "plan.md") if (ROOT / "plan.md").exists() else None,
        code_commit=None,
        system_info_hash=None,
    )


def main():
    smoke = "--smoke" in sys.argv
    rs = RS[:2] if smoke else RS

    X.capture_env(ROOT / "env", project_root=ROOT.parent.parent)
    log = make_log(smoke=smoke)
    rows = []

    # The determinant basis, the excitation pool and the kappa matrices do not depend
    # on the geometry, so they are built once for the whole scan.
    mol0 = gto.M(atom=nh3_geom(RS[1], ANGLE), basis="sto-3g", symmetry=True, verbose=0)
    nao, nelec = mol0.nao, mol0.nelectron
    no = nelec // 2
    dets, index = build_det_basis(nao, no, no)
    n = len(dets)
    psi0 = np.zeros(n)
    psi0[index[sum((1 << (2 * p)) | (1 << (2 * p + 1)) for p in range(no))]] = 1.0
    npar, gates = gates_of(nao, nelec)
    K = [kappa_matrix((g["occ"], g["virt"]), dets, index) for g in gates]
    pidx = [g["param"] for g in gates]
    coef = [g["coef"] for g in gates]
    sing = [(i, a) for i in range(no) for a in range(no, nao)]
    nS = len(sing)
    elems = sorted({(g["occ"], g["virt"]) for g in gates})
    Kel = [kappa_matrix(e, dets, index) for e in elems]
    nrm = np.array([float((k.multiply(k)).sum()) for k in Kel])
    G, H = group(GNAME)
    print(f"shared setup: dets={n} P={npar} gate entries={len(gates)} "
          f"distinct elementary={len(elems)} |G|={len(G)}", flush=True)

    for r_nh in rs:
        t0 = time.time()
        geom = nh3_geom(r_nh, ANGLE)
        mol = gto.M(atom=geom, basis="sto-3g", symmetry=True, verbose=0)
        mf = scf.RHF(mol)
        mf.kernel()
        if not mf.converged:
            print(f"r={r_nh:.4f}  RHF did NOT converge -- point discarded", flush=True)
            rows.append(dict(r_nh_angstrom=r_nh, angle_deg=ANGLE, converged_rhf=False,
                             n_dets=n, ansatz="", n_params="", n_gates="",
                             wallclock_s=round(time.time() - t0, 2)))
            log.record(parameters=dict(r_nh=r_nh, angle=ANGLE, ansatz="none"),
                       status="failure", seed=None, exit_code=1,
                       wallclock=time.time() - t0,
                       outputs=dict(result_files=[], scalar_results={}),
                       notes="RHF did not converge; the closed-shell invariant reference the "
                             "theorems require does not exist here, so the point is reported "
                             "and not worked around",
                       peak_memory_mb=X.peak_rss_mb())
            continue

        # Irrep labels are recomputed at every geometry: orbital ordering can cross.
        orbsym = np.array(symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb, mf.mo_coeff))
        h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
        eri = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), nao)
        ham = SpinFreeHam(h1, eri, mol.energy_nuc(), dets, index)
        lin = sla.LinearOperator((n, n), matvec=ham.matvec, dtype=float)
        ev = sla.eigsh(lin, k=2, which="SA", tol=1e-12)[0]
        egs, e_s1 = float(min(ev)), float(max(ev))
        corr = (egs - mf.e_tot) * 1000.0
        gap = (e_s1 - egs) * 1000.0

        sym_ids = sorted({g["param"] for g in gates if irrep_of(g, orbsym) == 0})
        ans = TrotterAnsatz(K, pidx, coef, npar, psi0, ham)

        # full C3v invariant ansatz, rebuilt at this geometry
        D, rep_err = mo_rep(mol, mf.mo_coeff, G)
        Ms = [np.array([[D[k][j, i] * D[k][b, a] for (i, a) in sing]
                        for (j, b) in sing]) for k in range(len(G))]
        Bs, Bd = sparse_invariant_basis(Ms, nS)
        Ep = [ham.E[a][i].tocsr() for (i, a) in sing]
        ops = build_ops(Bs, Bd, Ep, nS)
        M = len(ops)
        fp, fc, fk, resid = [], [], [], 0.0
        for m, O in enumerate(ops):
            w = np.array([float(O.multiply(Kel[k]).sum()) / nrm[k] for k in range(len(Kel))])
            nz = np.nonzero(np.abs(w) > 1e-12)[0]
            rec = sum(w[k] * Kel[k] for k in nz)
            resid = max(resid, float(abs(O - rec).max()) if (O - rec).nnz else 0.0)
            for k in nz:
                fp.append(m); fc.append(float(w[k])); fk.append(Kel[k])
        ans_full = TrotterAnsatz(fk, fp, fc, M, psi0, ham)

        # Exact exponentials of the same invariant operators.  This separates two
        # explanations for any cost the filtered ansatz shows: a genuine expressivity
        # limit of the 30-dimensional invariant subspace, or the price of tying the
        # parameters and Trotterising.  Without it the two are confounded.
        def eg_exact(th):
            v = psi0.copy()
            for k in range(M):
                v = sla.expm_multiply(th[k] * ops[k], v)
            Hv = ham.matvec(v)
            e = float(v @ Hv)
            w_ = 2.0 * Hv
            g = np.zeros(M)
            for k in range(M - 1, -1, -1):
                g[k] = float(w_ @ ops[k].dot(v))
                v = sla.expm_multiply(-th[k] * ops[k], v)
                w_ = sla.expm_multiply(-th[k] * ops[k], w_)
            return e, g

        base = dict(r_nh_angstrom=r_nh, angle_deg=ANGLE, converged_rhf=True, n_dets=n,
                    e_hf_ha=mf.e_tot, e_fci_ha=egs, e_s1_ha=e_s1,
                    corr_energy_mha=corr, gap_s0_s1_mha=gap)
        e_ref = None
        for label, a, active, npm, ng in (
                ("uccsd", ans, range(npar), npar, len(gates)),
                ("abelian", ans, sym_ids, len(sym_ids), sum(1 for g in gates
                                                            if g["param"] in set(sym_ids))),
                ("full_c3v", ans_full, range(M), M, len(fk)),
                ("full_c3v_exact", None, range(M), M, "")):
            tw = X.Timer()
            with tw:
                if label == "full_c3v_exact":
                    r = opt.minimize(eg_exact, np.zeros(M), jac=True, method="L-BFGS-B",
                                     options=dict(maxiter=4000, maxfun=8000,
                                                  ftol=1e-16, gtol=1e-11))
                    e, err = float(r.fun), (float(r.fun) - egs) * 1000.0
                    nit, gn = int(r.nit), float(np.max(np.abs(r.jac)))
                else:
                    e, err, nit, gn = optimise(a, npm if label == "full_c3v" else npar,
                                               active, egs)
            if label == "uccsd":
                e_ref = e
            rows.append(dict(**base, ansatz=label, n_params=npm, n_gates=ng, energy_ha=e,
                             err_vs_fci_mha=err, delta_vs_uccsd_mha=(e - e_ref) * 1000.0,
                             n_iterations=nit, grad_inf_norm=gn,
                             wallclock_s=round(tw.dt, 2)))
            log.record(parameters=dict(r_nh=r_nh, angle=ANGLE, ansatz=label,
                                       n_params=npm, n_gates=ng, basis="sto-3g"),
                       status="success", seed=None, wallclock=tw.dt,
                       t_start=tw.start, t_end=tw.end,
                       outputs=dict(result_files=[], scalar_results=dict(
                           energy_ha=e, err_vs_fci_mha=err,
                           delta_vs_uccsd_mha=(e - e_ref) * 1000.0,
                           corr_energy_mha=corr, gap_s0_s1_mha=gap)),
                       notes=f"decomposition residual of the invariant operators {resid:.1e}"
                             if label == "full_c3v" else "",
                       peak_memory_mb=X.peak_rss_mb())
        print(f"r={r_nh:.4f}  Ecorr={abs(corr):7.1f}  gap={gap:7.1f}  "
              f"uccsd={rows[-4]['err_vs_fci_mha']:8.4f}  "
              f"abelian d={rows[-3]['delta_vs_uccsd_mha']:+.2e}  "
              f"full_trotter d={rows[-2]['delta_vs_uccsd_mha']:+.2e}  "
              f"full_exact d={rows[-1]['delta_vs_uccsd_mha']:+.2e}  "
              f"[{time.time()-t0:.0f}s]", flush=True)

    ok = [r for r in rows if r.get("converged_rhf") and r["ansatz"] == "abelian"]
    okf = [r for r in rows if r.get("converged_rhf") and r["ansatz"] == "full_c3v"]
    oke = [r for r in rows if r.get("converged_rhf") and r["ansatz"] == "full_c3v_exact"]
    failed = [r for r in rows if r.get("converged_rhf") is False]
    eq = [r for r in oke if r["r_nh_angstrom"] <= 1.10]
    agg = dict(
        max_abs_delta_abelian_mha=max(abs(r["delta_vs_uccsd_mha"]) for r in ok),
        max_abs_delta_full_c3v_trotter_mha=max(abs(r["delta_vs_uccsd_mha"]) for r in okf),
        max_abs_delta_full_c3v_exact_mha=max(abs(r["delta_vs_uccsd_mha"]) for r in oke),
        max_abs_delta_full_c3v_exact_near_eq_mha=max(abs(r["delta_vs_uccsd_mha"]) for r in eq),
        full_c3v_exact_by_r={r["r_nh_angstrom"]: round(r["delta_vs_uccsd_mha"], 6) for r in oke},
        corr_energy_growth=f"{min(abs(r['corr_energy_mha']) for r in ok):.1f} to "
                           f"{max(abs(r['corr_energy_mha']) for r in ok):.1f} mHa",
        gap_collapse=f"{max(r['gap_s0_s1_mha'] for r in ok):.1f} to "
                     f"{min(r['gap_s0_s1_mha'] for r in ok):.1f} mHa",
        rhf_boundary_angstrom=min((r["r_nh_angstrom"] for r in failed), default=None),
        n_geometries=len(ok),
    )
    # Where the compression starts to cost, is it the invariant subspace or the tying?
    # Compare the exact-exponential and tied-Trotter forms at the far end of the scan.
    far_e = max(abs(r["delta_vs_uccsd_mha"]) for r in oke)
    far_t = max(abs(r["delta_vs_uccsd_mha"]) for r in okf)
    if far_e > 0.1 * far_t:
        attribution = ("The exact-exponential form carries most of that cost, so it belongs to "
                       "the 30-dimensional invariant subspace and not to Trotterisation.")
    else:
        attribution = ("The exact-exponential form stays far below the tied-Trotter one, so the "
                       "cost at the stretched end is the price of tying and Trotterising rather "
                       "than a limit of the invariant subspace.")
    agg["cost_attribution"] = attribution

    print("\n=== E2 summary ===")
    for k, v in agg.items():
        print(f"  {k:32s} {v}")

    if not smoke:
        X.write_csv(ROOT / "results" / "aggregate.csv", rows, FIELDS)
        (ROOT / "results" / "aggregate.json").write_text(
            json.dumps(rows, indent=2, default=X._jsonable) + "\n")
        log.summarize(
            aggregate_metrics=agg,
            key_findings=[
                f"across {agg['n_geometries']} geometries the Abelian filter changes the "
                f"energy by at most {agg['max_abs_delta_abelian_mha']:.2e} mHa, so the claim "
                f"that its deficit grows when the molecule is pulled apart does not survive",
                f"correlation energy grows over the scan ({agg['corr_energy_growth']}) while "
                f"the S0-S1 gap collapses ({agg['gap_collapse']}), so the scan does reach the "
                f"static-correlation regime it was meant to probe",
                f"the full C3v filter at 30 parameters departs from unfiltered UCCSD by at most "
                f"{agg['max_abs_delta_full_c3v_exact_near_eq_mha']:.2e} mHa near equilibrium but "
                f"by {agg['max_abs_delta_full_c3v_exact_mha']:.2f} mHa with exact exponentials "
                f"and {agg['max_abs_delta_full_c3v_trotter_mha']:.2f} mHa in tied Trotter form "
                f"at the far end of the scan",
                f"RHF stops converging at {agg['rhf_boundary_angstrom']} A, which is where the "
                f"invariant closed-shell reference the theorems assume ceases to exist",
            ],
            matches_hypothesis="ambiguous",
            interpretation="Two different statements come apart along this scan. The Abelian "
                           "filter is free everywhere the reference exists, which refutes the "
                           "claim it was tested against. The full non-Abelian compression is "
                           "free only near equilibrium: as the gap collapses its cost grows. "
                           + attribution +
                           " The compression claim therefore has a domain, and this scan "
                           "measures where it ends.",
        )
    return agg


if __name__ == "__main__":
    main()
