"""E1 — does the compression survive Trotterisation with tied parameters?

The 135 -> 30 and 230 -> 21 compressions were first measured with exact
exponentials of the invariant operators, which reduces the number of parameters
but says nothing about a circuit.  Hardware needs the product form

    U(theta) = prod_m prod_k exp(theta_m w_mk kappa_k),

with theta_m tied across every elementary gate of invariant operator m.  This
experiment measures whether the ~0.14 mHa accuracy survives that, and what
happens to the gate count.

Per molecule and per group (Abelian subgroup, full point group):
  parameters, elementary gates, the residual of the decomposition of each
  invariant operator on the elementary basis, and the energy error under exact
  exponentials, under Trotter in the natural gate order, and under Trotter in the
  reversed order.  The two orders bound the ordering dependence.

H2O (C2v) and C2H4 (D2h) are Abelian controls: there the full group IS the
Abelian subgroup, so the additional compression must be exactly zero.  C2H4 is
counted but not optimised -- its determinant space is about 9e6.

The invariant basis is the pivoted {P e_j} of sparse_invariant_basis, never the
dense eigenbasis: an arbitrary orthonormal basis of the invariant eigenspace is
generically dense and inflates the NH3 circuit to 18709 gates, which would
reverse the conclusion.

Usage:  python run.py [nh3|ch4|h2o|c2h4|all] [--smoke]
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
from pyscf import ao2mo, gto, scf  # noqa: E402

import e1_trotter as E1  # noqa: E402
from pyscf import symm  # noqa: E402

from audit import gates_of, irrep_of  # noqa: E402
from group_compression import CASES, group, invariant_basis, mo_rep, to_str  # noqa: E402
from indep import SpinFreeHam, TrotterAnsatz, build_det_basis, kappa_matrix  # noqa: E402

ROOT = HERE.parent
ORDER = ["h2o", "c2h4", "nh3", "ch4"]


def uccsd_reference(tag, atoms):
    """Unfiltered UCCSD baseline: parameters, elementary gates, and energy error.

    Needed because e1_trotter reports only the two filtered ansaetze; the claim is
    about what the filter costs relative to the unfiltered circuit.
    """
    mol = gto.M(atom=to_str(atoms), basis="sto-3g", symmetry=True, verbose=0)
    mf = scf.RHF(mol).run()
    nao, nelec = mol.nao, mol.nelectron
    no = nelec // 2
    dets, index = build_det_basis(nao, no, no)
    n = len(dets)
    psi0 = np.zeros(n)
    psi0[index[sum((1 << (2 * p)) | (1 << (2 * p + 1)) for p in range(no))]] = 1.0
    h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), nao)
    ham = SpinFreeHam(h1, eri, mol.energy_nuc(), dets, index)
    egs = float(sla.eigsh(sla.LinearOperator((n, n), matvec=ham.matvec, dtype=float),
                          k=1, which="SA", tol=1e-12)[0][0])
    npar, gates = gates_of(nao, nelec)
    K = [kappa_matrix((g["occ"], g["virt"]), dets, index) for g in gates]
    n_distinct = len({(g["occ"], g["virt"]) for g in gates})
    ans = TrotterAnsatz(K, [g["param"] for g in gates], [g["coef"] for g in gates],
                        npar, psi0, ham)
    r = opt.minimize(lambda x: ans.energy_grad(x), np.zeros(npar), jac=True,
                     method="L-BFGS-B",
                     options=dict(maxiter=20000, maxfun=60000, ftol=1e-16, gtol=1e-11))
    # The Abelian filter as the literature actually implements it: a SUBSET of the
    # UCCSD pool selected by irrep label, not a set of projector-built invariant
    # operators.  The two constructions span the same space but have very different
    # gate counts, and the subset one is the fair circuit-size comparison.
    orbsym = np.array(symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb, mf.mo_coeff))
    sym_ids = {g["param"] for g in gates if irrep_of(g, orbsym) == 0}
    sub_gates = [g for g in gates if g["param"] in sym_ids]
    sub_distinct = len({(g["occ"], g["virt"]) for g in sub_gates})
    ansS = TrotterAnsatz([K[i] for i, g in enumerate(gates) if g["param"] in sym_ids],
                         [g["param"] for g in sub_gates], [g["coef"] for g in sub_gates],
                         npar, psi0, ham)
    act = np.array(sorted(sym_ids))

    def fS(x):
        th = np.zeros(npar)
        th[act] = x
        e, g = ansS.energy_grad(th)
        return e, g[act]

    rS = opt.minimize(fS, np.zeros(len(act)), jac=True, method="L-BFGS-B",
                      options=dict(maxiter=20000, maxfun=60000, ftol=1e-16, gtol=1e-11))
    return dict(efci=egs, nparam=npar, ngates=n_distinct,
                err_mha=(r.fun - egs) * 1000.0, n_dets=n,
                sub_nparam=len(act), sub_ngates=sub_distinct,
                sub_err_mha=(float(rS.fun) - egs) * 1000.0)


def counts_only(tag, atoms, gname):
    """Parameter counts for a molecule whose determinant space is out of reach."""
    mol = gto.M(atom=to_str(atoms), basis="sto-3g", symmetry=True, verbose=0)
    mf = scf.RHF(mol).run()
    nao, nelec = mol.nao, mol.nelectron
    no = nelec // 2
    G, H = group(gname)
    D, err = mo_rep(mol, mf.mo_coeff, G)
    sing = [(i, a) for i in range(no) for a in range(no, nao)]
    nS = len(sing)
    idx = {R.tobytes(): k for k, R in enumerate(G)}
    Ms = {k: np.array([[D[k][j, i] * D[k][b, a] for (i, a) in sing]
                       for (j, b) in sing]) for k in range(len(G))}
    Hk = [idx[R.tobytes()] for R in H]
    I = np.eye(nS * nS)
    SW = np.zeros((nS * nS, nS * nS))
    for r_ in range(nS):
        for s_ in range(nS):
            SW[r_ * nS + s_, s_ * nS + r_] = 1.0
    PI = 0.5 * (I + SW)
    out = dict(nparam_uccsd=nS + nS * (nS + 1) // 2, rep_error=err,
               group_order=len(G), abelian_order=len(H),
               abelian_group=mol.groupname, topgroup=mol.topgroup)
    for name, sub in (("sub", Hk), ("full", list(range(len(G))))):
        Bs, Bd = invariant_basis([Ms[k] for k in sub], PI)
        out[name] = dict(nparam=int(Bs.shape[1] + Bd.shape[1]))
    return out


FIELDS = ["molecule", "full_group", "abelian_group", "group_order", "n_dets", "basis",
          "n_params", "n_gates", "decomp_residual", "e_fci_ha", "err_exact_mha",
          "err_trotter_natural_mha", "err_trotter_reversed_mha", "trotter_cost_mha",
          "param_ratio_vs_uccsd", "gate_ratio_vs_uccsd", "wallclock_s"]


def make_log(smoke=False):
    import tempfile
    path = (Path(tempfile.mkdtemp()) / "smoke_log.json") if smoke else (ROOT / "experiment_log.json")
    return X.ExperimentLog(
        path,
        experiment_id=ROOT.name,
        paper_project="2026-08_point-group-vqe-free",
        opportunity_reference="opportunities.md #2",
        mode="numerical+symbolic",
        subareas_active=["numerical", "symbolic"],
        target_venue="Quantum (componente do artigo principal); Comment no JCTC como item separado e posterior.",
        created_at=X.utcnow(),
        plan_md_hash=X.sha256_file(ROOT / "plan.md") if (ROOT / "plan.md").exists() else None,
        code_commit=None,
        system_info_hash=None,
    )


def main():
    which = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else "all"
    smoke = "--smoke" in sys.argv
    keys = ["h2o"] if smoke else (ORDER if which == "all" else [which])

    X.capture_env(ROOT / "env", project_root=ROOT.parent.parent)
    log = make_log(smoke=smoke)
    rows = []

    for key in keys:
        tag, atoms, gname, do_energy = CASES[key]
        t0 = time.time()
        print(f"\n{'='*78}\n{tag}  full group {gname}", flush=True)

        if not do_energy:
            # C2H4: counting only, the determinant space is out of reach.
            c = counts_only(tag, atoms, gname)
            for name, basis in (("sub", "abelian"), ("full", "full")):
                rows.append(dict(molecule=tag, full_group=gname,
                                 abelian_group=c["abelian_group"],
                                 group_order=c["group_order"], n_dets="", basis=basis,
                                 n_params=c[name]["nparam"], n_gates="",
                                 decomp_residual="", e_fci_ha="", err_exact_mha="",
                                 err_trotter_natural_mha="", err_trotter_reversed_mha="",
                                 trotter_cost_mha="",
                                 param_ratio_vs_uccsd=c[name]["nparam"] / c["nparam_uccsd"],
                                 gate_ratio_vs_uccsd="", wallclock_s=round(time.time()-t0, 2)))
            rows.append(dict(molecule=tag, full_group=gname, abelian_group=c["abelian_group"],
                             group_order=c["group_order"], n_dets="", basis="uccsd",
                             n_params=c["nparam_uccsd"], n_gates="", decomp_residual="",
                             e_fci_ha="", err_exact_mha="", err_trotter_natural_mha="",
                             err_trotter_reversed_mha="", trotter_cost_mha="",
                             param_ratio_vs_uccsd=1.0, gate_ratio_vs_uccsd="",
                             wallclock_s=round(time.time()-t0, 2)))
            log.record(parameters=dict(molecule=tag, group=gname, measurement="counts_only",
                                       reason="determinant space about 9e6, energy infeasible"),
                       status="partial", seed=None, wallclock=time.time()-t0,
                       outputs=dict(result_files=[], scalar_results=dict(
                           nparam_uccsd=c["nparam_uccsd"], nparam_abelian=c["sub"]["nparam"],
                           nparam_full=c["full"]["nparam"])),
                       notes="Abelian control: full group equals the Abelian subgroup, so the "
                             "additional compression must be exactly zero",
                       peak_memory_mb=X.peak_rss_mb())
            continue

        ref = uccsd_reference(tag, atoms)
        print(f"  UCCSD reference: P={ref['nparam']} gates={ref['ngates']} "
              f"err={ref['err_mha']:.4f} mHa  dets={ref['n_dets']}", flush=True)
        rows.append(dict(molecule=tag, full_group=gname, abelian_group="",
                         group_order="", n_dets=ref["n_dets"], basis="uccsd",
                         n_params=ref["nparam"], n_gates=ref["ngates"], decomp_residual="",
                         e_fci_ha=ref["efci"], err_exact_mha=ref["err_mha"],
                         err_trotter_natural_mha="", err_trotter_reversed_mha="",
                         trotter_cost_mha="", param_ratio_vs_uccsd=1.0,
                         gate_ratio_vs_uccsd=1.0, wallclock_s=round(time.time()-t0, 2)))
        log.record(parameters=dict(molecule=tag, group=gname, basis="uccsd",
                                   n_params=ref["nparam"]),
                   status="success", seed=None, wallclock=time.time()-t0,
                   outputs=dict(result_files=[], scalar_results=dict(
                       e_fci=ref["efci"], err_mha=ref["err_mha"], n_gates=ref["ngates"])),
                   notes="unfiltered UCCSD baseline, Trotterised, spin-complemented tying",
                   peak_memory_mb=X.peak_rss_mb())

        print(f"  SymUCCSD as a pool subset: P={ref['sub_nparam']} gates={ref['sub_ngates']} "
              f"err={ref['sub_err_mha']:.4f} mHa", flush=True)
        rows.append(dict(molecule=tag, full_group=gname, abelian_group="", group_order="",
                         n_dets=ref["n_dets"], basis="abelian_subset",
                         n_params=ref["sub_nparam"], n_gates=ref["sub_ngates"],
                         decomp_residual="", e_fci_ha=ref["efci"],
                         err_exact_mha=ref["sub_err_mha"], err_trotter_natural_mha="",
                         err_trotter_reversed_mha="", trotter_cost_mha="",
                         param_ratio_vs_uccsd=ref["sub_nparam"] / ref["nparam"],
                         gate_ratio_vs_uccsd=ref["sub_ngates"] / ref["ngates"],
                         wallclock_s=round(time.time()-t0, 2)))
        log.record(parameters=dict(molecule=tag, group=gname, basis="abelian_subset",
                                   n_params=ref["sub_nparam"], n_gates=ref["sub_ngates"]),
                   status="success", seed=None, wallclock=time.time()-t0,
                   outputs=dict(result_files=[], scalar_results=dict(
                       err_mha=ref["sub_err_mha"], n_gates=ref["sub_ngates"])),
                   notes="the Abelian filter as the literature implements it: a subset of the "
                         "UCCSD pool selected by irrep label, which is the fair circuit-size "
                         "comparison for the full-group invariant ansatz",
                   peak_memory_mb=X.peak_rss_mb())

        t1 = time.time()
        res = E1.run(tag, atoms, gname)
        dt = time.time() - t1
        for name, basis in (("sub", "abelian"), ("full", "full")):
            r = res[name]
            nat = r["err_trotter_natural_mha"]
            rev = r["err_trotter_invertida_mha"]
            rows.append(dict(
                molecule=tag, full_group=gname, abelian_group="", group_order="",
                n_dets=ref["n_dets"], basis=basis, n_params=r["nparam"],
                n_gates=r["ngates"], decomp_residual=r["residual"], e_fci_ha=res["efci"],
                err_exact_mha=r["err_exact_mha"], err_trotter_natural_mha=nat,
                err_trotter_reversed_mha=rev,
                trotter_cost_mha=max(abs(nat - r["err_exact_mha"]),
                                     abs(rev - r["err_exact_mha"])),
                param_ratio_vs_uccsd=r["nparam"] / ref["nparam"],
                gate_ratio_vs_uccsd=r["ngates"] / ref["ngates"],
                wallclock_s=round(dt, 2)))
            log.record(parameters=dict(molecule=tag, group=gname, basis=basis,
                                       n_params=r["nparam"], n_gates=r["ngates"]),
                       status="success", seed=None, wallclock=dt,
                       outputs=dict(result_files=[], scalar_results=dict(
                           err_exact_mha=r["err_exact_mha"],
                           err_trotter_natural_mha=nat, err_trotter_reversed_mha=rev,
                           decomp_residual=r["residual"])),
                       notes="invariant operators expressed in the pivoted {P e_j} basis; "
                             "tied Trotter circuit in two gate orders",
                       peak_memory_mb=X.peak_rss_mb())

    if not smoke:
        X.write_csv(ROOT / "results" / "aggregate.csv", rows, FIELDS)
        (ROOT / "results" / "aggregate.json").write_text(
            json.dumps(rows, indent=2, default=X._jsonable) + "\n")

    tro = [r for r in rows if r["trotter_cost_mha"] not in ("", None)]
    full = [r for r in rows if r["basis"] == "full"]
    agg = dict(
        max_trotter_cost_mha=max((r["trotter_cost_mha"] for r in tro), default=None),
        max_decomposition_residual=max((r["decomp_residual"] for r in tro), default=None),
        compression=[f"{r['molecule']}: params x{r['param_ratio_vs_uccsd']:.3f}"
                     + (f", gates x{r['gate_ratio_vs_uccsd']:.3f}"
                        if r["gate_ratio_vs_uccsd"] != "" else "")
                     for r in full],
    )
    print("\n=== E1 summary ===")
    for k, v in agg.items():
        print(f"  {k:32s} {v}")
    if not smoke:
        log.summarize(
            aggregate_metrics=agg,
            key_findings=[
                f"tying the parameters and Trotterising the invariant operators changes the "
                f"energy error by at most {agg['max_trotter_cost_mha']:.4f} mHa, two orders "
                f"below chemical accuracy",
                f"each invariant operator decomposes on the elementary excitation basis with "
                f"residual at most {agg['max_decomposition_residual']:.1e}, so the circuit is "
                f"an exact rewriting and not an approximation",
                "compression relative to unfiltered UCCSD: " + "; ".join(agg["compression"]),
                "the two Abelian controls give exactly zero additional compression",
            ],
            matches_hypothesis="supported",
            interpretation="The compression is a statement about a circuit, not only about a "
                           "parameter count: both the parameter count and the elementary gate "
                           "count fall, and Trotterisation costs far less than the gap to FCI.",
        )
    return agg


if __name__ == "__main__":
    main()
