"""E8 -- what does the full-group compression buy in two-qubit gates?

The compression 75 -> 30 (NH3) and 65 -> 21 (CH4) against the Abelian filter is a
statement about PARAMETERS.  A tied parameter driving twelve elementary excitations
still costs twelve gates, so a parameter count says nothing about a circuit on its
own.  This experiment prices every pool in CNOTs under the one ruler that can be
checked against a published figure.

The ruler is He et al.'s Qubit-Excitation-Based scheme:

    CNOTs = 2 * (spin-orbital singles) + 13 * (spin-orbital doubles)

counted over operators, not over parameters.

THE CURRENCY IS NOT A FREE CHOICE, and getting it wrong is the whole experiment.  There
are two ways to count a pool:

  distinct     the set of distinct (occ, virt) elementary excitations
  incidence    every (operator, elementary excitation) pair in the Trotterised product,
               so an excitation driven by three different parameters counts three times

They differ because the same-spin double (i_a j_a -> a_a b_a) appears in TWO singlet-UCCSD
parameters with coefficients +2 and -2 -- the antisymmetrisation t_ijab - t_ijba.  For NH3
the unfiltered pool has 375 incidences against 315 distinct operators.

Figure 5 of He et al. is in the DISTINCT currency, and this is settled by measurement, not
by taste: all four of their UCCSD and SymUCCSD bars are reproduced exactly by the distinct
count and by no other.  Both currencies are therefore reported for every pool, with the
distinct one as the headline, because it is the one the published gate fixes.

This matters beyond this experiment.  E1 recorded its `uccsd` and `abelian_subset` rows in
the distinct currency and its `abelian` and `full` rows in the incidence currency, so E1's
headline 163 -> 147 for NH3 compares one against the other.  In a single currency it is
163 -> 117, and CH4's apparent 146 -> 146 -- zero gate reduction -- is 146 -> 128.  The
regression check below asserts E1's numbers pool by pool in the currency E1 used, which is
what proves the mismatch is in E1's bookkeeping rather than in this run's sources.

Four pools per molecule:

  uccsd              unfiltered, the denominator of every ratio
  abelian_subset     the Abelian filter as the literature implements it: a subset of
                     the UCCSD pool selected by irrep label.  This is SymUCCSD, and
                     it is the headline baseline
  abelian_projector  the same subgroup through the symmetric projector.  Reported
                     because it is the like-for-like construction against `full`, and
                     because it is five times larger -- the favourable-looking
                     comparison, shown and declined
  full               the full point group through the symmetric projector

plus HiUCCSD in both bases, as a secondary check on the ruler.

Nothing here is an energy.  No VQE, no diagonalisation, no optimisation: the pools are
enumerated and the arithmetic is integer.

Usage:  python run.py [nh3|ch4|h2o|all] [--smoke]
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
from pyscf import ao2mo, gto, scf, symm  # noqa: E402

from audit import gates_of, irrep_of  # noqa: E402
from e1_trotter import sparse_invariant_basis  # noqa: E402
from group_compression import CASES, build_ops, group, mo_rep, to_str  # noqa: E402
from hi_count import ham_term_keys  # noqa: E402
from indep import SpinFreeHam, build_det_basis, kappa_matrix  # noqa: E402

ROOT = HERE.parent

# --- the ruler ---------------------------------------------------------------
# He et al., JCTC 2026, 22, 6008-6019: the QEB scheme costs 2 CNOTs for a
# spin-orbital single excitation and 13 for a spin-orbital double.
CNOT_SINGLE = 2
CNOT_DOUBLE = 13

# --- the validation gate -----------------------------------------------------
# Figure 5 of the same paper, read as (n_operators, n_CNOTs).  The first two rows of
# each molecule are ASSERTED: if the counting convention cannot reproduce them the
# experiment is measuring something else and nothing downstream is worth reporting.
# HiUCCSD is reported but NOT asserted -- their pool depends on a compress(1e-8)
# threshold and on the basis, so a divergence there is information, not error.
PUBLISHED = {
    ("nh3", "uccsd"): (315, 3765),
    ("nh3", "abelian_subset"): (163, 1921),
    ("nh3", "hiuccsd_nonadapted"): (197, 2407),
    ("ch4", "uccsd"): (560, 6840),
    ("ch4", "abelian_subset"): (146, 1788),
    ("ch4", "hiuccsd_nonadapted"): (410, 5088),
}
ASSERTED = {"uccsd", "abelian_subset"}

# E1's measured counts, as a regression check against a folder already validated.
# The third field names the currency E1 used for that row -- `uccsd` and `abelian_subset`
# came from `uccsd_reference`, which counts distinct operators, while `abelian` and `full`
# came from `E1.run`, which counts incidences.  Asserting each row in its own currency is
# what shows the sources have not drifted AND that the mismatch is E1's bookkeeping.
E1_COUNTS = {  # pool -> (n_params, n_gates, currency)
    "nh3": {"uccsd": (135, 315, "distinct"), "abelian_subset": (75, 163, "distinct"),
            "abelian_projector": (75, 779, "incidence"), "full": (30, 147, "incidence")},
    "ch4": {"uccsd": (230, 560, "distinct"), "abelian_subset": (65, 146, "distinct"),
            "abelian_projector": (65, 182, "incidence"), "full": (21, 146, "incidence")},
    "h2o": {"uccsd": (65, 140, "distinct"), "abelian_subset": (26, 48, "distinct"),
            "abelian_projector": (26, 54, "incidence"), "full": (26, 54, "incidence")},
}

ORDER = ["nh3", "ch4"]


def classify(items):
    """Split (occ, virt) excitations into spin-orbital singles and doubles.

    An excitation is a single when one spin-orbital is annihilated and one created, a
    double when two are.  Anything else has no price under the ruler and is a bug, so
    it raises rather than being silently dropped (verification criterion 3).
    """
    n_s = n_d = 0
    for occ, virt in items:
        k = len(occ)
        if k != len(virt):
            raise ValueError(f"unbalanced excitation: {occ} -> {virt}")
        if k == 1:
            n_s += 1
        elif k == 2:
            n_d += 1
        else:
            raise ValueError(f"rank-{k} excitation has no price under the QEB ruler")
    return n_s, n_d


def cnots(n_s, n_d):
    return CNOT_SINGLE * n_s + CNOT_DOUBLE * n_d


def price(incidences):
    """Price one pool in both currencies.  See the module docstring: `distinct` is the
    headline because it is the currency Figure 5 is written in; `incidence` is what a
    tied-parameter Trotter circuit executes, and is the currency E1's projector rows used.
    """
    out = {}
    for currency, items in (("distinct", sorted(set(incidences))),
                            ("incidence", incidences)):
        n_s, n_d = classify(items)
        out[currency] = dict(n_gates=len(items), n_singles=n_s, n_doubles=n_d,
                             cnots=cnots(n_s, n_d))
    return out


def projector_pool(Ms, sub, Ep, Kel, nrm, elems, nS):
    """Pool built by the symmetric projector over the subgroup indexed by `sub`.

    Each invariant operator is projected onto the elementary excitation basis by a
    Frobenius inner product, exact because distinct generators connect disjoint pairs
    of determinants.  Returns the incidence list, the parameter count, and the worst
    reconstruction residual over the operators.
    """
    Bs, Bd = sparse_invariant_basis([Ms[k] for k in sub], nS)
    ops = build_ops(Bs, Bd, Ep, nS)
    inc, resid = [], 0.0
    for O in ops:
        w = np.array([float(O.multiply(Kel[k]).sum()) / nrm[k] for k in range(len(Kel))])
        nz = np.nonzero(np.abs(w) > 1e-12)[0]
        rec = sum(w[k] * Kel[k] for k in nz)
        d = O - rec
        resid = max(resid, float(abs(d).max()) if d.nnz else 0.0)
        inc.extend(elems[k] for k in nz)
    return inc, len(ops), resid


def hiuccsd_pool(atoms, gates, symmetry):
    """HiUCCSD: keep an excitation iff its term survives in the compressed Hamiltonian.

    `symmetry=False` is He et al.'s own basis, where symmetry-forbidden integrals do
    not vanish numerically; `symmetry=True` is the adapted basis, where they do.
    """
    mol = gto.M(atom=to_str(atoms), basis="sto-3g", symmetry=symmetry, verbose=0)
    mf = scf.RHF(mol).run()
    keys = ham_term_keys(mol, mf.mo_coeff)
    keep = [g for g in gates if (g["virt"], g["occ"]) in keys]
    inc = [(g["occ"], g["virt"]) for g in keep]
    return inc, len({g["param"] for g in keep})


def measure(key, tag, atoms, gname):
    """Every pool for one molecule, priced in CNOTs."""
    t0 = time.time()
    mol = gto.M(atom=to_str(atoms), basis="sto-3g", symmetry=True, verbose=0)
    mf = scf.RHF(mol).run()
    nao, nelec = mol.nao, mol.nelectron
    no = nelec // 2

    # 1. representation validity
    G, H = group(gname)
    D, rep_err = mo_rep(mol, mf.mo_coeff, G)
    Hk = [i for i, R in enumerate(G) if any(np.allclose(R, h) for h in H)]

    sing = [(i, a) for i in range(no) for a in range(no, nao)]
    nS = len(sing)
    Ms = [np.array([[Dg[j, i] * Dg[b, a] for (i, a) in sing] for (j, b) in sing]) for Dg in D]

    dets, index = build_det_basis(nao, no, no)
    h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), nao)
    ham = SpinFreeHam(h1, eri, mol.energy_nuc(), dets, index)
    Ep = [ham.E[a][i].tocsr() for (i, a) in sing]

    npar, gates = gates_of(nao, nelec)
    elems = sorted({(g["occ"], g["virt"]) for g in gates})
    Kel = [kappa_matrix(e, dets, index) for e in elems]
    nrm = np.array([float(k.multiply(k).sum()) for k in Kel])

    # 2. incidence bookkeeping.  The two currencies diverge for the UCCSD pool because
    #    each same-spin double appears in two parameters with coefficients +2 and -2.
    #    That is a property of the singlet parameterisation, not a bug, so it is
    #    recorded rather than raised -- but it is the reason every pool below is priced
    #    in both currencies instead of one.
    inc_uccsd = [(g["occ"], g["virt"]) for g in gates]
    currency_gap = len(inc_uccsd) - len(elems)

    orbsym = np.array(symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb, mf.mo_coeff))
    sym_ids = {g["param"] for g in gates if irrep_of(g, orbsym) == 0}
    inc_sub = [(g["occ"], g["virt"]) for g in gates if g["param"] in sym_ids]

    pools = {
        "uccsd": dict(inc=inc_uccsd, nparam=npar, resid=None,
                      construction="unfiltered"),
        "abelian_subset": dict(inc=inc_sub, nparam=len(sym_ids), resid=None,
                               construction="irrep subset of the UCCSD pool"),
    }
    for name, sub in (("abelian_projector", Hk), ("full", list(range(len(G))))):
        inc, nparam, resid = projector_pool(Ms, sub, Ep, Kel, nrm, elems, nS)
        pools[name] = dict(inc=inc, nparam=nparam, resid=resid,
                           construction="symmetric projector")
    for name, sym in (("hiuccsd_nonadapted", False), ("hiuccsd_adapted", True)):
        inc, nparam = hiuccsd_pool(atoms, gates, sym)
        pools[name] = dict(inc=inc, nparam=nparam, resid=None,
                           construction=f"Hamiltonian term set, symmetry={sym}")

    out = dict(molecule=tag, key=key, full_group=gname, abelian_group=mol.groupname,
               topgroup=mol.topgroup, group_order=len(G), abelian_order=len(H),
               n_dets=len(dets), rep_error=float(rep_err), nao=nao, nelec=nelec,
               uccsd_currency_gap=currency_gap,
               wallclock_s=round(time.time() - t0, 2), pools={})
    for name, p in pools.items():
        out["pools"][name] = dict(       # 3. class split is exhaustive, inside price()
            n_params=p["nparam"], decomp_residual=p["resid"],
            construction=p["construction"], **price(p["inc"]))
    return out


def verify(res):
    """Verification criteria 4, 5 and 6, in dependency order.  Returns the gate report."""
    key = res["key"]
    report = []

    # 4. parameter and operator counts reproduce E1, each row in the currency E1 used
    for name, (p_e1, g_e1, cur) in E1_COUNTS.get(key, {}).items():
        p, g = res["pools"][name]["n_params"], res["pools"][name][cur]["n_gates"]
        ok = (p, g) == (p_e1, g_e1)
        report.append(dict(check="e1_regression", pool=name, currency=cur,
                           expected=[p_e1, g_e1], measured=[p, g], passed=ok))
        if not ok:
            raise AssertionError(
                f"{key}/{name}: E1 measured {p_e1} params / {g_e1} operators in the "
                f"{cur} currency, this run gives {p} / {g}; the vendored sources have drifted")

    # 5. decomposition residual at machine precision
    for name in ("abelian_projector", "full"):
        r = res["pools"][name]["decomp_residual"]
        ok = r is not None and r < 1e-12
        report.append(dict(check="decomposition_residual", pool=name,
                           expected="< 1e-12", measured=r, passed=ok))
        if not ok:
            raise AssertionError(f"{key}/{name}: decomposition residual {r} is not machine zero; "
                                 f"the projector pool is an approximation, not a rewriting")

    # 6. the published gate, in the distinct currency and in that one only
    for (k, name), (n_pub, c_pub) in PUBLISHED.items():
        if k != key:
            continue
        d = res["pools"][name]["distinct"]
        i = res["pools"][name]["incidence"]
        ok = (d["n_gates"], d["cnots"]) == (n_pub, c_pub)
        report.append(dict(check="published_figure_5", pool=name, currency="distinct",
                           expected=[n_pub, c_pub], measured=[d["n_gates"], d["cnots"]],
                           also_incidence=[i["n_gates"], i["cnots"]],
                           passed=ok, asserted=name in ASSERTED))
        if not ok and name in ASSERTED:
            raise AssertionError(
                f"{key}/{name}: Figure 5 reports {n_pub} operators / {c_pub} CNOTs, this "
                f"run gives {d['n_gates']} / {d['cnots']} distinct and {i['n_gates']} / "
                f"{i['cnots']} by incidence; the counting convention is wrong")
    return report


FIELDS = ["molecule", "key", "full_group", "abelian_group", "group_order", "n_dets",
          "pool", "currency", "construction", "n_params", "n_gates", "n_singles",
          "n_doubles", "cnots", "decomp_residual", "param_ratio_vs_uccsd",
          "gate_ratio_vs_uccsd", "cnot_ratio_vs_uccsd", "cnot_ratio_vs_abelian_subset",
          "published_n_gates", "published_cnots", "published_check", "wallclock_s"]

POOL_ORDER = ["uccsd", "abelian_subset", "abelian_projector", "full",
              "hiuccsd_nonadapted", "hiuccsd_adapted"]
CURRENCIES = ["distinct", "incidence"]


def rows_of(res):
    """One row per (pool, currency).  Ratios are always taken inside a currency -- mixing
    them is the defect this experiment found in E1, and the CSV must make it impossible.
    """
    rows = []
    for cur in CURRENCIES:
        ref = res["pools"]["uccsd"][cur]
        base = res["pools"]["abelian_subset"][cur]
        for name in POOL_ORDER:
            p = res["pools"][name]
            q = p[cur]
            pub = PUBLISHED.get((res["key"], name))
            if pub is None or cur != "distinct":
                check = ""
            elif (q["n_gates"], q["cnots"]) == tuple(pub):
                check = "exact"
            else:
                check = "differs"
            rows.append(dict(
                molecule=res["molecule"], key=res["key"], full_group=res["full_group"],
                abelian_group=res["abelian_group"], group_order=res["group_order"],
                n_dets=res["n_dets"], pool=name, currency=cur,
                construction=p["construction"], n_params=p["n_params"],
                n_gates=q["n_gates"], n_singles=q["n_singles"], n_doubles=q["n_doubles"],
                cnots=q["cnots"],
                decomp_residual="" if p["decomp_residual"] is None else p["decomp_residual"],
                param_ratio_vs_uccsd=p["n_params"] / res["pools"]["uccsd"]["n_params"],
                gate_ratio_vs_uccsd=q["n_gates"] / ref["n_gates"],
                cnot_ratio_vs_uccsd=q["cnots"] / ref["cnots"],
                cnot_ratio_vs_abelian_subset=q["cnots"] / base["cnots"],
                published_n_gates="" if (pub is None or cur != "distinct") else pub[0],
                published_cnots="" if (pub is None or cur != "distinct") else pub[1],
                published_check=check, wallclock_s=res["wallclock_s"]))
    return rows


def make_log(smoke=False):
    import tempfile
    path = (Path(tempfile.mkdtemp()) / "smoke_log.json") if smoke else (ROOT / "experiment_log.json")
    return X.ExperimentLog(
        path,
        experiment_id=ROOT.name,
        paper_project="2026-08_point-group-vqe-free",
        opportunity_reference="opportunities.md #2",
        specification="EXPERIMENT_gate_counts.md",
        calibrated_by="experiments/2026-08-11_e1_invariant_trotter_ansatz",
        mode="symbolic",
        subareas_active=["symbolic"],
        target_venue="Journal of Mathematical Chemistry",
        ruler=f"QEB: {CNOT_SINGLE} CNOTs per spin-orbital single, "
              f"{CNOT_DOUBLE} per spin-orbital double, over operator incidences",
        created_at=X.utcnow(),
        plan_md_hash=X.sha256_file(ROOT / "plan.md") if (ROOT / "plan.md").exists() else None,
        code_commit=None,
        system_info_hash=None,
    )


def main():
    which = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else "all"
    smoke = "--smoke" in sys.argv
    keys = ["h2o"] if smoke else (ORDER if which == "all" else [which])

    if not smoke:
        X.capture_env(ROOT / "env", project_root=ROOT.parent.parent)
    log = make_log(smoke=smoke)
    rows, all_res = [], []

    for key in keys:
        tag, atoms, gname, _ = CASES[key]
        print(f"\n{'=' * 78}\n{tag}   full group {gname}", flush=True)
        t0 = time.time()
        with X.Timer() as T:
            res = measure(key, tag, atoms, gname)
            report = verify(res)
        res["verification"] = report
        all_res.append(res)

        print(f"  representation error |U^T S U - S| = {res['rep_error']:.2e}   "
              f"determinants {res['n_dets']}   UCCSD currency gap "
              f"{res['uccsd_currency_gap']}")
        all_rows = rows_of(res)
        for cur in CURRENCIES:
            print(f"\n  --- {cur} currency"
                  + ("  [headline: this is the one Figure 5 is written in]"
                     if cur == "distinct" else
                     "  [what a tied-parameter Trotter circuit executes]") + " ---")
            print(f"  {'pool':22s} {'params':>7s} {'gates':>7s} {'sing':>6s} {'doub':>6s} "
                  f"{'CNOTs':>7s}  {'vs UCCSD':>9s}  {'vs Abel':>8s}  published")
            for r in [x for x in all_rows if x["currency"] == cur]:
                pub = (f"{r['published_cnots']} {r['published_check']}"
                       if r["published_cnots"] != "" else "")
                print(f"  {r['pool']:22s} {r['n_params']:7d} {r['n_gates']:7d} "
                      f"{r['n_singles']:6d} {r['n_doubles']:6d} {r['cnots']:7d}  "
                      f"{r['cnot_ratio_vs_uccsd']:9.4f}  "
                      f"{r['cnot_ratio_vs_abelian_subset']:8.4f}  {pub}")
        rows.extend(all_rows)

        for name in POOL_ORDER:
            for cur in CURRENCIES:
                p, q = res["pools"][name], res["pools"][name][cur]
                log.record(
                    parameters=dict(molecule=tag, key=key, full_group=gname,
                                    abelian_group=res["abelian_group"], pool=name,
                                    currency=cur, construction=p["construction"],
                                    basis="sto-3g", ruler_single=CNOT_SINGLE,
                                    ruler_double=CNOT_DOUBLE),
                    status="success", seed=None,
                    wallclock=T.dt / (len(POOL_ORDER) * len(CURRENCIES)),
                    t_start=T.start, t_end=T.end,
                    inputs=dict(data_files=[], data_hashes={
                        "geometry": X.sha256_text(to_str(atoms))}),
                    outputs=dict(result_files=["results/aggregate.csv"], scalar_results=dict(
                        n_params=p["n_params"], n_gates=q["n_gates"],
                        n_singles=q["n_singles"], n_doubles=q["n_doubles"],
                        cnots=q["cnots"], decomp_residual=p["decomp_residual"])),
                    notes=f"{p['construction']}; priced in the {cur} currency; counts are "
                          f"integer enumerations, no optimisation and no energy is computed "
                          f"in this experiment",
                    peak_memory_mb=X.peak_rss_mb())
        print(f"  [{time.time() - t0:.0f}s]", flush=True)

    if smoke:
        return rows, all_res

    X.write_csv(ROOT / "results" / "aggregate.csv", rows, FIELDS)
    (ROOT / "results" / "aggregate.json").write_text(
        json.dumps(all_res, indent=2, default=X._jsonable) + "\n")

    # --- summary -------------------------------------------------------------
    agg = {}
    for res in all_res:
        k = res["key"]
        P = res["pools"]
        entry = dict(param_reduction=f"{P['abelian_subset']['n_params']} -> "
                                     f"{P['full']['n_params']}",
                     param_saving_pct=round(100.0 * (1.0 - P["full"]["n_params"]
                                                     / P["abelian_subset"]["n_params"]), 2))
        for cur in CURRENCIES:
            base, full = P["abelian_subset"][cur], P["full"][cur]
            entry[cur] = dict(
                gate_reduction=f"{base['n_gates']} -> {full['n_gates']}",
                gate_saving_pct=round(100.0 * (1.0 - full["n_gates"] / base["n_gates"]), 2),
                cnot_reduction=f"{base['cnots']} -> {full['cnots']}",
                cnot_saving_pct=round(100.0 * (1.0 - full["cnots"] / base["cnots"]), 2),
                removed_singles=base["n_singles"] - full["n_singles"],
                removed_doubles=base["n_doubles"] - full["n_doubles"],
                projector_abelian_cnots=P["abelian_projector"][cur]["cnots"])
        agg[k] = entry

    print("\n=== E8 summary: full group against the Abelian filter ===")
    for k, v in agg.items():
        print(f"  {k}: params {v['param_reduction']} ({v['param_saving_pct']}%)")
        for cur in CURRENCIES:
            c = v[cur]
            print(f"      {cur:9s} gates {c['gate_reduction']:12s} "
                  f"({c['gate_saving_pct']:5.1f}%)   CNOTs {c['cnot_reduction']:14s} "
                  f"({c['cnot_saving_pct']:5.1f}%)   removed {c['removed_singles']}s "
                  f"{c['removed_doubles']}d")

    gate_rows = [c for res in all_res for c in res["verification"]
                 if c["check"] == "published_figure_5"]
    n_exact = sum(c["passed"] for c in gate_rows)
    D = "distinct"
    findings = [
        f"the QEB ruler reproduces {n_exact} of the {len(gate_rows)} published Figure 5 bars "
        f"exactly in the distinct-operator currency, with no fitted quantity; the incidence "
        f"currency reproduces none of them, which is what fixes the currency",
    ]
    for k, v in agg.items():
        c = v[D]
        findings.append(
            f"{k}: against the Abelian filter the full point group cuts parameters by "
            f"{v['param_saving_pct']}% and CNOTs by {c['cnot_saving_pct']}% "
            f"({c['cnot_reduction']}), removing {c['removed_singles']} singles and "
            f"{c['removed_doubles']} doubles")
    findings.append(
        "the CNOT saving tracks the flat operator saving rather than exceeding it ("
        + "; ".join(f"{k} {v[D]['cnot_saving_pct']}% against {v[D]['gate_saving_pct']}%"
                    for k, v in agg.items())
        + "), so the filter shows no useful preference for the expensive operator class")
    findings.append(
        "E1 priced its subset pools in the distinct currency and its projector pools by "
        "incidence, so its headline 163 -> 147 for NH3 and 146 -> 146 for CH4 compare one "
        "currency against the other; in the distinct currency they are "
        + "; ".join(f"{k} {v[D]['gate_reduction']}" for k, v in agg.items()))
    findings.append(
        "the projector-built Abelian pool costs "
        + ", ".join(f"{k} {v[D]['projector_abelian_cnots']} CNOTs" for k, v in agg.items())
        + ", against the irrep-subset pool the literature actually implements; the headline "
          "comparison uses the subset, which is the smaller and less favourable baseline")

    supported = all(v[D]["cnot_saving_pct"] > v[D]["gate_saving_pct"] + 5.0
                    for v in agg.values())
    log.summarize(
        aggregate_metrics=agg,
        key_findings=findings,
        matches_hypothesis="supported" if supported else "refuted",
        interpretation=(
            "The parameter compression transfers to circuit size only in part. The full-group "
            "filter does cut CNOTs, by a quarter for NH3 and an eighth for CH4, but the saving "
            "tracks the flat operator count instead of beating it: the filter removes no "
            "disproportionate share of the expensive doubles. The compression is therefore "
            "mostly a statement about optimiser dimension, and the honest circuit claim is the "
            "measured percentage rather than the parameter ratio."
            if not supported else
            "The full-group filter removes preferentially doubles, so the CNOT saving exceeds "
            "the flat operator saving and the compression carries a circuit consequence beyond "
            "the parameter count."),
    )
    return agg


if __name__ == "__main__":
    main()
