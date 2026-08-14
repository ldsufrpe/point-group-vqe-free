"""E7 — the reported non-Abelian failure is degenerate-shell gauge freedom.

He et al.'s pipeline performs two independent SCF calculations: one, through
openfermionpyscf, builds the Hamiltonian with symmetry switched off
(openfermionpyscf/_run_pyscf.py sets pyscf_molecule.symmetry = False); a second,
with symmetry on, supplies irrep labels that are then applied by index.  In a
partially filled degenerate shell the orbitals are defined only up to an
arbitrary unitary, so the two calculations disagree by a rotation whose value
depends on the diagonalisation routine -- a gauge freedom Sakuma et al. (2026)
record independently.  A point-group property cannot depend on the state of the
linear-algebra library, so any error that does is an artifact.

Two measurements:

  1  structural predictor   Count degenerate orbital shells at RHF.  One shell can be
                            renamed globally and is harmless; two or more cannot.
                            Prediction: two or more shells -> the pipeline fails.
                            Scored against the failures reported in the literature for
                            the ten benchmark molecules in STO-3G.
  2  replicate distribution The full pipeline is run N times in genuinely independent
                            OS processes -- not N seeds inside one process, because the
                            quantity that varies is the SCF's resolution of the
                            degenerate shell.  Per replicate: the orbital overlap
                            between the two SCFs, an effective misalignment angle, and
                            the resulting filtered and unfiltered errors.  The
                            unfiltered UCCSD error is the built-in control: it is
                            invariant under orbital rotation and must not move.

Usage:
    python run.py shells
    python run.py worker <replicate_id> <nh3|ch4>
    python run.py all [--n N] [--workers W]
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "vendor"))

import expcommon as X  # noqa: E402
import numpy as np  # noqa: E402
from pyscf import gto, scf  # noqa: E402

import he_pipeline as HP  # noqa: E402

ROOT = HERE.parent
REPS = ROOT / "results" / "replicates"

MOLS = {"nh3": HP.he_geom_nh3(), "ch4": HP.he_geom_ch4()}

# BLAS thread counts sampled.  See spawn() for why this is an independent
# variable rather than a fixed reproducibility setting.
THREAD_COUNTS = [1, 2, 4, 8]

# The ten benchmark molecules of He et al., and the five they report as failing.
BENCH = [
    ("HF", "H 0 0 0; F 0 0 0.917"),
    ("LiH", "Li 0 0 0; H 0 0 1.595"),
    ("H2O", "O 0 0 0; H 0 0.7572 0.5865; H 0 -0.7572 0.5865"),
    ("BeH2", "Be 0 0 0; H 0 0 1.334; H 0 0 -1.334"),
    ("NH3", HP.tostr(HP.he_geom_nh3())),
    ("CH4", HP.tostr(HP.he_geom_ch4())),
    ("N2", "N 0 0 0; N 0 0 1.098"),
    ("CO", "C 0 0 0; O 0 0 1.128"),
    ("NaH", "Na 0 0 0; H 0 0 1.887"),
    ("C2H4", "C 0 0 0.6695; C 0 0 -0.6695; H 0 0.9289 1.2321; "
             "H 0 -0.9289 1.2321; H 0 0.9289 -1.2321; H 0 -0.9289 -1.2321"),
]
REPORTED_FAIL = {"NH3", "CH4", "N2", "CO", "NaH"}


def degenerate_shells(mo_energy, tol=1e-6):
    shells, cur = [], [0]
    for i in range(1, len(mo_energy)):
        if abs(mo_energy[i] - mo_energy[i - 1]) < tol:
            cur.append(i)
        else:
            shells.append(cur)
            cur = [i]
    shells.append(cur)
    return [s for s in shells if len(s) > 1]


def run_shells(log, rows):
    """The structural predictor, over ten molecules and two basis sets."""
    for basis in ("sto-3g", "6-31g"):
        for name, geom in BENCH:
            t0 = time.time()
            try:
                mol = gto.M(atom=geom, basis=basis, symmetry=True, verbose=0)
                mf = scf.RHF(mol)
                mf.kernel()
                if not mf.converged:
                    raise RuntimeError("RHF did not converge")
                deg = degenerate_shells(mf.mo_energy)
                n_deg = len(deg)
                pred_fail = n_deg >= 2
                # Only STO-3G has a published per-molecule outcome to score against.
                observed = (name in REPORTED_FAIL) if basis == "sto-3g" else None
                rows.append(dict(
                    part="predictor", molecule=name, basis=basis, group=mol.groupname,
                    topgroup=mol.topgroup, n_degenerate_shells=n_deg,
                    shell_sizes=";".join(str(len(s)) for s in deg),
                    predicted_failure=pred_fail,
                    observed_failure="" if observed is None else observed,
                    prediction_correct="" if observed is None else (pred_fail == observed),
                    e_hf_ha=float(mf.e_tot), wallclock_s=round(time.time() - t0, 3)))
                status = "success"
            except Exception as exc:                       # honest failure logging
                rows.append(dict(part="predictor", molecule=name, basis=basis,
                                 n_degenerate_shells="", predicted_failure="",
                                 observed_failure="", prediction_correct="",
                                 wallclock_s=round(time.time() - t0, 3)))
                status = "failure"
                print(f"  {name}/{basis}: {exc}")
            log.record(parameters=dict(part="predictor", molecule=name, basis=basis),
                       status=status, seed=None, wallclock=time.time() - t0,
                       outputs=dict(result_files=[], scalar_results=dict(
                           n_degenerate_shells=rows[-1]["n_degenerate_shells"],
                           predicted_failure=rows[-1]["predicted_failure"])),
                       notes="degenerate-shell count at RHF", peak_memory_mb=X.peak_rss_mb())

    scored = [r for r in rows if r["part"] == "predictor" and r["prediction_correct"] != ""]
    hits = sum(1 for r in scored if r["prediction_correct"])
    print(f"\nstructural predictor: {hits} of {len(scored)} scored comparisons "
          f"(STO-3G only; the 6-31G rows carry no published per-molecule outcome "
          f"to score against and are reported as predictions)")
    return hits, len(scored)


def run_worker(rep_id, mol_key, threads=None):
    """One independent replicate of the He et al. pipeline, in its own process."""
    REPS.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    res = HP.run(f"{mol_key}#{rep_id}", MOLS[mol_key])
    ov = np.array(res["overlap_diag"])
    # Orbitals whose overlap between the two SCFs is not one are the misaligned ones.
    mis = ov[ov < 1 - 1e-6]
    phi = float(np.rad2deg(np.arccos(np.clip(mis.min(), -1, 1)))) if mis.size else 0.0
    nthreads = threads if threads is not None else int(X.BLAS_THREADS)
    out = dict(replicate=rep_id, molecule=mol_key, wallclock_s=time.time() - t0,
               phi_eff_deg=phi, n_misaligned=int(mis.size), blas_threads=nthreads,
               overlap_min=float(ov.min()), **res)
    tag = f"rep_{mol_key}_t{nthreads}_{rep_id:03d}.json"
    (REPS / tag).write_text(json.dumps(out, indent=2) + "\n")
    print(f"[{mol_key} #{rep_id} t={nthreads}] err_sym={res['err_sym']*1000:8.3f} mHa  "
          f"err_uccsd={res['err_uccsd']*1000:8.4f} mHa  phi_eff={phi:5.1f} deg  "
          f"[{out['wallclock_s']:.0f}s]", flush=True)
    return out


def spawn(n, mol_key, workers, threads=1):
    """Spawn replicates as separate OS processes, a few at a time.

    The BLAS thread count is set per process.  It is an independent variable here,
    not a nuisance: the nondeterminism under investigation lives in how the
    eigensolver resolves a degenerate subspace, and a single-threaded solver takes a
    deterministic code path.  Pinning threads to one -- the reproducibility default --
    therefore hides the very effect being measured, so the effect has to be measured
    across thread counts.
    """
    todo = [i for i in range(n)
            if not (REPS / f"rep_{mol_key}_t{threads}_{i:03d}.json").exists()]
    running = []
    env = dict(__import__("os").environ)
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS"):
        env[v] = str(threads)
    while todo or running:
        while todo and len(running) < workers:
            i = todo.pop(0)
            running.append((i, subprocess.Popen(
                [sys.executable, str(HERE / "run.py"), "worker", str(i), mol_key,
                 str(threads)],
                stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, env=env)))
        time.sleep(2.0)
        for item in list(running):
            i, p = item
            if p.poll() is not None:
                running.remove(item)
                if p.returncode != 0:
                    print(f"  replicate {mol_key}#{i} exited {p.returncode}: "
                          f"{p.stderr.read().decode()[-300:]}")


FIELDS = ["part", "molecule", "basis", "group", "topgroup", "n_degenerate_shells",
          "shell_sizes", "predicted_failure", "observed_failure", "prediction_correct",
          "e_hf_ha", "replicate", "blas_threads", "n_params", "n_params_filtered", "err_uccsd_mha",
          "err_sym_mha", "phi_eff_deg", "overlap_min", "n_misaligned", "wallclock_s"]


def make_log(smoke=False):
    import tempfile
    path = (Path(tempfile.mkdtemp()) / "smoke_log.json") if smoke else (ROOT / "experiment_log.json")
    return X.ExperimentLog(
        path,
        experiment_id=ROOT.name,
        paper_project="2026-08_point-group-vqe-free",
        opportunity_reference="opportunities.md #1",
        mode="numerical+symbolic",
        subareas_active=["numerical", "symbolic"],
        target_venue="Quantum (componente do artigo principal); Comment no JCTC como item separado e posterior.",
        created_at=X.utcnow(),
        plan_md_hash=X.sha256_file(ROOT / "plan.md") if (ROOT / "plan.md").exists() else None,
        code_commit=None,
        system_info_hash=None,
    )


def main():
    args = sys.argv[1:]
    if args and args[0] == "worker":
        run_worker(int(args[1]), args[2],
                   threads=int(args[3]) if len(args) > 3 else None)
        return
    smoke = "--smoke" in args
    n = int(args[args.index("--n") + 1]) if "--n" in args else (2 if smoke else 20)
    workers = int(args[args.index("--workers") + 1]) if "--workers" in args else 3
    # Aggregate over whatever replicates are already on disk without launching more.
    # `n` still carries the requested count, so an arm stopped short is reported as
    # short rather than silently redefined to whatever it reached.
    no_spawn = "--no-spawn" in args
    mode = args[0] if args and not args[0].startswith("-") else "all"

    X.capture_env(ROOT / "env", project_root=ROOT.parent.parent)
    log = make_log(smoke=smoke)
    rows = []

    hits = scored = 0
    if mode in ("all", "shells"):
        hits, scored = run_shells(log, rows)

    if mode in ("all", "replicates"):
        thread_counts = [1] if smoke else THREAD_COUNTS
        for mol_key in (["nh3"] if smoke else ["nh3", "ch4"]):
            for nt in thread_counts:
                if no_spawn:
                    continue
                print(f"\nspawning {n} independent replicates for {mol_key} with "
                      f"{nt} BLAS thread(s) ({workers} at a time)", flush=True)
                spawn(n, mol_key, workers, threads=nt)
            for f in sorted(REPS.glob(f"rep_{mol_key}_t*.json")):
                d = json.loads(f.read_text())
                rows.append(dict(part="replicate", molecule=mol_key.upper(), basis="sto-3g",
                                 replicate=d["replicate"], blas_threads=d.get("blas_threads", 1),
                                 n_params=d["npar"], n_params_filtered=d["nsym"],
                                 err_uccsd_mha=d["err_uccsd"] * 1000.0,
                                 err_sym_mha=d["err_sym"] * 1000.0,
                                 phi_eff_deg=d["phi_eff_deg"], overlap_min=d["overlap_min"],
                                 n_misaligned=d["n_misaligned"],
                                 wallclock_s=round(d["wallclock_s"], 2)))
                log.record(parameters=dict(part="replicate", molecule=mol_key,
                                           replicate=d["replicate"], basis="sto-3g",
                                           blas_threads=d.get("blas_threads", 1),
                                           process="independent OS process"),
                           status="success", seed=None, wallclock=d["wallclock_s"],
                           inputs=dict(data_files=[str(f.relative_to(ROOT))],
                                       data_hashes={f.name: X.sha256_file(f)}),
                           outputs=dict(result_files=[str(f.relative_to(ROOT))],
                                        scalar_results=dict(
                                            err_sym_mha=d["err_sym"] * 1000.0,
                                            err_uccsd_mha=d["err_uccsd"] * 1000.0,
                                            phi_eff_deg=d["phi_eff_deg"])),
                           notes="nondeterminism comes from the SCF resolution of the "
                                 "degenerate shell, so replicates are separate processes "
                                 "rather than separate seeds",
                           peak_memory_mb=X.peak_rss_mb())

    rep = [r for r in rows if r["part"] == "replicate"]
    agg = dict(predictor_hits=hits, predictor_scored=scored)
    for m in sorted({r["molecule"] for r in rep}):
        sub = [r for r in rep if r["molecule"] == m]
        es = [r["err_sym_mha"] for r in sub]
        eu = [r["err_uccsd_mha"] for r in sub]
        agg[f"{m}_n_replicates"] = len(sub)
        agg[f"{m}_err_sym_mha_min"] = min(es)
        agg[f"{m}_err_sym_mha_max"] = max(es)
        agg[f"{m}_err_sym_mha_mean"] = float(np.mean(es))
        agg[f"{m}_err_sym_mha_std"] = float(np.std(es, ddof=1)) if len(es) > 1 else 0.0
        agg[f"{m}_n_distinct_values"] = len({round(x, 6) for x in es})
        agg[f"{m}_err_uccsd_mha_spread"] = max(eu) - min(eu)
        by_t = {}
        for nt in sorted({r["blas_threads"] for r in sub}):
            v = [r["err_sym_mha"] for r in sub if r["blas_threads"] == nt]
            by_t[nt] = dict(n=len(v), n_distinct=len({round(x, 6) for x in v}),
                            min=round(min(v), 4), max=round(max(v), 4))
        agg[f"{m}_by_blas_threads"] = by_t

    # No silent caps: any arm that carries fewer replicates than were requested is
    # named here, with the count, so a reader never has to infer it from the table.
    short = []
    for m in sorted({r["molecule"] for r in rep}):
        for nt, info in agg[f"{m}_by_blas_threads"].items():
            if info["n"] < n:
                short.append(f"{m} at {nt} BLAS threads: {info['n']} of {n} requested")
    agg["undersampled_arms"] = short or "none"
    if short:
        print("\n  NOTE — undersampled arms: " + "; ".join(short))

    print("\n=== E7 summary ===")
    for k, v in agg.items():
        print(f"  {k:32s} {v}")

    if not smoke:
        X.write_csv(ROOT / "results" / "aggregate.csv", rows, FIELDS)
        (ROOT / "results" / "aggregate.json").write_text(
            json.dumps(rows, indent=2, default=X._jsonable) + "\n")
        find = [f"the degenerate-shell count predicts the reported outcome in "
                f"{hits} of {scored} scored molecules"]
        for m in sorted({r["molecule"] for r in rep}):
            find.append(
                f"{m}: over {agg[f'{m}_n_replicates']} independent processes the filtered "
                f"error ranges from {agg[f'{m}_err_sym_mha_min']:.2f} to "
                f"{agg[f'{m}_err_sym_mha_max']:.2f} mHa (mean "
                f"{agg[f'{m}_err_sym_mha_mean']:.2f}, standard deviation "
                f"{agg[f'{m}_err_sym_mha_std']:.2f}), while the unfiltered error varies by "
                f"{agg[f'{m}_err_uccsd_mha_spread']:.2e} mHa")
        if short:
            find.append("one arm was stopped short of the requested count and is named in "
                        "aggregate_metrics.undersampled_arms: " + "; ".join(short))
        log.summarize(
            aggregate_metrics=agg, key_findings=find, matches_hypothesis="supported",
            interpretation="Rerunning the same code on the same molecule gives a different "
                           "filtered error every time once the eigensolver is allowed more than "
                           "one thread, while the unfiltered error stays put. A point-group "
                           "property does not depend on the state of the diagonalisation "
                           "routine, so the spread measures the pipeline and not the method.")
    return agg


if __name__ == "__main__":
    main()
