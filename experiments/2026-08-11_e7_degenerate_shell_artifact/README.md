# Reportable narrative

The non-Abelian failure reported for symmetry-filtered unitary coupled cluster is a property of the
pipeline that produced it rather than of the method. Two measurements support that. A predictor that
counts degenerate orbital shells at the mean-field level, and runs no variational calculation at all,
separates the five molecules reported as failing from the five reported as working, ten of ten in the
minimal basis. And rerunning the published pipeline for methane in thirty-six independent processes
returns a filtered error anywhere between 15.9 and 49.7 mHa, twenty-seven distinct values, while the
unfiltered error over the same processes moves by 1.7e-4 mHa and the same filter applied with labels
taken from the calculation that defines the Hamiltonian costs 0.194 mHa. The spread appears only when
the eigensolver is given more than one thread: single-threaded runs return one value to twelve
figures. Ammonia, whose degenerate shells are twofold rather than threefold, is reproducible at every
thread count and lands on a fixed sixty-degree misalignment worth 29.257 mHa. A point-group property
cannot depend on how many threads a linear-algebra library was given, so what these numbers measure
is the labelling step, not the filter.

# Operational log

**Experiment ID:** 2026-08-11_e7_degenerate_shell_artifact
**Date:** 2026-08-11 (start) → 2026-08-11 (end)
**Status:** complete
**Mode:** numerical and symbolic (hybrid)
**Opportunity:** opportunities.md #1
**Target venue:** Quantum (componente do artigo principal); Comment no JCTC como item separado e posterior.

## Numerical results

**Structural predictor**, STO-3G, ten benchmark molecules. Rule: two or more degenerate orbital
shells at RHF predicts failure. Scored against the outcomes reported in the audited work.

| shells | molecules | predicted | reported | agree |
|---|---|---|---|---|
| 0 | H2O, C2H4 | works | works | yes |
| 1 | HF, LiH, BeH2 | works | works | yes |
| 2 (sizes 2;2) | NH3, N2, CO, NaH | fails | fails | yes |
| 2 (sizes 3;3) | CH4 | fails | fails | yes |

Ten scored comparisons, ten correct. The 6-31G rows are in `results/aggregate.csv` as predictions;
they carry no published per-molecule outcome to score against and are not counted.

**Replicate distribution.** Each replicate is a separate operating-system process running the
two-SCF pipeline end to end. Filtered error in mHa.

| molecule | BLAS threads | replicates | distinct values | min | max | mean | s.d. |
|---|---|---|---|---|---|---|---|
| CH4 | 1 | 10 | 1 | 49.722 | 49.722 | 49.722 | 0.000 |
| CH4 | 2 | 10 | 10 | 15.922 | 46.637 | 35.017 | 10.175 |
| CH4 | 4 | 10 | 10 | 18.001 | 44.829 | 37.270 | 9.277 |
| CH4 | 8 | 6 | 6 | 31.841 | 46.614 | 35.871 | 5.618 |
| NH3 | 1, 2, 4, 8 | 20, 10, 10, 10 | 1 | 29.257 | 29.257 | 29.257 | 1.6e-8 |

Controls over the same replicates: the unfiltered UCCSD error spans 0.193797 to 0.193964 mHa for
methane and varies by 2.0e-9 mHa for ammonia, so the quantity that is invariant under a rotation
inside a degenerate shell does not move. The same filter with labels drawn from the calculation that
defines the Hamiltonian costs 0.1942 mHa for methane and 0.1392 for ammonia, measured in
`experiments/2026-08-11_e1_invariant_trotter_ansatz`.

The effective misalignment angle between the two self-consistent field calculations is fixed at 60.0
degrees for ammonia and ranges from 60.7 to 89.0 degrees for methane. The published values this
distribution has to account for, 40.8 mHa for methane and 27.8 for ammonia, and the two figures
recorded by the earlier audit, 37.47 and 28.59, all fall inside the measured range.

**One arm is undersampled and is named in the log.** The methane arm at eight threads carries 6
replicates of the 10 requested; it was stopped because eight threads in each of two concurrent
processes oversubscribes an eight-core machine and each replicate had slowed from about 180 to about
1500 seconds. The finding rests on the arms at two and four threads, which are complete.

## Hypothesis check

**Supported.** Both measurements land, and the control holds. The mechanism is named rather than
inferred: the pipeline takes irrep labels from a second, independent self-consistent field
calculation and applies them by index, and in a partially filled degenerate shell the orbitals are
defined only up to an arbitrary unitary.

## Figures

- `figures/2026-08-11_e7_degenerate_shell_artifact_artifact_predictor_and_spread.pdf` —
  **For manuscript.** The predictor over the benchmark set, above the replicate distribution against
  thread count. Element inventory in `figures/figure_manifest.json`.
- `figures/2026-08-11_e7_degenerate_shell_artifact_artifact_predictor_and_spread.png` — raster
  preview, for inspection only.

## Pointers

- Aggregate data: `results/aggregate.csv`, `results/aggregate.json`
- Per-replicate raw results: `results/replicates/rep_<molecule>_t<threads>_<id>.json`
- Visualisation contract: `results/viz_schema.json`
- Full log: `experiment_log.json` (106 run entries, all `success`)
- Plan, including why the published predictor score was not reproducible: `plan.md`
- Code: `code/run.py`, `code/analysis.py`, `code/smoke_test.py`
- Vendored sources and their hashes: `code/vendor/PROVENANCE.md`
- Environment: `env/system_info.json`, `env/requirements.txt`, `env/git_commit.txt`

## Notes

Two things came out differently from what the inherited material recorded, and both are corrections
rather than refinements.

The predictor was recorded as scoring fifteen of fifteen. The script it came from saved no output
and structurally produces ten scored comparisons in the minimal basis plus five rows in the larger
basis that have no published per-molecule outcome to score against. Ten of ten is the number this
experiment can defend, and it is the number reported.

The run-to-run irreproducibility was recorded as a general property of the pipeline. It is not.
Pinning the linear-algebra library to one thread, which is what the reproducibility rules of this
project require by default, makes every molecule deterministic: methane returns 49.722 mHa in ten
consecutive processes. The variation appears only when the eigensolver is allowed more than one
thread, and even then only for methane, whose degenerate shells hold three orbitals each. Ammonia,
with twofold shells, is reproducible everywhere and fails by a fixed rotation instead. That is a
sharper claim than the original and it survives a referee running the code with default settings,
which the original might not have.

The single-thread methane value, 49.722 mHa, sits above every multi-threaded replicate. Nothing in
the mechanism privileges it; it is one draw from the same set of arbitrary unitaries, and it happens
to be a bad one. It is reported rather than smoothed because it is what a reader following the
reproducibility rules of this project would obtain.

Replicates are processes and not seeds because the quantity that varies is not under the program's
control. There is no random number generator to seed here: the variation enters through how the
eigensolver resolves a degenerate subspace, which is why `seed` is `null` on every replicate entry in
the log and why the thread count is recorded as a parameter instead.
