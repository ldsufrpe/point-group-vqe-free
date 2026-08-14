# Plan — 2026-08-11_e7_degenerate_shell_artifact

**Date:** 2026-08-11
**Scale:** canonical
**Mode:** hybrid — `numerical` governs the variational errors and their distribution; `symbolic`
governs the irrep labelling and the degenerate-shell structure that drives the predictor.
**Opportunity:** `opportunities.md` #1
**Target venue:** Quantum

---

## Stage 1A — conceptual plan

### Question

Is the published non-Abelian failure a property of the method, or an artifact of the pipeline that
produced it?

### Hypothesis

**Supports "artifact":** rerunning the identical pipeline on the identical molecule in independent
processes gives a different filtered error each time, spanning the range reported in the literature,
while the unfiltered error stays fixed. And a predictor that counts degenerate orbital shells, using
no variational calculation at all, reproduces which molecules were reported as failing.

**Refutes it:** a filtered error that is stable across processes. That would mean the deviation is a
reproducible property of the filtered ansatz, and the published interpretation would stand.

### Variables

- **Independent:** the replicate index — realised as a separate operating-system process, since what
  varies is the SCF's resolution of the degenerate shell, not any seed the code controls — and the
  molecule.
- **Dependent:** the filtered and unfiltered errors, the orbital overlap between the two SCF
  calculations, and an effective misalignment angle derived from it. For the predictor: the number of
  degenerate shells and whether the prediction matches the reported outcome.
- **Controlled:** code, geometry, basis, optimiser, thresholds, machine. Nothing that the experiment
  can set is allowed to differ between replicates.

### Method

The pipeline under audit builds the Hamiltonian through `openfermionpyscf`, which fixes
`symmetry = False`, and takes irrep labels from a second, independent symmetric SCF, applying them by
index. In a partially filled degenerate shell the orbitals are defined only up to an arbitrary
unitary, so the two calculations disagree by a rotation whose value depends on the diagonalisation
routine — a gauge freedom that Sakuma et al. (2026) record independently. Each replicate reports the
overlap between the two orbital sets, which measures that rotation directly.

The unfiltered UCCSD error is the internal control: it is invariant under rotations within a
degenerate shell, so it must not move between replicates.

### Stopping criterion

Twenty independent replicates per molecule for NH3 and CH4, and the predictor over the ten benchmark
molecules in two basis sets.

### Acceptance preview

For Quantum the claim is about someone else's published result, so the standard is higher than usual:
the mechanism must be named, measured per replicate, and accompanied by a control that stays fixed.
A spread on its own would only show instability.

### Anticipated risks

- The audited third-party sources are browser transcriptions rather than a repository checkout. This
  is the experiment where a transcription divergence would matter most.
- The claim must not be extended to linear molecules, whose groups are continuous and where the
  averaging projector does not apply.

---

## Stage 1B — execution plan

### Methodological critique of the plan-mini

**The predictor's published score is not reproducible, and the honest number is smaller.** The
handoff records the structural predictor as scoring 15 of 15. The script it came from has no saved
output, and structurally it produces ten scored comparisons in STO-3G plus five rows in 6-31G that
have no per-molecule published outcome to score against. This experiment therefore reports ten of
ten scored comparisons in STO-3G and reports the 6-31G rows as predictions, with the reason stated.
Carrying 15 of 15 into a manuscript would be an unsupported number, and a referee who reran the
script would find it.

**A histogram alone is weak; a mechanism per replicate is not.** The plan-mini asked for the
distribution of the error. Recording the two-SCF orbital overlap and an effective misalignment angle
in the same replicate turns a spread into a measurement of a named cause, and lets the distribution
be checked against the independently measured relation between misalignment and error.

**Replicates must be processes, not seeds.** The nondeterminism lives in how the eigensolver resolves
a degenerate subspace. Looping inside one process would reuse the same library state and could
produce a spuriously narrow spread. Each replicate is spawned as a separate process, and the
replicate files are written individually so an interrupted run resumes rather than restarting.

**The control is what makes it an argument.** Without the unfiltered error held fixed, a spread in
the filtered error would be consistent with a badly conditioned optimisation. With it, the only
thing that changed between replicates is the gauge of the degenerate shell.

### Code reuse

| what | where | used for |
|---|---|---|
| `run(label, geom)` | `vendor/he_pipeline.py:51` | the end-to-end two-SCF pipeline, including the orbital overlap diagnostic, called unchanged |
| `he_geom_nh3`, `he_geom_ch4`, `tostr` | `vendor/he_pipeline.py:34,42,47` | the geometries of the audited work |
| `gates_of`, `irrep_of` | `vendor/audit.py:53,119` | pool and filter, through `he_pipeline` |
| `build_det_basis`, `SpinFreeHam`, `kappa_matrix`, `TrotterAnsatz` | `vendor/indep.py` | numerical core |
| benchmark molecule list and the reported-failure set | `vendor/shells.py:17,22` | the predictor's scoring set |

The predictor is rewritten here rather than imported, because `vendor/shells.py` executes at import
time and has no saved output; its molecule list and reported-failure set are reused.

### New functions

- `degenerate_shells(mo_energy, tol)` — the shell decomposition the predictor counts.
- `run_shells(log, rows)` — the predictor over ten molecules and two basis sets, scoring only where a
  published per-molecule outcome exists.
- `run_worker(rep_id, mol_key)` — one replicate, in its own process, writing its own result file.
- `spawn(n, mol_key, workers)` — process pool that skips replicates already on disk.

### Verification order

1. Each replicate's two SCF calculations agree on the total energy and the orbital energies. If they
   did not, the comparison would be between two different problems rather than two gauges of one.
2. The orbital overlap identifies which shell is misaligned, and the effective angle is derived only
   from orbitals whose overlap is not one.
3. The unfiltered error is stable across replicates. This is the control; if it moved, the spread
   would not be attributable to the labelling.
4. Only then is the spread of the filtered error interpretable as an artifact.
5. The predictor is scored only where a published outcome exists.

### File inventory

**Create:** `code/run.py`, `code/expcommon.py`, `experiment_log.json`,
`results/replicates/rep_*.json`, `results/aggregate.csv`, `results/aggregate.json`,
`results/viz_schema.json`, `figures/`, `README.md`, `env/`.
**Reuse read-only:** five vendored modules.
**Update outside this folder:** none.

### Time budget

About 33 s per NH3 replicate and a few minutes per CH4 replicate, two processes at a time. Under an
hour for both molecules.

### Known limitation carried forward

The audited third-party sources are browser transcriptions, not a repository checkout. The agreement
on eight integers and the exact reproduction of one published pool count make a relevant divergence
unlikely, but a clean checkout should be obtained before submission — it is cheap and it closes the
objection precisely where this experiment is most exposed.
