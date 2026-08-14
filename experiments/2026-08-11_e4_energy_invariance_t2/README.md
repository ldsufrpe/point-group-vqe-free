# Reportable narrative

The energy of a tied product ansatz is unchanged when the characters of the symmetry group flip the
signs of its parameters, and that identity had never been measured. Twenty parameter vectors per
molecule were drawn uniformly on the interval from minus one half to one half, deliberately away
from the symmetric submanifold where the sign flip acts trivially and the test would be empty, and
both energies were evaluated through the same Trotterised circuit on the determinant basis of
ammonia in its reflection subgroup and methane in its dihedral subgroup. Across the eighty samples
the largest difference was 5.3e-12 Ha, and the stronger state-level identity held to exactly zero;
at the same time the gradient along the removed directions stayed below 2.2e-10 Ha everywhere on the
submanifold. Rotating one degenerate orbital pair by thirty degrees while keeping its irrep labels,
which is the misalignment a second independent self-consistent field introduces, raised the same
difference to 8.1e-3 Ha, nine orders of magnitude higher. The invariance therefore holds to machine
precision wherever its hypotheses hold, and the measurement can detect their failure, which is what
separates the identity from an arithmetic tautology.

# Operational log

**Experiment ID:** 2026-08-11_e4_energy_invariance_t2
**Date:** 2026-08-11 (start) → 2026-08-11 (end)
**Status:** complete
**Mode:** numerical and symbolic (hybrid)
**Opportunity:** opportunities.md #3
**Target venue:** Quantum (componente do artigo principal); Comment no JCTC como item separado e posterior.

## Numerical results

| quantity | NH3 / Cs | CH4 / D2 |
|---|---|---|
| determinants | 3136 | 15876 |
| parameters, total / symmetric / removed | 135 / 75 / 60 | 230 / 65 / 165 |
| nontrivial group elements | 1 | 3 |
| parameter vectors sampled | 20 | 20 |
| max abs energy difference (Ha) | 5.3e-12 | 3.6e-15 |
| max relative energy difference | 9.5e-14 | 9.1e-17 |
| max state-level residual | 0.0 | 0.0 |
| max gradient in removed directions (Ha) | 2.2e-10 | 3.1e-15 |
| control, degenerate pair rotated 30 deg (Ha) | 5.4e-3 | 8.1e-3 |

The reference determinant is a fixed point of every element's determinant-space sign operator, to
exactly zero, so hypothesis (A2) is satisfied rather than assumed. No parameter mixes irreps, which
is hypothesis (A3), asserted inside the pool builder. Every sampled vector had a component outside
the symmetric submanifold of magnitude at least 0.49, so the sign flip acted nontrivially in all 80
samples.

The gradient figure of 2.2e-10 Ha for ammonia reproduces the 2.151e-10 recorded independently by the
earlier audit at the Hartree-Fock point, which is a cross-check on the vendored numerical core
rather than a new measurement.

## Hypothesis check

**Supported.** The invariance holds at machine precision, uniformly over the sampled parameter
magnitudes and over every group element, and the state-level identity holds exactly. The negative
control separates from the adapted-basis measurements by nine orders of magnitude, so a systematic
violation of the size the withdrawn literature reports would have been visible.

## Figures

- `figures/2026-08-11_e4_energy_invariance_t2_invariance_violation.pdf` — **For manuscript.**
  Measured violation of the invariance identity in the adapted basis against the rotated-label
  control, with chemical accuracy marked. Element inventory in `figures/figure_manifest.json`.
- `figures/2026-08-11_e4_energy_invariance_t2_invariance_violation.png` — raster preview of the same
  figure, for inspection only.

## Pointers

- Aggregate data: `results/aggregate.csv`, `results/aggregate.json`
- Visualisation contract: `results/viz_schema.json`
- Full log: `experiment_log.json` (90 run entries, all `success`)
- Plan, including the methodological critique that changed the design: `plan.md`
- Code: `code/run.py`, `code/analysis.py`, `code/smoke_test.py`
- Vendored sources and their hashes: `code/vendor/PROVENANCE.md`
- Environment: `env/system_info.json`, `env/requirements.txt`, `env/git_commit.txt`

## Notes

The design changed in Stage 1B for a reason worth recording. Sampling parameter vectors inside the
symmetric submanifold, as the earlier evidence for the neighbouring corollary had done, makes the
sign flip act as the identity: the measured difference is then exactly zero for a reason unrelated
to the statement under test. Forcing the sampled vectors off the submanifold, and recording how far
off each one sits, is what turns the check into evidence.

Two additions beyond the original specification cost almost nothing and changed what the experiment
establishes. The state-level residual is stronger than the energy difference and closes the
equivariance lemma, which also lacked a direct check. The rotated-label control makes the result
discriminating; without it a reader cannot tell a theorem from a tautology.

Methane's differences sit three orders of magnitude below ammonia's. The likely cause is that its
three nontrivial elements each flip a larger share of the 230 parameters, so cancellation in the
accumulated floating-point error is more complete; this was not investigated further because both
figures are far below any threshold that matters.
