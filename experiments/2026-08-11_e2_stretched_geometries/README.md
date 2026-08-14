# Reportable narrative

Pulling ammonia apart separates two claims that had been treated as one. Over nine nitrogen-hydrogen
distances from 0.90 to 2.20 angstrom at fixed angle, the correlation energy grows from 50.7 to 439.5
mHa while the gap to the first excited state collapses from 628.1 to 20.3 mHa, so the scan does reach
the static-correlation regime it was built to probe. The Abelian filter costs nothing anywhere in it:
the largest difference from unfiltered unitary coupled cluster over the whole scan is 7.1e-10 mHa,
while the unfiltered error itself grows a hundredfold. That refutes the claim the scan was designed
to test. The full non-Abelian filter behaves differently. Its cost stays below 2.9e-4 mHa out to 1.10
angstrom, then rises through 9.3e-3 at 1.40 and 0.12 at 1.60 to 6.94 mHa at 2.00, crossing chemical
accuracy between 1.60 and 1.80. Exact exponentials and the tied Trotter product grow together, within
a factor of two, so the cost belongs to the thirty-dimensional invariant subspace and not to
Trotterisation. At 2.20 angstrom the restricted reference stops converging, which is where the
hypothesis behind every statement here ceases to hold.

# Operational log

**Experiment ID:** 2026-08-11_e2_stretched_geometries
**Date:** 2026-08-11 (start) → 2026-08-11 (end)
**Status:** complete
**Mode:** numerical and symbolic (hybrid)
**Opportunity:** opportunities.md #3
**Target venue:** Quantum (componente do artigo principal); Comment no JCTC como item separado e posterior.

## Numerical results

NH3, STO-3G, H–N–H angle fixed at 106.7 degrees, 3136 determinants, L-BFGS-B with `ftol = 1e-16` and
`gtol = 1e-11`. The filtered pool size was recomputed at every geometry from the current irrep
labels and stayed at 75 throughout, so no orbital crossing changed the Abelian filter.

Difference from unfiltered UCCSD, in mHa. Negative means the filtered ansatz reached a marginally
lower energy.

| r (Å) | E_corr (mHa) | S0–S1 gap (mHa) | UCCSD error (mHa) | Abelian, 75 par. | full C3v, 30 par., exact | full C3v, 30 par., Trotter |
|---|---|---|---|---|---|---|
| 0.90 | 50.7 | 628.1 | 0.0992 | −1.8e-10 | −4.0e-05 | +5.1e-06 |
| 1.0124 | 65.1 | 480.4 | 0.1757 | −1.5e-10 | −5.9e-05 | +2.1e-03 |
| 1.10 | 79.1 | 385.9 | 0.2831 | −1.4e-11 | −2.9e-04 | −2.0e-05 |
| 1.25 | 110.4 | 260.4 | 0.6033 | +7.1e-12 | −3.6e-03 | +8.8e-03 |
| 1.40 | 153.2 | 171.8 | 1.1049 | −7.1e-10 | +9.3e-03 | +8.6e-02 |
| 1.60 | 230.5 | 94.0 | 2.0740 | +9.9e-11 | +1.2e-01 | +7.3e-01 |
| 1.80 | 328.6 | 46.4 | 4.1418 | 0.0 | +3.8e+00 | +6.2e+00 |
| 2.00 | 439.5 | 20.3 | 9.5875 | +4.3e-10 | +6.9e+00 | +1.3e+01 |
| 2.20 | — | — | — | RHF did not converge; point logged as a failure and discarded | | |

Largest absolute deviation over the scan: 7.1e-10 mHa for the Abelian filter, 6.94 mHa for the full
C3v filter with exact exponentials, 12.59 mHa for the same filter in tied Trotter form.

## Hypothesis check

**Ambiguous**, and the ambiguity is the finding. For the Abelian filter the hypothesis under test is
refuted: the deficit does not become more pronounced at stretched geometries; it does not appear at
all. For the full non-Abelian filter the freeness holds near equilibrium and then fails, so the
compression claim has a domain and this scan measures where it ends.

## Figures

- `figures/2026-08-11_e2_stretched_geometries_filter_cost_vs_stretch.pdf` — **For manuscript.**
  Cost of each filter against distance, above the two energy scales that define the regime. Element
  inventory in `figures/figure_manifest.json`.
- `figures/2026-08-11_e2_stretched_geometries_filter_cost_vs_stretch.png` — raster preview, for
  inspection only.

## Pointers

- Aggregate data: `results/aggregate.csv`, `results/aggregate.json`
- Visualisation contract: `results/viz_schema.json`
- Full log: `experiment_log.json` (33 run entries; 32 `success`, 1 `failure` at 2.20 Å)
- Plan, including why the exact-exponential control was added: `plan.md`
- Code: `code/run.py`, `code/analysis.py`
- Vendored sources and their hashes: `code/vendor/PROVENANCE.md`
- Environment: `env/system_info.json`, `env/requirements.txt`, `env/git_commit.txt`

## Notes

The scan as originally designed compared unfiltered UCCSD against the Abelian filter only. That
answers the question the withdrawn manuscript raised but leaves the filter this project proposes
untested in exactly the regime where it might fail. Adding the full C3v ansatz at every geometry is
what produced the result, and it changed the conclusion from a clean refutation to a refutation plus
a boundary.

Adding the exact-exponential form was needed to make the number interpretable. A cost appearing in
the tied Trotter ansatz alone could mean either that the thirty-dimensional invariant subspace cannot
represent the stretched wavefunction, or that tying and Trotterising is what costs. The two forms
grow together within a factor of two across the stretched end of the scan, which attributes the cost
to the subspace.

The failure at 2.20 Å is logged rather than worked around. Every statement in this project assumes a
symmetric closed-shell reference, and that is the distance at which one stops existing for this
molecule in this basis. Retrying with a different initial guess would have hidden the boundary rather
than found it.

Two small negative deviations, at 0.90 and 1.10 Å, have the filtered ansatz reaching a marginally
lower energy than the unfiltered one. Both are at the 1e-4 mHa level and reflect the optimiser
stopping at slightly different points on a flat surface, not a violation of the variational ordering,
which does not apply between two different ansaetze in any case.
