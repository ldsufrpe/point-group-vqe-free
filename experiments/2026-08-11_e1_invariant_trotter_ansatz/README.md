# Reportable narrative

The compression the full point group buys survives the product form a device would execute. Each
invariant operator was decomposed on the basis of distinct elementary excitations and re-exponentiated
as a tied Trotter product, in both the natural and the reversed gate order, and compared against exact
exponentials of the same operators and against unfiltered unitary coupled cluster. The decomposition
is a rewriting rather than an approximation: its largest residual over every molecule and every group
is 8.9e-16. Tying and Trotterising shifts the energy error by at most 0.0075 mHa, which is two orders
below chemical accuracy and comparable to the spread between the two gate orders. For ammonia the
full group carries 30 parameters and 147 elementary gates against 135 and 315 unfiltered, and for
methane 21 and 146 against 230 and 560. Measured instead against the Abelian filter as the literature
builds it, the parameters fall by factors of 2.5 and 3.1 while the gate count falls only from 163 to
147 for ammonia and not at all for methane. The compression is therefore a statement about parameters
first and about circuit size second, and water and ethylene, whose groups are already Abelian, return
exactly zero additional compression.

# Operational log

**Experiment ID:** 2026-08-11_e1_invariant_trotter_ansatz
**Date:** 2026-08-11 (start) → 2026-08-11 (end)
**Status:** complete
**Mode:** numerical and symbolic (hybrid)
**Opportunity:** opportunities.md #2
**Target venue:** Quantum (componente do artigo principal); Comment no JCTC como item separado e posterior.

## Numerical results

Basis STO-3G, geometries of the audited work, L-BFGS-B with `ftol = 1e-16` and `gtol = 1e-11`.
Errors are against the exact ground state of the same Hamiltonian, in mHa. "Gates" counts distinct
elementary excitation operators, not two-qubit gates.

| molecule | ansatz | parameters | gates | exact exp. | Trotter, natural | Trotter, reversed |
|---|---|---|---|---|---|---|
| H2O / C2v | unfiltered UCCSD | 65 | 140 | 0.1203 | — | — |
| H2O / C2v | Abelian, pool subset | 26 | 48 | 0.1203 | — | — |
| H2O / C2v | Abelian, invariant op. | 26 | 54 | 0.1199 | 0.1199 | 0.1144 |
| H2O / C2v | full group, invariant op. | 26 | 54 | 0.1199 | 0.1199 | 0.1144 |
| C2H4 / D2h | unfiltered UCCSD | 1224 | — | — | — | — |
| C2H4 / D2h | Abelian, invariant op. | 231 | — | — | — | — |
| C2H4 / D2h | full group, invariant op. | 231 | — | — | — | — |
| NH3 / C3v | unfiltered UCCSD | 135 | 315 | 0.1392 | — | — |
| NH3 / C3v | Abelian, pool subset | 75 | 163 | 0.1392 | — | — |
| NH3 / C3v | Abelian, invariant op. | 75 | 779 | 0.1388 | 0.1463 | 0.1434 |
| NH3 / C3v | full group, invariant op. | 30 | 147 | 0.1406 | 0.1476 | 0.1415 |
| CH4 / Td | unfiltered UCCSD | 230 | 560 | 0.1942 | — | — |
| CH4 / Td | Abelian, pool subset | 65 | 146 | 0.1942 | — | — |
| CH4 / Td | Abelian, invariant op. | 65 | 182 | 0.1939 | 0.1941 | 0.1937 |
| CH4 / Td | full group, invariant op. | 21 | 146 | 0.1940 | 0.1952 | 0.1944 |

Largest Trotterisation penalty over all rows: 0.0075 mHa, for ammonia under the Abelian filter.
Largest residual of the decomposition of an invariant operator on the elementary basis: 8.9e-16.
Determinant space: 441 for water, 3136 for ammonia, 15876 for methane; about nine million for
ethylene, which is why that row carries counts only.

The two Abelian controls return the full-group count equal to the Abelian count, 26 against 26 and
231 against 231, as a group with no degenerate irreps must.

## Hypothesis check

**Supported**, with one qualification the manuscript has to carry. The accuracy survives tying and
Trotterisation with room to spare, and the parameter count falls by the expected factors. The gate
count falls sharply against unfiltered unitary coupled cluster but only marginally against the
Abelian filter implemented as a pool subset, and for methane not at all. A claim that the full group
shortens the circuit relative to existing symmetry-adapted practice is not supported by these
numbers; a claim about parameters is.

## Figures

- `figures/2026-08-11_e1_invariant_trotter_ansatz_compression_params_gates.pdf` — **For manuscript.**
  Parameters and elementary gates for four ansaetze across four molecules. Element inventory in
  `figures/figure_manifest.json`.
- `figures/2026-08-11_e1_invariant_trotter_ansatz_compression_params_gates.png` — raster preview,
  for inspection only.

## Pointers

- Aggregate data: `results/aggregate.csv`, `results/aggregate.json`
- Visualisation contract: `results/viz_schema.json`
- Full log: `experiment_log.json` (13 run entries; 12 `success`, 1 `partial` for the ethylene
  counts-only measurement)
- Plan, including the basis-choice trap that decides the result: `plan.md`
- Code: `code/run.py`, `code/analysis.py`
- Vendored sources and their hashes: `code/vendor/PROVENANCE.md`
- Environment: `env/system_info.json`, `env/requirements.txt`, `env/git_commit.txt`

## Notes

The choice of basis for the invariant subspace decides the outcome and is not a matter of taste. An
eigensolver returns an arbitrary orthonormal basis of the eigenvalue-one subspace, and that basis is
generically dense: every invariant operator becomes a combination of hundreds of elementary
excitations. The first attempt at this measurement produced 18709 gates for ammonia against 315
unfiltered, inverting the conclusion. The pivoted basis obtained by applying the projector to
coordinate vectors and selecting independent columns by QR with column pivoting gives operators with
at most the group order in terms each, and produces the 147 in the table.

The row that surprised us is the Abelian filter built from invariant operators for ammonia: 779
gates, more than twice the unfiltered 315. Building the same 75-parameter space as a subset of the
pool instead gives 163. Both span the same space and reach the same energy to four decimal places,
so the difference is entirely in how the operators are packaged. That row is kept in the table
because it is the reason the comparison in the manuscript has to be made against the pool-subset
construction and not against the projector construction of the Abelian filter, which would flatter
the result.

The reversed gate order lowers the water error by 0.005 mHa relative to the natural order while
raising the ammonia error by a comparable amount. The two orders bound the ordering dependence and
neither is preferred; reporting only one would understate the uncertainty.

Circuit depth is not measured here. Mapping elementary excitation operators onto a device gate set
was not performed, so these counts bound circuit size in this generator basis and say nothing about
depth.
