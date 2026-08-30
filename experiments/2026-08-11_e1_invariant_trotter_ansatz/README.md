# Reportable narrative

The compression the full point group buys survives the product form a device would execute. Each
invariant operator was decomposed on the basis of distinct elementary excitations and re-exponentiated
as a tied Trotter product, in both the natural and the reversed gate order, and compared against exact
exponentials of the same operators and against unfiltered unitary coupled cluster. The decomposition
is a rewriting rather than an approximation: its largest residual over every molecule and every group
is 8.9e-16. Tying and Trotterising shifts the energy error by at most 0.0075 mHa, which is two orders
below chemical accuracy and comparable to the spread between the two gate orders. Counted over the
distinct elementary excitation operators a pool touches, ammonia under the full group carries 30
parameters and 117 operators against 135 and 315 unfiltered, and methane 21 and 128 against 230 and
560. Measured instead against the Abelian filter as the literature builds it, the parameters fall by
factors of 2.5 and 3.1 and the operator count falls from 163 to 117 for ammonia and from 146 to 128
for methane. The compression is therefore a statement about parameters first and about circuit size
second, and water and ethylene, whose groups are already Abelian, return exactly zero additional
compression.

# Operational log

**Experiment ID:** 2026-08-11_e1_invariant_trotter_ansatz
**Date:** 2026-08-11 (start) → 2026-08-11 (end)
**Status:** complete
**Mode:** numerical and symbolic (hybrid)
**Opportunity:** opportunities.md #2
**Target venue:** Quantum (componente do artigo principal); Comment no JCTC como item separado e posterior.

## Numerical results

Basis STO-3G, geometries of the audited work, L-BFGS-B with `ftol = 1e-16` and `gtol = 1e-11`.
Errors are against the exact ground state of the same Hamiltonian, in mHa. No count here is a
two-qubit gate count.

**The operator column of this run is not in one currency, and the `currency` column says which one
each row is in.** A pool admits two counts: the *distinct* elementary spin-orbital excitation
operators it touches, and the *incidences* of those operators in the Trotterised product. They
differ because each same-spin double appears in two singlet parameters with coefficients of opposite
sign. This run took its `uccsd` and pool-subset rows from a counter that returns the first and its
projector rows from one that returns the second. The single-currency numbers are in the block below
this table and they are what the manuscript quotes.

| molecule | ansatz | parameters | operators | currency | exact exp. | Trotter, natural | Trotter, reversed |
|---|---|---|---|---|---|---|---|
| H2O / C2v | unfiltered UCCSD | 65 | 140 | distinct | 0.1203 | — | — |
| H2O / C2v | Abelian, pool subset | 26 | 48 | distinct | 0.1203 | — | — |
| H2O / C2v | Abelian, invariant op. | 26 | 54 | incidence | 0.1199 | 0.1199 | 0.1144 |
| H2O / C2v | full group, invariant op. | 26 | 54 | incidence | 0.1199 | 0.1199 | 0.1144 |
| C2H4 / D2h | unfiltered UCCSD | 1224 | — | — | — | — | — |
| C2H4 / D2h | Abelian, invariant op. | 231 | — | — | — | — | — |
| C2H4 / D2h | full group, invariant op. | 231 | — | — | — | — | — |
| NH3 / C3v | unfiltered UCCSD | 135 | 315 | distinct | 0.1392 | — | — |
| NH3 / C3v | Abelian, pool subset | 75 | 163 | distinct | 0.1392 | — | — |
| NH3 / C3v | Abelian, invariant op. | 75 | 779 | incidence | 0.1388 | 0.1463 | 0.1434 |
| NH3 / C3v | full group, invariant op. | 30 | 147 | incidence | 0.1406 | 0.1476 | 0.1415 |
| CH4 / Td | unfiltered UCCSD | 230 | 560 | distinct | 0.1942 | — | — |
| CH4 / Td | Abelian, pool subset | 65 | 146 | distinct | 0.1942 | — | — |
| CH4 / Td | Abelian, invariant op. | 65 | 182 | incidence | 0.1939 | 0.1941 | 0.1937 |
| CH4 / Td | full group, invariant op. | 21 | 146 | incidence | 0.1940 | 0.1952 | 0.1944 |

### The same pools counted in one currency

Every pool of this experiment, over distinct elementary excitation operators. The ammonia and
methane rows are E8's re-count of these pools
(`experiments/2026-08-14_e8_qeb_cnot_cost/results/aggregate.csv`, rows with `currency=distinct`).
E8 did not run water; its projector pool is a basis of the invariant subspace, and the union of the
supports of any basis of a subspace is the support of the subspace, so its count is
`|supp P_G| = 48`, which is the 48 the `abelian_subset` row of this experiment already reports.

| molecule | unfiltered | Abelian, pool subset | Abelian, invariant op. | full group |
|---|---|---|---|---|
| H2O / C2v | 140 | 48 | 48 | 48 |
| NH3 / C3v | 315 | 163 | 309 | 117 |
| CH4 / Td | 560 | 146 | 146 | 128 |

Against the Abelian filter as the literature builds it, the full group therefore removes 163 → 117
operators for ammonia and 146 → 128 for methane, while the parameters fall 75 → 30 and 65 → 21. E8
prices the same pools in controlled-NOT gates under the scheme of the audited work: 1921 → 1411 and
1788 → 1554.

Largest Trotterisation penalty over all rows: 0.0075 mHa, for ammonia under the Abelian filter.
Largest residual of the decomposition of an invariant operator on the elementary basis: 8.9e-16.
Determinant space: 441 for water, 3136 for ammonia, 15876 for methane; about nine million for
ethylene, which is why that row carries counts only.

The two Abelian controls return the full-group count equal to the Abelian count, 26 against 26 and
231 against 231, as a group with no degenerate irreps must.

## Hypothesis check

**Supported**, with one qualification the manuscript has to carry. The accuracy survives tying and
Trotterisation with room to spare, and the parameter count falls by the expected factors. The
operator count falls sharply against unfiltered unitary coupled cluster and by a smaller factor
against the Abelian filter implemented as a pool subset, 28.2 per cent for ammonia and 12.3 per cent
for methane. A claim that the full group shortens the circuit relative to existing symmetry-adapted
practice is supported at that size and no larger, and the parameter claim is the stronger of the two.

The first version of this section read the operator saving as marginal for ammonia and absent for
methane. That reading came from the mixed currency and not from the pools: it compared a distinct
count with an incidence count, and methane's apparent zero was the coincidence of the two numbers at
146. Circuit depth is still not measured here.

## Figures

- `figures/2026-08-11_e1_invariant_trotter_ansatz_compression_params_gates.pdf` — **For manuscript.**
  Parameters and distinct elementary excitation operators for four ansaetze across four molecules.
  Regenerated on 2026-08-14: the lower panel now draws the single-currency counts of the block
  above, which `code/analysis.py` carries as a table with its provenance, because this experiment's
  own `n_gates` column mixes the two currencies. The upper panel is unchanged. Element inventory in
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
at most the group order in terms each, and produces the 147 incidences, 117 distinct operators, in
the tables above.

The row that surprised us is the Abelian filter built from invariant operators for ammonia: 779
incidences, more than twice the unfiltered 315, and 309 distinct operators once the currency is
fixed. Building the 75-parameter space as a subset of the pool instead gives 163 distinct operators.

**That pair of numbers is an open item, and the corrected currency is what exposes it.** Only 163 of
the 315 elementary excitations carry the trivial character of the labelled Abelian subgroup, so a
pool invariant under *that* subgroup cannot touch 309 of them. The projector row therefore is not
the same pool as the subset row, whatever the two energies suggest, and 0.1388 against 0.1392 says
they are not the same space either. Two readings fit the arithmetic and this experiment does not
separate them: the group matrices may implement a reflection plane other than the one PySCF's orbital
labels use, which is the misalignment mechanism E7 measures elsewhere in this project, or the
coefficient cut of 1e-12 may be admitting representation noise. Methane and water show no such gap,
their projector and subset counts agreeing exactly, which is what a resolution has to explain as
well. Until it is resolved the row stands as the reason the comparison in the manuscript is made
against the pool-subset construction, which is what the literature implements.

The reversed gate order lowers the water error by 0.005 mHa relative to the natural order while
raising the ammonia error by a comparable amount. The two orders bound the ordering dependence and
neither is preferred; reporting only one would understate the uncertainty.

Circuit depth is not measured here. Mapping elementary excitation operators onto a device gate set
was not performed, so these counts bound circuit size in this generator basis and say nothing about
depth.
