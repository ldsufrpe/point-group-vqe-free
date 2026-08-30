# Reportable narrative

The parameter compression the full point group buys does not transfer to circuit size, and
the gap is now measured rather than asserted. Every pool was resolved into elementary
spin-orbital excitations and priced under the qubit-excitation scheme of the audited work,
at two controlled-NOT gates for a single excitation and thirteen for a double. The ruler
reproduces all six published bars of their Figure 5 exactly, ammonia 3765 and 1921,
methane 6840 and 1788 and both Hamiltonian-informed pools, with no fitted quantity, and that is what licenses it to price the pool they never report. Against the Abelian filter
the full group then cuts ammonia from 1921 to 1411 controlled-NOT gates and methane from
1788 to 1554, while the parameter counts fall from 75 to 30 and from 65 to 21. The gate
saving tracks the flat operator saving instead of exceeding it, so the filter shows no
useful preference for the expensive doubles. One bookkeeping correction came out of the
same measurement. The counts have to be taken over distinct operators, which is the
currency the published figure is written in; mixing that with the incidence count is what
made an earlier reading of methane show no gate reduction at all.

# Operational log

**Experiment ID:** 2026-08-14_e8_qeb_cnot_cost
**Date:** 2026-08-14 (start) → 2026-08-14 (end)
**Status:** complete
**Mode:** symbolic
**Opportunity:** opportunities.md #2
**Specification:** `EXPERIMENT_gate_counts.md`
**Calibrated by:** `experiments/2026-08-11_e1_invariant_trotter_ansatz`
**Target venue:** Journal of Mathematical Chemistry. The per-opportunity field of
`opportunities.md` still carries the superseded string, quoted verbatim on the next line so
that the mechanical venue check compares like with like. The same file supersedes it in its
own change note of 2026-08-12, and `CLAUDE.md` records the decision.
**Superseded venue field:** Quantum (componente do artigo principal); Comment no JCTC como item separado e posterior.

## Numerical results

Basis STO-3G, geometries of the audited work, integer enumeration throughout. No energy is
computed in this experiment: no variational optimisation, no diagonalisation, no floating
point in any reported quantity. The ruler is two controlled-NOT gates per spin-orbital
single and thirteen per spin-orbital double.

Counts are over **distinct** elementary spin-orbital excitation operators, which is the
currency the published figure uses. The incidence currency is in `results/aggregate.csv`
and in the notes below.

| molecule | pool | construction | params | operators | singles | doubles | CNOTs | published |
|---|---|---|---|---|---|---|---|---|
| NH3 / C3v | unfiltered UCCSD | — | 135 | 315 | 30 | 285 | 3765 | 3765 exact |
| NH3 / C3v | Abelian filter | irrep subset | 75 | 163 | 18 | 145 | 1921 | 1921 exact |
| NH3 / C3v | Abelian filter | projector | 75 | 309 | 30 | 279 | 3687 | not reported |
| NH3 / C3v | full point group | projector | 30 | 117 | 10 | 107 | 1411 | not reported |
| NH3 / C3v | HiUCCSD | their basis | 93 | 197 | 14 | 183 | 2407 | 2407 exact |
| NH3 / C3v | HiUCCSD | adapted basis | 59 | 117 | 10 | 107 | 1411 | not reported |
| CH4 / Td | unfiltered UCCSD | — | 230 | 560 | 40 | 520 | 6840 | 6840 exact |
| CH4 / Td | Abelian filter | irrep subset | 65 | 146 | 10 | 136 | 1788 | 1788 exact |
| CH4 / Td | Abelian filter | projector | 65 | 146 | 10 | 136 | 1788 | not reported |
| CH4 / Td | full point group | projector | 21 | 128 | 10 | 118 | 1554 | not reported |
| CH4 / Td | HiUCCSD | their basis | 188 | 410 | 22 | 388 | 5088 | 5088 exact |
| CH4 / Td | HiUCCSD | adapted basis | 65 | 128 | 10 | 118 | 1554 | not reported |

What the full point group removes, measured against the Abelian filter as the literature
builds it:

| molecule | parameters | operators | CNOTs | removed singles | removed doubles |
|---|---|---|---|---|---|
| NH3 / C3v | 75 → 30 (60.0 %) | 163 → 117 (28.2 %) | 1921 → 1411 (26.5 %) | 8 | 38 |
| CH4 / Td | 65 → 21 (67.7 %) | 146 → 128 (12.3 %) | 1788 → 1554 (13.1 %) | 0 | 18 |

Representation error, the largest departure of the point-group matrices from orthogonality
in the overlap metric: 3.0e-13 for ammonia, 0.0 for methane. Largest residual of the
decomposition of an invariant operator on the elementary basis: 3.3e-16. Determinant space
3136 for ammonia, 15876 for methane, 441 for the water smoke case.

## Hypothesis check

**Refuted**, and the refutation is the useful outcome. The hypothesis was that the
full-group filter removes preferentially doubles, so that its saving in controlled-NOT
gates would materially exceed its saving in flat operator count. It does not: the two
agree to within two points in ammonia and one in methane, both directions included. What
the filter removes is optimiser dimension, and only incidentally circuit size.

The measured circuit saving is nevertheless real and is larger than the project has been
claiming. Ammonia loses about a quarter of its controlled-NOT gates and methane about an
eighth. That is the honest hardware statement available from these numbers, and it is
weaker than the parameter ratio and stronger than nothing.

Two secondary observations survive the same run. HiUCCSD in a symmetry-adapted basis lands
on exactly the full-group operator set in both molecules, 117 operators and 1411 gates
for ammonia and 128 and 1554 for methane, while carrying more parameters, 59 against 30 and
65 against 21. The full-group filter is therefore a reparameterisation of a set the
Hamiltonian already selects, once the basis is adapted. And the projector-built Abelian
pool costs 3687 gates for ammonia against 1921 for the subset construction, which is why
the comparison in this folder is made against the subset.

## Figures

- `figures/2026-08-14_e8_qeb_cnot_cost_cnot_cost.pdf` — **For manuscript.** Upper panel,
  the controlled-NOT cost of all six pools with the published values marked on the bars
  that have one; lower panel, what the full group removes relative to the Abelian filter,
  in parameters, in operators and in gates. Element inventory in
  `figures/figure_manifest.json`.
- `figures/2026-08-14_e8_qeb_cnot_cost_cnot_cost.png` — raster preview, for inspection only.

## Pointers

- Aggregate data: `results/aggregate.csv` (24 rows, one per pool and currency),
  `results/aggregate.json`
- Visualisation contract: `results/viz_schema.json`
- Full log: `experiment_log.json` (24 run entries, all `success`)
- Console transcript: `results/e8_run.log`
- Plan, including the two currencies and why the headline baseline is the subset: `plan.md`
- Code: `code/run.py`, `code/analysis.py`, `code/smoke_test.py`
- Vendored sources and the three patches applied to them: `code/vendor/PROVENANCE.md`
- Inputs and why `data/` is empty: `data/README.md`
- Environment: `env/system_info.json`, `env/requirements.txt`, `env/git_commit.txt`

## Notes

**The counting currency is not a matter of taste, and getting it wrong inverts a
conclusion.** A pool can be counted as the set of distinct elementary excitations, or as
every operator-excitation incidence in the Trotterised product. They differ because each
same-spin double appears in two singlet parameters with coefficients of opposite sign, the
antisymmetrisation of the amplitude. For ammonia the unfiltered pool has 375 incidences
against 315 distinct operators. The published figure is in the distinct currency, and this
is settled by measurement: all six of its bars are reproduced by the distinct count and
none by the incidence count.

**E1 mixed the two, and the project inherited the mixture.** E1 recorded its unfiltered
and subset rows in the distinct currency and its projector rows by incidence, so its
ammonia figure of 163 against 147 compares one currency with the other. In a single
currency the same measurement reads 163 against 117. The methane figure of 146 against 146,
no gate reduction at all, which `CLAUDE.md` carries as a house-rule number and the
manuscript carries in its scope section, is the same artefact. the distinct count of the
Abelian pool coincides numerically with the incidence count of the full-group pool. In the
distinct currency methane reads 146 against 128, and in the incidence currency 182 against
146. Neither is zero. This experiment does not edit the manuscript; the correction is
recorded here and belongs to the writer.

**The geometry is the audited work's, deliberately.** `geo_nh3` takes the tetrahedral
angle as the angle to the threefold axis, which is the convention this project identifies
elsewhere as an error in that code. Reproducing their published bars requires their
geometry, and the counts do not depend on the angle in any case, because the point group
is the same for every value of it. This is not to be repaired.

**What this does not measure.** Transpiled depth. There is no device gate set, no
connectivity graph and no scheduling here, so the numbers bound circuit size in the
qubit-excitation scheme and say nothing about the depth of a compiled circuit. The scheme
itself is the audited work's choice, adopted here for comparability with their figure and
for no other reason.
