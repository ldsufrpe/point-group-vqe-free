# Plan — 2026-08-14_e8_qeb_cnot_cost

**Date:** 2026-08-14
**Scale:** canonical
**Mode:** symbolic — the measured quantities are exact integers and correctness is exact, not
statistical. A numerical dependency (RHF, the point-group representation on the MO basis) supplies the
orbital basis the enumeration runs on, but no quantity reported here is a floating-point result.
**Opportunity:** `opportunities.md` #2
**Target venue:** Journal of Mathematical Chemistry
**Specification source:** `EXPERIMENT_gate_counts.md` (project root)
**Calibrated by:** `experiments/2026-08-11_e1_invariant_trotter_ansatz/` — this experiment re-uses E1's
validated pool construction and adds the spin-orbital class split and the CNOT ruler.

---

## Stage 1A — conceptual plan

### Question

The v3 records in Scope that "circuit depth was never measured, and no device-level resource claim
follows from" the operator counts. The headline compression is in **parameters** — 75 → 30 for NH₃ and
65 → 21 for CH₄ against the Abelian filter. What does that compression buy in **two-qubit gates**,
counted in the one scheme where the comparison can be made exactly against a published figure?

### Hypothesis

**Supports a circuit consequence:** the CNOT reduction from the Abelian filter to the full-group filter
is materially larger than the flat operator-count reduction, because the full-group filter removes
preferentially *doubles*, which cost 13 CNOTs against 2 for a single.

**Refutes it:** the CNOT reduction tracks the flat operator count, or vanishes. Then the compression is
a statement about parameters and optimiser dimension with no claim available about circuit size, and the
Scope caveat stands as the final word rather than as a placeholder.

**The refuting outcome is the more likely one, and half of it is already fixed** — see the critique.

### Variables

- **Independent:** molecule (NH₃, CH₄); pool (unfiltered UCCSD, Abelian filter, full-group filter,
  HiUCCSD); and, for the Abelian filter, the *construction* — irrep subset versus symmetric projector.
- **Dependent:** spin-orbital singles count, spin-orbital doubles count, CNOT count under the QEB ruler,
  and the reduction ratios in parameters and in CNOTs side by side.
- **Controlled:** basis (STO-3G), geometry, reference determinant, the elementary generator set, and the
  counting convention — every pool is resolved on the same set of elementary spin-orbital excitations,
  built once per molecule.

### Method

The ruler is He *et al.*'s Qubit-Excitation-Based scheme: **2 CNOTs per spin-orbital single and 13 per
spin-orbital double**, counted over operator *incidences* in the Trotterised product, not over
parameters. Parameter sharing does not reduce gates — a tied parameter driving twelve elementary
excitations still costs twelve gates.

Each pool is resolved into elementary spin-orbital excitations by machinery already validated in E1:
`gates_of` for the subset-type pools, and the Frobenius projection of each invariant operator onto the
elementary basis for the projector-built pools. Each incidence is classified as a single or a double by
the number of annihilated spin-orbitals, and the ruler is applied.

### Stopping criterion

Two molecules, four pools each, plus the validation gate against the six published bars of Figure 5.
Done when the UCCSD and SymUCCSD CNOT counts reproduce the published values **exactly** and the
full-group row is measured.

### Acceptance preview

For the JMC the claim is mathematical, so what has to hold is the exactness: the counting convention
reproduces published integers with no fitted parameter, and the new row follows by the same convention.
A CNOT reduction that is merely *reported* is not enough — it has to be attributed to which operator
class was removed, because that is the part that generalises beyond two molecules.

This experiment can close nothing and still be worth having. If the full-group filter buys no CNOTs,
that is the honest bound on the paper's own claim, and it converts the Scope caveat from an admission of
ignorance into a measured statement.

### Anticipated risks

- **The CH₄ row is decided before the run.** Risk of writing the experiment to look like discovery. The
  README must state that E1 already fixed it.
- **The QEB ruler is theirs, not a universal.** Right for comparability with Figure 5, wrong for a
  hardware claim. Transpiled depth stays out of scope by decision, and the README has to say so.

---

## Stage 1B — execution plan

### Methodological critique of the plan-mini

> **Post-run note, 2026-08-14.** Point 1 below is pre-registration and the run overturned it. It
> reasons from E1's `163 → 147` for NH₃ and `146 → 146` for CH₄, and those two comparisons each
> cross the distinct operator count with the incidence count. Re-counted in one currency the same
> pools read **163 → 117** and **146 → 128**; methane's "exactly zero" was the coincidence of a
> distinct count with an incidence count at 146. The headline decision recorded in point 2 stands
> — the baseline is still the irrep subset, not the projector construction — but its numbers are
> the corrected ones. The text below is left as written, because it is what was decided before the
> measurement. `README.md` carries the outcome.

Reading E1's `results/aggregate.csv` before writing any code surfaced two things that change what the
experiment can honestly claim. Both are settled decisions, not open questions.

**1. Half the result is already determined, and the experiment must not pretend otherwise.** E1 records
flat operator counts of 163 → 147 for NH₃ and **146 → 146 for CH₄**, Abelian irrep-subset against
full-group projector. CH₄ shows a 3.1× parameter compression (65 → 21) and *exactly zero* operator
reduction. Under any positive weights, zero removed operators is zero CNOTs saved, so the CH₄ answer is
fixed before the run and this experiment confirms rather than discovers it. The live question is NH₃,
where 16 operators are removed and the saving depends entirely on their class: between **1.7%** (all
singles) and **10.8%** (all doubles), a six-fold spread the flat count cannot resolve. That spread is
the reason the experiment is worth running.

**2. There are two Abelian baselines and they differ by a factor of five.** E1's `abelian_subset` row is
163 operators; its `abelian` row — same subgroup, built through the symmetric projector instead of by
selecting symmetric parameters — is **779**. Comparing 163 to 147 crosses constructions. The
like-for-like comparison is 779 → 147, which would look spectacular and would be misleading, because
nobody implements the Abelian filter that way.

**Decision, confirmed by the researcher on 2026-08-14: the headline stays 163 → 147.** 163 is what the
state of the art actually costs. The projector-built Abelian row is measured and reported beside it,
with the construction named, so the reader can see that the favourable-looking number was available and
declined. This is the project's house rule — no figure without its baseline and its variant — applied to
a new quantity.

**3. The ruler is confirmed before any code runs.** Solving 2s + 13d = C against s + d = N for the four
published UCCSD/SymUCCSD bars returns integer, physically correct singles counts in every case:

| bar | N (published) | C (published) | s implied | check |
|---|---|---|---|---|
| NH₃ UCCSD | 315 | 3765 | 30 | = 2 × 15 spatial singles |
| NH₃ SymUCCSD | 163 | 1921 | 18 | even, as spin complementation requires |
| CH₄ UCCSD | 560 | 6840 | 40 | = 2 × 20 spatial singles |
| CH₄ SymUCCSD | 146 | 1788 | 10 | even |

Four independent integer solutions from a two-equation system is not something a wrong convention
produces. The experiment reproduces these by construction rather than by back-solving, which is the
actual gate.

**4. The geometry is He *et al.*'s, deliberately.** `geo_nh3` in `group_compression.py` uses the same
`theta/2`-to-the-C₃-axis convention as `he_geom_nh3`, which is the convention this project elsewhere
identifies as their error. That is correct here and must not be "fixed": reproducing their published
counts requires their geometry, and the counts are invariant to the angle anyway, because the point
group is C₃ᵥ for any value of it. Worth stating in the README so a later reader does not repair it.

### Code inspection and reuse

All reuse is vendored into `code/vendor/` with a `PROVENANCE.md`, following E1's pattern, so the folder
reproduces on its own. Source: `/home/leon-denis/PycharmProjects/Experimentos/auditoria_independente/`.

| what | where | used for |
|---|---|---|
| `gates_of(norb, nelec)` | `vendor/audit.py:53` | the UCCSD pool resolved into elementary spin-orbital excitations, tagged by parameter — the load-bearing piece |
| `irrep_of(g, orbsym)` | `vendor/audit.py:119` | the irrep-subset Abelian filter, i.e. SymUCCSD as the literature implements it |
| `sparse_invariant_basis(Ms, nS)` | `vendor/e1_trotter.py:34` | the pivoted `{P e_j}` basis; never the dense eigenbasis, which inflated NH₃ to 18709 gates on E1's first attempt |
| `build_ops(Bs, Bd, Ep, nS)` | `vendor/group_compression.py:152` | anti-Hermitian invariant operators |
| `group`, `mo_rep`, `to_str`, `CASES` | `vendor/group_compression.py:76,108,66,253` | representation, projector, geometries |
| `build_det_basis`, `SpinFreeHam`, `kappa_matrix` | `vendor/indep.py:24,82,180` | determinant basis, the `E_pq` operators, and the elementary generators |
| `ham_term_keys(mol, C, tol)` | `vendor/hi_count.py:23` | the HiUCCSD pool, for the two secondary bars |
| `expcommon.py` | `../2026-08-11_e1_invariant_trotter_ansatz/code/expcommon.py` | env capture, append-only log, CSV writer, figure style already calibrated for `svjour3 [smallextended]` |

**Two patches are required and both are mechanical.** `hi_count.py:18` hard-codes `HERE` to a scratchpad
that no longer exists, and `he_pipeline.py:28` hard-codes an absolute path; both become
`os.path.dirname(os.path.abspath(__file__))`, which is the same patch E1 already applied to its four
vendored modules.

**One deliberate departure from E1's path.** E1's `run` calls `eigsh` for the FCI energy and then runs
three optimisations per pool. None of that is needed here — no quantity reported is an energy — so the
decomposition block is re-used and the optimisation block is dropped. `SpinFreeHam` is still constructed
in full rather than extracting only the `E_pq` I need, because using the validated constructor verbatim
costs seconds and removes any chance of the operators drifting from E1's.

### New functions

- `pool_uccsd(nao, nelec)` — the unfiltered pool as a list of `(occ, virt)` incidences.
- `pool_irrep_subset(gates, orbsym)` — the Abelian filter by irrep label.
- `pool_projector(Ms, sub, Ep, Kel, nrm, nS)` — the projector-built pool for a given subgroup, returning
  incidences plus the maximum decomposition residual.
- `pool_hiuccsd(mol, C, gates)` — the HiUCCSD pool, in the adapted and non-adapted bases.
- `classify(incidences)` → `(n_singles, n_doubles)` by `len(occ)`.
- `cnots(n_singles, n_doubles)` → `2 * n_singles + 13 * n_doubles`.
- `check_published(rows)` — the validation gate; asserts the four UCCSD/SymUCCSD bars match exactly and
  reports the HiUCCSD bars without asserting.

### Verification order

Dependency-ordered. A failure at step *n* makes every later step meaningless.

1. **Representation validity** — `mo_rep` returns `|UᵀSU − S|` at machine precision. If the group
   representation is wrong, every pool below it is wrong.
2. **Incidence bookkeeping** — for the UCCSD pool, the number of incidences equals the number of
   *distinct* `(occ, virt)` pairs. E1 reports 315 for NH₃ from both counts; if they diverge, the two
   pools are being counted in different currencies and the ruler cannot be applied uniformly.
3. **Class split is exhaustive** — every incidence has `len(occ)` equal to 1 or 2. Anything else means
   the pool contains an operator the ruler has no price for.
4. **Parameter counts reproduce E1** — 135/75/30 for NH₃ and 230/65/21 for CH₄, and the flat operator
   counts 315/163/147 and 560/146/146. This is a regression check against a folder already validated.
5. **Decomposition residual** — each invariant operator reproduced on the elementary basis to machine
   precision, as E1 reports (≤ 6.7e-16). Only then is the projector pool an exact rewriting.
6. **The published gate** — UCCSD and SymUCCSD CNOTs equal 3765/1921 (NH₃) and 6840/1788 (CH₄)
   **exactly**. If these fail, the counting convention is wrong and nothing downstream is worth
   reporting. HiUCCSD (2407, 5088) is reported but not asserted: their pool depends on a
   `compress(1e-8)` threshold and the basis, so divergence there is information rather than error.
7. Only then does the full-group CNOT count mean anything.

### File inventory

**Create:** `code/run.py`, `code/analysis.py`, `code/smoke_test.py`, `code/expcommon.py`,
`code/vendor/{audit,e1_trotter,group_compression,indep,hi_count,he_pipeline}.py`,
`code/vendor/PROVENANCE.md`, `data/README.md`, `experiment_log.json`, `results/aggregate.csv`,
`results/aggregate.json`, `results/viz_schema.json`, `results/e8_run.log`, `figures/` with its
`figure_manifest.json` and `README.md`, `env/`, `README.md`, `report_validation.md`.

**Reuse read-only:** the six vendored modules at their source, `expcommon.py` from E1.

**Update outside this folder:** none. `EXPERIMENT_gate_counts.md` at the project root is the
specification and stays read-only; recording the measured numbers into it is the writer's job, not this
skill's.

### Time budget

No VQE, no diagonalisation, no optimisation. The costs are the determinant basis and the elementary
generators: NH₃ 3136 determinants and 315 generators, CH₄ 15876 and 560. E1 spent 207 s and 220 s on
these molecules *including* six optimisations each; dropping those should leave well under a minute per
molecule. Budget: **under five minutes total**, and if it exceeds fifteen the design is wrong and I
should stop and say so.
