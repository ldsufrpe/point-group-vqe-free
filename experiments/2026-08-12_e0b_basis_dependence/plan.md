# Plan — E0b: basis-set dependence of the point-group filter

**Experiment ID:** 2026-08-12_e0b_basis_dependence
**Parent:** `experiments/2026-08-11_e0_classical_ccsd_fullgroup/` (extends), opportunities.md #1 and #2
**Mode:** numerical and symbolic (hybrid), canonical scale
**Target venue:** Quantum

## Why this experiment exists

Two independent reasons, one defensive and one offensive.

**Defensive.** `acceptance_report.md` records, under Notes: *"Every energy claim rests on STO-3G
alone. Rerunning E0 in 6-31G for NH3 and CH4 costs 0.1-0.3 s per row and is the cheapest available
strengthening."* A single-basis result is the cheapest objection a referee can raise.

**Offensive, and it is the stronger reason.** He *et al.* (JCTC 2026, 22, 6008), Table 3, is a
**6-31G** table, and the accompanying text states: *"under the 6-31G basis set, SymUCCSD retains
high accuracy for Abelian molecules but fails systematically for non-Abelian molecules, even for
those molecules (HF, LiH, BeH2, NH3) that performed successfully with the STO-3G basis set."*

So the failure they report is **basis-dependent and worsens in 6-31G** — precisely the regime E0 did
not test. Measuring the classical filter at 6-31G lands on the claim where it is strongest, not
where it is weakest.

## Hypothesis

Confining coupled-cluster amplitudes to the totally symmetric subspace of the full molecular point
group leaves the correlation energy unchanged **independently of the one-electron basis**. If the
freeness were an artifact of the minimal basis, the deficit would appear on going to 6-31G.

Falsifiable and pre-registered: a deficit above chemical accuracy (1.6 mHa) for NH₃ or CH₄ in 6-31G
would refute the hypothesis and would corroborate He *et al.*

## Design

Six molecules and the same code path as E0, with the basis lifted to a parameter. The first three
are E0's set; the last three complete the list He *et al.* name.

| molecule | full group | Abelian subgroup PySCF uses | role |
|---|---|---|---|
| H₂O | C₂ᵥ | C₂ᵥ | control: group already Abelian, counts must coincide |
| NH₃ | C₃ᵥ | C_s | named by He *et al.* as failing at 6-31G; full group finite |
| CH₄ | T_d | D₂ | largest finite group, strongest compression |
| HF | C∞ᵥ | Coov | named by He *et al.*; **continuous** group — Abelian filter only |
| LiH | C∞ᵥ | Coov | named by He *et al.*; **continuous** group — Abelian filter only |
| BeH₂ | D∞ₕ | Dooh | named by He *et al.*; **continuous** group — Abelian filter only |

The three linear molecules complete the comparison against the exact list He *et al.* give (HF, LiH,
BeH₂, NH₃). Their full point groups are continuous, so the group-average projector does not apply
and the full-group compression claim is out of scope for them — that boundary is declared in
`PROVAS.md`. What they do test, and what matters for the refutation, is the **Abelian** filter,
which is what SymUCCSD implements.

Bases: `sto-3g` (reproduces E0, continuity check) and `6-31g`. Three filters per row: none
(unconstrained CCSD), Abelian, full group. RHF reference in the adapted basis, `conv_tol = 1e-10`,
`conv_tol_normt = 1e-8`.

## Measured quantities

- amplitude classes surviving each filter (singles, doubles, total)
- correlation energy per filter, and the delta against unconstrained CCSD in mHa
- full-group invariance residual of the **converged unconstrained** amplitudes — this is Čársky's
  vanishing rule measured directly, not assumed
- independent combinatorial oracle for the Abelian counts, as a check on the character machinery
- representation validation |UᵀSU − S| and the max off-shell element of D

## Why a new experiment directory rather than editing E0

E0 is complete, was validated by the Phase 5 acceptance run, and carries provenance pointing at
commit `649d5d79`. Its `experiment_log.json` is append-only and its `report_validation.md` was
generated against the STO-3G state. Re-running inside E0 would either duplicate its log entries or
overwrite its `aggregate.csv`, and would make its README stale. E0b is self-contained: it runs both
bases, so the comparison lives inside one experiment and needs no cross-experiment join.

## Code reuse

`code/` is copied from E0 and the basis is lifted to a parameter in four places: the `gto.M` call,
the two `log.record` parameter dicts, and a new `basis` column in `FIELDS` and in the row dicts.
`main()` gains a `BASES` sweep. `code/vendor/` carries the audit sources unchanged, with hashes in
`vendor/PROVENANCE.md`. No scientific logic was touched — the STO-3G rows must reproduce E0 exactly,
and that reproduction is itself a check.
