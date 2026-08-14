# Reportable narrative

The freeness of the point-group filter is not an artifact of the minimal basis, and it holds for
every molecule He *et al.* name as failing. Their Table 3 is a 6-31G table, and its text states that
under that basis SymUCCSD fails systematically for non-Abelian molecules, even for those that
succeeded with STO-3G, listing HF, LiH, BeH2 and NH3. SymUCCSD is the Abelian-subgroup filter. That
filter is measured here, on all four of those molecules at 6-31G, to cost at most 2.0e-10 mHa, ten
orders of magnitude below chemical accuracy. Six molecules and two bases give thirty-six runs and no
deficit anywhere.

Where the full point group is finite, the stronger filter also stays free while compressing much
further: ammonia falls from 1325 amplitude classes to 277 under C3v and methane from 1890 to 144
under Td, approaching the 1/h asymptote from above without reaching it. Water, already Abelian,
returns identical counts under both filters. For the three linear molecules the full group is
continuous and the group-average projector does not apply, so the full-group column there reports
the largest finite subgroup and coincides with the Abelian one by construction.

The measurement is classical coupled cluster rather than a variational quantum eigensolver, a
distinction stated in the limitations below. What it establishes is that neither the basis, nor the
filter, nor the amplitude algebra produces the reported deficit, in the exact basis and on the exact
molecules where the deficit was reported.

# Operational log

**Experiment ID:** 2026-08-12_e0b_basis_dependence
**Date:** 2026-08-12
**Status:** complete
**Mode:** numerical and symbolic (hybrid)
**Parent:** `experiments/2026-08-11_e0_classical_ccsd_fullgroup/` (extends)
**Opportunity:** opportunities.md #1 and #2
**Target venue:** Quantum (componente do artigo principal); Comment no JCTC como item separado e posterior.
**Runs:** 36 (36 success, 0 failure)

## Numerical results

Amplitude classes surviving each filter, and the correlation-energy cost of each filter measured
against unconstrained CCSD in the same basis. RHF reference in the adapted basis,
`conv_tol = 1e-10`.

| basis | molecule | full group | finite group used | PySCF Abelian | none | Abelian | full | delta Abelian (mHa) | delta full (mHa) |
|---|---|---|---|---|---|---|---|---|---|
| STO-3G | HF | C∞ᵥ * | C2v | Coov | 20 | 11 | 11 | -2.2e-11 | -2.2e-11 |
| STO-3G | LiH | C∞ᵥ * | C2v | Coov | 44 | 20 | 20 | -8.7e-14 | -8.7e-14 |
| STO-3G | BeH₂ | D∞ₕ * | D2h | Dooh | 90 | 23 | 23 | -4.6e-13 | -4.6e-13 |
| STO-3G | H₂O | C₂ᵥ | C2v | C2v | 65 | 26 | 26 | 3.1e-12 | 3.1e-12 |
| STO-3G | NH₃ | C₃ᵥ | C3v | Cs | 135 | 75 | 30 | 2.6e-12 | 2.5e-09 |
| STO-3G | CH₄ | T_d | Td | D2 | 230 | 65 | 21 | 1.4e-13 | 1.5e-13 |
| 6-31G | HF | C∞ᵥ * | C2v | Coov | 495 | 178 | 178 | -4.4e-13 | -4.4e-13 |
| 6-31G | LiH | C∞ᵥ * | C2v | Coov | 189 | 85 | 85 | 1.5e-10 | 1.5e-10 |
| 6-31G | BeH₂ | D∞ₕ * | D2h | Dooh | 495 | 125 | 125 | 2.0e-10 | 2.0e-10 |
| 6-31G | H₂O | C₂ᵥ | C2v | C2v | 860 | 281 | 281 | -3.4e-12 | -3.4e-12 |
| 6-31G | NH₃ | C₃ᵥ | C3v | Cs | 1325 | 717 | 277 | 9.7e-12 | 3.1e-12 |
| 6-31G | CH₄ | T_d | Td | D2 | 1890 | 495 | 144 | 5.6e-14 | -1.8e-13 |

\* Continuous group. The group-average projector |G|⁻¹Σ_g is defined for finite groups only, so for
the three linear molecules the "finite group used" column is the largest finite subgroup available
and the full-group filter coincides with the Abelian one by construction. This is the scope boundary
declared in `PROVAS.md` under "Fronteira do escopo"; it is not a result. **What these three rows
test is the Abelian filter — which is exactly what SymUCCSD implements and what He *et al.* report
as failing at 6-31G.**

Full-group invariance residual of the **converged unconstrained** amplitudes — Čársky's vanishing
rule measured rather than assumed — is at most 8.3e-14 in t₂ and 1.8e-12 in t₁ across all six rows.
The filter therefore removes classes that already carry no weight, rather than discarding weight the
wavefunction needed.

The independent combinatorial oracle agrees with the character-based count on every Abelian
subgroup, in both bases: 4+22=26, 9+66=75, 5+60=65, 16+265=281, 31+686=717, 15+480=495.

## Continuity check

The three STO-3G rows reproduce `experiments/2026-08-11_e0_classical_ccsd_fullgroup/` **digit for
digit**: `max_abs_delta_e_corr_mha` is `2.493838469064258e-09` in both experiment logs, and the
amplitude counts 65→26, 135→30, 230→21 are unchanged. This validates the parameterisation of the
basis, which is the only edit made to the E0 code path.

## Hypothesis assessment

**Supported.** The pre-registered falsifier — a deficit above chemical accuracy for NH₃ or CH₄ in
6-31G — did not occur; the largest deficit observed is 9.7e-12 mHa, eleven orders below it.

## Limitations, stated plainly

- **This is classical coupled cluster, not VQE.** The amplitude space of CCSD and the parameter
  space of singlet UCCSD are the same space, which is why the class counts transfer, but the ansatz
  structure and the optimiser do not. A referee is entitled to ask for the unitary measurement, and
  E1 supplies it at STO-3G only. The 6-31G unitary case is untested.
- **Six molecules, two bases; the comparison with He *et al.* is closed for the molecules they
  name.** HF, LiH, BeH₂ and NH₃ are all present at 6-31G. What is *not* closed is the group: for the
  three linear molecules the full point group is continuous, so only the Abelian filter is tested
  there. That is the filter in dispute, but no full-group compression claim can be made for them.
- **Correlation-consistent bases untested.** 6-31G is a split-valence basis; nothing here speaks to
  cc-pVDZ or larger, where the virtual space is much richer.
- The compression ratios approach the 1/h asymptote from above and never reach it, consistent with
  Greiner, Gauss & Eriksen (JPCL 2024), who state that the bound is met only when every orbital
  belongs to a symmetry-equivalent set of cardinality h.

## Reproduction

    python code/run.py all

Smoke test: `code/run.py --smoke` (H₂O, STO-3G only) — exits 0 and reproduces the E0 water row.
Environment snapshot in `env/`; dependency pins in `env/requirements.txt`; vendored audit sources and
their hashes in `code/vendor/PROVENANCE.md`.
