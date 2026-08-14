# Reportable narrative

Confining coupled-cluster amplitudes to the totally symmetric subspace of the full molecular point
group leaves the correlation energy where it was. PySCF cannot be asked to do this: it detects the
full group and then works in the largest Abelian subgroup it supports, reducing ammonia from its
threefold group to a single reflection and methane from its tetrahedral group to three orthogonal
twofold axes. The projector was therefore built here from the Cartesian action of each group
element on the molecular orbitals and applied to the singles and doubles amplitudes at every
iteration. Under the full group the amplitude classes fall from 65 to 26 for water, from 135 to 30
for ammonia and from 230 to 21 for methane, and the correlation energy moves by at most 2.5e-9 mHa,
nine orders of magnitude below chemical accuracy. The converged unfiltered amplitudes already
satisfy the full-group invariance to 1.1e-12, so the filter removes classes that carry no weight
rather than discarding weight the wavefunction needed. Water, whose group is already Abelian, gives
identical counts under both filters, as it must. These are the same integers the unitary ansatz
produces, because the amplitude space of coupled-cluster singles and doubles and the parameter space
of singlet unitary coupled cluster are the same space.

# Operational log

**Experiment ID:** 2026-08-11_e0_classical_ccsd_fullgroup
**Date:** 2026-08-11 (start) → 2026-08-11 (end)
**Status:** complete
**Mode:** numerical and symbolic (hybrid)
**Opportunity:** opportunities.md #1 and #2
**Target venue:** Quantum (componente do artigo principal); Comment no JCTC como item separado e posterior.

## Numerical results

Amplitude classes surviving each filter, and the correlation energy each filtered calculation
returns. Basis STO-3G, RHF reference in the adapted basis, `conv_tol = 1e-10`.

| molecule | full group | subgroup PySCF uses | no filter | Abelian | full group | ΔE Abelian (mHa) | ΔE full (mHa) |
|---|---|---|---|---|---|---|---|
| H2O | C2v | C2v | 65 = 10+55 | 26 = 4+22 | 26 = 4+22 | 3.1e-12 | 3.1e-12 |
| NH3 | C3v | Cs | 135 = 15+120 | 75 = 9+66 | 30 = 4+26 | 2.6e-12 | 2.5e-9 |
| CH4 | Td | D2 | 230 = 20+210 | 65 = 5+60 | 21 = 3+18 | 1.4e-13 | 1.5e-13 |

Each cell is written as singles plus doubles. Unfiltered correlation energies:
−0.057903425686 Ha for water, −0.086320628246 for ammonia, −0.079100918245 for methane.

Invariance of the converged unfiltered amplitudes under the full group, which is the vanishing rule
of Čársky et al. stated as a measurement rather than inferred from an energy:

| molecule | max abs residual, singles | max abs residual, doubles | max abs amplitude, singles | max abs amplitude, doubles |
|---|---|---|---|---|
| H2O | 7.3e-16 | 4.2e-17 | 1.7e-2 | 8.7e-2 |
| NH3 | 1.1e-12 | 7.7e-14 | 1.2e-2 | 9.2e-2 |
| CH4 | 9.7e-16 | 1.1e-16 | 7.1e-4 | 5.0e-2 |

Two independent routes to the Abelian counts agree exactly for all three molecules: character
summation over the representation matrices, and a combinatorial count over irrep labels using only
the direct-product rule. The agreement is asserted in code, so a mismatch would stop the run.

Representation validation: the overlap matrix is preserved to 3.0e-13 for ammonia and exactly for
the other two, and no element of the representation connects orbitals of different energy above
2.8e-11.

## Hypothesis check

**Supported.** The non-Abelian filter is free in classical coupled cluster, by both available
measures, and the Abelian control behaves as required.

## Figures

- `figures/2026-08-11_e0_classical_ccsd_fullgroup_amplitude_classes.pdf` — **For manuscript.**
  Surviving amplitude classes under three filters, with the correlation-energy change below.
  Element inventory in `figures/figure_manifest.json`.
- `figures/2026-08-11_e0_classical_ccsd_fullgroup_amplitude_classes.png` — raster preview, for
  inspection only.

## Pointers

- Aggregate data: `results/aggregate.csv`, `results/aggregate.json`
- Visualisation contract: `results/viz_schema.json`
- Full log: `experiment_log.json` (9 run entries, all `success`)
- Plan, including why the original design was not executable: `plan.md`
- Code: `code/run.py`, `code/analysis.py`
- Vendored sources and their hashes: `code/vendor/PROVENANCE.md`
- Environment: `env/system_info.json`, `env/requirements.txt`, `env/git_commit.txt`

## Notes

The experiment was specified on the assumption that PySCF would do the calculation directly, since
it does coupled cluster with symmetry. Measured on this machine with PySCF 2.12.1, `mol.topgroup`
returns `C3v` for ammonia while `mol.groupname` returns `Cs`, and `Td` becomes `D2` for methane;
`pyscf.symm.geom.symm_ops` returns the eight operations of the orthorhombic table whatever group was
detected. Switching symmetry on therefore exercises the Abelian filter and nothing beyond it, which
is the ceiling this project is about. The full-group machinery had to be written here.

That turned out to improve the experiment. Measuring the invariance residual of the converged
unfiltered amplitudes tests the vanishing rule directly, whereas an energy difference only tests its
consequence, and the residual is three to four orders of magnitude smaller than the smallest
quantity anyone would call converged.

One guard was needed and is worth recording. A projector built with a transposed index convention
is still a projector, so a null result could come from projecting onto the wrong subspace without
announcing itself. Two checks close that: the representation is validated against the overlap matrix
and the orbital energies, and the character count is required to reproduce the combinatorial count
exactly wherever both apply.

The ammonia full-group energy change, 2.5e-9 mHa, is three orders of magnitude larger than the other
five entries. It tracks the representation error for that molecule, whose group elements involve
irrational rotation matrices while the other two groups are built from signed permutations and are
exact in floating point. Nothing here depends on the distinction, but it explains the pattern.
