# Plan — 2026-08-11_e0_classical_ccsd_fullgroup

**Date:** 2026-08-11
**Scale:** canonical
**Mode:** hybrid — `numerical` governs the coupled-cluster energies and residuals; `symbolic`
governs the character counting and the independent combinatorial oracle.
**Opportunity:** `opportunities.md` #1 and #2
**Target venue:** Quantum

---

## Stage 1A — conceptual plan

### Question

In classical coupled cluster, does confining the amplitudes to the totally symmetric subspace of the
**full** point group — C3v for NH3, Td for CH4 — change the correlation energy?

### Hypothesis

**Supports:** the correlation energy is unchanged to the convergence threshold, and the converged
unconstrained amplitudes already satisfy the full-group invariance to numerical precision. Cársky
et al. proved exactly this in 1987 and ran NH3 in C3v with the symmetry inside the coupled-cluster
step; reproducing it with our own data makes the point without relying on the audit of anyone's code.

**Refutes:** a correlation-energy shift above the convergence threshold, or a converged amplitude
tensor with a significant component outside the symmetric subspace. Either would mean the vanishing
rule fails for degenerate irreps, and the claim that non-Abelian filtering is harmful would have a
classical foothold.

### Variables

- **Independent:** the molecule (H2O/C2v, NH3/C3v, CH4/Td) and the filter (none, Abelian subgroup,
  full point group).
- **Dependent:** the correlation energy, the number of surviving amplitude classes, and the
  invariance residual of the converged unconstrained amplitudes.
- **Controlled:** basis (STO-3G), reference (RHF in the adapted basis), convergence thresholds,
  and the amplitude equations themselves — only the projection differs.

### Method

The group acts on molecular orbitals by orthogonal matrices `D(g)` built from the Cartesian action
of each element. Occupied and virtual orbitals sit at different energies and never mix, so the
action factorises on `t1[i,a]` and `t2[i,j,a,b]`. The totally symmetric projector `(1/|G|) sum_g g`
— Eq. (6) of Häser, Almlöf and Feyereisen (1991), not a new construction — is applied to the
amplitudes at every coupled-cluster iteration by overriding `update_amps`.

Amplitude-class counts come from characters: the singles space is `V = occ (x) virt` and the doubles
space is `Sym^2 V`, so the invariant dimensions are `(1/|G|) sum_g chi_V(g)` and
`(1/|G|) sum_g [chi_V(g)^2 + chi_V(g^2)]/2`.

### Stopping criterion

Three molecules, three filters each, all converged to `conv_tol = 1e-10`.

### Acceptance preview

For Quantum this is the classical anchor of the argument, so it must be (i) reproducible without any
third-party code, (ii) counted by two independent routes, and (iii) accompanied by an Abelian control
where the additional compression is required to be exactly zero.

### Anticipated risks

- PySCF does not support non-Abelian point groups, so the full-group machinery has to be built and
  validated here rather than switched on.
- A wrong orbital-transformation convention would produce a projector that is still idempotent and
  would therefore not announce itself. It has to be caught by a post-condition, not by inspection.

---

## Stage 1B — execution plan

### Methodological critique of the plan-mini

**The experiment as originally specified cannot be run.** The plan said "PySCF does CCSD with
symmetry", so this would be cheap. It is not: PySCF detects the full point group and then works in
the largest Abelian subgroup it supports. Measured directly on this machine with PySCF 2.12.1,
`mol.topgroup` is `C3v` for NH3 while `mol.groupname` is `Cs`, and `Td` becomes `D2` for CH4;
`pyscf.symm.geom.symm_ops` returns only the eight D2h-type operations regardless of the group
detected. Switching `symmetry=True` on therefore tests the Abelian filter and nothing else, which is
precisely the ceiling the paper is about.

The fix, agreed before building: implement the full-group amplitude projector here and measure three
things instead of one — the class counts, the invariance residual of the converged amplitudes, and
the constrained energy. This is a better experiment than the original, because the residual measures
Cársky's vanishing rule directly rather than inferring it from an energy difference.

**The counting turns out to be the same problem as the quantum one.** The CCSD amplitude space is
`V (+) Sym^2 V` with the pair-exchange symmetry, which is exactly the singlet-UCCSD parameter space.
The classical and unitary counts are therefore directly comparable integers rather than analogous
quantities, which strengthens the comparison in the manuscript and costs nothing.

**Fair-comparison check.** All three filters use the same amplitude equations, the same reference and
the same thresholds; only the projection differs. The unconstrained run is the baseline for both.

**Alternative interpretation that must be ruled out.** A projector built with the wrong index
convention is still a projector, so a null result could come from projecting onto something other
than the intended subspace. Two guards: the representation is validated (`|U^T S U - S|` and the
largest off-shell element of `D`), and the character count is checked against an independent
combinatorial count over irrep labels for every Abelian subgroup, where the two routes must agree
exactly.

### Code reuse

| what | where | used for |
|---|---|---|
| `group(name)` | `vendor/group_compression.py:76` | C3v as rotations and reflections; Td as the 24 signed permutations with positive sign product; the Abelian subgroup |
| `mo_rep(mol, mo, G)` | `vendor/group_compression.py:108` | `D(g)` in the MO basis, with validation |
| `geo_h2o`, `geo_nh3`, `geo_ch4`, `to_str` | `vendor/group_compression.py:53,37,46,66` | geometries |

### New functions

- `AmplitudeProjector` — the totally symmetric projector on `t1` and `t2`.
- `counts_from_rep(D_sub, nocc, nmo)` — invariant dimensions of `V` and `Sym^2 V` by characters.
  Using characters rather than building the `|V|^2 x |V|^2` projector keeps the cost negligible.
- `abelian_oracle(orbsym, nocc, nmo)` — the independent combinatorial count.
- `SymCCSD` — `pyscf.cc.ccsd.CCSD` with the projector applied in `update_amps` and `init_amps`.

### Verification order

1. Representation validity: `|U^T S U - S|` at machine precision and no element of `D` connecting
   orbitals of different energy.
2. Counting agreement: the character count equals the combinatorial count for every Abelian
   subgroup. Asserted, so a mismatch stops the run.
3. Abelian control: for H2O in C2v and any group already Abelian, the full-group count must equal
   the Abelian count exactly.
4. Invariance residual of the converged unconstrained amplitudes.
5. Only then is the constrained-versus-unconstrained energy difference interpretable.

### File inventory

**Create:** `code/run.py`, `code/smoke_test.py`, `code/expcommon.py`, `experiment_log.json`,
`results/aggregate.csv`, `results/aggregate.json`, `results/viz_schema.json`, `figures/`, `README.md`,
`env/`.
**Reuse read-only:** three vendored modules.
**Update outside this folder:** none.

### Time budget

Under two seconds per molecule. The whole experiment runs in about one second of coupled-cluster
time; the cost is entirely in building and validating the representation.
