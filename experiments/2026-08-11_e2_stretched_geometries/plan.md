# Plan — 2026-08-11_e2_stretched_geometries

**Date:** 2026-08-11
**Scale:** canonical
**Mode:** hybrid — `numerical` governs the geometry scan and the variational energies; `symbolic`
governs the per-geometry irrep labelling and the invariant subspace.
**Opportunity:** `opportunities.md` #3
**Target venue:** Quantum

---

## Stage 1A — conceptual plan

### Question

Does the cost of point-group filtering grow when the molecule is pulled apart and static correlation
takes over?

### Hypothesis

**Supports the filter being free:** the difference between filtered and unfiltered energies stays at
the level of numerical noise across the whole scan, even as the unfiltered error itself grows by two
orders of magnitude.

**Refutes it:** a difference that grows with the N–H distance, tracking the collapse of the S0–S1
gap. That would vindicate the withdrawn manuscript's claim that the deficit "becomes more pronounced
at stretched geometries where full orbital relaxation is required".

### Variables

- **Independent:** the N–H distance, from 0.90 to 2.20 A at fixed H–N–H angle; the ansatz.
- **Dependent:** the energy error against FCI, the difference against unfiltered UCCSD, the
  correlation energy, and the S0–S1 gap.
- **Controlled:** basis, angle, determinant basis, generator set, optimiser, thresholds. The
  determinant basis and the generators do not depend on geometry and are built once for the scan.

### Method

RHF at each geometry, spin-free Hamiltonian, two lowest eigenvalues by sparse diagonalisation so the
gap can be reported alongside the energies. Irrep labels are recomputed at every geometry because
orbital ordering can cross. Three ansaetze are compared under an identical Hamiltonian and optimiser.

### Stopping criterion

The nine distances of the scan. A geometry where RHF does not converge is reported as a failure, not
retried with a different guess: the hypothesis of the theorems is a symmetric closed-shell reference,
so its absence is the answer, not an obstacle.

### Acceptance preview

For Quantum the scan has to be severe enough to matter — the correlation energy and the gap have to
move by a large factor — and the boundary of the claim has to be stated rather than avoided.

### Anticipated risks

- Orbital crossing could change the irrep labels mid-scan and silently change the filtered pool size.
- RHF will stop converging somewhere in the stretched region.

---

## Stage 1B — execution plan

### Methodological critique of the plan-mini

**The plan-mini tests the wrong filter.** As specified — and as the earlier scan was actually run —
the comparison is unfiltered UCCSD against the *Abelian* SymUCCSD. That answers the question the
withdrawn manuscript raised, but it leaves the filter this paper proposes, the full non-Abelian one,
untested under static correlation. The full C3v invariant ansatz is added at every geometry. This is
the single most consequential change to the experiment, and it is the one that produced the result.

**A confound has to be separated.** If the full-group ansatz does start to cost something, two very
different things could be responsible: the 30-dimensional invariant subspace may genuinely fail to
represent the stretched wavefunction, or the price may be tying the parameters and Trotterising. The
tied Trotter form alone cannot distinguish them. The same invariant operators are therefore also run
with exact exponentials, and the two are reported side by side. Without this control the experiment
would produce a number that cannot be interpreted.

**Fair comparison.** One Hamiltonian, one reference, one optimiser and one set of thresholds per
geometry; the ansatz is the only thing that changes.

**Honest reporting.** The scan is expected to end at an RHF convergence failure. That point is logged
with `status: failure` and reported, because it is where the hypothesis of every theorem in the paper
stops holding, and saying so converts a caveat into a delimited result.

### Code reuse

| what | where | used for |
|---|---|---|
| `nh3_geom(r, ang_deg)` | `vendor/audit.py:32` | the geometry family |
| `gates_of`, `irrep_of` | `vendor/audit.py:53,119` | pool and per-geometry filter |
| `sparse_invariant_basis` | `vendor/e1_trotter.py:34` | the pivoted invariant basis |
| `group`, `mo_rep`, `build_ops` | `vendor/group_compression.py:76,108,152` | representation and invariant operators |
| `build_det_basis`, `SpinFreeHam`, `kappa_matrix`, `TrotterAnsatz` | `vendor/indep.py` | numerical core |

### New functions

- `optimise(ans, npar, active, egs)` — masked L-BFGS-B driver returning energy, error, iteration
  count and gradient norm.
- `eg_exact(th)` — energy and adjoint gradient under exact exponentials of the invariant operators,
  the control that separates subspace expressivity from Trotterisation.

### Verification order

1. RHF converged at this geometry. If not, the point is a logged failure and nothing else is
   computed there.
2. The filtered pool size is recomputed from the current labels, so an orbital crossing shows up as a
   change in the count rather than silently corrupting the filter.
3. The invariant operators decompose on the elementary basis to machine precision.
4. Only then are the four energies comparable.
5. The exact-exponential control is what licenses any statement about *why* a cost appears.

### File inventory

**Create:** `code/run.py`, `code/expcommon.py`, `experiment_log.json`, `results/aggregate.csv`,
`results/aggregate.json`, `results/viz_schema.json`, `figures/`, `README.md`, `env/`.
**Reuse read-only:** four vendored modules.
**Update outside this folder:** none.

### Time budget

Roughly 20 to 50 s per geometry for the three variational ansaetze, plus the exact-exponential
control. Under an hour for the scan.
