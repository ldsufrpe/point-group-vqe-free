# Plan — 2026-08-11_e4_energy_invariance_t2

**Date:** 2026-08-11
**Scale:** canonical
**Mode:** hybrid — `numerical` governs the energy and gradient measurements; `symbolic` governs
the character algebra, the (A1)–(A3) hypothesis checks and the exactness claims (an exact
identity is asserted, not a converged approximation).
**Opportunity:** `opportunities.md` #3
**Target venue:** Quantum

---

## Stage 1A — conceptual plan

### Question

Does the energy of the tied product ansatz actually satisfy `E(sigma_n theta) = E(theta)` for every
element `n` of the Abelian symmetry group, at parameter vectors that are not already symmetric?

### Hypothesis

**Supports:** the difference sits at machine precision, of order `1e-12` absolute, with no
dependence on the sampled parameter vector, and it survives at parameter magnitudes comparable to
converged amplitudes.

**Refutes:** any systematic deviation that grows with the magnitude of `theta`, or that depends on
which group element is applied. That would mean the invariance holds only near `theta = 0`, and the
critical-submanifold corollary and the `theta = 0` protocol theorem would both lose their support.

### Variables

- **Independent:** the parameter vector `theta` (20 draws per molecule), the group element `n`
  (1 nontrivial element for NH3 in Cs, 3 for CH4 in D2), the molecule, and the orbital basis
  (symmetry-adapted, or rotated with stale labels).
- **Dependent:** `|E(theta) - E(sigma_n theta)|`, its relative version, the state-level residual
  `|| U(sigma_n theta)|phi0> - S_n U(theta)|phi0> ||`, and the largest gradient component in the
  removed directions.
- **Controlled:** the Hamiltonian, the reference determinant, the gate ordering, the tying map, and
  the code path — both energies are evaluated by the same `TrotterAnsatz` instance, so nothing but
  the sign pattern differs between the two evaluations.

### Method

The ansatz is `U(theta) = prod_k exp(theta_{p(k)} c_k kappa_k)` with `kappa_k = tau_k - tau_k^dagger`
built on the determinant basis. The Hamiltonian is the spin-free form assembled from `E_pq`
generators. The group is the Abelian subgroup of real characters that PySCF exposes — Cs for NH3,
D2 for CH4 — whose irreps are labelled by bit patterns, so the character of irrep `j` at the element
indexed by `i` is `(-1)^popcount(i & j)`. The sign flip `sigma_n` is applied to the parameter vector
through those characters.

### Stopping criterion

20 parameter vectors per molecule per group element, plus 20 gradient evaluations on the
submanifold, plus 5 negative-control vectors per molecule. Both molecules complete.

### Acceptance preview

For Quantum, the theorem must be shown to hold at machine precision, and the measurement must be
shown to be capable of failing. A test that returns zero for a trivial reason is not evidence, so
the negative control is part of the acceptance bar, not an extra.

### Anticipated risks

- The test is vacuous if `theta` is drawn inside the symmetric submanifold, because `sigma_n theta`
  then equals `theta` identically.
- Taking the Hamiltonian from `openfermionpyscf` would silently reintroduce the very label mismatch
  under audit, since that package fixes `symmetry = False`.

---

## Stage 1B — execution plan

### Methodological critique of the plan-mini

**One real problem, fixed before building.** The plan as originally stated — "draw 20 theta, apply
sigma_n, check the difference" — is vacuous as written. The existing evidence for the neighbouring
corollary was collected *on* the symmetric submanifold, where `theta_A = 0`; there `sigma_n theta`
and `theta` are the same vector and the measured difference is exactly zero for a reason that has
nothing to do with the theorem. The fix is to force `theta_A` nonzero and to record how far off the
submanifold each sample sits, which the code asserts and the aggregate reports
(`theta_offmanifold_max`).

**Second, the magnitude matters.** A draw at `theta ~ 1e-6` is nearly vacuous too, because every
term is then linear and the invariance is trivial. Parameters are drawn uniformly on `[-0.5, 0.5]`,
which is well above converged amplitude magnitudes (`max |t1| ~ 1e-2`, `max |t2| ~ 9e-2`), so the
test probes the nonlinear regime.

**Third, the original scope was too narrow to be convincing.** Two additions cost almost nothing and
change what the experiment establishes. The state-level identity is strictly stronger than the
energy identity and closes the equivariance lemma, which also had no direct check. And the negative
control makes the result discriminating: in a basis where a degenerate pair has been rotated while
the labels are kept, the hypothesis of the theorem fails, so the measured quantity must blow up.
Without that, a reader cannot distinguish the theorem from an arithmetic tautology.

### Code reuse

Everything load-bearing already exists and is vendored verbatim into `code/vendor/`.

| what | where | used for |
|---|---|---|
| `build_det_basis(norb, nalpha, nbeta)` | `vendor/indep.py:24` | determinant basis |
| `SpinFreeHam.__init__/matvec` | `vendor/indep.py:82,99` | spin-free Hamiltonian via `E_pq` |
| `kappa_matrix(elem, dets, index)` | `vendor/indep.py:180` | anti-Hermitian generators |
| `TrotterAnsatz` | `vendor/indep.py:189` | product ansatz, tied parameters, adjoint gradient |
| `gates_of(norb, nelec)` | `vendor/audit.py:53` | the canonical singlet-UCCSD pool |
| `irrep_of(g, orbsym)` | `vendor/audit.py:119` | irrep label of an excitation |
| `prep(name, geom, conv, symmetry)` | `vendor/audit.py:126` | one-stop builder; also asserts (A3) |
| `nh3_geom`, `ch4_geom` | `vendor/audit.py:32,46` | geometries |

The only patch applied to the vendored sources is the removal of a hard-coded absolute path; see
`code/vendor/PROVENANCE.md`.

### New functions

- `characters(par_irrep, elem)` — the character vector of the parameter classes at one group element.
- `det_sign_operator(dets, orbsym, elem)` — the diagonal representation of the element on the
  determinant basis, used for the state-level test and for checking (A2).
- `rotate_degenerate_pair(M, phi_deg)` — the negative control.

### Verification order

1. The representation is sound: `(A2)` holds exactly, i.e. the reference is a fixed point of every
   element's determinant-space sign operator. Asserted, not merely reported.
2. `(A3)` holds: no parameter mixes irreps. Asserted inside `prep`.
3. The sampled parameter vectors are genuinely off the submanifold. Asserted.
4. Only then are the energy and state differences meaningful.
5. The negative control breaks the invariance by orders of magnitude; if it does not, criteria 1–4
   were satisfied for the wrong reason.

### File inventory

**Create:** `code/run.py`, `code/smoke_test.py`, `code/expcommon.py`, `experiment_log.json`,
`results/aggregate.csv`, `results/aggregate.json`, `results/viz_schema.json`, `figures/`, `README.md`,
`env/`.
**Reuse read-only:** the four vendored modules.
**Update outside this folder:** none.

### Time budget

Setup dominates: 2.6 s for NH3, 16 s for CH4. The measured runs total about 20 s. Well inside the
budget; the experiment is cheap enough that re-running it is never a reason to import.
