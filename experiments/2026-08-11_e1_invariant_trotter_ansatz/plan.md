# Plan — 2026-08-11_e1_invariant_trotter_ansatz

**Date:** 2026-08-11
**Scale:** canonical
**Mode:** hybrid — `numerical` governs the variational energies; `symbolic` governs the invariant
subspace construction and the exactness of the operator decomposition.
**Opportunity:** `opportunities.md` #2
**Target venue:** Quantum

---

## Stage 1A — conceptual plan

### Question

The compression from 135 to 30 parameters for NH3 and from 230 to 21 for CH4 was measured with exact
exponentials of the invariant operators. Does the accuracy survive the product form a device would
actually execute — Trotterised, with one parameter tied across every elementary gate of an invariant
operator — and what happens to the gate count?

### Hypothesis

**Supports:** the energy error stays at the same order as unfiltered UCCSD, with the Trotterisation
penalty far below chemical accuracy, and the elementary gate count falls together with the parameter
count.

**Refutes:** a Trotterisation penalty comparable to the error itself, or a gate count that rises when
the parameter count falls. Either would mean the compression is a bookkeeping statement about
parameters with no consequence for a circuit, and no hardware claim could be made.

### Variables

- **Independent:** the molecule, the group used to build the invariant subspace (Abelian subgroup or
  full point group), the exponentiation scheme (exact, Trotter in natural gate order, Trotter in
  reversed order).
- **Dependent:** parameter count, elementary gate count, the residual of the decomposition of each
  invariant operator on the elementary basis, and the energy error against FCI.
- **Controlled:** basis, geometry, reference, optimiser, convergence thresholds, and the elementary
  generator set — all three schemes act on the same generators.

### Method

Invariant operators are built from the totally symmetric projector on the singles space `V` and the
doubles space `Sym^2 V`, then each is decomposed on the basis of distinct elementary excitations by
a Frobenius inner product, which is exact because distinct generators connect disjoint pairs of
determinants. The two Trotter orders bound the ordering dependence, which is the only free choice in
the product form.

H2O in C2v and C2H4 in D2h are Abelian controls: there the full group *is* the Abelian subgroup, so
the additional compression must be exactly zero. C2H4 is counted but not optimised, because its
determinant space is around nine million.

### Stopping criterion

Four molecules; two of them with the full energy comparison, two as controls; both Trotter orders.

### Acceptance preview

For Quantum the claim has to be about a circuit, so parameter counts alone are not enough: the gate
count and the Trotterisation penalty are both required, and the decomposition has to be shown exact
rather than approximate.

### Anticipated risks

- The invariant subspace has many orthonormal bases and the natural one from an eigensolver is
  generically dense, which would inflate the circuit rather than shrink it.
- Real circuit *depth* depends on a two-qubit gate decomposition that this experiment does not
  perform, so depth must not be claimed.

---

## Stage 1B — execution plan

### Methodological critique of the plan-mini

**The central trap is the choice of basis, and it is not a detail.** `numpy.linalg.eigh` returns an
arbitrary orthonormal basis of the eigenvalue-one subspace, generically dense: every invariant
operator becomes a combination of hundreds of elementary excitations, and the first attempt at this
measurement produced 18709 gates for NH3 against 315 for unfiltered UCCSD — the conclusion inverted.
The correct basis is the pivoted `{P e_j}`: the projector applied to coordinate vectors, each with at
most `|G|` terms, with independent columns selected by QR with column pivoting. This is the
difference between a method that is implementable and one that is not, so it is stated in the plan
rather than left to the code.

**The plan-mini omitted the baseline.** Comparing the two filtered ansaetze to each other says nothing
about what the filter costs; the unfiltered UCCSD circuit has to be measured under identical
conditions. It is added here as a separate run with its own parameter and gate counts.

**Fair comparison.** All ansaetze use the same generators, the same reference, the same optimiser and
the same thresholds; only the operator set and the tying differ.

**Alternative interpretation to rule out.** A small energy error could come from a lucky optimiser
path rather than from the ansatz. The decomposition residual guards the other half: if each invariant
operator is reproduced on the elementary basis to machine precision, the Trotter circuit is an exact
rewriting of the operator, and any energy difference is due to operator ordering alone.

### Code reuse

| what | where | used for |
|---|---|---|
| `sparse_invariant_basis(Ms, nS)` | `vendor/e1_trotter.py:34` | the pivoted `{P e_j}` basis — the load-bearing piece |
| `run(tag, atoms, gname)` | `vendor/e1_trotter.py:71` | the validated per-molecule measurement, called unchanged |
| `build_ops(Bs, Bd, Ep, nS)` | `vendor/group_compression.py:152` | anti-Hermitian invariant operators |
| `group`, `mo_rep`, `invariant_basis`, `CASES`, `to_str` | `vendor/group_compression.py:76,108,140,253,66` | representation, projector, control counts |
| `gates_of` | `vendor/audit.py:53` | the unfiltered pool for the baseline |
| `build_det_basis`, `SpinFreeHam`, `kappa_matrix`, `TrotterAnsatz` | `vendor/indep.py` | the numerical core |

The validated measurement is called rather than reimplemented, so the physics is bit-for-bit the
audited one; this experiment adds the baseline, the controls, the logging and the aggregation around
it.

### New functions

- `uccsd_reference(tag, atoms)` — the unfiltered UCCSD baseline: parameters, distinct elementary
  gates, and energy error.
- `counts_only(tag, atoms, gname)` — parameter counts for a molecule whose determinant space is out
  of reach.

### Verification order

1. Representation validity, inside `mo_rep`.
2. Abelian controls: for H2O and C2H4 the full-group count must equal the Abelian count exactly. If
   it does not, the representation or the projector is wrong and nothing downstream matters.
3. Decomposition residual per invariant operator at machine precision.
4. Only then are the three energies comparable, and only then does the gate count mean anything.

### File inventory

**Create:** `code/run.py`, `code/expcommon.py`, `experiment_log.json`, `results/aggregate.csv`,
`results/aggregate.json`, `results/viz_schema.json`, `figures/`, `README.md`, `env/`.
**Reuse read-only:** four vendored modules.
**Update outside this folder:** none.

### Time budget

About 18 s for H2O, seconds for the C2H4 counts, and roughly three minutes each for NH3 and CH4
including the baseline. Under ten minutes in total.

### Known limitation carried forward

Gate count here means elementary excitation operators, not two-qubit gates. The mapping to a device
gate set is not performed, so the result bounds circuit *size* in this generator basis and says
nothing about depth.
