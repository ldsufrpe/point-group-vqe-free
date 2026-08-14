# Provenance of the vendored sources

These files are copied verbatim from the independent audit that produced the evidence this
project is built on, so that this experiment folder reproduces without reaching outside
itself.

**Source directory:** `/home/leon-denis/PycharmProjects/Experimentos/auditoria_independente/`

**Copied:** 2026-08-11. The source directory is not a git repository, so the sha256 of each
file at copy time is the identifier.

| file | sha256 of the source | what it provides |
|---|---|---|
| `audit.py` | `65b6912feefc5fe96be55d112a1f580581025f21da1efa69c192805977c31d52` | canonical singlet-UCCSD pool (gates_of), irrep filter (irrep_of), one-stop builder (prep), L-BFGS-B driver (optimize), geometries |
| `e1_trotter.py` | `0e723d8cda59a54c3cdf4a0ff275419a858baa38ab0ca889d4df0b2199e7a484` | pivoted {P e_j} sparse invariant basis, and the tied Trotter circuit built from it |
| `group_compression.py` | `a19da2a8f3d17b247cdd2c65d9c068e05ca52efa4a9c5f2c386a47068a8fb8c5` | point-group representation for C3v/Td/C2v/D2h, MO-basis representation with validation, totally symmetric projector, invariant operator construction |
| `indep.py` | `90ad483c2bf769662c93ad1c439362d71eabb12566feecc8b0a5fa777d69bca5` | core library: determinant basis, spin-free Hamiltonian via E_pq, Trotterised product ansatz with tied parameters, adjoint gradient |

## The one modification

Each file declared a module-level `HERE` holding an absolute path, used both to locate its
sibling modules and to write its output. Two of those paths pointed at a scratch directory
that no longer holds the modules, so outputs were being written somewhere the author did not
intend. The line was replaced by

```python
HERE = os.path.dirname(os.path.abspath(__file__))
```

and nothing else was touched. Diff the files against the sha256 above to confirm.
