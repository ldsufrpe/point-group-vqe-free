# Provenance of the vendored modules

Origin: `/home/leon-denis/PycharmProjects/Experimentos/auditoria_independente/`, the independent
audit code described in `EXPERIMENTS_HANDOFF.md` §1. Vendored here so that this experiment folder
reproduces on its own, without a path into a directory outside the project tree.

`audit.py`, `e1_trotter.py`, `group_compression.py` and `indep.py` were copied from
`experiments/2026-08-11_e1_invariant_trotter_ansatz/code/vendor/`, which is where they were first
vendored and where their patches were first applied. Copying from E1 rather than from the origin
guarantees this experiment runs the same bytes E1 ran — the parameter and operator counts of E1's
`results/aggregate.csv` are a regression check for this folder, and that check is only meaningful
against identical sources.

`hi_count.py` and `he_pipeline.py` were copied from the origin, having no E1 counterpart.

| module | what it provides here |
|---|---|
| `indep.py` | determinant basis, the `E_pq` operators, the elementary generators `kappa_matrix` |
| `audit.py` | `gates_of` — the UCCSD pool resolved into elementary spin-orbital excitations; `irrep_of` — the Abelian filter by irrep label |
| `group_compression.py` | point-group representation on the MO basis, the symmetric projector, `build_ops`, the geometries |
| `e1_trotter.py` | `sparse_invariant_basis` — the pivoted `{P e_j}` basis |
| `hi_count.py` | `ham_term_keys` — the HiUCCSD pool, from the compressed Hamiltonian's term set |
| `he_pipeline.py` | He *et al.*'s geometries, imported by `hi_count.py` |

## Patches applied, and nothing else

Three, all mechanical. No numerical or algorithmic line was touched.

1. `audit.py`, `e1_trotter.py`, `group_compression.py` — `HERE` was a hard-coded absolute path;
   now `os.path.dirname(os.path.abspath(__file__))`. Applied by E1, inherited here.
2. `hi_count.py:18` and `he_pipeline.py:28` — the same `HERE` patch, applied here for the first
   time. `hi_count.py`'s original value pointed into an agent scratchpad that no longer exists, so
   the module could not be imported at all before this.
3. `hi_count.py` end of file — the three driver statements (`count(...)`, `count(...)`, `print`)
   ran at **import** time, so `from hi_count import ham_term_keys` executed the whole script. Now
   guarded by `if __name__ == '__main__':`. The statements themselves are unchanged.

## Third-party code

None. `indep.py` is a from-scratch reimplementation; its header records that no code was inherited
from He *et al.* (MindQuantum) or from any prior local script. The audited third-party sources live
in `/home/leon-denis/PycharmProjects/Experimentos/He etal code/` and are read, never imported.
