# Reproduction Report

What was actually executed before this repository was published. Written so that the README
never promises a command that was not run.

| | |
|---|---|
| Date | 2026-08-13 |
| Machine | Linux 7.0.0-28-generic, x86_64, AMD Ryzen 7 3700U, 8 logical cores |
| Environment | fresh `python3 -m venv`, `pip install -r requirements.txt`, nothing else |
| Python | 3.12.3 |

## 1. Environment installs from `requirements.txt` alone — PASS

A clean virtual environment was created outside this repository and populated only from
`requirements.txt`. Installation exited 0. The resolved versions match the interpreter that
produced the evidence, exactly:

```
pyscf 2.12.1   openfermion 1.7.1   numpy 2.4.3   scipy 1.17.1
```

## 2. Smoke tests — 6 of 6 PASS

Each experiment's `code/smoke_test.py` runs a reduced configuration and asserts the invariants
that experiment exists to establish. All six were run in the clean environment above, from a
scratch copy of this tree so that no output could contaminate the published artifacts.

| Experiment | Result | Representative assertion |
|---|---|---|
| `2026-08-11_e0_classical_ccsd_fullgroup` | OK | `oracle_agreement: True`; invariance residual t2 = 4.2e-17 |
| `2026-08-11_e1_invariant_trotter_ansatz` | OK | decomposition residual 3.3e-16; H₂O/STO-3G params ×0.400 |
| `2026-08-11_e2_stretched_geometries` | OK | gap collapse 628.1 → 480.4 mHa over the reduced scan |
| `2026-08-11_e4_energy_invariance_t2` | OK | max gradient along removed directions 2.2e-10 |
| `2026-08-11_e7_degenerate_shell_artifact` | OK | NH₃ spread 2.0e-09 mHa; no undersampled arm |
| `2026-08-12_e0b_basis_dependence` | OK | `oracle_agreement: True` |

## 3. Artifact integrity — 17 of 17 identical

Every exported artifact was compared by md5 against the authors' working copy, and the five
figures were additionally compared against the PDFs compiled into the manuscript.

- 6 / 6 `results/aggregate.csv` identical
- 6 / 6 `experiment_log.json` identical (287 runs in total)
- 5 / 5 figures identical to the manuscript's

## What was NOT run

**The full canonical runs were not re-executed for this release.** They take hours, and the
published `aggregate.csv` files are the originals from the canonical execution of 2026-08-11 and
2026-08-12, byte-verified above rather than regenerated. Anyone wishing to regenerate them from
scratch should use the `python code/run.py all` commands in the README's table.

The E7 thread-sensitivity arms were exercised only at the smoke scale. Reproducing the full
spread requires the BLAS thread sweep described in that experiment's `plan.md`.
