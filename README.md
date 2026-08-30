# Point-group filtering in unitary coupled-cluster ansätze

Reproducibility repository for the paper *"Point-group filtering in unitary coupled-cluster
ansätze: exact criticality under the Abelian filter, freeness measured near equilibrium under
the full group"* — version 3 of [arXiv:2603.21009](https://arxiv.org/abs/2603.21009), under
review at the *Journal of Mathematical Chemistry*.

## What this is

A curated export: the code, aggregated results, environment snapshots, run logs and figures
behind the paper's seven canonical experiments (**311 runs**). It is not the authors' working
directory — drafts, review notes and third-party PDFs are deliberately absent.

The paper itself is not mirrored here; the public text is the arXiv version, which is what this
repository's numbers should be read against.

### Relation to `nh3-symuccsd-confinement`

An earlier repository, [`ldsufrpe/nh3-symuccsd-confinement`](https://github.com/ldsufrpe/nh3-symuccsd-confinement),
accompanies versions 1 and 2 of the same arXiv entry. **Its central result is superseded.**
Those versions reported a 21.8 mHa energy deviation for NH₃ in STO-3G and read it as a
variational cost of the symmetry filter; it was an artifact of an incomplete excitation pool in
that code, and the pool arithmetic is diagnosed in
`experiments/2026-08-11_e7_degenerate_shell_artifact/`. That repository is archived and kept
online because the superseded arXiv versions link to it. This repository is the current one.

## Repository structure

```
experiments/
  2026-08-11_e0_classical_ccsd_fullgroup/   classical CCSD under the full point group (STO-3G)
  2026-08-11_e1_invariant_trotter_ansatz/   unitary ansatz: parameter and operator counts
  2026-08-11_e2_stretched_geometries/       filter cost along a bond-stretching scan
  2026-08-11_e4_energy_invariance_t2/       energy invariance and gradients along removed directions
  2026-08-11_e7_degenerate_shell_artifact/  incomplete-pool artifact and its predictor
  2026-08-12_e0b_basis_dependence/          the classical result in 6-31G, six molecules
  2026-08-14_e8_qeb_cnot_cost/              operator and CNOT counts in a single counting currency
    each containing:
      plan.md               the pre-registered plan
      README.md             protocol, hypotheses, what the run showed
      code/                 run.py, analysis.py, smoke_test.py, expcommon.py, vendor/
      results/              aggregate.csv, aggregate.json, run logs
      figures/              PDF/PNG plus figure_manifest.json
      env/                  requirements.txt, system_info.json, git_commit.txt
      experiment_log.json   append-only per-run log
      report_validation.md  post-run validation
requirements.txt            project-level pinned environment
```

## Requirements

Python 3.12 and the pinned versions in `requirements.txt` — PySCF 2.12.1, OpenFermion 1.7.1,
openfermionpyscf 0.5, NumPy 2.4.3, SciPy 1.17.1, matplotlib 3.10.8, SciencePlots 2.2.2.

```sh
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
```

Tested on a clean virtual environment on 2026-08-13; see `reproduction_report.md`.
`experiments/<id>/env/requirements.txt` holds the full transitive freeze and is authoritative
for reproducing one specific experiment.

## Reproducing

Every experiment ships a smoke test that runs a reduced configuration in under a minute and
asserts the experiment's invariants. **These are the commands that were executed for this
release**, all six passing:

```sh
cd experiments/2026-08-11_e0_classical_ccsd_fullgroup && python code/smoke_test.py
```

…and likewise for the other six directories.

The full runs are the ones that produce the paper's numbers. They are **not** executed by the
smoke tests and take hours on a laptop; the commands are:

| Paper item | Command (run from the experiment directory) | Produces |
|---|---|---|
| Fig. 1, Table 2 (STO-3G half) | `python code/run.py all` in `…_e0_classical_ccsd_fullgroup` | `results/aggregate.csv`, `figures/*_amplitude_classes.pdf` |
| Table 2 (6-31G half) | `python code/run.py all` in `…_e0b_basis_dependence` | `results/aggregate.csv` |
| Fig. 2, Table 3 | `python code/run.py all` in `…_e1_invariant_trotter_ansatz` | `results/aggregate.csv`, `figures/*_compression_params_gates.pdf` |
| Fig. 3 | `python code/run.py all` in `…_e4_energy_invariance_t2` | `results/aggregate.csv`, `figures/*_invariance_violation.pdf` |
| Fig. 4, Table 4 | `python code/run.py all` in `…_e2_stretched_geometries` | `results/aggregate.csv`, `figures/*_filter_cost_vs_stretch.pdf` |
| Fig. 5, Table 5 | `python code/run.py all` in `…_e7_degenerate_shell_artifact` | `results/aggregate.csv`, `figures/*_artifact_predictor_and_spread.pdf` |
| Operator counts behind Fig. 2 and Table 3 | `python code/run.py all` in `…_e8_qeb_cnot_cost` | `results/aggregate.csv`, `figures/*_cnot_cost.pdf` |

`run.py` accepts a molecule key (`h2o`, `nh3`, `ch4`, `all`) and `--smoke`. The five figures
shipped in `experiments/*/figures/` that the paper prints are byte-identical to the printed ones
(md5 verified); E8's figure is evidence for the counts, not a figure of the paper.
Table numbering follows the arXiv version; check the captions if it has been revised.

## What the numbers say, and against which baseline

Compression figures are meaningless without their baseline, their ansatz variant and their
domain, so each is stated with all three. All of these read directly off the `aggregate.csv`
files.

- **Unitary ansatz, parameters (E1, STO-3G only).** Against *unfiltered* UCCSD, NH₃ falls 135 →
  30 and CH₄ 230 → 21. Against the *Abelian* filter — which is what current practice uses, and
  the honest comparison — NH₃ falls 75 → 30 and CH₄ 65 → 21.
- **Unitary ansatz, operator counts (E1, re-counted in E8).** Against the Abelian baseline, NH₃
  goes 163 → 117 and CH₄ 146 → 128, in *distinct elementary excitation operators*. A pool admits
  two counts — the distinct operators it touches, and their incidences in the Trotterised product —
  and E1's own `n_gates` column mixes them: its `uccsd` and pool-subset rows are distinct counts,
  its projector rows incidences. E8 re-counts every pool in both currencies, and its
  `results/aggregate.csv` carries a `currency` column so the two are never crossed again. Read the
  `currency=distinct` rows for the numbers above.
- **Unitary ansatz, CNOT counts (E8).** Priced under the qubit-excitation scheme of the audited
  work — two CNOTs per spin-orbital single, thirteen per double — the full group takes NH₃ from
  1921 → 1411 and CH₄ from 1788 → 1554 against the Abelian filter. The ruler is fixed by
  reproducing all six published bars of that work's Figure 5 exactly, with no fitted quantity. The
  saving tracks the flat operator count rather than beating it, so the filter shows no useful
  preference for the expensive doubles. **Transpiled depth was never measured** — no device gate
  set, no connectivity, no scheduling — so this bounds circuit *size* in that scheme and supports
  no device-level resource claim.
- **Classical CCSD (E0, E0b).** The full-group filter changes the correlation energy by at most
  3.4 × 10⁻¹² mHa in 6-31G for the molecules whose point group is finite (H₂O, NH₃, CH₄), against
  2.5 × 10⁻⁹ mHa in STO-3G. The Abelian filter costs at most 2.0 × 10⁻¹⁰ mHa in 6-31G across all
  six molecules. Amplitude classes compress 1325 → 277 for NH₃ and 1890 → 144 for CH₄. Three of
  the six molecules are linear, so their full point group is continuous and only the Abelian
  filter is tested there.
- **Geometric domain (E2).** The Abelian filter is unrestricted over the whole 0.90–2.00 Å scan.
  The full-group filter is free only near equilibrium, crossing chemical accuracy between 1.60
  and 1.80 Å. Degradation figures differ by roughly a factor of two between the exact-exponential
  and tied-Trotter variants; `aggregate.csv` names the variant per row.
- **Basis.** The classical result is measured in STO-3G *and* 6-31G. The unitary result (E1) is
  STO-3G only, and no correlation-consistent basis was tested.

## Data & experiments

`results/aggregate.csv` is the aggregated, analysis-ready table per experiment; `aggregate.json`
carries the same content with the summary block. `experiment_log.json` is the append-only
per-run record, and the seven of them sum to 311 runs. `env/system_info.json` records the machine,
BLAS implementation and thread count — E7 depends on the last of these, since the artifact it
diagnoses is thread-sensitive. `env/git_commit.txt` refers to the authors' private working
repository and is kept as a provenance token only.

## Citation

See `CITATION.cff`. Cite the paper; this repository is its evidence.

## License

MIT for code, CC BY 4.0 for data, results and figures. See `LICENSE` and `LICENSE-DATA.md`.
