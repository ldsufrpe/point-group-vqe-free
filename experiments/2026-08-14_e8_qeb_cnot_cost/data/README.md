# data/

Empty by design, and that is a statement about the experiment rather than an omission.

E8 consumes no external data. Every input is generated in `code/run.py` from two
closed-form geometries and a basis set name:

| input | where it comes from |
|---|---|
| NH3 and CH4 geometries | `vendor/group_compression.py`, `geo_nh3(d=1.09, ang=107.3)` and `geo_ch4(d=1.09)` — He et al.'s convention, deliberately, so that their published bars are reproducible |
| basis set | STO-3G, from PySCF's internal library |
| the QEB ruler | 2 CNOTs per spin-orbital single, 13 per double, from He et al., JCTC 2026, 22, 6008-6019 |
| the published Figure 5 values | transcribed into the `PUBLISHED` table of `code/run.py`, six integers |

The geometry strings are hashed into every `experiment_log.json` entry under
`inputs.data_hashes.geometry`, so the inputs are recoverable from the log even though no
file carries them.

The comparison targets — E1's parameter and operator counts — are read from
`../2026-08-11_e1_invariant_trotter_ansatz/results/aggregate.csv` by a human and
transcribed into the `E1_COUNTS` table of `code/run.py`, where they are asserted at run
time. That transcription is the one place a copy exists, and the assertion is what keeps
it honest.
