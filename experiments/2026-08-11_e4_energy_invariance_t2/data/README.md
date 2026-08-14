# data/

This experiment stores no input files. Every input is generated deterministically by
`code/run.py` from values written in the source: molecular geometries, basis-set name,
group elements, and the sampled parameter vectors, which come from `numpy.random.default_rng`
with the seeds recorded in `experiment_log.json`.

Keeping the inputs in code rather than on disk is deliberate here, because the inputs are
short and the generating expressions are what a reader needs in order to check them. The
provenance of the reused numerical core is in `code/vendor/PROVENANCE.md`.
