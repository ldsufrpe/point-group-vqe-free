"""Smoke test for E4 — runs in well under a minute on NH3 with three theta samples.

Checks the pipeline end to end and asserts the qualitative outcome the full run
must reproduce: invariance at machine precision in the adapted basis, and a
failure many orders of magnitude larger in the rotated-label control.
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "vendor"))

import run as R  # noqa: E402


def main():
    sys.argv = ["run.py", "nh3", "--smoke"]
    agg = R.main()
    assert agg["max_delta_energy_ha"] < 1e-9, agg
    assert agg["max_state_residual"] < 1e-9, agg
    assert agg["max_grad_removed"] < 1e-7, agg
    assert agg["control_max_delta_energy_ha"] > 1e-5, \
        "the negative control must break T2; if it does not, the test measures nothing"
    print("\nsmoke_test: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
