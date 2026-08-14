"""Smoke test for 2026-08-11_e7_degenerate_shell_artifact — runs a reduced configuration in under a minute.

Exercises the whole pipeline end to end before anything expensive is started.
The smoke run writes its log to a temporary directory, so it never contaminates
the canonical append-only experiment_log.json.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "vendor"))

import run as R  # noqa: E402


def main():
    sys.argv = ["run.py", "all", "--n", "2", "--workers", "2", "--smoke"]
    agg = R.main()
    assert agg["predictor_scored"] == 10, agg
    assert agg["predictor_hits"] == agg["predictor_scored"], agg
    print("\nsmoke_test: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
