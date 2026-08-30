"""Smoke test -- H2O / C2v, the Abelian control, in well under a minute.

H2O is the right smoke case because its full point group IS its Abelian subgroup, so
three independent things must hold at once and any one of them failing localises the
bug:

  1. the projector pool for the full group and for the Abelian subgroup are identical,
     since the two groups are the same set.  A difference means `group()` or the
     subgroup indexing is wrong;
  2. the counts reproduce E1's H2O row exactly, which catches vendored-source drift;
  3. every incidence is a single or a double, and the CNOT arithmetic is consistent
     with the class split.

It writes nothing into the experiment folder -- the log goes to a temporary directory
-- so a failed smoke test cannot leave a half-written aggregate behind.

Usage:  python smoke_test.py
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import run as R  # noqa: E402


def main():
    sys.argv = ["smoke_test.py", "h2o", "--smoke"]
    rows, all_res = R.main()
    res = all_res[0]
    P = res["pools"]
    fails = []

    def check(name, cond, detail):
        print(f"  [{'ok  ' if cond else 'FAIL'}] {name}: {detail}")
        if not cond:
            fails.append(name)

    print("\n=== smoke assertions ===")
    for cur in R.CURRENCIES:
        check(f"abelian control [{cur}]",
              (P["full"][cur]["n_gates"], P["full"][cur]["cnots"]) ==
              (P["abelian_projector"][cur]["n_gates"], P["abelian_projector"][cur]["cnots"]),
              f"full {P['full'][cur]['n_gates']} gates / {P['full'][cur]['cnots']} CNOTs "
              f"against projector-Abelian {P['abelian_projector'][cur]['n_gates']} / "
              f"{P['abelian_projector'][cur]['cnots']}; C2v is its own Abelian subgroup, "
              f"so these must be identical")

    check("representation validity", res["rep_error"] < 1e-10,
          f"|U^T S U - S| = {res['rep_error']:.2e}")

    for name, (p_e1, g_e1, cur) in R.E1_COUNTS["h2o"].items():
        check(f"E1 regression [{name}, {cur}]",
              (P[name]["n_params"], P[name][cur]["n_gates"]) == (p_e1, g_e1),
              f"{P[name]['n_params']} params / {P[name][cur]['n_gates']} gates "
              f"against E1's {p_e1} / {g_e1}")

    for name in ("abelian_projector", "full"):
        check(f"decomposition residual [{name}]", P[name]["decomp_residual"] < 1e-12,
              f"{P[name]['decomp_residual']:.2e}")

    for name, p in P.items():
        for cur in R.CURRENCIES:
            q = p[cur]
            check(f"class split exhaustive [{name}, {cur}]",
                  q["n_singles"] + q["n_doubles"] == q["n_gates"],
                  f"{q['n_singles']} singles + {q['n_doubles']} doubles = {q['n_gates']} gates")
            check(f"ruler arithmetic [{name}, {cur}]",
                  q["cnots"] == R.CNOT_SINGLE * q["n_singles"] + R.CNOT_DOUBLE * q["n_doubles"],
                  f"{q['cnots']} = 2*{q['n_singles']} + 13*{q['n_doubles']}")
        check(f"distinct never exceeds incidence [{name}]",
              p["distinct"]["n_gates"] <= p["incidence"]["n_gates"],
              f"{p['distinct']['n_gates']} <= {p['incidence']['n_gates']}")

    # The singles count is forced: every spatial single is spin-complemented into
    # exactly two spin-orbital singles, so the UCCSD pool must carry 2 * n_occ * n_virt,
    # and singles are never duplicated, so the two currencies must agree on them.
    no = res["nelec"] // 2
    expected_singles = 2 * no * (res["nao"] - no)
    for cur in R.CURRENCIES:
        check(f"spin complementation [{cur}]",
              P["uccsd"][cur]["n_singles"] == expected_singles,
              f"{P['uccsd'][cur]['n_singles']} spin-orbital singles against 2 x {no} x "
              f"{res['nao'] - no} = {expected_singles}")

    # The currency gap must be doubles-only: it comes from the antisymmetrisation of the
    # same-spin doubles, so no single may be duplicated anywhere.
    check("currency gap is doubles-only",
          all(p["distinct"]["n_singles"] == p["incidence"]["n_singles"] for p in P.values()),
          "every pool has the same singles count in both currencies")

    print(f"\n{len(fails)} failure(s)" + (": " + ", ".join(fails) if fails else ""))
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
