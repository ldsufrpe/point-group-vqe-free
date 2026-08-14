"""E7 figures and the visualisation contract.

One manuscript figure: the structural predictor over the benchmark set, and the
filtered error measured over independent processes against the value the same filter
gives when the irrep labels describe the orbitals they are attached to.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import expcommon as X  # noqa: E402
import numpy as np  # noqa: E402

ROOT = HERE.parent
EID = ROOT.name

# Measured in experiments/2026-08-11_e1_invariant_trotter_ansatz on the same geometry
# and basis: the filtered error when the labels come from the same calculation that
# defines the Hamiltonian.
CORRECT_MHA = {"NH3": 0.1392, "CH4": 0.1942}


def main():
    with open(ROOT / "results" / "aggregate.csv") as fh:
        rows = list(csv.DictReader(fh))
    plt = X.use_style()

    pred = [r for r in rows if r["part"] == "predictor" and r["basis"] == "sto-3g"
            and r["n_degenerate_shells"] != ""]
    rep = [r for r in rows if r["part"] == "replicate"]

    fig, (ax, bx) = plt.subplots(2, 1, figsize=X.FIGSIZE_WIDE,
                                 gridspec_kw=dict(height_ratios=[1.0, 1.25], hspace=0.16))

    # --- (a) structural predictor
    pred = sorted(pred, key=lambda r: (int(r["n_degenerate_shells"]), r["molecule"]))
    names = [r["molecule"] for r in pred]
    nsh = [int(r["n_degenerate_shells"]) for r in pred]
    obs = [r["observed_failure"] == "True" for r in pred]
    for i, (n, o) in enumerate(zip(nsh, obs)):
        ax.scatter([i], [n], s=42, marker=X.MARKERS[0] if o else X.MARKERS[1],
                   color=X.OKABE_ITO[6] if o else X.OKABE_ITO[5],
                   edgecolors="black", linewidths=0.4, zorder=3)
    ax.axhline(1.5, color=X.OKABE_ITO[0], ls="--", lw=0.9, zorder=2)
    ax.text(-0.4, 1.60, "predictor threshold: two or more shells", ha="left", va="bottom",
            fontsize=7, color=X.OKABE_ITO[0])
    ax.set_xticks(range(len(names)))
    sub = {"C2H4": r"$\mathrm{C_2H_4}$", "H2O": r"$\mathrm{H_2O}$", "BeH2": r"$\mathrm{BeH_2}$",
           "HF": "HF", "LiH": "LiH", "CH4": r"$\mathrm{CH_4}$", "CO": "CO",
           "N2": r"$\mathrm{N_2}$", "NH3": r"$\mathrm{NH_3}$", "NaH": "NaH"}
    ax.set_xticklabels([sub.get(x, x) for x in names], fontsize=7)
    ax.set_ylabel("degenerate orbital shells")
    # Integer count: half-integer ticks advertise values the quantity cannot take.
    ax.set_yticks([0, 1, 2, 3])
    ax.set_ylim(-0.45, 3.2)
    ax.set_xlim(-0.5, len(names) - 0.5)
    hit = sum(1 for r in pred if r["prediction_correct"] == "True")
    ax.set_title(f"reported as failing: circles ({hit} of {len(pred)} predicted)",
                 fontsize=7, loc="left")

    # --- (b) replicate distribution
    mols = sorted({r["molecule"] for r in rep})
    threads = sorted({int(r["blas_threads"]) for r in rep})
    # Panel (a) encodes the REPORTED OUTCOME in colour (#D55E00 failing,
    # #0072B2 working); panel (b) encodes the MOLECULE.  Reusing either colour
    # here made CH4 vermillion above and blue below, and made a square mean
    # "working" above and "NH3" below.  Pick two colours (a) does not use, far
    # apart in grayscale luminance: #009E73 ~106 against #E69F00 ~168.
    colours = {m: c for m, c in zip(mols, [X.OKABE_ITO[3], X.OKABE_ITO[1]])}
    rng = np.random.default_rng(11)
    for mi, m in enumerate(mols):
        for ti, t in enumerate(threads):
            v = [float(r["err_sym_mha"]) for r in rep
                 if r["molecule"] == m and int(r["blas_threads"]) == t]
            if not v:
                continue
            xpos = ti + (mi - 0.5) * 0.3
            # s=14 printed a 2.57pt marker, at which circle and square are
            # indistinguishable and the molecules separated by colour alone.
            bx.scatter(np.full(len(v), xpos) + rng.uniform(-0.05, 0.05, len(v)), v,
                       s=30, color=colours[m], marker=X.MARKERS[mi], alpha=0.75,
                       linewidths=0, zorder=3,
                       label=sub.get(m, m) if ti == 0 else None)
        # The unfiltered error and the correctly labelled filtered error coincide to
        # three decimals, because the filter is free when the labels describe the
        # orbitals.  One line carries both.
        u = [float(r["err_uccsd_mha"]) for r in rep if r["molecule"] == m]
        bx.axhline(CORRECT_MHA[m], color=colours[m], ls=":", lw=1.1, zorder=2)
        # The two lines are only 0.15 decades apart, which at the 7pt floor is
        # exactly one line of text, so the labels cannot both sit above.  Put
        # CH4's above its line and NH3's below, anchored on the edge that faces
        # away from the line: va="bottom" at *0.80 previously drew the NH3 line
        # straight through the cap height of its own label.
        above = (m == "CH4")
        bx.text(-0.45, CORRECT_MHA[m] * (1.15 if above else 0.87),
                f"{sub.get(m, m)}: correct labels, {CORRECT_MHA[m]:.3f} mHa "
                f"(unfiltered {np.mean(u):.3f})",
                ha="left", va="bottom" if above else "top",
                # Amber at 2.25:1 on white is the faintest ink in the figure, so
                # the LABEL is darkened while the marker keeps the colour the
                # caption names.
                fontsize=7, color="#B37B00" if m == "NH3" else colours[m])

    bx.set_yscale("log")
    # Headroom below for the NH3 reference label, which now hangs under its line.
    bx.set_ylim(0.055, 90)
    bx.set_xticks(range(len(threads)))
    bx.set_xticklabels([str(t) for t in threads])
    bx.set_xlabel("BLAS threads in the replicate process")
    bx.set_ylabel("filtered error (mHa)")
    bx.set_xlim(-0.5, len(threads) - 0.5)
    bx.legend(loc="center left", frameon=False, fontsize=7, handletextpad=0.4,
              labelspacing=0.3, borderaxespad=0.4)

    fig.savefig(ROOT / "figures" / f"{EID}_artifact_predictor_and_spread.pdf")
    fig.savefig(ROOT / "figures" / f"{EID}_artifact_predictor_and_spread.png", format="png")
    plt.close(fig)

    manifest = {
        "artifact_predictor_and_spread": {
            "panels": ["(a) top: number of degenerate orbital shells for the ten benchmark "
                       "molecules in STO-3G, against the outcome reported in the literature",
                       "(b) bottom: filtered error measured in independent processes, one point "
                       "per replicate, grouped by the BLAS thread count of the process"],
            "elements": [
                {"artist": "filled circles", "encodes":
                 "molecules the audited work reports as failing",
                 "visual": "vermillion #D55E00 circles, 42 pt, black edge", "panel": "(a)"},
                {"artist": "filled squares", "encodes":
                 "molecules reported as working",
                 "visual": "blue #0072B2 squares, 42 pt, black edge", "panel": "(a)"},
                {"artist": "horizontal dashed line at 1.5", "encodes":
                 "the predictor threshold: two or more degenerate shells predicts failure",
                 "visual": "black dashed, 0.9 pt, labelled", "panel": "(a)"},
                {"artist": "panel title", "encodes":
                 "how many of the scored molecules the predictor gets right",
                 "visual": "left-aligned 7 pt text above the axes", "panel": "(a)"},
                {"artist": "scatter, CH4 replicates", "encodes":
                 "filtered error of one independent process; horizontal jitter is cosmetic",
                 "visual": "bluish green #009E73 circles, alpha 0.75, size 30. Panel (b) deliberately uses colours panel (a) does not, so no colour reading carries between panels", "panel": "(b)"},
                {"artist": "scatter, NH3 replicates", "encodes": "same, for ammonia",
                 "visual": "orange #E69F00 squares, alpha 0.75, size 30", "panel": "(b)"},
                {"artist": "dotted horizontal line per molecule", "encodes":
                 "two coinciding quantities: the cost of the same filter when the irrep labels "
                 "come from the calculation that defines the Hamiltonian, measured in the "
                 "companion experiment E1, and the mean unfiltered UCCSD error over these "
                 "replicates, which is the rotation-invariant control. They agree to three "
                 "decimals because the filter is free when the labels are correct, so one line "
                 "carries both and its label prints both numbers",
                 "visual": "dotted in the molecule's panel-(b) colour, 1.1 pt, labelled at left in 7 pt; the CH4 label sits above its line and the NH3 label below, because at 7 pt the two lines are one line of text apart",
                 "panel": "(b)"},
                {"artist": "legend of (b)", "encodes": "the two molecules",
                 "visual": "centre left, no frame, 7 pt", "panel": "(b)"},
            ],
            "axes": {"x": "(a) molecule, none; (b) BLAS thread count of the replicate process",
                     "y": "(a) count of degenerate orbital shells; (b) filtered error in mHa",
                     "xscale": "linear", "yscale": "(a) linear; (b) log"},
            "params": {"basis_set": "STO-3G", "geometries": "those of the audited work",
                       "replicates_per_thread_count": "10, except NH3 at one thread which has 20",
                       "degeneracy_tolerance_ha": "1e-6",
                       "reference_values_source": "experiments/2026-08-11_e1_invariant_trotter_ansatz",
                       "horizontal_jitter": "uniform +-0.05, seed 11, cosmetic only"},
            "reference_lines": [
                "predictor threshold at 1.5 shells in (a) (dashed)",
                "filter cost with correct labels per molecule in (b) (dotted)",
                "unfiltered UCCSD error per molecule, coinciding with the line above"],
            "normalization": "none; absolute counts and absolute mHa",
        }
    }
    (ROOT / "figures" / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    viz = {"figures": [
        {"title": "Filtered error across independent processes",
         "plot_type": "scatter", "x_axis": "replicate", "y_axis": "err_sym_mha",
         "z_axis": None, "series": "molecule", "filters": {"part": "replicate"},
         "log_scale_x": False, "log_scale_y": True},
        {"title": "Unfiltered error across the same processes (control)",
         "plot_type": "scatter", "x_axis": "replicate", "y_axis": "err_uccsd_mha",
         "z_axis": None, "series": "molecule", "filters": {"part": "replicate"},
         "log_scale_x": False, "log_scale_y": False},
        {"title": "Degenerate shells against reported outcome",
         "plot_type": "bar", "x_axis": "molecule", "y_axis": "n_degenerate_shells",
         "z_axis": None, "series": "basis", "filters": {"part": "predictor"},
         "log_scale_x": False, "log_scale_y": False},
    ]}
    (ROOT / "results" / "viz_schema.json").write_text(json.dumps(viz, indent=2) + "\n")
    print("figures:", sorted(p.name for p in (ROOT / "figures").glob("*")))


if __name__ == "__main__":
    main()
