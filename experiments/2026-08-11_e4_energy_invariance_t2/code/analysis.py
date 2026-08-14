"""E4 figures and the visualisation contract.

One manuscript figure: the measured violation of the energy-invariance theorem in a
symmetry-adapted basis against the same quantity in a basis where a degenerate pair
has been rotated and the labels left stale.
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
CHEM_ACC_MHA = 1.6


def load():
    with open(ROOT / "results" / "aggregate.csv") as fh:
        return list(csv.DictReader(fh))


def main():
    rows = load()
    plt = X.use_style(X.FIGSIZE_SINGLE)

    inv = [r for r in rows if r["test"] == "T2_energy_invariance"]
    ctl = [r for r in rows if r["test"] == "T2_negative_control"]
    grad = [r for r in rows if r["test"] == "T3_removed_gradient"]

    pretty = {"NH3": r"$\mathrm{NH_3}$", "CH4": r"$\mathrm{CH_4}$"}
    groups = []
    for mol in ("NH3", "CH4"):
        groups.append((f"{pretty[mol]}" + "\nadapted",
                       [float(r["delta_energy_ha"]) for r in inv if r["molecule"] == mol]))
    for mol in ("NH3", "CH4"):
        groups.append((f"{pretty[mol]}" + "\nrotated",
                       [float(r["delta_energy_ha"]) for r in ctl if r["molecule"] == mol]))

    fig, ax = plt.subplots()
    colours = [X.OKABE_ITO[5], X.OKABE_ITO[2], X.OKABE_ITO[6], X.OKABE_ITO[1]]
    markers = [X.MARKERS[0], X.MARKERS[1], X.MARKERS[2], X.MARKERS[3]]
    rng = np.random.default_rng(7)
    for i, ((label, vals), c, m) in enumerate(zip(groups, colours, markers)):
        v = np.array([x for x in vals if x > 0])
        floor = 1e-17
        v = np.maximum(v, floor)
        # Wider jitter and smaller, semi-transparent markers: at s=7 with
        # +/-0.16 the 60 CH4-adapted points fused into two solid rectangles and
        # read as a broken line rather than as the individual points the
        # caption calls squares.
        jitter = rng.uniform(-0.26, 0.26, len(v))
        ax.scatter(np.full(len(v), i) + jitter, v, s=4.5, color=c, marker=m,
                   alpha=0.65, linewidths=0, zorder=3)

    ax.axhline(CHEM_ACC_MHA * 1e-3, color=X.OKABE_ITO[0], ls="--", lw=0.9, zorder=2)
    ax.text(-0.42, CHEM_ACC_MHA * 1e-3 * 2.2, "chemical accuracy", ha="left", va="bottom",
            fontsize=7, color=X.OKABE_ITO[0])
    # Shifted right and down: the wider jitter brought a slot-0 circle to within
    # 0.24 pt of the old anchor, which read as a smudge on the first glyph.
    ax.annotate("machine precision", xy=(0.72, 6e-14), ha="center", va="center",
                fontsize=7, color="0.35")
    ax.set_yscale("log")
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g[0] for g in groups])
    ax.set_ylabel(r"$|E(\theta)-E(\sigma_n\theta)|$  (Ha)")
    ax.set_xlim(-0.5, len(groups) - 0.5)
    ax.set_ylim(1e-17, 1e-1)
    ax.axvline(1.5, color="0.75", lw=0.7, zorder=1)
    fig.savefig(ROOT / "figures" / f"{EID}_invariance_violation.pdf")
    fig.savefig(ROOT / "figures" / f"{EID}_invariance_violation.png", format="png")
    plt.close(fig)

    manifest = {
        "invariance_violation": {
            "panels": ["single panel: measured violation of the energy-invariance identity, "
                       "one point per sampled parameter vector and group element"],
            "elements": [
                {"artist": "scatter, NH3 adapted basis", "encodes":
                 "|E(theta) - E(sigma_n theta)| for 20 parameter vectors, 1 nontrivial Cs element",
                 "visual": "blue #0072B2 circles, alpha 0.65, size 4.5", "panel": "all"},
                {"artist": "scatter, CH4 adapted basis", "encodes":
                 "same quantity for 20 parameter vectors, 3 nontrivial D2 elements",
                 "visual": "sky blue #56B4E9 squares, alpha 0.65, size 4.5", "panel": "all"},
                {"artist": "scatter, NH3 rotated basis", "encodes":
                 "same quantity with a degenerate pair rotated by 30 degrees and the irrep "
                 "labels left stale", "visual": "vermillion #D55E00 triangles", "panel": "all"},
                {"artist": "scatter, CH4 rotated basis", "encodes":
                 "same quantity, rotated-label control",
                 "visual": "orange #E69F00 diamonds", "panel": "all"},
                {"artist": "horizontal dashed line", "encodes": "chemical accuracy, 1.6 mHa",
                 "visual": "black dashed, 0.9 pt, labelled", "panel": "all"},
                {"artist": "vertical rule between category 2 and 3", "encodes":
                 "separates the adapted-basis categories from the rotated-label controls",
                 "visual": "grey solid, 0.7 pt", "panel": "all"},
                {"artist": "x tick labels", "encodes":
                 "molecule and basis of each category", "visual": "two-line text", "panel": "all"},
                {"artist": "annotation 'machine precision'", "encodes":
                 "the band the adapted-basis points occupy",
                 "visual": "grey text, 7 pt, no leader", "panel": "all"},
            ],
            "axes": {"x": "category (molecule and orbital basis), none",
                     "y": "absolute energy difference (Ha)",
                     "xscale": "linear", "yscale": "log"},
            "params": {"n_theta_per_molecule": "20", "theta_distribution": "uniform on [-0.5, 0.5]",
                       "control_rotation_deg": "30", "basis_set": "STO-3G",
                       "horizontal_jitter": "uniform +-0.26, seed 7, cosmetic only",
                       "zero_differences_dropped": "the y axis is logarithmic, so only "
                       "strictly positive differences can be drawn. Of the CH4 adapted "
                       "series 54 of 60 differences are EXACTLY zero and are not plotted, "
                       "leaving 6 markers; NH3 adapted draws all 20. In the controls, CH4 "
                       "draws 10 of 15 and NH3 all 5. The absences mean the identity held "
                       "to the last bit, which is a stronger statement than the plotted "
                       "points make"},
            "reference_lines": ["chemical accuracy 1.6 mHa (black dashed, labelled at left)"],
            "normalization": "none; absolute Hartree",
        }
    }
    (ROOT / "figures" / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n")

    viz = {"figures": [
        {"title": "Energy-invariance violation by molecule and basis",
         "plot_type": "scatter", "x_axis": "seed", "y_axis": "delta_energy_ha",
         "z_axis": None, "series": "molecule",
         "filters": {"test": "T2_energy_invariance"},
         "log_scale_x": False, "log_scale_y": True},
        {"title": "Negative control: same quantity with stale labels",
         "plot_type": "scatter", "x_axis": "seed", "y_axis": "delta_energy_ha",
         "z_axis": None, "series": "molecule",
         "filters": {"test": "T2_negative_control"},
         "log_scale_x": False, "log_scale_y": True},
        {"title": "Gradient in the removed directions on the symmetric submanifold",
         "plot_type": "scatter", "x_axis": "seed", "y_axis": "grad_max_removed",
         "z_axis": None, "series": "molecule",
         "filters": {"test": "T3_removed_gradient"},
         "log_scale_x": False, "log_scale_y": True},
    ]}
    (ROOT / "results" / "viz_schema.json").write_text(json.dumps(viz, indent=2) + "\n")

    print("figures:", sorted(p.name for p in (ROOT / "figures").glob("*")))
    print("max adapted   :", max(float(r["delta_energy_ha"]) for r in inv))
    print("max control   :", max(float(r["delta_energy_ha"]) for r in ctl))
    print("max grad removed:", max(float(r["grad_max_removed"]) for r in grad))


if __name__ == "__main__":
    main()
