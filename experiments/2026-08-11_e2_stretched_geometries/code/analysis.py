"""E2 figures and the visualisation contract.

One manuscript figure: what each filter costs as NH3 is pulled apart, against the
collapsing S0-S1 gap that marks the onset of static correlation.
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
FLOOR = 1e-11


def main():
    with open(ROOT / "results" / "aggregate.csv") as fh:
        rows = [r for r in csv.DictReader(fh) if r["converged_rhf"] == "True"]
    plt = X.use_style()
    import matplotlib.ticker as mticker  # after use_style pins the backend

    # The exact and Trotter series are 13/255 apart in grayscale luminance
    # (#009E73 = 106, #D55E00 = 119) and cross four times, so line style, not
    # colour, has to carry the distinction.
    series = [("abelian", r"Abelian $C_s$ filter, 75 par.", X.OKABE_ITO[5], X.MARKERS[0], "-"),
              ("full_c3v_exact", r"full $C_{3v}$, 30 par., exact", X.OKABE_ITO[3], X.MARKERS[1], "-"),
              ("full_c3v", r"full $C_{3v}$, 30 par., Trotter", X.OKABE_ITO[6], X.MARKERS[2], "-.")]

    fig, (ax, bx) = plt.subplots(2, 1, figsize=X.FIGSIZE_WIDE, sharex=True,
                                 gridspec_kw=dict(height_ratios=[1.6, 1.0], hspace=0.12))
    for key, label, colour, marker, ls in series:
        sub = sorted([r for r in rows if r["ansatz"] == key],
                     key=lambda r: float(r["r_nh_angstrom"]))
        x = [float(r["r_nh_angstrom"]) for r in sub]
        y = [max(abs(float(r["delta_vs_uccsd_mha"])), FLOOR) for r in sub]
        ax.plot(x, y, marker=marker, ms=3.2, lw=1.0, ls=ls, color=colour,
                label=label, zorder=3)

    ax.axhline(1.6, color=X.OKABE_ITO[0], ls="--", lw=0.9, zorder=2)
    ax.text(0.92, 3.5, "chemical accuracy", ha="left", va="bottom", fontsize=7,
            color=X.OKABE_ITO[0])
    ax.set_yscale("log")
    ax.set_ylim(FLOOR / 3, 2000)
    ax.set_ylabel("filter cost vs.\nunfiltered UCCSD (mHa)")
    # handlelength 3.0, not 2.0: at 2.0 the triangle marker covers the dot of the
    # dash-dot cycle, so the Trotter key read "dashed" while the plotted line was
    # correctly dash-dot.
    ax.legend(loc="center left", bbox_to_anchor=(0.01, 0.30), frameon=False,
              fontsize=7, handlelength=3.0, handletextpad=0.5, labelspacing=0.3)

    unc = sorted([r for r in rows if r["ansatz"] == "uccsd"],
                 key=lambda r: float(r["r_nh_angstrom"]))
    x = [float(r["r_nh_angstrom"]) for r in unc]
    # EN DASH, not "--": the scienceplots "no-latex" style renders text through
    # mathtext/DejaVu, where "--" stays two literal hyphens on the page.
    bx.plot(x, [float(r["gap_s0_s1_mha"]) for r in unc], marker=X.MARKERS[3], ms=3.2, lw=1.0,
            color=X.OKABE_ITO[0], label="$S_0$–$S_1$ gap", zorder=3)
    bx.plot(x, [abs(float(r["corr_energy_mha"])) for r in unc], marker=X.MARKERS[4], ms=3.2,
            lw=1.0, ls="--", color=X.OKABE_ITO[1], label="correlation energy", zorder=3)
    bx.set_yscale("log")
    # The data span 20-628 mHa, i.e. 1.5 decades, so the decade locator found a
    # single major tick ("10^2") and the panel could not be read.  Tick it
    # explicitly and print plain numbers.
    bx.set_ylim(15, 900)
    bx.set_yticks([20, 50, 100, 200, 500])
    bx.yaxis.set_major_formatter(mticker.ScalarFormatter())
    bx.yaxis.set_minor_formatter(mticker.NullFormatter())
    bx.set_ylabel("energy scale (mHa)")
    bx.set_xlabel("N–H distance ($\\mathrm{\\AA}$)")
    bx.legend(loc="lower left", frameon=False, fontsize=7, handlelength=2.0,
              handletextpad=0.5, labelspacing=0.3)
    bx.set_xlim(0.85, 2.28)
    for a in (ax, bx):
        a.axvline(2.2, color="0.7", lw=0.8, ls=":", zorder=1)
    bx.text(2.17, 22, "RHF fails", rotation=90, ha="right", va="bottom",
            fontsize=7, color="0.45")

    fig.savefig(ROOT / "figures" / f"{EID}_filter_cost_vs_stretch.pdf")
    fig.savefig(ROOT / "figures" / f"{EID}_filter_cost_vs_stretch.png", format="png")
    plt.close(fig)

    manifest = {
        "filter_cost_vs_stretch": {
            "panels": ["upper (no letter drawn): cost of each symmetry filter relative to unfiltered UCCSD, "
                       "against N-H distance",
                       "lower (no letter drawn): the two energy scales that define the regime, shared x axis"],
            "elements": [
                {"artist": "line with circles", "encodes":
                 "absolute energy difference between the Abelian Cs filter (75 parameters) and "
                 "unfiltered UCCSD; values below 1e-11 mHa are drawn at 1e-11",
                 "visual": "blue #0072B2, 1.0 pt, circles 3 pt", "panel": "(a)"},
                {"artist": "line with squares", "encodes":
                 "same difference for the full C3v filter (30 parameters) with exact "
                 "exponentials of the invariant operators",
                 "visual": "bluish green #009E73, squares 3 pt", "panel": "(a)"},
                {"artist": "line with triangles", "encodes":
                 "same difference for the full C3v filter in tied Trotter product form",
                 "visual": "vermillion #D55E00, DASH-DOTTED line, triangles 3.2 pt", "panel": "(a)"},
                {"artist": "horizontal dashed line", "encodes": "chemical accuracy, 1.6 mHa",
                 "visual": "black dashed, 0.9 pt, labelled", "panel": "(a)"},
                {"artist": "line with diamonds", "encodes":
                 "S0-S1 gap of the exact ground state, in mHa",
                 "visual": "black #000000, diamonds 3 pt", "panel": "(b)"},
                {"artist": "line with inverted triangles", "encodes":
                 "absolute correlation energy, in mHa",
                 "visual": "orange #E69F00, DASHED line, inverted triangles 3.2 pt", "panel": "(b)"},
                {"artist": "vertical dotted rule at 2.2 A", "encodes":
                 "the distance at which RHF stops converging, so the invariant closed-shell "
                 "reference the theorems assume ceases to exist",
                 "visual": "grey dotted, 0.8 pt in both panels, labelled 'RHF fails' in (b)",
                 "panel": "all"},
                {"artist": "legend of (a)", "encodes": "the three filters",
                 "visual": "centre left, no frame, 7 pt", "panel": "(a)"},
                {"artist": "legend of (b)", "encodes": "the two energy scales",
                 "visual": "lower left, no frame, 7 pt", "panel": "(b)"},
            ],
            "axes": {"x": "N-H distance in angstrom (shared)",
                     "y": "(a) absolute energy difference from unfiltered UCCSD in mHa; "
                          "(b) energy in mHa",
                     "xscale": "linear", "yscale": "log in both panels"},
            "params": {"basis_set": "STO-3G", "hnh_angle_deg": "106.7 (fixed)",
                       "n_determinants": "3136", "uccsd_parameters": "135",
                       "abelian_parameters": "75", "full_c3v_parameters": "30",
                       "optimiser": "L-BFGS-B, ftol 1e-16, gtol 1e-11"},
            "reference_lines": ["chemical accuracy 1.6 mHa in (a) (dashed)",
                                "RHF convergence boundary at 2.2 A (dotted, both panels)"],
            "normalization": "none; absolute mHa. Panel (a) floors values below 1e-11 mHa at "
                             "1e-11 so the Abelian series stays visible on the log axis.",
        }
    }
    (ROOT / "figures" / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    viz = {"figures": [
        {"title": "Filter cost against N-H distance",
         "plot_type": "line", "x_axis": "r_nh_angstrom", "y_axis": "delta_vs_uccsd_mha",
         "z_axis": None, "series": "ansatz",
         "filters": {"ansatz": ["abelian", "full_c3v", "full_c3v_exact"]},
         "log_scale_x": False, "log_scale_y": False},
        {"title": "Error against FCI for every ansatz",
         "plot_type": "line", "x_axis": "r_nh_angstrom", "y_axis": "err_vs_fci_mha",
         "z_axis": None, "series": "ansatz", "filters": {},
         "log_scale_x": False, "log_scale_y": True},
        {"title": "Onset of static correlation along the scan",
         "plot_type": "line", "x_axis": "r_nh_angstrom", "y_axis": "gap_s0_s1_mha",
         "z_axis": None, "series": None, "filters": {"ansatz": "uccsd"},
         "log_scale_x": False, "log_scale_y": True},
    ]}
    (ROOT / "results" / "viz_schema.json").write_text(json.dumps(viz, indent=2) + "\n")
    print("figures:", sorted(p.name for p in (ROOT / "figures").glob("*")))


if __name__ == "__main__":
    main()
