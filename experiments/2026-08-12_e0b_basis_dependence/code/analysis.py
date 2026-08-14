"""E0 figures and the visualisation contract.

One manuscript figure: the number of coupled-cluster amplitude classes surviving each
filter, with the correlation-energy change the filter causes written on the bars.
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


def main():
    with open(ROOT / "results" / "aggregate.csv") as fh:
        rows = list(csv.DictReader(fh))
    plt = X.use_style(X.FIGSIZE_SINGLE)

    mols = ["H2O", "NH3", "CH4"]
    pretty = {"H2O": r"$\mathrm{H_2O}$" + "\n" + r"$C_{2v}$",
              "NH3": r"$\mathrm{NH_3}$" + "\n" + r"$C_{3v}$",
              "CH4": r"$\mathrm{CH_4}$" + "\n" + r"$T_d$"}
    filters = [("none", "no filter", X.OKABE_ITO[0]),
               ("abelian", "Abelian subgroup", X.OKABE_ITO[2]),
               ("full", "full point group", X.OKABE_ITO[3])]

    def get(m, f, key):
        for r in rows:
            if r["molecule"] == m and r["filter"] == f:
                return r[key]
        return None

    fig, (ax, bx) = plt.subplots(2, 1, figsize=X.FIGSIZE_WIDE, sharex=True,
                                 gridspec_kw=dict(height_ratios=[1.5, 1.0], hspace=0.12))
    w = 0.26
    xs = np.arange(len(mols))
    for k, (fkey, flabel, colour) in enumerate(filters):
        vals = [int(get(m, fkey, "n_amplitudes")) for m in mols]
        pos = xs + (k - 1) * w
        bars = ax.bar(pos, vals, w, label=flabel, color=colour,
                      edgecolor="black", linewidth=0.4, zorder=3)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v * 1.15, str(v),
                    ha="center", va="bottom", fontsize=5.5)
        if fkey == "none":
            continue
        dv = [max(abs(float(get(m, fkey, "delta_e_corr_mha"))), 1e-13) for m in mols]
        bx.scatter(pos, dv, s=16, color=colour, marker=X.MARKERS[k],
                   edgecolors="black", linewidths=0.3, zorder=3)

    ax.set_yscale("log")
    ax.set_ylim(1, 1400)
    ax.set_ylabel("CCSD amplitude\nclasses")
    ax.legend(loc="upper left", frameon=False, fontsize=6, ncol=1)

    bx.axhline(1.6, color=X.OKABE_ITO[0], ls="--", lw=0.9, zorder=2)
    bx.text(-0.45, 2.4, "chemical accuracy", ha="left", va="bottom", fontsize=6,
            color=X.OKABE_ITO[0])
    bx.set_yscale("log")
    bx.set_ylim(1e-13, 60)
    bx.set_ylabel(r"$|\Delta E_{\mathrm{corr}}|$ (mHa)")
    bx.set_xticks(xs)
    bx.set_xticklabels([pretty[m] for m in mols])
    bx.set_xlim(-0.5, len(mols) - 0.5)
    fig.savefig(ROOT / "figures" / f"{EID}_amplitude_classes.pdf")
    fig.savefig(ROOT / "figures" / f"{EID}_amplitude_classes.png", format="png")
    plt.close(fig)

    manifest = {
        "amplitude_classes": {
            "panels": ["(a) top: surviving CCSD amplitude classes per molecule under three filters",
                       "(b) bottom: absolute change in CCSD correlation energy caused by each filter, "
                       "shared x axis with (a)"],
            "elements": [
                {"artist": "bar group 1 (left of each tick)", "encodes":
                 "amplitude classes with no symmetry filter",
                 "visual": "black #000000 bars, 0.4 pt black edge", "panel": "all"},
                {"artist": "bar group 2 (centre)", "encodes":
                 "amplitude classes surviving the Abelian subgroup PySCF actually works in "
                 "(Cs for NH3, D2 for CH4, C2v for H2O)",
                 "visual": "sky blue #56B4E9 bars", "panel": "all"},
                {"artist": "bar group 3 (right)", "encodes":
                 "amplitude classes surviving the full point group",
                 "visual": "bluish green #009E73 bars", "panel": "all"},
                {"artist": "numeric label above every bar", "encodes":
                 "the exact class count", "visual": "black text, 5.5 pt", "panel": "(a)"},
                {"artist": "scatter, Abelian", "encodes":
                 "absolute correlation-energy change caused by the Abelian filter, in mHa; "
                 "values below 1e-13 are drawn at 1e-13",
                 "visual": "sky blue #56B4E9 squares, black edge 0.3 pt", "panel": "(b)"},
                {"artist": "scatter, full group", "encodes":
                 "absolute correlation-energy change caused by the full-group filter, in mHa",
                 "visual": "bluish green #009E73 triangles, black edge 0.3 pt", "panel": "(b)"},
                {"artist": "horizontal dashed line", "encodes": "chemical accuracy, 1.6 mHa",
                 "visual": "black dashed, 0.9 pt, labelled at left", "panel": "(b)"},
                {"artist": "legend", "encodes": "the three filters",
                 "visual": "upper left of (a), no frame, 6 pt", "panel": "(a)"},
            ],
            "axes": {"x": "molecule and its full point group, none (shared)",
                     "y": "(a) number of surviving CCSD amplitude classes (singles plus doubles); "
                          "(b) absolute correlation-energy change in mHa",
                     "xscale": "linear", "yscale": "log in both panels"},
            "params": {"basis_set": "STO-3G", "reference": "RHF in the adapted basis",
                       "cc_conv_tol": "1e-10", "cc_conv_tol_normt": "1e-8",
                       "h2o_role": "Abelian control: its full group is already Abelian, so the "
                                   "second and third bars must be equal"},
            "reference_lines": ["chemical accuracy 1.6 mHa in panel (b) (dashed)"],
            "normalization": "none; absolute counts and absolute mHa. Panel (b) floors values "
                             "below 1e-13 mHa at 1e-13 so they remain visible on the log axis.",
        }
    }
    (ROOT / "figures" / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    viz = {"figures": [
        {"title": "Surviving CCSD amplitude classes by filter",
         "plot_type": "bar", "x_axis": "molecule", "y_axis": "n_amplitudes",
         "z_axis": None, "series": "filter", "filters": {},
         "log_scale_x": False, "log_scale_y": True},
        {"title": "Correlation-energy change caused by the filter",
         "plot_type": "bar", "x_axis": "molecule", "y_axis": "delta_e_corr_mha",
         "z_axis": None, "series": "filter",
         "filters": {"filter": ["abelian", "full"]},
         "log_scale_x": False, "log_scale_y": False},
    ]}
    (ROOT / "results" / "viz_schema.json").write_text(json.dumps(viz, indent=2) + "\n")
    print("figures:", sorted(p.name for p in (ROOT / "figures").glob("*")))


if __name__ == "__main__":
    main()
