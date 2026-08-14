"""E1 figures and the visualisation contract.

One manuscript figure: what the full point group buys in parameters and in elementary
gates, against the unfiltered ansatz and against the Abelian filter in the two forms it
can be built.
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
    plt = X.use_style()

    mols = ["H2O / STO-3G", "C2H4 / STO-3G", "NH3 / STO-3G", "CH4 / STO-3G"]
    pretty = [r"$\mathrm{H_2O}$" + "\n" + r"$C_{2v}$",
              r"$\mathrm{C_2H_4}$" + "\n" + r"$D_{2h}$",
              r"$\mathrm{NH_3}$" + "\n" + r"$C_{3v}$",
              r"$\mathrm{CH_4}$" + "\n" + r"$T_d$"]
    series = [("uccsd", "unfiltered UCCSD", X.OKABE_ITO[0]),
              ("abelian_subset", "Abelian filter, pool subset", X.OKABE_ITO[5]),
              ("abelian", "Abelian filter, invariant operators", X.OKABE_ITO[2]),
              ("full", "full point group, invariant operators", X.OKABE_ITO[3])]

    def get(m, b, key):
        for r in rows:
            if r["molecule"] == m and r["basis"] == b:
                return r[key]
        return ""

    fig, (ax, bx) = plt.subplots(2, 1, figsize=X.FIGSIZE_WIDE, sharex=True,
                                 gridspec_kw=dict(height_ratios=[1, 1], hspace=0.10))
    w = 0.2
    xs = np.arange(len(mols))
    for k, (bkey, label, colour) in enumerate(series):
        pos = xs + (k - 1.5) * w
        for axis, field in ((ax, "n_params"), (bx, "n_gates")):
            vals, px = [], []
            for i, m in enumerate(mols):
                v = get(m, bkey, field)
                if v not in ("", None):
                    vals.append(int(v))
                    px.append(pos[i])
            # The full-group green (#009E73, luminance 106) and the pool-subset
            # blue (#0072B2, luminance 87) are 19/255 apart and merge into one
            # tone in grayscale; a hatch on the green separates them without
            # touching the palette.
            bars = axis.bar(px, vals, w, label=label if axis is ax else None,
                            color=colour, edgecolor="black", linewidth=0.4,
                            hatch="///" if bkey == "full" else None, zorder=3)
            for b, v in zip(bars, vals):
                # Rotated: at 4.677in the bars are ~0.2in apart and a horizontal
                # 4-digit label at the 7pt floor would collide with its neighbour.
                axis.text(b.get_x() + b.get_width() / 2, v * 1.25, str(v),
                          ha="center", va="bottom", fontsize=7, rotation=90)

    # On a log axis floored at 10^0 an absent bar is pixel-identical to a bar of
    # value 1, so the vacant C2H4 pool-subset slot has to say so in the figure
    # and not only in the caption.
    ax.text(xs[1] + (1 - 1.5) * w, 1.5, "n/a", ha="center", va="bottom",
            fontsize=7, rotation=90, color="0.45")

    ax.set_yscale("log")
    ax.set_ylim(1, 20000)
    ax.set_ylabel("variational parameters")
    # Legend above the axes: inside the panel it collided with the 1224 label,
    # which cleared it by 4px (0.99pt at print size).
    # Matplotlib fills a multi-column legend column-major, which would present the
    # series 1,3,2,4 to a reader scanning by row; reorder so the legend reads in
    # the same left-to-right order as the bars.
    h, l = ax.get_legend_handles_labels()
    order = [0, 2, 1, 3]
    ax.legend([h[i] for i in order], [l[i] for i in order],
              loc="lower left", bbox_to_anchor=(0.0, 1.01), frameon=False,
              fontsize=7, ncol=2, handlelength=1.2, columnspacing=1.0,
              handletextpad=0.4, borderaxespad=0.0)
    bx.set_yscale("log")
    bx.set_ylim(1, 8000)
    # Not "gates": Tab:unitary explicitly disclaims a two-qubit gate count, and
    # the qualifier "elementary" is kept so the axis matches caption and table.
    bx.set_ylabel("distinct elementary\nexcitation operators")
    bx.set_xticks(xs)
    bx.set_xticklabels(pretty)
    bx.set_xlim(-0.5, len(mols) - 0.5)
    bx.text(1.0, 2.2, "energy not computed:\n" + r"$\sim 9\times 10^{6}$ determinants",
            ha="center", va="bottom", fontsize=7, color="0.35")

    fig.savefig(ROOT / "figures" / f"{EID}_compression_params_gates.pdf")
    fig.savefig(ROOT / "figures" / f"{EID}_compression_params_gates.png", format="png")
    plt.close(fig)

    manifest = {
        "compression_params_gates": {
            "panels": ["upper (no letter drawn): number of variational parameters of each ansatz",
                       "lower (no letter drawn): number of distinct elementary excitation operators of the "
                       "same ansaetze; shared x axis with (a). NOT a two-qubit gate count: "
                       "Tab:unitary disclaims that reading explicitly"],
            "elements": [
                {"artist": "bar group 1 (leftmost of each tick)", "encodes":
                 "unfiltered singlet UCCSD",
                 "visual": "black #000000 bars, 0.4 pt black edge", "panel": "all"},
                {"artist": "bar group 2", "encodes":
                 "the Abelian filter built as a subset of the UCCSD pool selected by irrep "
                 "label, which is how the literature implements SymUCCSD",
                 "visual": "blue #0072B2 bars", "panel": "all"},
                {"artist": "bar group 3", "encodes":
                 "the Abelian filter built instead from projector-derived invariant operators",
                 "visual": "sky blue #56B4E9 bars", "panel": "all"},
                {"artist": "bar group 4 (rightmost)", "encodes":
                 "the full point group, built from projector-derived invariant operators in the "
                 "pivoted {P e_j} basis", "visual": "bluish green #009E73 bars with a /// hatch, which separates them from the pool-subset blue in grayscale", "panel": "all"},
                {"artist": "numeric label above every bar", "encodes":
                 "the exact count", "visual": "black text, 7 pt, rotated 90 degrees", "panel": "all"},
                {"artist": "missing pool-subset bar at the C2H4 tick, upper panel", "encodes":
                 "the pool-subset count was not computed for C2H4, which was measured by "
                 "character counting only",
                 "visual": "an empty slot carrying a grey 'n/a' label, 7 pt, rotated 90 "
                 "degrees. The label is required because the log axis is floored at 10^0, "
                 "where an absent bar is otherwise pixel-identical to a bar of value 1",
                 "panel": "upper"},
                {"artist": "missing bars at the C2H4 tick in (b)", "encodes":
                 "gate counts and energies were not computed for C2H4",
                 "visual": "absence, explained by the grey note below", "panel": "(b)"},
                {"artist": "grey note under the C2H4 tick", "encodes":
                 "why C2H4 carries counts only", "visual": "grey text, 7 pt", "panel": "(b)"},
                {"artist": "legend", "encodes": "the four ansaetze",
                 "visual": "above the (a) axes, two columns, no frame, 7 pt", "panel": "(a)"},
            ],
            "axes": {"x": "molecule and its full point group, none (shared)",
                     "y": "(a) count of variational parameters; (b) count of distinct elementary "
                          "excitation operators", "xscale": "linear",
                     "yscale": "log in both panels"},
            "params": {"basis_set": "STO-3G", "optimiser": "L-BFGS-B, ftol 1e-16, gtol 1e-11",
                       "h2o_and_c2h4_role": "Abelian controls: their full group is already "
                                            "Abelian, so groups 3 and 4 must be equal",
                       "gate_definition": "distinct elementary excitation operators, not "
                                          "two-qubit gates; no device decomposition was done"},
            "normalization": "none; absolute counts",
        }
    }
    (ROOT / "figures" / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    viz = {"figures": [
        {"title": "Variational parameters by ansatz",
         "plot_type": "bar", "x_axis": "molecule", "y_axis": "n_params",
         "z_axis": None, "series": "basis", "filters": {},
         "log_scale_x": False, "log_scale_y": True},
        {"title": "Distinct elementary gates by ansatz",
         "plot_type": "bar", "x_axis": "molecule", "y_axis": "n_gates",
         "z_axis": None, "series": "basis", "filters": {},
         "log_scale_x": False, "log_scale_y": True},
        {"title": "Energy error under exact exponentials and under Trotter",
         "plot_type": "bar", "x_axis": "molecule", "y_axis": "err_trotter_natural_mha",
         "z_axis": None, "series": "basis",
         "filters": {"basis": ["abelian", "full"]},
         "log_scale_x": False, "log_scale_y": False},
    ]}
    (ROOT / "results" / "viz_schema.json").write_text(json.dumps(viz, indent=2) + "\n")
    print("figures:", sorted(p.name for p in (ROOT / "figures").glob("*")))


if __name__ == "__main__":
    main()
