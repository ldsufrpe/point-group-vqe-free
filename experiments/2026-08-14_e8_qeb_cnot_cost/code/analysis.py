"""E8 figure and the visualisation contract.

One manuscript figure, two panels answering the two questions the experiment asks.

(a) is the validation: the CNOT cost of every pool under the QEB ruler, with He et al.'s
    six published Figure 5 values marked on the bars that have one.  All six coincide,
    which is what licenses the ruler to price the pool they never report.

(b) is the finding: going from the Abelian filter to the full point group, how much of
    each quantity is removed.  Parameters fall by two thirds; CNOTs fall by an eighth to
    a quarter, and by almost exactly as much as the flat operator count.  The gap between
    the first bar and the last is the result.

Everything is in the distinct-operator currency, the one Figure 5 is written in.
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

POOLS = [("uccsd", "UCCSD\nunfilt."),
         ("abelian_subset", "Abelian\nsubset"),
         ("abelian_projector", "Abelian\nproj."),
         ("full", "full grp\nproj."),
         ("hiuccsd_nonadapted", "HiUCCSD\ntheirs"),
         ("hiuccsd_adapted", "HiUCCSD\nadapted")]
MOLS = [("nh3", r"$\mathrm{NH_3}$ / $C_{3v}$", X.OKABE_ITO[5]),
        ("ch4", r"$\mathrm{CH_4}$ / $T_d$", X.OKABE_ITO[6])]


def main():
    with open(ROOT / "results" / "aggregate.csv") as fh:
        rows = [r for r in csv.DictReader(fh) if r["currency"] == "distinct"]

    def get(key, pool, field):
        for r in rows:
            if r["key"] == key and r["pool"] == pool:
                return r[field]
        raise KeyError((key, pool, field))

    plt = X.use_style()
    # No explicit hspace: constrained_layout is on and an hspace here fights it, which
    # is what produced a half-inch of dead band between the panels on the first draw.
    fig, (ax, bx) = plt.subplots(2, 1, figsize=X.FIGSIZE_WIDE,
                                 gridspec_kw=dict(height_ratios=[1.3, 1]))

    # --- (a) CNOT cost per pool, with the published values marked -------------
    w = 0.38
    xs = np.arange(len(POOLS))
    for k, (key, label, colour) in enumerate(MOLS):
        vals = [int(get(key, p, "cnots")) for p, _ in POOLS]
        bars = ax.bar(xs + (k - 0.5) * w, vals, w, label=label, color=colour,
                      edgecolor="black", linewidth=0.4,
                      hatch="///" if key == "ch4" else None, zorder=3)
        for b, v in zip(bars, vals):
            # Cleared to 1.6x: at 1.15x the rotated digits ran straight through the
            # published-value marker, which sits ON the bar top by design.
            ax.text(b.get_x() + b.get_width() / 2, v * 1.6, str(v), ha="center",
                    va="bottom", fontsize=7, rotation=90)
        # Published values, where Figure 5 reports one.  Drawn ON the bar top: a
        # marker beside the bar would read as a separate quantity rather than as the
        # same number measured twice.
        px = [xs[i] + (k - 0.5) * w for i, (p, _) in enumerate(POOLS)
              if get(key, p, "published_cnots") != ""]
        pv = [int(get(key, p, "published_cnots")) for p, _ in POOLS
              if get(key, p, "published_cnots") != ""]
        ax.plot(px, pv, "x", color="#000000", markersize=4.5, markeredgewidth=0.9,
                zorder=5, label="He et al., Fig. 5" if k == 0 else None)

    ax.set_yscale("log")
    ax.set_ylim(900, 30000)
    ax.set_ylabel("CNOTs (QEB scheme)")
    ax.set_xticks(xs)
    ax.set_xticklabels([lab for _, lab in POOLS], fontsize=7)
    ax.set_xlim(-0.6, len(POOLS) - 0.4)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.01), frameon=False, fontsize=7,
              ncol=3, handlelength=1.2, columnspacing=1.0, handletextpad=0.4,
              borderaxespad=0.0)

    # --- (b) what the full group removes, relative to the Abelian filter ------
    quantities = [("n_params", "parameters", X.OKABE_ITO[1], None),
                  ("n_gates", "operators", X.OKABE_ITO[2], None),
                  ("cnots", "CNOTs", X.OKABE_ITO[3], "///")]
    w2 = 0.22
    xs2 = np.arange(len(MOLS)) * 0.75
    for k, (field, label, colour, hatch) in enumerate(quantities):
        vals = []
        for key, _, _ in MOLS:
            base = int(get(key, "abelian_subset", field))
            full = int(get(key, "full", field))
            vals.append(100.0 * (1.0 - full / base))
        bars = bx.bar(xs2 + (k - 1) * w2, vals, w2, label=label, color=colour,
                      edgecolor="black", linewidth=0.4, hatch=hatch, zorder=3)
        for b, v in zip(bars, vals):
            bx.text(b.get_x() + b.get_width() / 2, v + 1.5, f"{v:.0f}", ha="center",
                    va="bottom", fontsize=7)

    # Ceiling at 108 rather than 82: the legend sits inside the axes and at 82 it
    # overlapped the 68 label on the CH4 parameter bar.
    bx.set_ylim(0, 108)
    bx.set_yticks([0, 25, 50, 75])
    bx.set_ylabel("removed vs the\nAbelian filter (%)")
    bx.set_xticks(xs2)
    bx.set_xticklabels([lab for _, lab, _ in MOLS])
    bx.set_xlim(-0.45, xs2[-1] + 0.45)
    bx.legend(loc="upper center", frameon=False, fontsize=7, ncol=3, handlelength=1.2,
              columnspacing=0.9, handletextpad=0.4, borderaxespad=0.2)

    fig.savefig(ROOT / "figures" / f"{EID}_cnot_cost.pdf")
    fig.savefig(ROOT / "figures" / f"{EID}_cnot_cost.png", format="png")
    plt.close(fig)

    manifest = {
        "cnot_cost": {
            "panels": [
                "upper (no letter drawn): CNOT cost of each of the six pools under the "
                "QEB ruler, log y axis, one bar pair per pool. Absolute counts, not ratios",
                "lower (no letter drawn): the percentage of each quantity that the full "
                "point group removes relative to the Abelian filter. Its x axis is the two "
                "molecules and is NOT shared with the upper panel, whose x axis is the pool",
            ],
            "elements": [
                {"artist": "left bar of every pair, both panels",
                 "encodes": "NH3 in C3v",
                 "visual": "blue #0072B2, 0.4 pt black edge, no hatch", "panel": "all"},
                {"artist": "right bar of every pair, upper panel",
                 "encodes": "CH4 in Td",
                 "visual": "vermillion #D55E00 with a /// hatch, which separates it from "
                           "the blue in grayscale", "panel": "upper"},
                {"artist": "black x marker on a bar top, upper panel",
                 "encodes": "the value He et al. report in Figure 5 for that same pool, "
                            "drawn on the bar rather than beside it because it is the same "
                            "quantity measured twice, not a second quantity",
                 "visual": "black x, 4.5 pt, drawn above the bars",
                 "panel": "upper"},
                {"artist": "bars with no x marker, upper panel",
                 "encodes": "the three pools Figure 5 does not report: the projector-built "
                            "Abelian pool, the full-group pool, and HiUCCSD in the adapted "
                            "basis. The full-group bar is the quantity this experiment adds",
                 "visual": "absence of the marker", "panel": "upper"},
                {"artist": "numeric label above every bar, upper panel",
                 "encodes": "the exact CNOT count",
                 "visual": "black text, 7 pt, rotated 90 degrees, cleared to 1.6x the bar "
                           "height so it does not run through the published-value marker",
                 "panel": "upper"},
                {"artist": "bar group 1 of 3 at each tick, lower panel",
                 "encodes": "percentage of variational parameters removed",
                 "visual": "orange #E69F00, no hatch", "panel": "lower"},
                {"artist": "bar group 2 of 3 at each tick, lower panel",
                 "encodes": "percentage of distinct elementary excitation operators removed",
                 "visual": "sky blue #56B4E9, no hatch", "panel": "lower"},
                {"artist": "bar group 3 of 3 at each tick, lower panel",
                 "encodes": "percentage of CNOTs removed. Its near-equality with group 2 "
                            "and its distance from group 1 is the finding of the experiment",
                 "visual": "bluish green #009E73 with a /// hatch", "panel": "lower"},
                {"artist": "numeric label above every bar, lower panel",
                 "encodes": "the percentage, rounded to the unit",
                 "visual": "black text, 7 pt, horizontal", "panel": "lower"},
                {"artist": "legend, upper panel",
                 "encodes": "the two molecules and the published-value marker",
                 "visual": "above the axes, three columns, no frame, 7 pt", "panel": "upper"},
                {"artist": "legend, lower panel", "encodes": "the three quantities",
                 "visual": "upper right inside the axes, three columns, no frame, 7 pt",
                 "panel": "lower"},
            ],
            "axes": {"x": "upper: the pool; lower: the molecule and its full point group. "
                          "The two x axes are different and are deliberately not shared",
                     "y": "upper: CNOT count under the QEB ruler; lower: percentage removed",
                     "xscale": "linear in both", "yscale": "upper log, lower linear"},
            "params": {"basis_set": "STO-3G",
                       "geometry": "He et al.'s, so that their published bars are "
                                   "reproducible; the counts do not depend on the angle",
                       "ruler": "QEB, 2 CNOTs per spin-orbital single and 13 per double",
                       "currency": "distinct elementary spin-orbital excitation operators, "
                                   "which is the currency Figure 5 is written in. The "
                                   "incidence currency is in aggregate.csv and is NOT plotted",
                       "no_transpilation": "no device gate set, no connectivity, no depth: "
                                           "this is circuit SIZE in the QEB scheme only"},
            "normalization": "upper panel none, absolute counts; lower panel percentage of "
                             "the Abelian-subset value",
        }
    }
    (ROOT / "figures" / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    viz = {"figures": [
        {"title": "CNOT cost by pool, distinct-operator currency",
         "plot_type": "bar", "x_axis": "pool", "y_axis": "cnots", "z_axis": None,
         "series": "molecule", "filters": {"currency": "distinct"},
         "log_scale_x": False, "log_scale_y": True},
        {"title": "CNOT cost by pool, incidence currency",
         "plot_type": "bar", "x_axis": "pool", "y_axis": "cnots", "z_axis": None,
         "series": "molecule", "filters": {"currency": "incidence"},
         "log_scale_x": False, "log_scale_y": True},
        {"title": "Variational parameters by pool",
         "plot_type": "bar", "x_axis": "pool", "y_axis": "n_params", "z_axis": None,
         "series": "molecule", "filters": {"currency": "distinct"},
         "log_scale_x": False, "log_scale_y": True},
        {"title": "CNOT cost relative to the Abelian filter",
         "plot_type": "bar", "x_axis": "pool", "y_axis": "cnot_ratio_vs_abelian_subset",
         "z_axis": None, "series": "molecule", "filters": {"currency": "distinct"},
         "log_scale_x": False, "log_scale_y": False},
        {"title": "Doubles against singles, the two classes the ruler prices differently",
         "plot_type": "scatter", "x_axis": "n_singles", "y_axis": "n_doubles",
         "z_axis": None, "series": "pool", "filters": {"currency": "distinct"},
         "log_scale_x": False, "log_scale_y": False},
    ]}
    (ROOT / "results" / "viz_schema.json").write_text(json.dumps(viz, indent=2) + "\n")
    print("figures:", sorted(p.name for p in (ROOT / "figures").glob("*")))


if __name__ == "__main__":
    main()
