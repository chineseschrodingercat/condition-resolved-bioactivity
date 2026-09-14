"""Five source-traceable, data-only Communications Chemistry candidate figures.

All panels are re-plotted from preserved project CSVs. Published observations
remain labelled as such; these figures do not imply new wet-lab experiments.
"""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
V4 = "--v4" in sys.argv
INK = "#243445"
MUTED = "#82929a"
TEAL = "#007e87"
ORANGE = "#bd6d2d"
PURPLE = "#765e99"
GREEN = "#5b7e64"
PAIR_IDS = ["JAK3_2019_ATP", "JAK3_2017_ATP", "renin_plasma_trypsin", "KDR_HTRF_cell"]
PAIR_NAMES = ["JAK3 '19\nn=20", "JAK3 '17\nn=6", "Renin\nn=27", "KDR\nn=16"]
METHODS = ["average_pIC50", "discard_second_assay", "context_global_offset", "context_structure_interaction"]
METHOD_NAMES = ["Average", "Keep first", "Global offset", "Structure × condition"]
METHOD_COLORS = [MUTED, ORANGE, TEAL, PURPLE]


def csv(name: str) -> pd.DataFrame:
    return pd.read_csv(OUT / name)


def style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11.1,
        "axes.labelsize": 11.2, "axes.labelweight": "medium",
        "xtick.labelsize": 10.5, "ytick.labelsize": 10.5,
        "axes.labelcolor": INK, "text.color": INK,
        "axes.spines.top": True, "axes.spines.right": True,
        "axes.linewidth": 1.55, "axes.edgecolor": "#25313d", "xtick.color": INK,
        "ytick.color": INK, "xtick.major.width": 1.25, "ytick.major.width": 1.25,
        "xtick.major.size": 5.0, "ytick.major.size": 5.0,
        "lines.linewidth": 2.35, "lines.markersize": 7.5,
        "savefig.facecolor": "white",
    })
    if V4:
        plt.rcParams.update({
            "font.size": 9.6, "axes.labelsize": 9.6,
            "xtick.labelsize": 8.8, "ytick.labelsize": 8.8,
            "axes.linewidth": 1.55, "lines.linewidth": 2.15,
            "lines.markersize": 6.8,
        })


def mark(ax, letter: str, title: str, *, grid: bool = True) -> None:
    """Use a separate prominent uppercase panel index and concise data title."""
    ax.text(0, 1.035, letter.upper(), transform=ax.transAxes, fontsize=15 if V4 else 17.5,
            fontweight="black", ha="left", va="bottom", clip_on=False)
    title_x = (.15 if ax.get_position().width * ax.figure.get_figwidth() < 5 else .075) if V4 else .095
    ax.text(title_x, 1.044, title, transform=ax.transAxes, fontsize=9.8 if V4 else 11.3,
            fontweight="bold", ha="left", va="bottom", clip_on=False)
    if grid and not V4:
        ax.grid(axis="y", color="#d7dfe2", lw=.8, alpha=.8)
    ax.set_axisbelow(True)


def save(fig, stem: str) -> None:
    panel_inches = [(round(ax.get_position().width * fig.get_figwidth(), 2),
                     round(ax.get_position().height * fig.get_figheight(), 2))
                    for ax in fig.axes]
    version = "v4" if V4 else "v3"
    fig.savefig(OUT / f"commchem_{stem}_{version}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"commchem_{stem}_{version}.pdf", bbox_inches="tight")
    print(f"{stem}: panel plotting rectangles (in) {panel_inches}")
    plt.close(fig)


def fig1() -> None:
    sources = [
        csv("jak3_atp_current_20_pairs.csv"),
        csv("jak3_2017_4uM_vs_1mM_matched_pairs.csv"),
        csv("renin_plasma_vs_trypsin_matched_pairs.csv"),
        csv("kdr_htrf_vs_cell_matched_pairs.csv"),
    ]
    metrics = csv("phase1_four_arm_metrics.csv")
    transfer = csv("phase1_source_transfer_metrics.csv")
    if V4:
        fig, axs = plt.subplots(2, 2, figsize=(9.5, 6.6))
        axs = axs.ravel()
        fig.subplots_adjust(left=.105, right=.98, top=.92, bottom=.13, hspace=.76, wspace=.32)
    else:
        fig, axs = plt.subplots(4, 1, figsize=(6.7, 12.5))
        fig.subplots_adjust(left=.16, right=.985, top=.97, bottom=.045, hspace=.55)
    for i, frame in enumerate(sources):
        d = frame.delta_log10_second_over_first.to_numpy(dtype=float)
        jitter = np.linspace(-.13, .13, len(d))
        axs[0].scatter(np.full(len(d), i) + jitter, d, color=[TEAL, ORANGE, GREEN, PURPLE][i],
                          s=38, alpha=.77)
        axs[0].plot([i-.23, i+.23], [np.median(d)]*2, color=INK, lw=2.8)
        axs[1].scatter(np.full(len(d), i)+jitter, d-np.median(d),
                          color=[TEAL, ORANGE, GREEN, PURPLE][i], s=38, alpha=.77)
    for ax in axs[:2]:
        ax.set_xticks(range(4), PAIR_NAMES)
        ax.axhline(0, color="#555f65", lw=1.35, ls="--")
        ax.set_ylabel("log10(second / first IC50)")
    mark(axs[0], "a", "Measured contrast and source shifts" if V4 else "Measured contrasts include source-level shifts")
    mark(axs[1], "b", "Residual molecular variation" if V4 else "Median-centred residuals retain molecular variation")
    ax = axs[2]
    x = np.arange(4)
    for j, method in enumerate(METHODS):
        vals = metrics.loc[(metrics.pair.isin(PAIR_IDS)) & (metrics.method.eq(method))]
        vals = vals.set_index("pair").loc[PAIR_IDS].mae_both_conditions_pIC50.to_numpy()
        ax.bar(x+(j-1.5)*.19, vals, width=.19, color=METHOD_COLORS[j], label=METHOD_NAMES[j])
    ax.set_xticks(x, PAIR_NAMES)
    ax.set_ylabel("Both-readout LOOM MAE (pIC50)")
    ax.legend(frameon=False, ncol=2, fontsize=7.4 if V4 else 9, loc="upper right")
    if V4:
        ax.set_ylim(0, 2.1)
    mark(ax, "c", "Same-molecule policy baselines" if V4 else "Same-molecule record-policy baselines")
    ax = axs[3]
    directions = list(transfer.train_source.drop_duplicates())
    method_order = ["zero_delta", "source_mean_delta", "source_structure_delta"]
    for j, (meth, color, name) in enumerate(zip(method_order, [MUTED, TEAL, PURPLE],
                                                 ["Zero contrast", "Source offset", "Structure × condition"])):
        vals = [float(transfer.loc[(transfer.train_source.eq(src)) & (transfer.method.eq(meth)),
                                   "mae_delta_log10"].iloc[0]) for src in directions]
        ax.bar(np.arange(2)+(j-1)*.24, vals, .24, color=color, label=name)
    ax.set_xticks([0, 1], ["2019 → 2017\n5 test compounds", "2017 → 2019\n19 test compounds"])
    ax.set_ylabel("New-source contrast MAE (log10)")
    ax.legend(frameon=False, fontsize=7.4 if V4 else 9)
    ax.set_ylim(0, 3.1 if V4 else 2.7)
    mark(ax, "d", "JAK3 cross-paper calibration" if V4 else "JAK3 paper transfer tests calibration")
    save(fig, "fig1_contrast_baselines")


def fig2() -> None:
    equal = csv("phase5_alk_equal_budget_10_condition_comparison.csv")
    srpk = csv("phase5_alk_srpk_9_coverage_equal_budget.csv")
    conditions = ["WT", "C1156Y", "F1174L", "L1196M", "L1152R", "1151Tins",
                  "G1202R", "G1269A", "S1206Y", "parental_BaF3"]
    strategies = ["contrast_matched_v0", "first_potency", "geometric_mean_potency"]
    names = ["Contrast", "WT potency", "Mean potency"]
    fig = plt.figure(figsize=(9.5, 6.7) if V4 else (6.7, 11.6))
    if V4:
        gs = fig.add_gridspec(2, 2, height_ratios=[1.06, 1], width_ratios=[1.16, 1],
                              left=.095, right=.98, top=.92, bottom=.13,
                              hspace=.79, wspace=.35)
        slots = [gs[0, :], gs[1, 0], gs[1, 1]]
    else:
        gs = fig.add_gridspec(3, 1, height_ratios=[1.5, 1, 1],
                              left=.16, right=.985, top=.97, bottom=.055, hspace=.63)
        slots = [gs[0], gs[1], gs[2]]
    ax = fig.add_subplot(slots[0])
    chosen = equal[equal.strategy.eq(strategies[0])].set_index("condition").loc[conditions]
    x = np.arange(len(conditions))
    for field, rel, label, color, marker in [
        ("a_ic50_nM", "a_relation", "Compound 17", TEAL, "o"),
        ("b_ic50_nM", "b_relation", "Compound 19", ORANGE, "s")]:
        y = chosen[field].to_numpy(float)
        ax.plot(x, y, color=color, lw=2.55, marker=marker, ms=8, label=label)
        for xx, yy, rr in zip(x, y, chosen[rel]):
            if rr != "=":
                ax.scatter(xx, yy, facecolors="white", edgecolors=color, s=95, marker="^", zorder=5)
    ax.axvspan(8.5, 9.5, color="#eceff0", zorder=-1)
    ax.set_yscale("log")
    ax.set_xticks(x, [s.replace("parental_BaF3", "Parental") for s in conditions],
                  rotation=32, ha="right")
    ax.set_ylabel("Ba/F3 72-h MTS IC50 (nM)")
    ax.legend(frameon=False, ncol=2)
    mark(ax, "a", "ALK mutant and parental-cell responses")
    ax = fig.add_subplot(slots[1])
    mat = np.full((3, 10), np.nan)
    for i, strat in enumerate(strategies):
        s = equal[equal.strategy.eq(strat)].set_index("condition")
        for j, cond in enumerate(conditions):
            value = s.loc[cond, "within_pair_abs_log10_gap"]
            if pd.notna(value):
                mat[i, j] = float(value)
    cmap = plt.get_cmap("YlOrBr").copy()
    cmap.set_bad("#e6e9e9")
    ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=2.1)
    for i in range(3):
        for j in range(10):
            ax.text(j, i, "bound" if np.isnan(mat[i,j]) else f"{10**mat[i,j]:.1f}×",
                    ha="center", va="center", fontsize=7.2 if V4 else 8.5,
                    color="white" if np.isfinite(mat[i,j]) and mat[i,j] > 1.3 else INK)
    ax.set_xticks(x, [s.replace("parental_BaF3", "Parental") for s in conditions], rotation=28, ha="right")
    ax.set_yticks(range(3), names)
    mark(ax, "b", "ALK state-specific pair gaps" if V4 else "Two-compound separation across ALK states", grid=False)
    ax = fig.add_subplot(slots[2])
    vals = []
    for st in strategies:
        row = srpk[srpk.strategy.eq(st)].iloc[0]
        gap = row.srpk1_pair_fold_gap_lower_bound_if_valid
        if pd.notna(gap):
            vals.append((float(gap), True))
        else:
            vals.append((10**float(row.srpk1_pair_abs_log10_gap_if_exact), False))
    ax.barh(np.arange(3), [v[0] for v in vals], color=[TEAL, MUTED, MUTED])
    ax.set_yticks(range(3), names)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(.8, 1500)
    ax.set_xlabel("SRPK1 within-pair separation (fold)")
    for i, (v, bounded) in enumerate(vals):
        ax.text(v*1.08, i, f">{v:.0f}×" if bounded else f"{v:.1f}×", va="center", fontsize=8.5 if V4 else 9.5)
    mark(ax, "c", "SRPK1 cross-target gap" if V4 else "SRPK1 cross-target separation")
    save(fig, "fig2_alk_spectra")


def fig3() -> None:
    long = csv("phase5_idh1_all_orthogonal_endpoints.csv")
    reveal = csv("phase5_idh1_equal_budget_reveal.csv")
    pivot = long.pivot_table(index="molecule_id", columns="assay_axis", values="value", aggfunc="first")
    chosen = ["CHEMBL2180728", "CHEMBL2180746"]
    if V4:
        fig = plt.figure(figsize=(9.5, 6.7))
        gs = fig.add_gridspec(2, 2, height_ratios=[1.06, 1], width_ratios=[1.12, 1],
                              left=.095, right=.98, top=.92, bottom=.13,
                              hspace=.82, wspace=.35)
        axs = [fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    else:
        fig, axs = plt.subplots(3, 1, figsize=(6.7, 10.8),
                                gridspec_kw={"height_ratios": [1.12, 1.12, .9]})
        fig.subplots_adjust(left=.16, right=.985, top=.97, bottom=.065, hspace=.65)
    ax = axs[0]
    for mid, row in pivot.iterrows():
        x, y = row["R132C_enzyme"], row["HT1080_2HG"]
        selected = mid in chosen
        ax.scatter(x, y, s=120 if selected else 60, color=ORANGE if selected else MUTED,
                   edgecolors=INK if selected else "none", zorder=3)
        offset = (-28, -12) if mid.endswith("0727") else (4, 5)
        ax.annotate(mid[-4:], (x, y), xytext=offset, textcoords="offset points", fontsize=9)
    ax.set(xscale="log", yscale="log", xlabel="R132C enzyme IC50 (nM)",
           ylabel="HT1080 2-HG IC50 (nM)")
    mark(ax, "a", "Eight matched compounds")
    ax = axs[1]
    cols = ["R132C_enzyme", "HT1080_2HG", "R132H_enzyme", "U87_2HG"]
    x = np.arange(4)
    for mid, color, marker in zip(chosen, [TEAL, ORANGE], ["o", "s"]):
        vals = [pivot.loc[mid, col] for col in cols]
        ax.plot(x, vals, lw=2.65, marker=marker, ms=8, color=color, label=mid[-4:])
    ax.set_yscale("log")
    ax.set_xticks(x, ["R132C\nenzyme", "HT1080\n2-HG", "R132H\nenzyme", "U87\n2-HG"])
    ax.set_ylabel("IC50 (nM)")
    ax.legend(title="CHEMBL218…", frameon=False, fontsize=9)
    mark(ax, "b", "Pair across four endpoints" if V4 else "Selected pair across four endpoints")
    ax = axs[2]
    labels = {"contrast_matched_v0": "Contrast", "first_potency": "Enzyme",
              "geometric_mean_potency": "Mean"}
    vals = 10**reveal.U87_2HG_within_pair_abs_log10_gap.to_numpy(float)
    ax.barh(range(3), vals, color=[TEAL, MUTED, MUTED])
    ax.set_yticks(range(3), [labels[s] for s in reveal.strategy])
    ax.invert_yaxis()
    ax.set_xlim(0, 3.1)
    ax.set_xlabel("U87 2-HG pair gap (fold)")
    for i, v in enumerate(vals):
        ax.text(v+.05, i, f"{v:.2f}×", va="center", fontsize=8.5 if V4 else 9.5)
    mark(ax, "c", "Withheld R132H-cell readout")
    save(fig, "fig3_idh1_context")


def fig4() -> None:
    jak = csv("jak3_2017_ATP_PBMC_blood_complete_cases.csv")
    kdr = csv("kdr_2012_orthogonal_endpoints_joined.csv")
    if V4:
        fig, axs = plt.subplots(2, 2, figsize=(9.5, 6.6))
        axs = axs.ravel()
        fig.subplots_adjust(left=.105, right=.98, top=.92, bottom=.13, hspace=.76, wspace=.32)
    else:
        fig, axs = plt.subplots(4, 1, figsize=(6.7, 12.5))
        fig.subplots_adjust(left=.17, right=.985, top=.97, bottom=.05, hspace=.55)
    colors = {"CHEMBL4085457": TEAL, "CHEMBL4085582": ORANGE}
    ax = axs[0]
    for _, r in jak.iterrows():
        mid = r.molecule_chembl_id
        ax.scatter(r.ATP_ratio_high_over_low, r.whole_blood_over_PBMC,
                   s=125 if mid in colors else 75, color=colors.get(mid, MUTED),
                   edgecolors=INK if mid in colors else "none", label=mid if mid in colors else None)
    ax.set(xscale="log", yscale="log", xlabel="High / low ATP IC50 ratio",
           ylabel="Whole blood / PBMC IC50 ratio")
    ax.set_xlim(77, 124)
    ax.set_ylim(2.0 if V4 else 2.5, 100)
    for mid, label, offset in [("CHEMBL4085457", "3.86×", (-5, 6) if V4 else (-5, -16)),
                               ("CHEMBL4085582", "71.3×", (6, 2))]:
        r = jak[jak.molecule_chembl_id.eq(mid)].iloc[0]
        ax.annotate(label, (r.ATP_ratio_high_over_low, r.whole_blood_over_PBMC),
                    xytext=offset, textcoords="offset points", fontsize=9)
    ax.legend(frameon=False, fontsize=7.5 if V4 else 9,
              title="Highlighted ChEMBL IDs", title_fontsize=7.5 if V4 else 9)
    mark(ax, "a", "JAK3: two context axes" if V4 else "JAK3: two measured context axes")
    ax = axs[1]
    for mid, color, label in [("CHEMBL4085457", TEAL, "PF-06651600"),
                               ("CHEMBL4085582", ORANGE, "CHEMBL4085582")]:
        r = jak[jak.molecule_chembl_id.eq(mid)].iloc[0]
        ax.plot(range(4), [r.nM_first, r.nM_second, r.nM_PBMC, r.nM_whole_blood],
                marker="o", lw=2.65, ms=8, color=color, label=label)
    ax.set(yscale="log", xticks=range(4), xticklabels=["4 µM ATP", "1 mM ATP", "PBMC\n75 min", "Blood\n45 min"],
           ylabel="IC50 (nM)")
    ax.legend(frameon=False, fontsize=7.5 if V4 else 9)
    mark(ax, "b", "JAK3: four assay records" if V4 else "JAK3: two molecules, four assay records")
    focus = {"CHEMBL2071202": TEAL, "CHEMBL2071272": ORANGE, "CHEMBL2071273": PURPLE}
    ax = axs[2]
    xx = np.logspace(0, 2.5, 50)
    ax.plot(xx, xx, color=MUTED, ls="--", lw=1.35, label="equal potency")
    for _, r in kdr.iterrows():
        mid = r.molecule_chembl_id
        ax.scatter(r.kdr_biochem_ic50_nM, r.kdr_cell_ic50_nM,
                   s=110 if mid in focus else 60, color=focus.get(mid, MUTED),
                   edgecolors=INK if mid in focus else "none", zorder=3)
    ax.set(xscale="log", yscale="log", xlabel="KDR biochemical IC50 (nM)",
           ylabel="KDR cell IC50 (nM)")
    mark(ax, "c", "KDR: 16 matched molecules" if V4 else "KDR: 16 enzyme–cell matched molecules")
    ax = axs[3]
    for _, r in kdr.iterrows():
        if pd.isna(r.polyploidy_cell_value):
            continue
        mid = r.molecule_chembl_id
        ratio = r.kdr_cell_ic50_nM / r.kdr_biochem_ic50_nM
        bounded = r.polyploidy_cell_relation == "<"
        ax.scatter(ratio, r.polyploidy_cell_value,
                   s=115 if mid in focus else 65,
                   marker="v" if bounded else "o", color=focus.get(mid, MUTED),
                   edgecolors=INK if mid in focus else "none", zorder=3)
    ax.set(xscale="log", yscale="log", xlabel="KDR cell / biochemical IC50",
           ylabel="Polyploidy EC15 (nM or upper bound)")
    ax.text(.02, .97, "▼  reported < bound", transform=ax.transAxes, va="top", fontsize=8 if V4 else 9)
    mark(ax, "d", "Polyploidy-cell readout" if V4 else "Distinct polyploidy-cell readout")
    save(fig, "fig4_jak3_kdr")


def fig5() -> None:
    btk = csv("btk_2025_corrected_15_binding_pairs.csv")
    kinetic = csv("btk_crosspaper_8_kinetic_rank.csv")
    parp = csv("parp_2012_cell_catalytic_vs_trapping.csv")
    if V4:
        fig = plt.figure(figsize=(9.5, 6.7))
        gs = fig.add_gridspec(2, 2, height_ratios=[1.06, 1], width_ratios=[1.12, 1],
                              left=.095, right=.98, top=.92, bottom=.13,
                              hspace=.82, wspace=.35)
        axs = [fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    else:
        fig, axs = plt.subplots(3, 1, figsize=(6.7, 9.7))
        fig.subplots_adjust(left=.17, right=.985, top=.97, bottom=.07, hspace=.66)
    ax = axs[0]
    for mech, color in [("covalent", TEAL), ("noncovalent", ORANGE)]:
        s = btk[btk.mechanism.eq(mech)]
        for rel, marker in [("=", "o"), (">", "^")]:
            sub = s[s.C481S_Kd_relation.eq(rel)]
            if not sub.empty:
                ax.scatter(sub.WT_Kd_nM, sub.C481S_Kd_nM_boundary, s=85,
                           color=color, marker=marker, label=f"{mech} ({rel})")
    span = np.logspace(-1.5, 3.5, 100)
    ax.plot(span, span, color=MUTED, ls="--", lw=1.35)
    ax.set(xscale="log", yscale="log", xlabel="WT BTK Kd (nM)",
           ylabel="C481S Kd (nM or lower bound)")
    ax.set_xlim(.07, 40)
    ax.set_ylim(.05, 2300)
    ax.legend(frameon=False, fontsize=8 if V4 else 8.8, ncol=2, loc="lower right")
    mark(ax, "a", "Corrected BTK mutation contrasts")
    ax = axs[1]
    for _, r in kinetic.iterrows():
        bound = r.ratio_relation == ">"
        ax.scatter(r.C481S_over_WT_ratio_boundary, r.kinact_over_KI_rank_high_to_low,
                   marker="^" if bound else "o", color=TEAL, s=85)
        offset, align = ((-4, 0) if V4 else (-5, -11), "right") if r.drug == "Branebrutinib" else ((4, 3), "left")
        ax.annotate(r.drug, (r.C481S_over_WT_ratio_boundary, r.kinact_over_KI_rank_high_to_low),
                    xytext=offset, textcoords="offset points", fontsize=7.7 if V4 else 8.8, ha=align)
    ax.set(xscale="log", ylim=(8.6,.4), yticks=range(1,9),
           xlabel="C481S / WT Kd ratio (bounds marked ▲)",
           ylabel="Independent kinact/KI rank (1 = highest)")
    ax.set_xlim(7, 6500)
    mark(ax, "b", "Affinity shift vs kinetic rank" if V4 else "Affinity shift is not a kinetic rank")
    ax = axs[2]
    color_map = {"olaparib": TEAL, "veliparib": ORANGE, "niraparib": PURPLE}
    for _, r in parp.iterrows():
        yy = [r.catalytic_rank_1_strongest, r.published_trapping_rank_1_strongest]
        ax.plot([0,1], yy, marker="o", lw=2.65, ms=8, color=color_map[r.drug], label=r.drug)
    ax.set(xlim=(-.1,1.1), ylim=(3.3,.7), yticks=[1,2,3], xticks=[0,1],
           xticklabels=["Cellular PAR\ninhibition", "PARP trapping\n(source ordinal)"],
           ylabel="Rank (1 = strongest)")
    if V4:
        for _, r in parp.iterrows():
            ax.annotate(r.drug, (0, r.catalytic_rank_1_strongest),
                        xytext=(8, -4), textcoords="offset points", ha="left", fontsize=7.3,
                        color=color_map[r.drug],
                        bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85, "pad": .2})
    else:
        ax.legend(frameon=False, fontsize=9, loc="lower center",
                  bbox_to_anchor=(.5, .01), ncol=3)
    mark(ax, "c", "Published PARP rank inversion" if V4 else "Published PARP source-ordinal rank inversion")
    save(fig, "fig5_btk_parp")


if __name__ == "__main__":
    style()
    for make in [fig1, fig2, fig3, fig4, fig5]:
        make()
    print(f"Generated 5 data-only figures as commchem_fig[1-5]_..._{'v4' if V4 else 'v3'} PNG/PDF")
