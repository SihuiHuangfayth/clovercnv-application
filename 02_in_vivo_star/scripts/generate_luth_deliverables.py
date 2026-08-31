#!/usr/bin/env python
"""
generate_luth_deliverables_v2.py

Fixed version of Luth correspondence figure generation:
  - Figure 1: Widened + increased bar spacing to prevent x-label collision
  - Figure 2: Rebalanced layout (narrower gene labels, wider heatmap strip,
              colorbar closer)
  - Table 1: Widened keep_suggested column so header doesn't truncate

Run:
    conda activate CLOverCNV
    python generate_luth_deliverables_v2.py
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ============================================================
# CONFIG (same paths)
# ============================================================

RESULTS_DIR = "/scratch/gaz6bs/clovercnv/analysis"
LUTH_FILE = "/scratch/gaz6bs/clovercnv-application/02_in_vivo_star/luth_data_s4_filtered.tsv"
PER_CNV_FILE = f"{RESULTS_DIR}/luth_correspondence_per_cnv.tsv"
BY_COMPOUND_FILE = f"{RESULTS_DIR}/luth_correspondence_by_compound.tsv"

OUT_DIR = f"{RESULTS_DIR}/deliverables_luth_crossref"
os.makedirs(OUT_DIR, exist_ok=True)

COLORS = {
    "ratio_detected": "#2E7D32",
    "ratio_missed": "#C62828",
    "kept": "#1565C0",
    "not_kept": "#E0E0E0",
    "MMV006901": "#7B1FA2",
    "Suloctidil": "#EF6C00",
    "MMV019662": "#00838F",
    "DSM267": "#455A64",
    "axis": "#212121",
    "text": "#212121",
}

GAIN_CMAP = LinearSegmentedColormap.from_list(
    "gain_cmap", ["#FFFFFF", "#C8E6C9", "#66BB6A", "#2E7D32", "#1B5E20"]
)


def guler_style_axes(ax, fontsize=10, tick_fontsize=None, remove_top_right=True):
    if tick_fontsize is None:
        tick_fontsize = fontsize - 1
    ax.set_facecolor("white")
    for spine_name, spine in ax.spines.items():
        if remove_top_right and spine_name in ("top", "right"):
            spine.set_visible(False)
        else:
            spine.set_linewidth(1.5)
            spine.set_color(COLORS["axis"])
    ax.tick_params(axis="both", which="major", labelsize=tick_fontsize,
                   width=1.2, length=4, color=COLORS["axis"])
    ax.grid(False)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_color(COLORS["text"])


def savefig_both(fig, basename):
    svg_path = f"{OUT_DIR}/{basename}.svg"
    png_path = f"{OUT_DIR}/{basename}.png"
    fig.savefig(svg_path, format="svg", bbox_inches="tight", facecolor="white")
    fig.savefig(png_path, format="png", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {basename}.svg + {basename}.png")


def load_data():
    luth = pd.read_csv(LUTH_FILE, sep="\t")
    per_cnv = pd.read_csv(PER_CNV_FILE, sep="\t")
    by_compound = pd.read_csv(BY_COMPOUND_FILE, sep="\t")
    print(f"Loaded {len(per_cnv)} per-CNV rows across {per_cnv['compound'].nunique()} compounds.")
    return luth, per_cnv, by_compound


# ============================================================
# FIGURE 1: Bar chart — WIDENED with more spacing
# ============================================================

def figure_1(by_compound, style):
    compounds_ordered = ["Suloctidil", "MMV006901", "MMV019662", "DSM267"]
    df = by_compound.set_index("compound").reindex(compounds_ordered).reset_index()
    df = df.dropna(subset=["n_luth_cnvs_checked"])  # drop compounds absent from data

    ratio_pct = df["ratio_recall_rate"] * 100
    keep_pct = df["keep_suggested_recall_rate"] * 100
    n_luth = df["n_luth_cnvs_checked"]

    if style == "detailed":
        fig, ax = plt.subplots(figsize=(9, 5.5))       # was 7,5 — widened
        base_fontsize = 10
        title_fontsize = 12
        label_fontsize = 10
        annot_fontsize = 9
    else:  # slide
        fig, ax = plt.subplots(figsize=(12, 6.5))     # was 9,6 — much wider
        base_fontsize = 14
        title_fontsize = 16
        label_fontsize = 14
        annot_fontsize = 12

    x = np.arange(len(df)) * 1.4                      # was default 1.0 — more spacing
    width = 0.45

    bars1 = ax.bar(x - width / 2, ratio_pct, width,
                   label="Ratio-level detection\n(coverage signal)",
                   color=COLORS["ratio_detected"], edgecolor=COLORS["axis"], linewidth=1)
    bars2 = ax.bar(x + width / 2, keep_pct, width,
                   label="keep_suggested = TRUE\n(default confidence pass)",
                   color=COLORS["kept"], edgecolor=COLORS["axis"], linewidth=1)

    # Value annotations — stacked vertically to save horizontal space
    for i, (r, k, n) in enumerate(zip(ratio_pct, keep_pct, n_luth)):
        r_val = int(df.iloc[i]["n_ratio_detected"])
        k_val = int(df.iloc[i]["n_keep_suggested_any"])
        # Green bar annotation
        ax.text(x[i] - width / 2, r + 3,
                f"{r_val}/{int(n)}",
                ha="center", va="bottom", fontsize=annot_fontsize,
                fontweight="bold", color=COLORS["ratio_detected"])
        ax.text(x[i] - width / 2, r + 8,
                f"({r:.0f}%)",
                ha="center", va="bottom", fontsize=annot_fontsize - 1,
                color=COLORS["ratio_detected"])
        # Blue bar annotation
        ax.text(x[i] + width / 2, k + 3,
                f"{k_val}/{int(n)}",
                ha="center", va="bottom", fontsize=annot_fontsize,
                fontweight="bold", color=COLORS["kept"])
        ax.text(x[i] + width / 2, k + 8,
                f"({k:.0f}%)",
                ha="center", va="bottom", fontsize=annot_fontsize - 1,
                color=COLORS["kept"])

    # X labels — much clearer separation now
    xtick_labels = [f"{c}\n(n={int(n)})" for c, n in zip(df["compound"], n_luth)]
    ax.set_xticks(x)
    ax.set_xticklabels(xtick_labels, fontsize=label_fontsize,
                        fontweight="bold" if style == "slide" else "normal")

    ax.set_ylabel("CLOverCNV Recall of Luth 2024 CNVs (%)",
                  fontsize=label_fontsize, fontweight="bold")
    ax.set_ylim(0, 122)
    ax.set_yticks([0, 25, 50, 75, 100])

    title_text = ("Figure 1. CLOverCNV recovers 88% of Luth 2024 CNVs at coverage level;\n"
                  "default confidence thresholds filter many real signals"
                  if style == "detailed" else
                  "CLOverCNV recovers 88% of Luth CNVs at coverage level")
    ax.set_title(title_text, fontsize=title_fontsize, fontweight="bold", pad=15)

    leg = ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
                    fontsize=base_fontsize, frameon=True, framealpha=1,
                    edgecolor=COLORS["axis"], title="Detection level",
                    title_fontsize=base_fontsize)
    leg.get_title().set_fontweight("bold")

    ax.axhline(y=100, color=COLORS["axis"], linestyle=":", linewidth=1, alpha=0.5)

    guler_style_axes(ax, fontsize=base_fontsize)

    plt.tight_layout()
    savefig_both(fig, f"figure_1_recall_by_compound_{style}")
    plt.close()


# ============================================================
# FIGURE 2: Heatmap — REBALANCED LAYOUT
# ============================================================

def figure_2(luth, per_cnv, style):
    per_cnv = per_cnv.copy()
    compound_order = {"Suloctidil": 0, "MMV006901": 1, "MMV019662": 2, "DSM267": 3}
    per_cnv["comp_order"] = per_cnv["compound"].map(compound_order)
    per_cnv = per_cnv.sort_values(["comp_order", "chrom", "cnv_start"]).reset_index(drop=True)

    ratios = per_cnv["max_ratio_in_region"].fillna(0).values
    kept = per_cnv["keep_suggested_any"].values
    n_rows = len(per_cnv)

    if style == "detailed":
        # Compact layout for close-up review
        fig = plt.figure(figsize=(8, max(9, n_rows * 0.24)))
        # Use gridspec for tight control over ratios
        # 3 columns: labels | heatmap strip | ratio annotations
        gs = fig.add_gridspec(1, 4, width_ratios=[3.5, 0.6, 0.8, 0.4], wspace=0.05)
        base_fontsize = 8
        title_fontsize = 12
        label_fontsize = 9
        annot_fontsize = 8
        marker_size = 11
    else:
        fig = plt.figure(figsize=(10, max(10, n_rows * 0.28)))
        gs = fig.add_gridspec(1, 4, width_ratios=[3.5, 0.7, 0.9, 0.4], wspace=0.05)
        base_fontsize = 10
        title_fontsize = 14
        label_fontsize = 11
        annot_fontsize = 10
        marker_size = 14

    # Panel 1: Y-axis labels
    ax_labels = fig.add_subplot(gs[0, 0])
    ax_labels.set_xlim(0, 1)
    ax_labels.set_ylim(n_rows - 0.5, -0.5)
    ax_labels.axis("off")

    ylabels = []
    for _, row in per_cnv.iterrows():
        chrom_short = row["chrom"].replace("Pf3D7_", "").replace("_v3", "")
        gene = str(row["first_gene"]).replace("PF3D7_", "PF")
        ylabels.append((row["compound"], f"chr{chrom_short}", f"{int(row['cnv_start'])/1e6:.2f}M", gene))

    # Compound-color mapping for label prefix
    compound_colors = {"Suloctidil": COLORS["Suloctidil"], "MMV006901": COLORS["MMV006901"],
                       "MMV019662": COLORS["MMV019662"], "DSM267": COLORS["DSM267"]}
    for i, (comp, chrom, pos, gene) in enumerate(ylabels):
        ax_labels.text(0.02, i, comp, fontsize=base_fontsize,
                       fontweight="bold", color=compound_colors.get(comp, COLORS["text"]),
                       va="center", ha="left", family="sans-serif")
        ax_labels.text(0.36, i, f"{chrom}:{pos}", fontsize=base_fontsize,
                       va="center", ha="left", family="monospace", color=COLORS["text"])
        ax_labels.text(0.72, i, gene, fontsize=base_fontsize,
                       va="center", ha="left", family="monospace", color=COLORS["text"])

    # Panel 2: Heatmap strip
    ax_heat = fig.add_subplot(gs[0, 1])
    matrix = ratios.reshape(-1, 1)
    im = ax_heat.imshow(matrix, aspect="auto", cmap=GAIN_CMAP,
                        vmin=0.8, vmax=3.0)
    ax_heat.set_xticks([])
    ax_heat.set_yticks([])
    for spine in ax_heat.spines.values():
        spine.set_linewidth(1.2)

    # Overlay markers on heatmap
    for i, (r, k) in enumerate(zip(ratios, kept)):
        if r < 1.5:
            ax_heat.plot(0, i, "x", color=COLORS["ratio_missed"], markersize=marker_size,
                         markeredgewidth=2.5, zorder=5)
        if k:
            ax_heat.plot(0, i, "o", markerfacecolor="none",
                         markeredgecolor=COLORS["kept"],
                         markersize=marker_size + 3, markeredgewidth=2, zorder=4)

    # Compound-transition lines
    prev_compound = None
    for i, row in enumerate(per_cnv.itertuples()):
        if prev_compound is not None and row.compound != prev_compound:
            ax_heat.axhline(y=i - 0.5, color=COLORS["axis"], linewidth=1.5, zorder=3)
            # Also draw in labels panel
            ax_labels.axhline(y=i - 0.5, color=COLORS["axis"], linewidth=1.5, zorder=3,
                              xmin=0, xmax=1)
        prev_compound = row.compound

    # Panel 3: Ratio annotations
    ax_annot = fig.add_subplot(gs[0, 2])
    ax_annot.set_xlim(0, 1)
    ax_annot.set_ylim(n_rows - 0.5, -0.5)
    ax_annot.axis("off")
    for i, r in enumerate(ratios):
        if r > 0:
            color = COLORS["text"] if r >= 1.5 else COLORS["ratio_missed"]
            ax_annot.text(0.1, i, f"{r:.2f}×", va="center", ha="left",
                          fontsize=annot_fontsize, color=color, family="monospace",
                          fontweight="bold" if r >= 1.5 else "normal")
        else:
            ax_annot.text(0.1, i, "no overlap", va="center", ha="left",
                          fontsize=annot_fontsize - 1, color=COLORS["ratio_missed"],
                          family="monospace", style="italic")

    # Panel 4: Colorbar
    ax_cbar = fig.add_subplot(gs[0, 3])
    cbar = plt.colorbar(im, cax=ax_cbar, aspect=30)
    cbar.set_label("max_ratio\n(1.0 = baseline)",
                   fontsize=base_fontsize, fontweight="bold")
    cbar.ax.tick_params(labelsize=base_fontsize - 1)

    # Title spanning full figure
    title_text = ("Figure 2. Per-CNV correspondence: CLOverCNV signal at each Luth 2024 locus\n"
                  "○ = keep_suggested=TRUE  |  × = ratio below 1.5 (missed)"
                  if style == "detailed" else
                  "Per-CNV correspondence with Luth 2024\n"
                  "○ = keep_suggested passes  |  × = missed at ratio level")
    fig.suptitle(title_text, fontsize=title_fontsize, fontweight="bold", y=0.98)

    plt.subplots_adjust(top=0.93)
    savefig_both(fig, f"figure_2_percnv_heatmap_{style}")
    plt.close()


# ============================================================
# TABLE 1: WIDENED so header doesn't truncate
# ============================================================

def table_1(by_compound):
    compounds_ordered = ["Suloctidil", "MMV006901", "MMV019662", "DSM267"]
    df = by_compound.set_index("compound").reindex(compounds_ordered).reset_index()
    df = df.dropna(subset=["n_luth_cnvs_checked"])

    df.to_csv(f"{OUT_DIR}/table_1_per_compound_summary.tsv", sep="\t", index=False)
    df.to_csv(f"{OUT_DIR}/table_1_per_compound_summary.csv", index=False)

    display_df = pd.DataFrame({
        "Compound": df["compound"],
        "N Luth CNVs": df["n_luth_cnvs_checked"].astype(int),
        "Ratio recall": [f"{int(n)}/{int(t)} ({r*100:.0f}%)"
                         for n, t, r in zip(df["n_ratio_detected"],
                                             df["n_luth_cnvs_checked"],
                                             df["ratio_recall_rate"])],
        "keep_suggested recall": [f"{int(n)}/{int(t)} ({r*100:.0f}%)"
                                   for n, t, r in zip(df["n_keep_suggested_any"],
                                                       df["n_luth_cnvs_checked"],
                                                       df["keep_suggested_recall_rate"])],
        "Median max_ratio": df["median_max_ratio_in_luth_region"].round(2),
    })

    total_n = df["n_luth_cnvs_checked"].sum()
    total_ratio = df["n_ratio_detected"].sum()
    total_kept = df["n_keep_suggested_any"].sum()
    total_row = pd.DataFrame({
        "Compound": ["TOTAL"],
        "N Luth CNVs": [int(total_n)],
        "Ratio recall": [f"{int(total_ratio)}/{int(total_n)} ({total_ratio/total_n*100:.1f}%)"],
        "keep_suggested recall": [f"{int(total_kept)}/{int(total_n)} ({total_kept/total_n*100:.1f}%)"],
        "Median max_ratio": [df["median_max_ratio_in_luth_region"].median().round(2)],
    })
    display_df = pd.concat([display_df, total_row], ignore_index=True)

    for style in ["detailed", "slide"]:
        if style == "detailed":
            fig_w, fig_h = 11, 3.2                    # was 9,3.2 — widened for header fit
            fontsize = 10
            header_fontsize = 11
        else:
            fig_w, fig_h = 14, 3.8                    # was 11,3.5 — widened
            fontsize = 13
            header_fontsize = 14

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax.axis("off")

        cell_text = display_df.values.tolist()
        col_labels = display_df.columns.tolist()

        # Explicit column widths (make keep_suggested column widest)
        col_widths = [0.13, 0.10, 0.20, 0.28, 0.15]   # 5 cols, sum ~0.86

        table = ax.table(cellText=cell_text,
                         colLabels=col_labels,
                         loc="center",
                         cellLoc="center",
                         colColours=[COLORS["kept"]] * len(col_labels),
                         colWidths=col_widths)

        table.auto_set_font_size(False)
        table.set_fontsize(fontsize)
        table.scale(1, 2)

        for i in range(len(col_labels)):
            cell = table[(0, i)]
            cell.set_text_props(color="white", fontweight="bold", fontsize=header_fontsize)
            cell.set_height(0.14)
            cell.set_edgecolor(COLORS["axis"])
            cell.set_linewidth(1.5)

        for r in range(1, len(cell_text) + 1):
            for c in range(len(col_labels)):
                cell = table[(r, c)]
                cell.set_edgecolor(COLORS["axis"])
                cell.set_linewidth(0.8)
                cell.set_height(0.11)
                if r == len(cell_text):
                    cell.set_text_props(fontweight="bold", fontsize=fontsize + 1)
                    cell.set_facecolor("#FFF9C4")
                elif r % 2 == 0:
                    cell.set_facecolor("#F5F5F5")

        title_text = ("Table 1. CLOverCNV recall of Luth 2024 CNVs per compound cohort\n"
                      "(41 CNVs total, excluding 1 for excluded sample SRR5161093)"
                      if style == "detailed" else
                      "Table 1. Per-compound recall of Luth 2024 CNVs")
        ax.set_title(title_text, fontsize=header_fontsize + 1,
                     fontweight="bold", pad=15)

        plt.tight_layout()
        savefig_both(fig, f"table_1_per_compound_summary_{style}")
        plt.close()


def main():
    luth, per_cnv, by_compound = load_data()
    print(f"\nOutput directory: {OUT_DIR}\n")

    print("Generating Figure 1 (recall bar chart) - v2...")
    figure_1(by_compound, "detailed")
    figure_1(by_compound, "slide")

    print("\nGenerating Figure 2 (per-locus heatmap) - v2...")
    figure_2(luth, per_cnv, "detailed")
    figure_2(luth, per_cnv, "slide")

    print("\nGenerating Table 1 (per-compound summary) - v2...")
    table_1(by_compound)

    print(f"\n{'='*60}")
    print(f"All deliverables re-saved to {OUT_DIR}")
    print(f"(Previous files were overwritten)")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
