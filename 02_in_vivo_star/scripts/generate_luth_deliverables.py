#!/usr/bin/env python
"""
generate_luth_deliverables.py

Generate deliverables for Jenny 1-on-1 review from Luth correspondence outputs:
  - Figure 1: Recall bar chart (ratio vs keep_suggested by compound)
  - Figure 2: Per-locus heatmap (Luth CNV × sample, cell = observed max_ratio)
  - Table 1: Per-compound summary CSV/TSV + rendered PNG for slides

Output for each figure: TWO versions
  - _detailed  : sit-down 1-on-1 review (dense labels, all annotations)
  - _slide     : presentation-mode (big fonts, minimal labels, clean)

Formats: SVG (for Illustrator/Inkscape post-edit) + PNG @ 300 DPI (drop into slides)

Guler style: no gridlines, bold axes, rotated x-labels, value annotations,
             legends in titled box outside plot, square aspect for curves.

Run:
    conda activate CLOverCNV
    python /scratch/gaz6bs/clovercnv-application/02_in_vivo_star/scripts/generate_luth_deliverables.py
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

# ============================================================
# CONFIG
# ============================================================

RESULTS_DIR = "/scratch/gaz6bs/clovercnv/analysis"
LUTH_FILE = "/scratch/gaz6bs/clovercnv-application/02_in_vivo_star/luth_data_s4_filtered.tsv"
PER_CNV_FILE = f"{RESULTS_DIR}/luth_correspondence_per_cnv.tsv"
BY_COMPOUND_FILE = f"{RESULTS_DIR}/luth_correspondence_by_compound.tsv"

OUT_DIR = f"{RESULTS_DIR}/deliverables_luth_crossref"
os.makedirs(OUT_DIR, exist_ok=True)

# Guler lab palette (professional academic)
COLORS = {
    "ratio_detected": "#2E7D32",       # deep green
    "ratio_missed": "#C62828",         # deep red
    "kept": "#1565C0",                 # deep blue
    "not_kept": "#E0E0E0",             # neutral grey
    "MMV006901": "#7B1FA2",            # purple
    "Suloctidil": "#EF6C00",           # orange
    "MMV019662": "#00838F",            # teal
    "DSM267": "#455A64",               # slate
    "grid": "#F5F5F5",
    "axis": "#212121",
    "text": "#212121",
}

# Heatmap colormap: white (baseline=1) -> deep green (gain), red (loss)
GAIN_CMAP = LinearSegmentedColormap.from_list(
    "gain_cmap", ["#FFFFFF", "#C8E6C9", "#66BB6A", "#2E7D32", "#1B5E20"]
)


# ============================================================
# STYLE HELPERS
# ============================================================

def guler_style_axes(ax, fontsize=10, tick_fontsize=None, remove_top_right=True):
    """Apply Guler lab style: bold axes, no gridlines, clean spines."""
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
    """Save both SVG and PNG."""
    svg_path = f"{OUT_DIR}/{basename}.svg"
    png_path = f"{OUT_DIR}/{basename}.png"
    fig.savefig(svg_path, format="svg", bbox_inches="tight", facecolor="white")
    fig.savefig(png_path, format="png", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {basename}.svg + {basename}.png")


# ============================================================
# LOAD DATA
# ============================================================

def load_data():
    luth = pd.read_csv(LUTH_FILE, sep="\t")
    per_cnv = pd.read_csv(PER_CNV_FILE, sep="\t")
    by_compound = pd.read_csv(BY_COMPOUND_FILE, sep="\t")
    print(f"Loaded {len(per_cnv)} per-CNV rows across {per_cnv['compound'].nunique()} compounds.")
    return luth, per_cnv, by_compound


# ============================================================
# FIGURE 1: Recall bar chart
# ============================================================

def figure_1(by_compound, style):
    """
    Grouped bar chart: per-compound recall at two levels.
    style: 'detailed' or 'slide'
    """
    compounds = by_compound["compound"].tolist()
    compounds_ordered = ["Suloctidil", "MMV006901", "MMV019662", "DSM267"]
    compounds_ordered = [c for c in compounds_ordered if c in compounds]
    df = by_compound.set_index("compound").loc[compounds_ordered].reset_index()

    ratio_pct = df["ratio_recall_rate"] * 100
    keep_pct = df["keep_suggested_recall_rate"] * 100
    n_luth = df["n_luth_cnvs_checked"]

    if style == "detailed":
        fig, ax = plt.subplots(figsize=(7, 5))
        base_fontsize = 10
        title_fontsize = 12
        label_fontsize = 10
        annot_fontsize = 8
        rotation = 0
    else:  # slide
        fig, ax = plt.subplots(figsize=(9, 6))
        base_fontsize = 14
        title_fontsize = 16
        label_fontsize = 14
        annot_fontsize = 11
        rotation = 0

    x = np.arange(len(compounds_ordered))
    width = 0.38

    bars1 = ax.bar(x - width / 2, ratio_pct, width,
                   label="Ratio-level detection\n(coverage signal)",
                   color=COLORS["ratio_detected"], edgecolor=COLORS["axis"], linewidth=1)
    bars2 = ax.bar(x + width / 2, keep_pct, width,
                   label="keep_suggested = TRUE\n(default confidence pass)",
                   color=COLORS["kept"], edgecolor=COLORS["axis"], linewidth=1)

    # Value annotations
    for i, (r, k, n) in enumerate(zip(ratio_pct, keep_pct, n_luth)):
        r_val = df.iloc[i]["n_ratio_detected"]
        k_val = df.iloc[i]["n_keep_suggested_any"]
        ax.text(i - width / 2, r + 2, f"{int(r_val)}/{int(n)}\n({r:.0f}%)",
                ha="center", va="bottom", fontsize=annot_fontsize,
                fontweight="bold", color=COLORS["ratio_detected"])
        ax.text(i + width / 2, k + 2, f"{int(k_val)}/{int(n)}\n({k:.0f}%)",
                ha="center", va="bottom", fontsize=annot_fontsize,
                fontweight="bold", color=COLORS["kept"])

    # Cohort n annotation under x-axis
    xtick_labels = [f"{c}\n(n={int(n)})" for c, n in zip(compounds_ordered, n_luth)]
    ax.set_xticks(x)
    ax.set_xticklabels(xtick_labels, fontsize=label_fontsize, rotation=rotation,
                        fontweight="bold" if style == "slide" else "normal")

    ax.set_ylabel("CLOverCNV Recall of Luth 2024 CNVs (%)",
                  fontsize=label_fontsize, fontweight="bold")
    ax.set_ylim(0, 115)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(["0", "25", "50", "75", "100"], fontsize=base_fontsize)

    title_text = ("Figure 1. CLOverCNV recovers 88% of Luth 2024 CNVs at coverage level;\n"
                  "default confidence thresholds filter many real signals"
                  if style == "detailed" else
                  "CLOverCNV recovers 88% of Luth CNVs at coverage level")
    ax.set_title(title_text, fontsize=title_fontsize, fontweight="bold", pad=15)

    # Legend outside on the right with titled box
    leg = ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
                    fontsize=base_fontsize, frameon=True, framealpha=1,
                    edgecolor=COLORS["axis"], title="Detection level",
                    title_fontsize=base_fontsize)
    leg.get_title().set_fontweight("bold")

    # Line at 100%
    ax.axhline(y=100, color=COLORS["axis"], linestyle=":", linewidth=1, alpha=0.5)

    guler_style_axes(ax, fontsize=base_fontsize)

    plt.tight_layout()
    savefig_both(fig, f"figure_1_recall_by_compound_{style}")
    plt.close()


# ============================================================
# FIGURE 2: Per-locus heatmap
# ============================================================

def figure_2(luth, per_cnv, style):
    """
    Heatmap of Luth CNV × sample; cell = max_ratio observed in Luth's region.
    Rows: Luth CNVs (ordered by compound, then by chromosome).
    Cols: Cohort samples.
    """
    # Build matrix: rows = Luth CNV (uniquely id'd by clone+chrom+start),
    # cols = SRR of the same-compound sample
    per_cnv = per_cnv.copy()
    per_cnv["luth_key"] = (per_cnv["compound"] + " | " +
                           per_cnv["chrom"].str.replace("Pf3D7_", "chr").str.replace("_v3", "") +
                           ":" + per_cnv["cnv_start"].astype(str).str[:4] + "kb" +
                           " | " + per_cnv["first_gene"])

    # Only 1 sample per Luth CNV (per-cnv table is per-clone), so matrix is diagonal-ish
    # But we want to show, for each Luth locus, what the matched sample looks like
    # Sort by compound (Suloctidil, MMV006901, MMV019662, DSM267) then chromosome/start
    compound_order = {"Suloctidil": 0, "MMV006901": 1, "MMV019662": 2, "DSM267": 3}
    per_cnv["comp_order"] = per_cnv["compound"].map(compound_order)
    per_cnv = per_cnv.sort_values(["comp_order", "chrom", "cnv_start"])

    # For display: y-axis = Luth CNV (labeled with compound + locus)
    # Values = max_ratio; NaN or below threshold shown differently
    ratios = per_cnv["max_ratio_in_region"].fillna(0).values
    kept = per_cnv["keep_suggested_any"].values
    n_rows = len(per_cnv)

    if style == "detailed":
        fig, ax = plt.subplots(figsize=(6, max(7, n_rows * 0.20)))
        base_fontsize = 8
        title_fontsize = 11
        label_fontsize = 9
        annot_fontsize = 7
        show_all_labels = True
    else:
        fig, ax = plt.subplots(figsize=(8, max(9, n_rows * 0.24)))
        base_fontsize = 10
        title_fontsize = 14
        label_fontsize = 11
        annot_fontsize = 8
        show_all_labels = True  # Keep labels — this is the whole point

    # Draw heatmap as single column
    matrix = ratios.reshape(-1, 1)
    im = ax.imshow(matrix, aspect="auto", cmap=GAIN_CMAP,
                   vmin=0.8, vmax=3.0, extent=[-0.5, 0.5, n_rows - 0.5, -0.5])

    # Overlay: circle for keep_suggested=TRUE, red X for missed at ratio level
    for i, (r, k) in enumerate(zip(ratios, kept)):
        if r < 1.5:  # ratio not detected
            ax.plot(0, i, "x", color=COLORS["ratio_missed"], markersize=12,
                    markeredgewidth=2.5, zorder=5)
        if k:  # kept_suggested = TRUE
            ax.plot(0, i, "o", markerfacecolor="none",
                    markeredgecolor=COLORS["kept"],
                    markersize=15, markeredgewidth=2, zorder=4)

    # Left axis: Luth CNV labels
    ylabels = []
    for _, row in per_cnv.iterrows():
        chrom_short = row["chrom"].replace("Pf3D7_", "").replace("_v3", "")
        gene = row["first_gene"].replace("PF3D7_", "PF")
        ylabels.append(f"{row['compound']:10s} | chr{chrom_short}:{int(row['cnv_start'])/1e6:.2f}M | {gene}")
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(ylabels, fontsize=base_fontsize, family="monospace")

    # Right side: max_ratio value annotation
    for i, r in enumerate(ratios):
        if r > 0:
            ax.text(0.55, i, f"{r:.2f}×", va="center", ha="left",
                    fontsize=annot_fontsize, color=COLORS["text"],
                    family="monospace")
        else:
            ax.text(0.55, i, "no overlap", va="center", ha="left",
                    fontsize=annot_fontsize, color=COLORS["ratio_missed"],
                    family="monospace", style="italic")

    ax.set_xticks([])
    ax.set_xlabel("")

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.5, pad=0.35, aspect=15)
    cbar.set_label("max_ratio in Luth region\n(1.0 = baseline)",
                   fontsize=base_fontsize, fontweight="bold")
    cbar.ax.tick_params(labelsize=base_fontsize - 1)

    title_text = (
        "Figure 2. Per-CNV correspondence: CLOverCNV signal at each Luth 2024 locus\n"
        "Circle = keep_suggested=TRUE  |  Red × = ratio below 1.5 (miss)"
        if style == "detailed" else
        "Per-CNV correspondence with Luth 2024\n"
        "○ = keep_suggested  |  × = miss")
    ax.set_title(title_text, fontsize=title_fontsize, fontweight="bold", pad=15)

    # Bold compound-transition lines
    prev_compound = None
    for i, row in enumerate(per_cnv.itertuples()):
        if prev_compound is not None and row.compound != prev_compound:
            ax.axhline(y=i - 0.5, color=COLORS["axis"], linewidth=1.5, zorder=3)
        prev_compound = row.compound

    guler_style_axes(ax, fontsize=base_fontsize, remove_top_right=False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.2)

    plt.tight_layout()
    savefig_both(fig, f"figure_2_percnv_heatmap_{style}")
    plt.close()


# ============================================================
# TABLE 1: Per-compound summary
# ============================================================

def table_1(by_compound):
    """
    Format per-compound summary as CSV/TSV + PNG table for slides.
    """
    compounds_ordered = ["Suloctidil", "MMV006901", "MMV019662", "DSM267"]
    df = by_compound.set_index("compound")
    df = df.loc[[c for c in compounds_ordered if c in df.index]].reset_index()

    # Save TSV / CSV (data)
    df.to_csv(f"{OUT_DIR}/table_1_per_compound_summary.tsv", sep="\t", index=False)
    df.to_csv(f"{OUT_DIR}/table_1_per_compound_summary.csv", index=False)

    # Rendered pretty PNG for slides
    display_df = pd.DataFrame({
        "Compound": df["compound"],
        "N Luth CNVs": df["n_luth_cnvs_checked"],
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

    # Add TOTAL row
    total_n = df["n_luth_cnvs_checked"].sum()
    total_ratio = df["n_ratio_detected"].sum()
    total_kept = df["n_keep_suggested_any"].sum()
    total_row = pd.DataFrame({
        "Compound": ["TOTAL"],
        "N Luth CNVs": [total_n],
        "Ratio recall": [f"{int(total_ratio)}/{int(total_n)} ({total_ratio/total_n*100:.1f}%)"],
        "keep_suggested recall": [f"{int(total_kept)}/{int(total_n)} ({total_kept/total_n*100:.1f}%)"],
        "Median max_ratio": [df["median_max_ratio_in_luth_region"].median().round(2)],
    })
    display_df = pd.concat([display_df, total_row], ignore_index=True)

    # Render two styles
    for style in ["detailed", "slide"]:
        if style == "detailed":
            fig_w, fig_h = 9, 3.2
            fontsize = 10
            header_fontsize = 11
        else:
            fig_w, fig_h = 11, 3.5
            fontsize = 13
            header_fontsize = 14

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax.axis("off")

        cell_text = display_df.values.tolist()
        col_labels = display_df.columns.tolist()

        table = ax.table(cellText=cell_text,
                         colLabels=col_labels,
                         loc="center",
                         cellLoc="center",
                         colColours=[COLORS["kept"]] * len(col_labels))

        table.auto_set_font_size(False)
        table.set_fontsize(fontsize)
        table.scale(1, 2)

        # Style header row
        for i in range(len(col_labels)):
            cell = table[(0, i)]
            cell.set_text_props(color="white", fontweight="bold", fontsize=header_fontsize)
            cell.set_height(0.14)
            cell.set_edgecolor(COLORS["axis"])
            cell.set_linewidth(1.5)

        # Style data rows
        for r in range(1, len(cell_text) + 1):
            for c in range(len(col_labels)):
                cell = table[(r, c)]
                cell.set_edgecolor(COLORS["axis"])
                cell.set_linewidth(0.8)
                cell.set_height(0.11)
                # Bold TOTAL row
                if r == len(cell_text):
                    cell.set_text_props(fontweight="bold", fontsize=fontsize + 1)
                    cell.set_facecolor("#FFF9C4")  # light yellow
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


# ============================================================
# MAIN
# ============================================================

def main():
    luth, per_cnv, by_compound = load_data()
    print(f"\nOutput directory: {OUT_DIR}\n")

    print("Generating Figure 1 (recall bar chart)...")
    figure_1(by_compound, "detailed")
    figure_1(by_compound, "slide")

    print("\nGenerating Figure 2 (per-locus heatmap)...")
    figure_2(luth, per_cnv, "detailed")
    figure_2(luth, per_cnv, "slide")

    print("\nGenerating Table 1 (per-compound summary)...")
    table_1(by_compound)

    print(f"\n{'='*60}")
    print(f"All deliverables saved to {OUT_DIR}")
    print(f"{'='*60}")
    for f in sorted(os.listdir(OUT_DIR)):
        size = os.path.getsize(f"{OUT_DIR}/{f}") / 1024
        print(f"  {f}  ({size:.1f} KB)")


if __name__ == "__main__":
    main()
