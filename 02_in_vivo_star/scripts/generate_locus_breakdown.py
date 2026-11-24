#!/usr/bin/env python
"""
Output:
  - Table S7 (per-locus breakdown TSV/CSV) 
  - Figure 3: per-compound per-locus stacked bar (detailed + slide versions)
  - Console summary for reference

Data source: 02_in_vivo_star/luth_data_s4_filtered.tsv (42 CNVs, filters exclusion)

Run:
    conda activate CLOverCNV
    python /scratch/gaz6bs/clovercnv-application/02_in_vivo_star/scripts/generate_locus_breakdown.py
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================
# CONFIG
# ============================================================
LUTH_FILE = "/scratch/gaz6bs/clovercnv-application/02_in_vivo_star/luth_data_s4_filtered.tsv"
OUT_DIR = "/scratch/gaz6bs/clovercnv/analysis/deliverables_luth_crossref"
os.makedirs(OUT_DIR, exist_ok=True)

# Target locus assignments (verified from PlasmoDB + literature)
TARGET_LOCI = {
    "Suloctidil":  {"main_chr": "Pf3D7_05_v3", "main_gene": "PfMDR1",   "main_gene_id": "PF3D7_0523000", "main_role": "ABC transporter (MDR)"},
    "MMV006901":   {"main_chr": "Pf3D7_03_v3", "main_gene": "PfABCI3",  "main_gene_id": "PF3D7_0319700", "main_role": "ABC transporter (multidrug)"},
    "MMV019662":   {"main_chr": "Pf3D7_01_v3", "main_gene": "PfNCR1",   "main_gene_id": "PF3D7_0107500", "main_role": "Cholesterol transporter (lipid homeostasis)"},
    "DSM267":      {"main_chr": "Pf3D7_06_v3", "main_gene": "PfDHODH",  "main_gene_id": "PF3D7_0603300", "main_role": "Dihydroorotate dehydrogenase (pyrimidine synthesis)"},
}

# Colors — Guler lab palette
COLORS = {
    "main": "#2E7D32",       # deep green (main target locus)
    "secondary": "#EF6C00",  # deep orange (secondary loci)
    "MMV006901": "#7B1FA2",
    "Suloctidil": "#EF6C00",
    "MMV019662": "#00838F",
    "DSM267": "#455A64",
    "axis": "#212121",
    "text": "#212121",
}


def guler_style_axes(ax, fontsize=10, remove_top_right=True):
    ax.set_facecolor("white")
    for spine_name, spine in ax.spines.items():
        if remove_top_right and spine_name in ("top", "right"):
            spine.set_visible(False)
        else:
            spine.set_linewidth(1.5)
            spine.set_color(COLORS["axis"])
    ax.tick_params(axis="both", which="major", labelsize=fontsize - 1,
                   width=1.2, length=4, color=COLORS["axis"])
    ax.grid(False)


def savefig_both(fig, basename):
    fig.savefig(f"{OUT_DIR}/{basename}.svg", format="svg", bbox_inches="tight", facecolor="white")
    fig.savefig(f"{OUT_DIR}/{basename}.png", format="png", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {basename}.svg + {basename}.png")


# ============================================================
# BUILD PER-LOCUS BREAKDOWN
# ============================================================

def build_locus_breakdown(luth_df):
    """
    Group CNVs by (compound, chromosome), assign as main-target-locus or secondary.
    Excludes CNVs from excluded samples.
    """
    # Exclude EXCLUDED samples
    luth_df = luth_df[~luth_df["srr"].astype(str).str.startswith("EXCLUDED")].copy()

    rows = []
    for compound, group in luth_df.groupby("compound"):
        target_chr = TARGET_LOCI.get(compound, {}).get("main_chr", "unknown")
        target_gene = TARGET_LOCI.get(compound, {}).get("main_gene", "unknown")
        target_id = TARGET_LOCI.get(compound, {}).get("main_gene_id", "")
        target_role = TARGET_LOCI.get(compound, {}).get("main_role", "")

        # Group by chromosome within compound
        for chrom, chr_group in group.groupby("chrom"):
            is_main = chrom == target_chr
            # Determine dominant gene per group (use first_gene mode)
            first_genes = chr_group["gene_list"].apply(lambda x: str(x).split(",")[0])
            dominant_first_gene = first_genes.mode().iloc[0] if len(first_genes) > 0 else ""

            # Position range
            min_start = chr_group["cnv_start"].min()
            max_start = chr_group["cnv_start"].max()
            n_cnvs = len(chr_group)
            n_unique_clones = chr_group["clone_name"].nunique()

            # If this is the main target locus, use annotated gene name
            if is_main:
                label = f"{target_gene} ({target_id})"
                locus_type = "Main target locus (published)"
                biology_note = f"Drug target: {target_role}"
            else:
                label = f"{dominant_first_gene}"
                locus_type = "Secondary locus"
                biology_note = "Requires further annotation; possible off-target amplification, secondary resistance mechanism, or artifact"

            chrom_short = chrom.replace("Pf3D7_", "chr").replace("_v3", "")
            rows.append({
                "compound": compound,
                "locus_type": locus_type,
                "chromosome": chrom_short,
                "region_start_Mb": round(min_start / 1e6, 3),
                "region_end_Mb": round(max_start / 1e6, 3),
                "target_gene_or_region": label,
                "n_cnvs": n_cnvs,
                "n_unique_clones": n_unique_clones,
                "biology_note": biology_note,
            })

    # Sort: within compound, main first, then by chromosome
    df = pd.DataFrame(rows)
    df["_sort_type"] = df["locus_type"].map({"Main target locus (published)": 0, "Secondary locus": 1})
    compound_order = {"Suloctidil": 0, "MMV006901": 1, "MMV019662": 2, "DSM267": 3}
    df["_sort_compound"] = df["compound"].map(compound_order)
    df = df.sort_values(["_sort_compound", "_sort_type", "chromosome"]).reset_index(drop=True)
    df = df.drop(columns=["_sort_type", "_sort_compound"])
    return df


# ============================================================
# TABLE OUTPUT
# ============================================================

def write_table_s7(df):
    """Write TSV + CSV for machine-readable use."""
    df.to_csv(f"{OUT_DIR}/table_S7_per_locus_breakdown.tsv", sep="\t", index=False)
    df.to_csv(f"{OUT_DIR}/table_S7_per_locus_breakdown.csv", index=False)
    print(f"  Saved: table_S7_per_locus_breakdown.tsv + .csv")

    # Console summary
    print("\n" + "=" * 80)
    print("Table S7 preview:")
    print("=" * 80)
    for compound in ["Suloctidil", "MMV006901", "MMV019662", "DSM267"]:
        print(f"\n--- {compound} ---")
        sub = df[df["compound"] == compound]
        total = sub["n_cnvs"].sum()
        for _, r in sub.iterrows():
            marker = "1" if "Main" in r["locus_type"] else "0"
            print(f"  {marker} {r['chromosome']:6} {r['region_start_Mb']:.2f}-{r['region_end_Mb']:.2f} Mb  |  {r['target_gene_or_region']:35}  |  {r['n_cnvs']} CNVs")
        print(f"  Total: {total} CNVs across {len(sub)} loci")


# ============================================================
# FIGURE 3: PER-LOCUS BREAKDOWN
# ============================================================

def figure_3_locus_breakdown(df, style="detailed"):
    """
    Stacked bar chart: per-compound CNV count with main vs secondary loci stacked.
    Includes explicit locus labels and CNV counts.
    """
    compounds = ["Suloctidil", "MMV006901", "MMV019662", "DSM267"]
    compounds = [c for c in compounds if c in df["compound"].unique()]

    # Aggregate main vs secondary per compound
    plot_data = []
    for compound in compounds:
        sub = df[df["compound"] == compound]
        main_row = sub[sub["locus_type"].str.contains("Main")]
        secondary_rows = sub[sub["locus_type"].str.contains("Secondary")]

        n_main = int(main_row["n_cnvs"].sum())
        n_secondary = int(secondary_rows["n_cnvs"].sum())
        main_label = main_row["target_gene_or_region"].iloc[0].split(" (")[0] if len(main_row) else "unknown"

        secondary_labels = []
        for _, r in secondary_rows.iterrows():
            gene_short = str(r["target_gene_or_region"]).replace("PF3D7_", "PF")
            secondary_labels.append(f"{r['chromosome']}:{gene_short} ({r['n_cnvs']})")
        secondary_label = ", ".join(secondary_labels) if secondary_labels else ""

        plot_data.append({
            "compound": compound,
            "n_main": n_main,
            "n_secondary": n_secondary,
            "main_label": main_label,
            "secondary_label": secondary_label,
        })

    plot_df = pd.DataFrame(plot_data)

    # Figure setup
    if style == "detailed":
        fig, ax = plt.subplots(figsize=(9, 6))
        title_fs = 12
        label_fs = 10
        annot_fs = 9
        secondary_fs = 8
    else:
        fig, ax = plt.subplots(figsize=(12, 7))
        title_fs = 16
        label_fs = 14
        annot_fs = 12
        secondary_fs = 10

    x = np.arange(len(plot_df)) * 1.3
    width = 0.55

    # Stacked bars
    ax.bar(x, plot_df["n_main"], width, label="Main target locus (published gene)",
           color=COLORS["main"], edgecolor=COLORS["axis"], linewidth=1)
    ax.bar(x, plot_df["n_secondary"], width, bottom=plot_df["n_main"],
           label="Secondary loci (other chromosomes)",
           color=COLORS["secondary"], edgecolor=COLORS["axis"], linewidth=1)

    # Annotations: main gene name INSIDE the main bar (if room)
    for i, row in plot_df.iterrows():
        # Main label centered on main portion
        y_main_center = row["n_main"] / 2
        if row["n_main"] >= 3:  # only annotate if bar is tall enough
            ax.text(x[i], y_main_center, row["main_label"], ha="center", va="center",
                    fontsize=annot_fs, fontweight="bold", color="white")
        # Count label at top of stack
        y_top = row["n_main"] + row["n_secondary"]
        ax.text(x[i], y_top + 0.3, f"n={y_top}", ha="center", va="bottom",
                fontsize=annot_fs, fontweight="bold", color=COLORS["axis"])
        # Secondary loci breakdown BELOW x-tick label
        if row["secondary_label"]:
            ax.text(x[i], -2.5, row["secondary_label"], ha="center", va="top",
                    fontsize=secondary_fs, color=COLORS["secondary"], style="italic")

    # X-axis
    xtick_labels = [f"{r['compound']}\n(→ {r['main_label']})" for _, r in plot_df.iterrows()]
    ax.set_xticks(x)
    ax.set_xticklabels(xtick_labels, fontsize=label_fs, fontweight="bold")

    ax.set_ylabel("Number of Luth 2024 CNVs (post-exclusion)",
                  fontsize=label_fs, fontweight="bold")
    ax.set_ylim(-4, plot_df[["n_main"]].values.max() + plot_df[["n_secondary"]].values.max() + 3)

    title = ("Figure 3. Per-locus breakdown of Luth 2024 CNVs by drug target\n"
             "Each cohort's CNVs distribute across a main published target + secondary loci"
             if style == "detailed" else
             "Per-locus distribution of Luth CNVs by drug target")
    ax.set_title(title, fontsize=title_fs, fontweight="bold", pad=15)

    leg = ax.legend(loc="upper right", fontsize=label_fs - 1, frameon=True,
                    framealpha=1, edgecolor=COLORS["axis"])
    guler_style_axes(ax, fontsize=label_fs)

    plt.tight_layout()
    savefig_both(fig, f"figure_3_locus_breakdown_{style}")
    plt.close()


# ============================================================
# MAIN
# ============================================================

def main():
    print("Loading Luth Data S4...")
    luth_df = pd.read_csv(LUTH_FILE, sep="\t")
    print(f"  Loaded {len(luth_df)} rows (pre-exclusion)")

    print("\nBuilding per-locus breakdown (Table S7)...")
    breakdown = build_locus_breakdown(luth_df)
    write_table_s7(breakdown)

    print("\nGenerating Figure 3 (per-locus visualization)...")
    figure_3_locus_breakdown(breakdown, "detailed")
    figure_3_locus_breakdown(breakdown, "slide")

    print(f"\n{'='*60}")
    print(f"Done. All deliverables in {OUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
