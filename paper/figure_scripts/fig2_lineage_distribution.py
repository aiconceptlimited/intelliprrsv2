#!/usr/bin/env python3

import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


# Resolve paths relative to this script
SCRIPT_DIR = Path(__file__).resolve().parent
PAPER_DIR = SCRIPT_DIR.parent

DATA_FILE = PAPER_DIR / "figure_data" / "Fig2_lineages.csv"
FIGURE_DIR = PAPER_DIR / "figures"

FIGURE_DIR.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# Load validated source data
# ------------------------------------------------------------
df = pd.read_csv(DATA_FILE)

required_columns = {"lineage", "count", "percentage"}

if not required_columns.issubset(df.columns):
    raise SystemExit(
        f"ERROR: required columns missing from {DATA_FILE}"
    )

if df["count"].sum() != 1866:
    raise SystemExit(
        f"ERROR: expected 1,866 isolates, found {df['count'].sum()}"
    )

if df["lineage"].nunique() != 4:
    raise SystemExit(
        f"ERROR: expected 4 lineages, found {df['lineage'].nunique()}"
    )


# ------------------------------------------------------------
# Create deterministic publication canvas
# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6), dpi=300)

bars = ax.bar(
    df["lineage"],
    df["count"],
)


# ------------------------------------------------------------
# Value labels
# ------------------------------------------------------------
for bar, count, pct in zip(
    bars,
    df["count"],
    df["percentage"],
):
    height = bar.get_height()

    ax.text(
        bar.get_x() + bar.get_width() / 2,
        height + 25,
        f"{count:,}\n({pct:.2f}%)",
        ha="center",
        va="bottom",
        fontsize=10,
    )


# ------------------------------------------------------------
# Labels
# ------------------------------------------------------------
ax.set_xlabel(
    "Lineage",
    fontsize=12,
)

ax.set_ylabel(
    "Number of Isolates",
    fontsize=12,
)

ax.set_title(
    "PRRSV Lineage Distribution",
    fontsize=14,
    fontweight="bold",
)


# Leave sufficient headroom for labels
ax.set_ylim(
    0,
    max(df["count"]) * 1.15,
)


# Y-axis formatting
ax.yaxis.set_major_formatter(
    plt.FuncFormatter(
        lambda x, p: format(int(x), ",")
    )
)


# Remove unnecessary spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)


# Deterministic margins; do not use bbox_inches='tight'
fig.subplots_adjust(
    left=0.14,
    right=0.97,
    bottom=0.14,
    top=0.88,
)


# ------------------------------------------------------------
# Save publication outputs
# ------------------------------------------------------------
PNG = FIGURE_DIR / "Fig2_v1.0.png"
PDF = FIGURE_DIR / "Fig2_v1.0.pdf"

fig.savefig(
    PNG,
    dpi=300,
)

fig.savefig(
    PDF,
)

plt.close(fig)

print(
    "PASS: Figure 2 generated:",
    PNG,
    "and",
    PDF,
)
