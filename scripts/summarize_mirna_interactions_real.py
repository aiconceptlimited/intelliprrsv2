#!/usr/bin/env python3
# ============================================================
# 🧠 IntelliPRRSV2 — miRNA Interaction Summary & Visualization Prep
# ============================================================

import os
import pandas as pd
from datetime import datetime

# ------------------------------------------------------------
# INPUT / OUTPUT
# ------------------------------------------------------------
IN_FILE = str(MIRNA_RESULT_FILE)
OUT_DIR = str(SUMMARY_RESULTS_DIR)
os.makedirs(OUT_DIR, exist_ok=True)

OUT_TOP_MIRNA_PER_VIRUS = os.path.join(OUT_DIR, "top_mirna_per_virus.csv")
OUT_TOP_VIRAL_TARGETS_PER_MIRNA = os.path.join(OUT_DIR, "top_viral_targets_per_mirna.csv")
OUT_GLOBAL_STATS = os.path.join(OUT_DIR, "mirna_summary_overview.csv")

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
df = pd.read_csv(IN_FILE)
print(f"📄 Loaded {len(df):,} miRNA–virus interactions from {IN_FILE}")

# ------------------------------------------------------------
# SUMMARY 1: Top miRNAs per PRRSV strain
# ------------------------------------------------------------
top_mirna = (
    df.sort_values(["accession", "binding_score"], ascending=[True, False])
    .groupby("accession")
    .head(10)
)
top_mirna.to_csv(OUT_TOP_MIRNA_PER_VIRUS, index=False)
print(f"✅ Saved top 10 miRNAs per PRRSV strain → {OUT_TOP_MIRNA_PER_VIRUS}")

# ------------------------------------------------------------
# SUMMARY 2: Top viral genomes targeted by each miRNA
# ------------------------------------------------------------
top_virus = (
    df.sort_values(["miRNA", "binding_score"], ascending=[True, False])
    .groupby("miRNA")
    .head(10)
)
top_virus.to_csv(OUT_TOP_VIRAL_TARGETS_PER_MIRNA, index=False)
print(f"✅ Saved top 10 viral targets per miRNA → {OUT_TOP_VIRAL_TARGETS_PER_MIRNA}")

# ------------------------------------------------------------
# SUMMARY 3: Global binding stats
# ------------------------------------------------------------
summary = (
    df.groupby("miRNA")["binding_score"]
    .agg(["count", "mean", "max"])
    .reset_index()
    .sort_values("mean", ascending=False)
)
summary["timestamp"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
summary.to_csv(OUT_GLOBAL_STATS, index=False)
print(f"✅ Global interaction summary saved → {OUT_GLOBAL_STATS}")

# ------------------------------------------------------------
# PRINT KEY INSIGHTS
# ------------------------------------------------------------
print("\n📊 Top 5 miRNAs with strongest mean binding scores:")
print(summary.head(5).to_string(index=False))
print("\n🏁 Summary complete — ready for dashboard visualization.")

