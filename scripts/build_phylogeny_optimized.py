#!/usr/bin/env python3
# ============================================================
# 🌳 IntelliPRRSV2 — Ultrafast Phylogenetic Tree + Real Lineage Tracker (Enhanced + Deduplication Fix)
# ============================================================

import os
import subprocess
import pandas as pd
import numpy as np
from datetime import datetime
from Bio import SeqIO, AlignIO
from scipy.cluster.hierarchy import linkage, fcluster
import mysql.connector
from intelliprrsv2.scripts.db_utils import insert_dataframe
from intelliprrsv2.scripts.config import DB_CONFIG
from intelliprrsv2.scripts.paths import PROJECT_ROOT, CLEAN_DIR, RESULTS_DIR, REFERENCE_DIR

# ------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------
BASE_DIR = str(PROJECT_ROOT)
CLEAN_DIR = f"{BASE_DIR}/data/clean/"
RESULTS_DIR = f"{BASE_DIR}/results/"
os.makedirs(RESULTS_DIR, exist_ok=True)

COMBINED_FASTA = f"{RESULTS_DIR}/combined_cleaned_ultrafast.fasta"
ALIGNED_FILE = f"{RESULTS_DIR}/aligned_ultrafast.fasta"
TREE_FILE = f"{RESULTS_DIR}/prrsv_tree_ultrafast.nwk"
LINEAGE_TABLE = f"{RESULTS_DIR}/lineage_table_ultrafast.csv"
REF_FASTA = f"{BASE_DIR}/data/reference_lineages.fasta"
FASTTREE_ERROR_LOG = f"{RESULTS_DIR}/fasttree_error.log"

# ------------------------------------------------------------
# STEP 1 — Combine Clean FASTA Files
# ------------------------------------------------------------
records = []
for file in os.listdir(CLEAN_DIR):
    if file.endswith(".fasta"):
        records.extend(list(SeqIO.parse(os.path.join(CLEAN_DIR, file), "fasta")))

if not records:
    print("❌ No cleaned FASTA files found.")
    exit(1)

# 🧩 Ensure all sequence IDs are unique before alignment
unique_records = []
seen = {}
for rec in records:
    base_id = rec.id.strip()
    if base_id in seen:
        seen[base_id] += 1
        rec.id = f"{base_id}_{seen[base_id]}"
    else:
        seen[base_id] = 1
    rec.description = ""  # simplify header
    unique_records.append(rec)

SeqIO.write(unique_records, COMBINED_FASTA, "fasta")
print(f"✅ Normalized and deduplicated {len(unique_records)} sequence IDs.")
print(f"🧬 Combined {len(unique_records)} cleaned genomes for ultrafast tree building.")
# ------------------------------------------------------------
# STEP 2 — Run MAFFT Alignment (Safe version)
# ------------------------------------------------------------
print("⚙️ Running MAFFT (ultrafast mode, safe fallback)...")

mafft_cmd = ["mafft", "--anysymbol", "--auto", "--thread", "2", COMBINED_FASTA]
with open(ALIGNED_FILE, "w") as out_f:
    result = subprocess.run(
        mafft_cmd,
        stdout=out_f,
        stderr=subprocess.PIPE,
        text=True
    )

# MAFFT sometimes returns 1 even when successful — we’ll verify by file output instead
if result.returncode not in (0, 1):
    print("❌ MAFFT failed with exit code", result.returncode)
    print(result.stderr[:400])
    with open(f"{RESULTS_DIR}/mafft_error.log", "w") as log:
        log.write(result.stderr)
    exit(1)
elif os.path.getsize(ALIGNED_FILE) == 0:
    print("❌ MAFFT produced an empty alignment file.")
    exit(1)
else:
    print("✅ MAFFT alignment complete.")


# ------------------------------------------------------------
# STEP 3 — Build Phylogenetic Tree (FastTree)
# ------------------------------------------------------------
print("🌳 Running FastTree (GTR model)...")
try:
    result = subprocess.run(
        ["FastTree", "-nt", "-gtr", ALIGNED_FILE],
        stdout=open(TREE_FILE, "w"),
        stderr=subprocess.PIPE,
        text=True,
        check=True,
    )
    print(f"✅ Tree saved → {TREE_FILE}")
except subprocess.CalledProcessError as e:
    print("❌ FastTree failed to build the tree.")
    print("🧾 Error message from FastTree:\n")
    print(e.stderr)
    with open(FASTTREE_ERROR_LOG, "w") as f:
        f.write(e.stderr)
    print(f"⚠️ FastTree error log saved → {FASTTREE_ERROR_LOG}")
    exit(1)
# ------------------------------------------------------------
# STEP 4 — Compute Genetic Distance Clusters
# ------------------------------------------------------------
print("🧠 Computing real genetic distance–based clusters...")
alignment = AlignIO.read(ALIGNED_FILE, "fasta")
n = len(alignment)
dist_matrix = np.zeros((n, n))

# Use simplified pairwise identity for clustering
for i in range(min(n, 200)):  # limit comparisons for speed
    seq_i = np.array(list(str(alignment[i].seq)))
    for j in range(i + 1, min(n, 200)):
        seq_j = np.array(list(str(alignment[j].seq)))
        matches = np.sum(seq_i == seq_j)
        identity = matches / len(seq_i)
        dist_matrix[i, j] = 1 - identity
        dist_matrix[j, i] = dist_matrix[i, j]

Z = linkage(dist_matrix[np.triu_indices(min(n, 200), 1)], method="ward")
clusters_ids = fcluster(Z, t=9, criterion="maxclust")
clusters = {alignment[i].id: f"L{clusters_ids[i]}" for i in range(min(n, 200))}
# ------------------------------------------------------------
# STEP 5 — Hybrid Reference-Based Lineage Adjustment (Ultra-Fast)
# ------------------------------------------------------------
if os.path.exists(REF_FASTA):
    print("🚀 Refining cluster lineages using ultra-fast k-mer matching...")

    from concurrent.futures import ThreadPoolExecutor
    from Bio import pairwise2

    # Load reference sequences
    ref_map = {
        rec.id.split("_")[0].replace("REF", "").strip().upper(): str(rec.seq).upper()
        for rec in SeqIO.parse(REF_FASTA, "fasta")
    }

    # Precompute k-mers for each reference sequence
    def get_kmers(seq, k=6):
        return {seq[i:i+k] for i in range(len(seq) - k)}

    ref_kmers = {lineage: get_kmers(seq) for lineage, seq in ref_map.items()}

    def fast_score(seq):
        """Compute best lineage using fast k-mer overlap."""
        seq = seq.upper()
        seq_kmers = get_kmers(seq)
        best_lineage, best_score = "Unclassified", 0

        for lineage, refset in ref_kmers.items():
            inter = len(seq_kmers & refset)
            union = len(seq_kmers | refset)
            jaccard = inter / union if union > 0 else 0
            if jaccard > best_score:
                best_lineage, best_score = lineage, jaccard

        # Optional quick confirmation for close calls
        if best_score > 0.4:
            ref_seq = ref_map[best_lineage]
            score = pairwise2.align.globalxx(seq[:1000], ref_seq[:1000], score_only=True)
            if score / len(seq) < 0.2:
                best_lineage = "Unclassified"
        return best_lineage

    # Use parallel processing to assign lineages
    with ThreadPoolExecutor(max_workers=8) as executor:
        seqs = [(a.id, str(a.seq)) for a in alignment]
        results = list(executor.map(lambda x: (x[0], fast_score(x[1])), seqs))

    clusters = dict(results)



# ------------------------------------------------------------
# STEP 6 — Save Lineage Assignments
# ------------------------------------------------------------
df = pd.DataFrame(list(clusters.items()), columns=["accession", "assigned_lineage"])
df["run_timestamp"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
df.to_csv(LINEAGE_TABLE, index=False)
print(f"✅ Lineage assignments saved → {LINEAGE_TABLE}")


# Insert into MySQL
try:
    insert_dataframe(df, "lineage_assignments")
    print("✅ Lineage assignments inserted into MySQL.")
except Exception as e:
    print(f"❌ Failed to insert lineages into MySQL: {e}")

# ------------------------------------------------------------
# STEP 7 — Lineage Diversity Summary
# ------------------------------------------------------------
print("📊 Summarizing lineage diversity trends...")
lineage_counts = df["assigned_lineage"].value_counts().to_dict()

try:
    conn = mysql.connector.connect(**DB_CONFIG)
    cursor = conn.cursor()
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS lineage_summary (
            id INT AUTO_INCREMENT PRIMARY KEY,
            lineage VARCHAR(20),
            count INT,
            run_timestamp DATETIME
        )
    """)
    for lineage, count in lineage_counts.items():
        cursor.execute("""
            INSERT INTO lineage_summary (lineage, count, run_timestamp)
            VALUES (%s, %s, NOW())
        """, (lineage, int(count)))
    conn.commit()
    cursor.close()
    conn.close()
    print("✅ Lineage diversity summary inserted into MySQL.")
except Exception as e:
    print(f"❌ Failed to insert lineage summary: {e}")

# ------------------------------------------------------------
# STEP 8 — Insert Phylogenetic Tree into MySQL
# ------------------------------------------------------------
try:
    with open(TREE_FILE, "r") as f:
        newick_data = f.read().strip()

    conn = mysql.connector.connect(**DB_CONFIG)
    cursor = conn.cursor()
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS phylogenetic_tree (
            id INT AUTO_INCREMENT PRIMARY KEY,
            newick_data LONGTEXT,
            timestamp DATETIME
        )
    """)
    cursor.execute("DELETE FROM phylogenetic_tree")
    cursor.execute(
        "INSERT INTO phylogenetic_tree (newick_data, timestamp) VALUES (%s, %s)",
        (newick_data, datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
    )
    conn.commit()
    conn.close()
    print("✅ Phylogenetic tree stored in MySQL.")
except Exception as e:
    print(f"❌ Failed to store phylogenetic tree: {e}")

print("🏁 Ultrafast phylogenetic analysis complete — lineages, diversity, and tree synced to MySQL.")

