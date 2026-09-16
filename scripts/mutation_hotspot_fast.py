#!/usr/bin/env python3
# ============================================================
# 🧬 IntelliPRRSV2 — Mutation Hotspot Analysis (Corrected Scientific Version)
# ============================================================
# Improvements:
# - Filters full genomes only (≥14500 bp)
# - Proper FASTA parsing using Biopython
# - Ignores gaps and ambiguous bases
# - Skips low coverage positions
# - Avoids reference gap artifacts
# - Maintains MySQL compatibility
# ============================================================

import os
import pandas as pd
import subprocess
import mysql.connector
from .config import DB_CONFIG
from .paths import CLEAN_DIR, TEMP_DIR, MUTATION_RESULT_FILE
from datetime import datetime
from tqdm import tqdm
from Bio import SeqIO
from intelliprrsv2.scripts.db_utils import insert_dataframe

# ============================================================
# PATH CONFIGURATION
# ============================================================
CLEAN_DIR = str(CLEAN_DIR)
TEMP_DIR = str(TEMP_DIR)
RESULT_FILE = str(MUTATION_RESULT_FILE)

os.makedirs(TEMP_DIR, exist_ok=True)
os.makedirs(os.path.dirname(RESULT_FILE), exist_ok=True)

# ============================================================
# ORF MAPPING (PRRSV-2 canonical reference)
# ============================================================
orf_map = [
    ("ORF1a", 1, 7080),
    ("ORF1b", 7081, 11830),
    ("ORF2", 11831, 12600),
    ("ORF3", 12601, 13320),
    ("ORF4", 13321, 14050),
    ("ORF5", 14051, 14780),
    ("ORF6", 14781, 15200),
    ("ORF7", 15201, 15500),
    ("3'UTR", 15501, 15600)
]

def get_orf(pos):
    for name, start, end in orf_map:
        if start <= pos <= end:
            return name
    return "Non-coding"

# ============================================================
# MYSQL CONNECTION
# ============================================================
def get_sequence_id_map():
    try:
        conn = mysql.connector.connect(**DB_CONFIG)
        cursor = conn.cursor()
        cursor.execute("SELECT accession, id FROM sequences")
        mapping = {a: i for a, i in cursor.fetchall()}
        cursor.close()
        conn.close()
        return mapping
    except Exception as e:
        print(f"⚠️ MySQL connection failed: {e}")
        return {}

# ============================================================
# STEP 1 — Combine FULL genomes only
# ============================================================
print("🧬 Preparing combined FASTA (full genomes only)...")

combined_fasta = os.path.join(TEMP_DIR, "all_sequences.fasta")

with open(combined_fasta, "w") as out:
    for f in os.listdir(CLEAN_DIR):
        if f.endswith(".fasta"):
            file_path = os.path.join(CLEAN_DIR, f)
            record = next(SeqIO.parse(file_path, "fasta"))

            # ✅ FILTER: keep only near-full genomes
            if len(record.seq) >= 14500:
                SeqIO.write(record, out, "fasta")

# ============================================================
# STEP 2 — Run MAFFT Alignment
# ============================================================
aligned_file = os.path.join(TEMP_DIR, "aligned.fasta")

print("⚙️ Running MAFFT alignment...")
subprocess.run(
    ["mafft", "--auto", "--quiet", combined_fasta],
    stdout=open(aligned_file, "w"),
)

# ============================================================
# STEP 3 — Mutation Frequency Calculation (Corrected)
# ============================================================
print("🔬 Analyzing mutation frequencies...")

alignment_records = list(SeqIO.parse(aligned_file, "fasta"))

if len(alignment_records) < 2:
    raise ValueError("Not enough aligned sequences for analysis.")

ref_seq = str(alignment_records[0].seq).upper()
other_seqs = [str(rec.seq).upper() for rec in alignment_records[1:]]

positions = []
run_ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

print(f"🔢 Reference length: {len(ref_seq)} bases")
print(f"🧬 Comparing against {len(other_seqs)} genomes")

for i in tqdm(range(len(ref_seq)), desc="Computing mutation frequencies"):

    # Skip if reference base is not valid nucleotide
    if ref_seq[i] not in ["A", "T", "G", "C"]:
        continue

    valid_bases = [
        seq[i] for seq in other_seqs
        if i < len(seq) and seq[i] in ["A", "T", "G", "C"]
    ]

    # Skip low coverage positions
    if len(valid_bases) < 100:
        continue

    freq = 1 - (valid_bases.count(ref_seq[i]) / len(valid_bases))

    if freq <= 0:
        continue

    variants = set(b for b in valid_bases if b != ref_seq[i])

    for var in variants:
        positions.append({
            "position": i + 1,
            "ref_base": ref_seq[i],
            "mutation_type": f"{ref_seq[i]}→{var}",
            "mutation_frequency": round(freq, 4),
            "orf": get_orf(i + 1),
            "timestamp": run_ts,
            "run_timestamp": run_ts
        })

# ============================================================
# STEP 4 — Save Results
# ============================================================
df = pd.DataFrame(positions)

if df.empty:
    print("⚠️ No valid mutation hotspots detected.")
    exit()

df.to_csv(RESULT_FILE, index=False)
print(f"✅ Mutation hotspot analysis complete → {RESULT_FILE}")

# ============================================================
# STEP 5 — Add sequence_id placeholder (safe)
# ============================================================
seq_map = get_sequence_id_map()
df["accession"] = None
df["sequence_id"] = None

# ============================================================
# STEP 6 — Insert into MySQL
# ============================================================
try:
    insert_dataframe(df, "mutation_hotspots")
    print(f"✅ Inserted {len(df)} mutation records into MySQL.")
except Exception as e:
    print(f"❌ MySQL insert failed: {e}")

print("🏁 Mutation hotspot layer complete.")
