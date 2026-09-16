#!/usr/bin/env python3
"""
🧬 IntelliPRRSV2 — Auto ORF Annotation Mapper (with MySQL Logging)
-----------------------------------------------------------------
Auto-detects PRRSV type (VR2332 / Lelystad) based on sequence metadata,
updates ORFs for mutation hotspots, and logs each run into MySQL.
-----------------------------------------------------------------
"""

import pandas as pd
from sqlalchemy import create_engine, text
from datetime import datetime
from .config import DB_CONFIG

# ============================================================
# DATABASE CONFIGURATION
# ============================================================
ENGINE_URL = (
    f"mysql+mysqlconnector://{DB_CONFIG['user']}:{DB_CONFIG['password']}"
    f"@{DB_CONFIG['host']}/{DB_CONFIG['database']}"
)
engine = create_engine(ENGINE_URL)

# ============================================================
# REFERENCE ORF MAPS
# ============================================================
ORF_MAPS = {
    "VR2332": {
        "5' UTR": (1, 189),
        "ORF1a": (190, 7260),
        "ORF1b": (7261, 11860),
        "GP2 (ORF2a)": (11861, 12588),
        "E (ORF2b)": (12589, 12812),
        "GP3 (ORF3)": (12813, 13558),
        "GP4 (ORF4)": (13559, 14245),
        "GP5 (ORF5)": (14246, 14913),
        "M (ORF6)": (14914, 15421),
        "N (ORF7)": (15422, 15922),
        "3' UTR": (15923, 16100),
    },
    "Lelystad": {
        "5' UTR": (1, 190),
        "ORF1a": (191, 7265),
        "ORF1b": (7266, 11760),
        "GP2 (ORF2a)": (11761, 12517),
        "E (ORF2b)": (12518, 12749),
        "GP3 (ORF3)": (12750, 13477),
        "GP4 (ORF4)": (13478, 14215),
        "GP5 (ORF5)": (14216, 14920),
        "M (ORF6)": (14921, 15480),
        "N (ORF7)": (15481, 16020),
        "3' UTR": (16021, 16100),
    },
}

# ============================================================
# AUTO-DETECTION OF REFERENCE TYPE
# ============================================================
def detect_reference(conn):
    """Determine which reference (VR2332 or Lelystad) best fits current data."""
    try:
        df_seq = pd.read_sql(
            "SELECT accession, LENGTH(sequence) AS genome_len FROM sequences LIMIT 100", conn
        )
        if df_seq.empty:
            print("⚠️ No sequence data found — defaulting to VR2332.")
            return "VR2332"

        avg_len = df_seq["genome_len"].mean()
        prefix = df_seq["accession"].astype(str).str[:2].mode().iloc[0]

        eu_prefixes = ("AM", "AY", "AF", "DQ", "EU", "FN", "AB")
        na_prefixes = ("KR", "KF", "KJ", "EF", "GU", "MN")

        if prefix in eu_prefixes or avg_len < 15500:
            ref = "Lelystad"
        elif prefix in na_prefixes or avg_len >= 15500:
            ref = "VR2332"
        else:
            ref = "VR2332"

        print(f"🧠 Auto-detected reference: {ref}  (avg_len ≈ {avg_len:.0f}, prefix={prefix})")
        return ref

    except Exception as e:
        print(f"⚠️ Detection failed → defaulting to VR2332: {e}")
        return "VR2332"

# ============================================================
# MAIN EXECUTION
# ============================================================
with engine.begin() as conn:
    # Detect reference
    ref = detect_reference(conn)
    orf_map = ORF_MAPS[ref]

    # Load mutation data
    df = pd.read_sql("SELECT id, position FROM mutation_hotspots", conn)
    print(f"✅ Loaded {len(df)} mutation records")

    # Assign ORF based on reference
    def assign_orf(pos):
        if pd.isna(pos):
            return "Unknown"
        for orf, (start, end) in orf_map.items():
            if start <= pos <= end:
                return orf
        return "Intergenic"

    df["position"] = pd.to_numeric(df["position"], errors="coerce")
    df["orf_real"] = df["position"].apply(assign_orf)

    # Update MySQL
    update_sql = text("UPDATE mutation_hotspots SET orf = :orf WHERE id = :id")
    for _, row in df.iterrows():
        conn.execute(update_sql, {"orf": row["orf_real"], "id": int(row["id"])})

    # Log run summary
    conn.execute(
        text("""
        CREATE TABLE IF NOT EXISTS orf_mapping_log (
            id INT AUTO_INCREMENT PRIMARY KEY,
            reference_used VARCHAR(20),
            total_records INT,
            run_timestamp DATETIME
        )
        """)
    )

    conn.execute(
        text("""
        INSERT INTO orf_mapping_log (reference_used, total_records, run_timestamp)
        VALUES (:ref, :total, :ts)
        """),
        {"ref": ref, "total": len(df), "ts": datetime.now()},
    )

print("✅ ORF mapping complete — auto-detected reference logged in MySQL.")

