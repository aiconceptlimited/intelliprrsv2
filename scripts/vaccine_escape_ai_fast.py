#!/usr/bin/env python3

"""
Layer 5 — AI Vaccine Escape Prediction (Real, Fast)

Controlled production implementation.

Reference selection is governed exclusively by the authoritative
vaccine reference manifest and validator.

The prediction mathematics is intentionally unchanged from the
pre-correction implementation:

    identity = globalxx identity over the first 1000 nt
    escape_probability =
        1 / (1 + exp((identity - 85) / -3))

Database insertion is disabled by default.

Use:
    python3 scripts/vaccine_escape_ai_fast.py --dry-run

The dry-run performs prediction and CSV generation but makes no
database changes.
"""

import argparse
import os
from datetime import datetime

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
from tqdm import tqdm

from intelliprrsv2.scripts.vaccine_reference_loader import (
    load_authoritative_references,
)


# -----------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------

from .paths import PROJECT_ROOT

BASE_DIR = str(PROJECT_ROOT)

DATA_DIR = os.path.join(
    BASE_DIR,
    "data/clean"
)

OUT_FILE = os.path.join(
    BASE_DIR,
    "results/vaccine_escape_ai_fast.csv"
)

os.makedirs(
    os.path.join(BASE_DIR, "results"),
    exist_ok=True
)


# -----------------------------------------------------
# COMMAND-LINE OPTIONS
# -----------------------------------------------------

parser = argparse.ArgumentParser(
    description="PRRSV vaccine escape prediction using "
                "authoritative vaccine references."
)

parser.add_argument(
    "--dry-run",
    action="store_true",
    help="Generate predictions without modifying MySQL."
)

parser.add_argument(
    "--write-db",
    action="store_true",
    help="Explicitly write generated predictions to MySQL."
)

args = parser.parse_args()

if args.dry_run and args.write_db:
    raise SystemExit(
        "ERROR: --dry-run and --write-db cannot be used together."
    )


# -----------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------

def sequence_identity(seq1, seq2):
    """
    Simple percentage identity based on globalxx alignment.

    This preserves the original predictor calculation.
    """
    alignments = pairwise2.align.globalxx(
        seq1,
        seq2,
        one_alignment_only=True
    )

    if not alignments:
        raise ValueError("No alignment was produced.")

    aln = alignments[0]

    matches = sum(
        a == b
        for a, b in zip(aln.seqA, aln.seqB)
    )

    if len(aln.seqA) == 0:
        raise ValueError("Alignment has zero length.")

    return 100 * matches / len(aln.seqA)


def predict_escape_probability(identity_score):
    """
    Convert identity percentage to escape probability.

    This is the original predictor equation and is deliberately
    unchanged during the reference-layer correction.
    """
    return 1 / (
        1 + np.exp((identity_score - 85) / -3)
    )


# -----------------------------------------------------
# LOAD AUTHORITATIVE VACCINE REFERENCES
# -----------------------------------------------------

print("=" * 100)
print("AUTHORITATIVE VACCINE ESCAPE PREDICTOR")
print("=" * 100)

print()
print("Loading authoritative vaccine references...")

passed_refs, failed_refs = load_authoritative_references(
    strict=False
)

if not passed_refs:
    raise SystemExit(
        "ERROR: No authoritative vaccine references passed QC."
    )

vaccine_refs = []

for ref in passed_refs:
    vaccine_refs.append(
        (
            f"{ref['vaccine_name']}_{ref['accession']}",
            str(ref["sequence"])
        )
    )

print()
print(
    f"Authoritative references accepted : {len(vaccine_refs)}"
)
print(
    f"References excluded               : {len(failed_refs)}"
)

print()
print("ACCEPTED REFERENCES")

for fasta, sequence in vaccine_refs:
    print(
        f"  PASS | {fasta:40} | "
        f"{len(sequence):6} bp"
    )

print()
print("EXCLUDED REFERENCES")

for ref in failed_refs:
    print(
        f"  FAIL | "
        f"{ref['vaccine_name']:25} | "
        f"{ref['accession']:12} | "
        f"{ref['reason']}"
    )


# -----------------------------------------------------
# LOAD PRRSV CLEANED GENOMES
# -----------------------------------------------------

if not os.path.isdir(DATA_DIR):
    raise SystemExit(
        f"ERROR: PRRSV input directory does not exist: {DATA_DIR}"
    )

prrsv_files = sorted(
    f
    for f in os.listdir(DATA_DIR)
    if f.endswith(".fasta")
)

if not prrsv_files:
    raise SystemExit(
        "ERROR: No cleaned PRRSV FASTA files found."
    )

print()
print(
    f"PRRSV sequences loaded : {len(prrsv_files)}"
)

expected_pairs = (
    len(prrsv_files) *
    len(vaccine_refs)
)

print(
    f"Expected prediction pairs : {expected_pairs}"
)


# -----------------------------------------------------
# MAIN LOOP — COMPUTE ESCAPE PROBABILITIES
# -----------------------------------------------------

results = []

for prrsv_file in tqdm(
    prrsv_files,
    desc="AI vaccine escape prediction"
):

    prrsv_record = next(
        SeqIO.parse(
            os.path.join(DATA_DIR, prrsv_file),
            "fasta"
        )
    )

    isolate_seq = str(
        prrsv_record.seq
    )

    for vaccine_name, vaccine_seq in vaccine_refs:

        # -------------------------------------------------
        # PRESERVED ORIGINAL CALCULATION
        # -------------------------------------------------

        identity = sequence_identity(
            isolate_seq[:1000],
            vaccine_seq[:1000]
        )

        escape_prob = predict_escape_probability(
            identity
        )

        # -------------------------------------------------
        # PRESERVED ORIGINAL RISK THRESHOLDS
        # -------------------------------------------------

        if escape_prob > 0.75:
            risk = "High"

        elif escape_prob > 0.45:
            risk = "Moderate"

        else:
            risk = "Low"

        # -------------------------------------------------
        # STORE RESULT
        # -------------------------------------------------

        results.append(
            {
                "isolate": prrsv_file,
                "vaccine": vaccine_name,
                "match_score_%": round(identity, 2),
                "escape_probability": round(
                    escape_prob,
                    3
                ),
                "predicted_escape_risk": risk,
            }
        )


# -----------------------------------------------------
# VALIDATE OUTPUT BEFORE ANY DATABASE ACTION
# -----------------------------------------------------

if not results:
    raise SystemExit(
        "ERROR: No predictions were generated."
    )

df = pd.DataFrame(results)

df["timestamp"] = datetime.now().strftime(
    "%Y-%m-%d %H:%M:%S"
)

actual_pairs = (
    df[["isolate", "vaccine"]]
    .drop_duplicates()
    .shape[0]
)

print()
print("=" * 100)
print("PREDICTION OUTPUT VALIDATION")
print("=" * 100)

print(
    f"Rows generated          : {len(df)}"
)

print(
    f"Unique isolates         : "
    f"{df['isolate'].nunique()}"
)

print(
    f"Unique vaccines         : "
    f"{df['vaccine'].nunique()}"
)

print(
    f"Unique isolate-vaccine  : "
    f"{actual_pairs}"
)

print(
    f"Expected isolate-vaccine: "
    f"{expected_pairs}"
)

if len(df) != expected_pairs:
    raise SystemExit(
        "ERROR: Generated row count does not match "
        "expected isolate × authoritative-reference count."
    )

if actual_pairs != expected_pairs:
    raise SystemExit(
        "ERROR: Unique isolate-vaccine pair count "
        "does not match expected count."
    )


# -----------------------------------------------------
# SAVE CSV
# -----------------------------------------------------

df.to_csv(
    OUT_FILE,
    index=False
)

print()
print(
    f"CSV saved: {OUT_FILE}"
)


# -----------------------------------------------------
# DATABASE CONTROL
# -----------------------------------------------------

if args.dry_run:

    print()
    print("=" * 100)
    print("DRY-RUN COMPLETE")
    print("=" * 100)
    print()
    print("NO DATABASE WRITE PERFORMED.")
    print()
    print(
        "The generated predictions have been written "
        "to CSV only."
    )

    raise SystemExit(0)


if not args.write_db:

    print()
    print("=" * 100)
    print("DATABASE WRITE NOT REQUESTED")
    print("=" * 100)
    print()
    print(
        "No database modification was performed."
    )
    print(
        "Use --write-db only after independent validation."
    )

    raise SystemExit(0)


# -----------------------------------------------------
# EXPLICIT DATABASE WRITE
# -----------------------------------------------------

from intelliprrsv2.scripts.db_utils import insert_dataframe

print()
print("=" * 100)
print("DATABASE WRITE REQUESTED")
print("=" * 100)

insert_dataframe(
    df,
    "vaccine_escape_predictions"
)

print()
print("Database insertion completed.")
