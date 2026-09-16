#!/usr/bin/env python3

"""
Authoritative PRRSV reference validator.

The manifest is the authoritative definition of eligible references.
A FASTA is accepted only when:

1. It is explicitly represented in vaccine_refs.csv.
2. The expected FASTA exists.
3. FASTA accession matches manifest accession.
4. Sequence contains only A/C/G/T/N.
5. Sequence has no ambiguous bases (N).
6. Sequence length is within the expected PRRSV complete-genome range.
7. The FASTA description identifies it as PRRSV.

This script is READ/VALIDATE ONLY.
It does not modify the database or reference files.
"""

from pathlib import Path
import csv
import hashlib
import re
import sys

from Bio import SeqIO


from .paths import PROJECT_ROOT

BASE_DIR = PROJECT_ROOT
REF_DIR = BASE_DIR / "reference" / "vaccines"
MANIFEST = REF_DIR / "vaccine_refs.csv"
OUTPUT_DIR = BASE_DIR / "results"
OUTPUT = OUTPUT_DIR / "vaccine_reference_qc.csv"

# Conservative complete-genome range for the PRRSV references currently
# present in this project.
MIN_GENOME_LENGTH = 14000
MAX_GENOME_LENGTH = 17000

VALID_NT = set("ACGTN")


def normalize_accession(value):
    """Normalize accession by removing version suffix."""
    return value.strip().split(".")[0]


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_manifest():
    if not MANIFEST.exists():
        raise FileNotFoundError(f"Manifest not found: {MANIFEST}")

    with open(MANIFEST, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    required = {
        "vaccine_name",
        "accession",
        "lineage",
        "type",
        "region",
    }

    if not rows:
        raise ValueError("Manifest contains no entries.")

    missing = required - set(rows[0])
    if missing:
        raise ValueError(
            f"Manifest missing required columns: {sorted(missing)}"
        )

    return rows


def validate_record(manifest_row):
    name = manifest_row["vaccine_name"].strip()
    accession = normalize_accession(manifest_row["accession"])

    expected_file = REF_DIR / f"{name}_{accession}.fasta"

    result = {
        "vaccine_name": name,
        "accession": accession,
        "lineage": manifest_row["lineage"].strip(),
        "type": manifest_row["type"].strip(),
        "region": manifest_row["region"].strip(),
        "fasta_path": str(expected_file),
        "fasta_accession": "",
        "length": "",
        "valid_nt": "",
        "n_count": "",
        "accession_match": "",
        "prrsv_description": "",
        "length_qc": "",
        "md5": "",
        "status": "",
        "reason": "",
    }

    if not expected_file.exists():
        result["status"] = "FAIL"
        result["reason"] = "FASTA_MISSING"
        return result

    try:
        records = list(SeqIO.parse(expected_file, "fasta"))
    except Exception as exc:
        result["status"] = "FAIL"
        result["reason"] = f"FASTA_PARSE_ERROR:{exc}"
        return result

    if len(records) != 1:
        result["status"] = "FAIL"
        result["reason"] = f"EXPECTED_ONE_RECORD_FOUND_{len(records)}"
        return result

    rec = records[0]

    fasta_accession = normalize_accession(rec.id)
    sequence = str(rec.seq).upper()
    description = rec.description

    result["fasta_accession"] = fasta_accession
    result["length"] = len(sequence)
    result["valid_nt"] = set(sequence) <= VALID_NT
    result["n_count"] = sequence.count("N")
    result["accession_match"] = fasta_accession == accession
    result["prrsv_description"] = bool(
        re.search(r"porcine reproductive and respiratory syndrome", description, re.I)
        or re.search(r"porcine respiratory and reproductive syndrome", description, re.I)
    )
    result["length_qc"] = (
        MIN_GENOME_LENGTH <= len(sequence) <= MAX_GENOME_LENGTH
    )
    result["md5"] = md5_file(expected_file)

    failures = []

    if not result["accession_match"]:
        failures.append("ACCESSION_MISMATCH")

    if not result["valid_nt"]:
        failures.append("INVALID_NUCLEOTIDE_ALPHABET")

    if result["n_count"] != 0:
        failures.append("AMBIGUOUS_BASES_PRESENT")

    if not result["length_qc"]:
        failures.append("LENGTH_OUTSIDE_COMPLETE_GENOME_RANGE")

    if not result["prrsv_description"]:
        failures.append("DESCRIPTION_NOT_IDENTIFIED_AS_PRRSV")

    if failures:
        result["status"] = "FAIL"
        result["reason"] = ";".join(failures)
    else:
        result["status"] = "PASS"
        result["reason"] = "ALL_REFERENCE_QC_PASSED"

    return result


def main():
    print("=" * 100)
    print("AUTHORITATIVE PRRSV REFERENCE VALIDATION")
    print("=" * 100)
    print(f"Manifest : {MANIFEST}")
    print(f"Reference: {REF_DIR}")
    print()

    manifest = load_manifest()

    print(f"Manifest entries: {len(manifest)}")
    print()

    results = []

    for row in manifest:
        result = validate_record(row)
        results.append(result)

        print(
            f"{result['status']:4} | "
            f"{result['vaccine_name']:25} | "
            f"{result['accession']:12} | "
            f"{str(result['length']):6} | "
            f"{result['reason']}"
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fieldnames = list(results[0].keys())

    with open(OUTPUT, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    passed = sum(r["status"] == "PASS" for r in results)
    failed = sum(r["status"] == "FAIL" for r in results)

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)
    print(f"Manifest references : {len(results)}")
    print(f"PASS                : {passed}")
    print(f"FAIL                : {failed}")
    print(f"QC report           : {OUTPUT}")

    if failed:
        print()
        print("WARNING: One or more manifest references failed QC.")
        print("No database changes were made.")

    print()
    print("Validation complete.")


if __name__ == "__main__":
    main()
