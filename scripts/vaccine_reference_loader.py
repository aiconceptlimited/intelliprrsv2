#!/usr/bin/env python3

"""
Authoritative PRRSV vaccine-reference loader.

The manifest (reference/vaccines/vaccine_refs.csv) is authoritative.
Only manifest-listed references that pass sequence QC are returned.

This module performs no database operations and does not modify files.
"""

from pathlib import Path
import csv
import re
from Bio import SeqIO


from .paths import PROJECT_ROOT

BASE_DIR = PROJECT_ROOT
REFERENCE_DIR = BASE_DIR / "reference" / "vaccines"
MANIFEST = REFERENCE_DIR / "vaccine_refs.csv"

# Complete PRRSV genome length guard.
# Deliberately broad enough to accommodate known PRRSV genome variation,
# while rejecting obviously unrelated short sequences.
MIN_GENOME_LENGTH = 13000
MAX_GENOME_LENGTH = 17000

VALID_NT = re.compile(r"^[ACGTN]+$", re.IGNORECASE)


class VaccineReferenceError(RuntimeError):
    """Raised when the authoritative reference configuration is invalid."""


def _load_manifest():
    if not MANIFEST.exists():
        raise VaccineReferenceError(
            f"Manifest not found: {MANIFEST}"
        )

    with MANIFEST.open(newline="") as fh:
        rows = list(csv.DictReader(fh))

    required = {
        "vaccine_name",
        "accession",
        "lineage",
        "type",
        "region",
    }

    if not rows:
        raise VaccineReferenceError("Manifest contains no reference entries.")

    missing_columns = required.difference(rows[0].keys())
    if missing_columns:
        raise VaccineReferenceError(
            f"Manifest missing required columns: "
            f"{sorted(missing_columns)}"
        )

    return rows


def _normalise_accession(value):
    return str(value).strip().split(".")[0]


def _validate_record(row):
    name = row["vaccine_name"].strip()
    accession = _normalise_accession(row["accession"])

    fasta = REFERENCE_DIR / f"{name}_{accession}.fasta"

    if not fasta.exists():
        return None, f"FASTA_MISSING:{fasta.name}"

    try:
        records = list(SeqIO.parse(fasta, "fasta"))
    except Exception as exc:
        return None, f"FASTA_PARSE_ERROR:{exc}"

    if len(records) != 1:
        return None, f"EXPECTED_ONE_FASTA_RECORD:found={len(records)}"

    record = records[0]
    fasta_accession = _normalise_accession(record.id)

    if fasta_accession != accession:
        return None, (
            f"ACCESSION_MISMATCH:"
            f"manifest={accession},fasta={fasta_accession}"
        )

    sequence = str(record.seq).upper()

    if not sequence:
        return None, "EMPTY_SEQUENCE"

    if not VALID_NT.fullmatch(sequence):
        return None, "INVALID_NUCLEOTIDE_ALPHABET"

    ambiguous = sequence.count("N")
    if ambiguous:
        return None, f"AMBIGUOUS_BASES:{ambiguous}"

    length = len(sequence)

    if not MIN_GENOME_LENGTH <= length <= MAX_GENOME_LENGTH:
        return None, f"INVALID_GENOME_LENGTH:{length}"

    if not re.search(
        r"porcine reproductive and respiratory syndrome virus|"
        r"porcine respiratory and reproductive syndrome virus|"
        r"PRRSV",
        record.description,
        re.IGNORECASE,
    ):
        return None, "DESCRIPTION_NOT_PRRSV"

    reference = {
        "vaccine_name": name,
        "accession": accession,
        "lineage": row["lineage"].strip(),
        "type": row["type"].strip(),
        "region": row["region"].strip(),
        "sequence": sequence,
        "fasta": str(fasta),
        "length": length,
        "description": record.description,
    }

    return reference, "PASS"


def load_authoritative_references(strict=True):
    """
    Return only manifest-listed references that pass QC.

    strict=True:
        Return PASS references but raise an exception if any manifest
        reference fails QC.

    strict=False:
        Return PASS references while reporting failures.
    """

    rows = _load_manifest()

    passed = []
    failed = []

    print("=" * 100)
    print("AUTHORITATIVE VACCINE REFERENCE LOADER")
    print("=" * 100)
    print(f"Manifest : {MANIFEST}")
    print(f"Entries  : {len(rows)}")
    print()

    for row in rows:
        reference, status = _validate_record(row)

        if reference is not None:
            passed.append(reference)

            print(
                f"PASS | "
                f"{reference['vaccine_name']:25} | "
                f"{reference['accession']:12} | "
                f"{reference['length']:5} bp"
            )
        else:
            failed.append({
                "vaccine_name": row["vaccine_name"].strip(),
                "accession": _normalise_accession(row["accession"]),
                "reason": status,
            })

            print(
                f"FAIL | "
                f"{row['vaccine_name'].strip():25} | "
                f"{_normalise_accession(row['accession']):12} | "
                f"{status}"
            )

    print()
    print("-" * 100)
    print(f"Manifest entries : {len(rows)}")
    print(f"PASS             : {len(passed)}")
    print(f"FAIL             : {len(failed)}")

    if failed:
        print()
        print("Excluded references:")
        for item in failed:
            print(
                f"  - {item['vaccine_name']} "
                f"({item['accession']}): {item['reason']}"
            )

    print("-" * 100)

    if strict and failed:
        raise VaccineReferenceError(
            f"{len(failed)} authoritative reference(s) failed QC."
        )

    return passed, failed


if __name__ == "__main__":
    load_authoritative_references(strict=False)
