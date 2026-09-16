import os
from Bio import SeqIO, pairwise2
import pandas as pd
from datetime import datetime

CLEAN_DIR = str(CLEAN_DIR)
REF_DIR = str(LINEAGE_REFERENCE_DIR)
OUT_TABLE = str(LINEAGE_RESULT_FILE)

def load_refs():
    refs = []
    for f in os.listdir(REF_DIR):
        if f.endswith(".fasta"):
            lineage = f.split("_")[1]  # e.g. Lineage_1A_KF501397.fasta → 1A
            seq = str(next(SeqIO.parse(os.path.join(REF_DIR, f), "fasta")).seq)
            refs.append((lineage, seq))
    return refs

def compute_identity(seq1, seq2):
    score = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True, score_only=True)
    return score / max(len(seq1), len(seq2))

def main():
    refs = load_refs()
    results = []

    print(f"🔍 Loaded {len(refs)} reference lineages.")

    for f in os.listdir(CLEAN_DIR):
        if not f.endswith(".fasta"):
            continue
        record = next(SeqIO.parse(os.path.join(CLEAN_DIR, f), "fasta"))
        seq = str(record.seq)
        best = ("Unknown", 0)
        for lineage, refseq in refs:
            identity = compute_identity(seq, refseq)
            if identity > best[1]:
                best = (lineage, identity)
        results.append({
            "accession": f.replace(".fasta", ""),
            "lineage": best[0],
            "similarity(%)": round(best[1] * 100, 2)
        })
        print(f"✅ {f}: {best[0]} ({round(best[1] * 100, 2)}%)")

    df = pd.DataFrame(results)
    df["timestamp"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    os.makedirs(os.path.dirname(OUT_TABLE), exist_ok=True)
    df.to_csv(OUT_TABLE, index=False)
    print(f"\n📁 Saved lineage assignments to {OUT_TABLE}")

if __name__ == "__main__":
    main()

