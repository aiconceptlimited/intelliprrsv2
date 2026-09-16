import os
from Bio import SeqIO
from joblib import Parallel, delayed
from tqdm import tqdm
from .paths import RAW_DIR, CLEAN_DIR

RAW_DIR = str(RAW_DIR)
CLEAN_DIR = str(CLEAN_DIR)
os.makedirs(CLEAN_DIR, exist_ok=True)

# ----------------------------------------
# Cleaning rules
# ----------------------------------------
def clean_file(file):
    if not file.endswith(".fasta"):
        return None

    file_path = os.path.join(RAW_DIR, file)
    try:
        seq_record = next(SeqIO.parse(file_path, "fasta"))
        seq = str(seq_record.seq).upper()

        # Rule 1: Minimum length
        if len(seq) < 1000:
            return f"❌ Too short: {file}"

        # Rule 2: Remove ambiguous bases
        seq = seq.replace("N", "")

        # Rule 3: Save the cleaned sequence
        seq_record.seq = type(seq_record.seq)(seq)
        clean_path = os.path.join(CLEAN_DIR, file)
        SeqIO.write(seq_record, clean_path, "fasta")
        return f"✅ Cleaned: {file}"
    except Exception as e:
        return f"⚠️ {file} failed: {e}"

# ----------------------------------------
# Run parallel cleaning
# ----------------------------------------
def main():
    files = [f for f in os.listdir(RAW_DIR) if f.endswith(".fasta")]
    print(f"🧹 Cleaning {len(files)} sequences...")

    results = Parallel(n_jobs=8)(
        delayed(clean_file)(f) for f in tqdm(files, desc="Processing")
    )

    for r in results:
        if r:
            print(r)

    print(f"✅ Cleaning complete — results saved in {CLEAN_DIR}/")


if __name__ == "__main__":
    main()

