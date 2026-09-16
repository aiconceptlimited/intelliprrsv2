
import os
import pandas as pd
from Bio import SeqIO
from joblib import Parallel, delayed
from tqdm import tqdm

# ----------------------------------------------------------
# PATH CONFIGURATION
# ----------------------------------------------------------
CLEAN_DIR = str(CLEAN_DIR)
META_DIR = str(METADATA_DIR)
os.makedirs(META_DIR, exist_ok=True)
OUTPUT_FILE = os.path.join(META_DIR, "metadata_summary.csv")

# ----------------------------------------------------------
# FUNCTION: Extract metadata from one FASTA file
# ----------------------------------------------------------
def extract_metadata(fasta_file):
    """Extracts accession, organism, date, country, and length from FASTA header."""
    try:
        path = os.path.join(CLEAN_DIR, fasta_file)
        rec = next(SeqIO.parse(path, "fasta"))
        desc = rec.description

        # Accession from filename
        accession = fasta_file.replace(".fasta", "")

        # Extract potential metadata fields
        organism = "Unknown"
        collection_date = "Unknown"
        country = "Unknown"

        # Parse from FASTA description (heuristic)
        desc_parts = desc.split()
        for i, token in enumerate(desc_parts):
            if token.lower() in ["usa", "china", "vietnam", "korea", "germany", "netherlands"]:
                country = token.capitalize()
            if token.isdigit() and len(token) == 4 and token.startswith("20"):
                collection_date = token
            if i == 0:
                organism = desc_parts[0]

        return {
            "accession": accession,
            "organism": organism,
            "collection_date": collection_date,
            "country": country,
            "length": len(rec.seq)
        }

    except Exception as e:
        return {"accession": fasta_file, "error": str(e)}

# ----------------------------------------------------------
# MAIN PIPELINE
# ----------------------------------------------------------
if __name__ == "__main__":
    files = [f for f in os.listdir(CLEAN_DIR) if f.endswith(".fasta")]
    print(f"📦 Extracting metadata from {len(files)} cleaned sequences...")

    # Run parallel extraction
    results = Parallel(n_jobs=8)(
        delayed(extract_metadata)(f) for f in tqdm(files)
    )

    df = pd.DataFrame(results)
    df.to_csv(OUTPUT_FILE, index=False)

    valid = df[df["error"].isna() if "error" in df.columns else [True] * len(df)]
    print(f"✅ Metadata extraction complete — {len(valid)} valid entries saved to {OUTPUT_FILE}")
    print(valid.head(10))

