import os
import json
import pandas as pd
from tqdm import tqdm

META_DIR = str(METADATA_DIR)
OUTPUT_FILE = str(CLEAN_METADATA_FILE)

def parse_metadata(file_path):
    """Extract key fields from one GenBank metadata file."""
    with open(file_path, "r") as f:
        data = json.load(f)

    try:
        seq = data["GBSet"]["GBSeq"]
    except KeyError:
        return None

    record = {
        "accession": seq.get("GBSeq_locus", ""),
        "definition": seq.get("GBSeq_definition", ""),
        "length": seq.get("GBSeq_length", ""),
        "organism": seq.get("GBSeq_organism", ""),
        "country": "",
        "collection_date": "",
        "host": ""
    }

    try:
        features = seq["GBSeq_feature-table"]["GBFeature"]
        if isinstance(features, dict):
            features = [features]
        for feature in features:
            if feature.get("GBFeature_key") == "source":
                quals = feature.get("GBFeature_quals", [])
                if isinstance(quals, dict):
                    quals = [quals]
                for q in quals:
                    name = q.get("GBQualifier_name")
                    value = q.get("GBQualifier_value", "")
                    if name in record:
                        record[name] = value
    except Exception:
        pass

    return record

def main():
    files = [f for f in os.listdir(META_DIR) if f.endswith(".json") and f != "accessions.json"]
    print(f"Found {len(files)} metadata files to process.")

    records = []
    for file in tqdm(files):
        path = os.path.join(META_DIR, file)
        result = parse_metadata(path)
        if result:
            records.append(result)

    df = pd.DataFrame(records)
    df["country"] = df["country"].str.strip().str.title()
    df["host"] = df["host"].str.strip().str.capitalize()
    df["collection_date"] = df["collection_date"].str.replace(" ", "-")

    df.to_csv(OUTPUT_FILE, index=False)
    print(f"✅ Clean metadata saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()

