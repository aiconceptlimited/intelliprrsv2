#!/usr/bin/env python3
"""
🧠 IntelliPRRSV2 — Smart PRRSV Data Fetcher (Full + Auto-Resume)
-----------------------------------------------------------------------
Fetches all PRRSV sequences (complete + fragments) from NCBI GenBank.

✨ Features:
- Backfill mode for small databases
- Incremental mode for ongoing updates
- Auto-resume after interruptions
- Dynamic batch fetching
- Duplicate avoidance (MySQL + JSON)
- Progress checkpointing for safety
-----------------------------------------------------------------------
"""

import os
import requests
import xmltodict
import time
import json
import mysql.connector
import pandas as pd
from tqdm import tqdm
from datetime import datetime, timedelta
from intelliprrsv2.scripts.db_utils import insert_dataframe
from intelliprrsv2.scripts.config import DB_CONFIG
from intelliprrsv2.scripts.paths import RAW_DIR, METADATA_DIR, FETCHED_SEQUENCES_DIR
# ============================================================
# CONFIGURATION
# ============================================================
BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
OUT_RAW = str(RAW_DIR)
OUT_META = str(METADATA_DIR)
OUT_JSON = str(FETCHED_SEQUENCES_DIR)  # 🆕 NEW
ACC_TRACKER = os.path.join(OUT_META, "accessions.json")
RESUME_TRACKER = os.path.join(OUT_META, "fetch_resume.json")

BATCH_SIZE = 500
REQUEST_DELAY = 2.0
SINCE_DAYS = 1825  # Default: 5 years
API_KEY = None  # Optional: add your NCBI API key

MYSQL_CONFIG = DB_CONFIG.copy()

# --- Dynamic fetch limits based on time of day ---
hour = datetime.now().hour
if 0 <= hour < 6:
    MAX_NEW_FETCH = 3000  # 🌙 Night (faster)
elif 6 <= hour < 18:
    MAX_NEW_FETCH = 1000  # 🌞 Daytime (slower)
else:
    MAX_NEW_FETCH = 2000  # 🌆 Evening (medium)
print(f"🕒 Time-based fetch limit → {MAX_NEW_FETCH} sequences this run.")

os.makedirs(OUT_RAW, exist_ok=True)
os.makedirs(OUT_META, exist_ok=True)
os.makedirs(OUT_JSON, exist_ok=True)  # 🆕 NEW

# ============================================================
# LOAD LOCAL TRACKERS
# ============================================================
if os.path.exists(ACC_TRACKER):
    with open(ACC_TRACKER, "r") as f:
        fetched = set(json.load(f))
else:
    fetched = set()

if os.path.exists(RESUME_TRACKER):
    with open(RESUME_TRACKER, "r") as f:
        resume_data = json.load(f)
        resume_offset = resume_data.get("offset", 0)
else:
    resume_offset = 0

# ============================================================
# MYSQL HELPERS
# ============================================================
def get_existing_accessions():
    """Get all existing sequence accessions from MySQL."""
    try:
        conn = mysql.connector.connect(**MYSQL_CONFIG)
        cursor = conn.cursor()
        cursor.execute("SELECT accession FROM sequences")
        existing = {row[0] for row in cursor.fetchall()}
        cursor.close()
        conn.close()
        print(f"✅ Retrieved {len(existing)} existing accessions from MySQL.")
        return existing
    except Exception as e:
        print(f"⚠️ Could not fetch existing accessions from DB: {e}")
        return set()

def get_last_update_date():
    """Auto-switch between backfill and incremental modes."""
    try:
        conn = mysql.connector.connect(**MYSQL_CONFIG)
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*), MAX(run_timestamp) FROM sequences")
        count, last_update = cursor.fetchone()
        cursor.close()
        conn.close()

        if count < 5000:
            forced = (datetime.now() - timedelta(days=1825)).strftime("%Y/%m/%d")
            print(f"🕒 Database small ({count} seqs) → Backfill mode since {forced}")
            return forced

        if last_update:
            last_date = last_update.strftime("%Y/%m/%d")
            print(f"🕒 Incremental mode: fetching updates since {last_date}")
            return last_date
        else:
            fallback = (datetime.now() - timedelta(days=SINCE_DAYS)).strftime("%Y/%m/%d")
            print(f"🕒 No timestamp found, fallback to {fallback}")
            return fallback
    except Exception as e:
        print(f"⚠️ Could not retrieve last update date: {e}")
        fallback = (datetime.now() - timedelta(days=SINCE_DAYS)).strftime("%Y/%m/%d")
        return fallback

# ============================================================
# STEP 1 — SEARCH ACCESSIONS
# ============================================================
def search_prrsv():
    """Search for all PRRSV accessions (complete + partial)."""
    mindate = get_last_update_date()
    all_ids = []
    retstart = 0
    retmax = 1000

    print(f"🔍 Searching NCBI for PRRSV sequences modified since {mindate}...")

    while True:
        params = {
            "db": "nuccore",
            "term": 'PRRSV[Organism] AND ("Porcine reproductive and respiratory syndrome virus"[All Fields])',
            "retmode": "json",
            "retstart": retstart,
            "retmax": retmax,
            "sort": "datemodified"
        }
        if API_KEY:
            params["api_key"] = API_KEY

        try:
            r = requests.get(BASE_URL + "esearch.fcgi", params=params, timeout=30)
            r.raise_for_status()
            data = r.json().get("esearchresult", {})
            ids = data.get("idlist", [])
            if not ids:
                break

            all_ids.extend(ids)
            retstart += retmax
            time.sleep(REQUEST_DELAY)
        except Exception as e:
            print(f"⚠️ Search error at offset {retstart}: {e}")
            break

    print(f"✅ Found {len(all_ids)} PRRSV accessions since {mindate}.")
    return all_ids

# ============================================================
# STEP 2 — FETCH BATCHES (with safe JSON saving)
# ============================================================
def fetch_batch(batch_ids):
    """Fetch a batch of FASTA and XML data, and save JSON copies."""
    ids_str = ",".join(batch_ids)
    fasta_url = BASE_URL + "efetch.fcgi"

    try:
        fasta_params = {"db": "nuccore", "id": ids_str, "rettype": "fasta", "retmode": "text"}
        xml_params = {"db": "nuccore", "id": ids_str, "rettype": "gb", "retmode": "xml"}
        if API_KEY:
            fasta_params["api_key"] = xml_params["api_key"] = API_KEY

        fasta_data = requests.get(fasta_url, params=fasta_params, timeout=60).text
        xml_data = requests.get(fasta_url, params=xml_params, timeout=60).text

        # Save locally (existing logic)
        fasta_file = os.path.join(OUT_RAW, f"batch_{time.time_ns()}.fasta")
        xml_file = os.path.join(OUT_META, f"batch_{time.time_ns()}.xml")
        with open(fasta_file, "w") as f1:
            f1.write(fasta_data)
        with open(xml_file, "w") as f2:
            f2.write(xml_data)

        # 🆕 NEW: Safe JSON saving for geo metadata
        try:
            if xml_data and xml_data.strip():  # ensure non-empty XML
                parsed = xmltodict.parse(xml_data)
                gbset = parsed.get("GBSet", {}).get("GBSeq", [])
                if isinstance(gbset, dict):  # handle single-record XMLs
                    gbset = [gbset]

                saved = 0
                for record in gbset:
                    uid = record.get("GBSeq_primary-accession", str(time.time_ns()))
                    json_path = os.path.join(OUT_JSON, f"{uid}.json")
                    with open(json_path, "w") as jf:
                        json.dump(record, jf)
                    saved += 1

                print(f"💾 Saved {saved} JSON records for geo metadata extraction.")
            else:
                print("⚠️ Empty XML response — skipping JSON save for this batch.")
        except Exception as e:
            print(f"⚠️ Could not parse/save JSON batch: {e}")

        return len(batch_ids)
    except Exception as e:
        print(f"⚠️ Batch fetch failed: {e}")
        return 0

# ============================================================
# MAIN EXECUTION
# ============================================================
def main():
    global resume_offset

    print("🧠 IntelliPRRSV2 — Batch Fetch Layer (Auto-Resume Enabled)")
    print("=" * 60)

    existing = get_existing_accessions()
    ids = search_prrsv()

    if resume_offset > 0:
        print(f"🔁 Resuming from offset {resume_offset}")
        ids = ids[resume_offset:]
    else:
        print("🚀 Starting fresh fetch session.")

    new_ids = [i for i in ids if i not in existing and i not in fetched][:MAX_NEW_FETCH]
    print(f"🧬 New sequences to fetch: {len(new_ids)}")

    if not new_ids:
        print("⚠️ No new PRRSV accessions found — database already up to date.")
        return

    total_fetched = 0
    for i in tqdm(range(0, len(new_ids), BATCH_SIZE), desc="Fetching PRRSV batches"):
        batch = new_ids[i:i + BATCH_SIZE]
        total_fetched += fetch_batch(batch)

        resume_offset += len(batch)
        with open(RESUME_TRACKER, "w") as f:
            json.dump({"offset": resume_offset, "timestamp": str(datetime.now())}, f, indent=2)

        time.sleep(REQUEST_DELAY)

    print(f"✅ Fetch complete — {total_fetched} sequences fetched this run.")

    fetched.update(new_ids)
    with open(ACC_TRACKER, "w") as f:
        json.dump(list(fetched), f, indent=2)

    df = pd.DataFrame({
        "accession": new_ids,
        "timestamp": [datetime.now().strftime("%Y-%m-%d %H:%M:%S")] * len(new_ids),
        "run_timestamp": [datetime.now().strftime("%Y-%m-%d %H:%M:%S")] * len(new_ids)
    })

    try:
        insert_dataframe(df, "sequences")
        print(f"✅ Inserted {len(df)} new records into MySQL.")
    except Exception as e:
        print(f"❌ MySQL insert failed: {e}")

    print("=" * 60)
    print(f"🏁 Run complete — {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

# ============================================================
# ENTRY POINT
# ============================================================
if __name__ == "__main__":
    main()

