import requests
import time
import mysql.connector
from intelliprrsv2.scripts.paths import UID_MAP_LOG
from intelliprrsv2.scripts.config import DB_CONFIG
from datetime import datetime
import logging

# -----------------------------
# 🧠 Logging Configuration
# -----------------------------
LOG_PATH = str(UID_MAP_LOG)
logging.basicConfig(
    filename=LOG_PATH,
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

logging.info("=== Starting UID → Accession Mapping ===")

# -----------------------------
# ⚙️ MySQL Connection
# -----------------------------
try:
    conn = mysql.connector.connect(**DB_CONFIG)
    cursor = conn.cursor()
    logging.info("✅ MySQL connection established successfully.")
except mysql.connector.Error as e:
    logging.error(f"❌ Database connection failed: {e}")
    raise SystemExit(f"Database connection failed: {e}")

# -----------------------------
# 🧬 Fetch distinct UIDs
# -----------------------------
cursor.execute("SELECT isolate FROM temp_uid_batch")
uids = [row[0] for row in cursor.fetchall()]
logging.info(f"Found {len(uids)} unique UIDs to map.")

# Create lookup table if it doesn't exist
cursor.execute("""
    CREATE TABLE IF NOT EXISTS sequence_lookup (
        gi_number VARCHAR(50) PRIMARY KEY,
        accession VARCHAR(50)
    )
""")
conn.commit()

insert_sql = "REPLACE INTO sequence_lookup (gi_number, accession) VALUES (%s, %s)"

# -----------------------------
# 🔄 Fetch accession for each UID
# -----------------------------
for i, uid in enumerate(uids, start=1):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id={uid}&retmode=json"
    success = False

    for attempt in range(3):  # retry 3 times
        try:
            r = requests.get(url, timeout=15)
            if r.status_code == 429:
                logging.warning(f"⚠️ 429 Too Many Requests for UID {uid}. Sleeping 5s...")
                time.sleep(5)
                continue
            r.raise_for_status()

            data = r.json()
            acc = data["result"][uid]["accessionversion"]

            cursor.execute(insert_sql, (uid, acc))
            conn.commit()

            logging.info(f"✅ {uid} → {acc}")
            success = True
            break

        except Exception as e:
            logging.warning(f"⚠️ Attempt {attempt+1}/3 failed for {uid}: {e}")
            time.sleep(2 ** attempt)  # exponential backoff

    if not success:
        logging.error(f"❌ Skipping {uid}: failed after 3 attempts")

    # Gentle pacing for NCBI (max ~3 requests/sec)
    time.sleep(0.34)

# -----------------------------
# ✅ Wrap Up
# -----------------------------
cursor.close()
conn.close()
logging.info("✅ UID → Accession mapping complete.")
logging.info(f"Total UIDs processed: {len(uids)}")
print("✅ UID → Accession mapping complete. Logs saved to", LOG_PATH)

