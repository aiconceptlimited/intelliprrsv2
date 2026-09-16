#!/usr/bin/env python3
# ============================================================
# IntelliPRRSV2 — Robust Geographic Metadata Extractor (v1.1)
# ============================================================

import re
import pandas as pd
from datetime import datetime
from sqlalchemy import create_engine, text
from .config import DB_CONFIG

# ------------------------------------------------------------
# DB CONFIG
# ------------------------------------------------------------
ENGINE_URL = (
    f"mysql+mysqlconnector://{DB_CONFIG['user']}:{DB_CONFIG['password']}"
    f"@{DB_CONFIG['host']}/{DB_CONFIG['database']}"
)

engine = create_engine(ENGINE_URL)

# ------------------------------------------------------------
# AUTHORITATIVE RUN TIMESTAMP (CRITICAL FIX)
# ------------------------------------------------------------
RUN_TIMESTAMP = datetime.now()

# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

COUNTRY_NORMALIZE = {
    "united states": "USA",
    "united states of america": "USA",
    "u s a": "USA",
    "u.s.a": "USA",
    "us": "USA",
    "uk": "UK",
    "england": "UK",
    "scotland": "UK",
    "wales": "UK",
    "russia": "Russia",
    "peoples republic of china": "China",
    "pr china": "China",
    "south korea": "South Korea",
    "republic of korea": "South Korea",
    "korea": "South Korea",
}

def _clean_text(x):
    if not x:
        return ""
    x = str(x).strip()
    x = re.sub(r"\s+", " ", x)
    return x

def normalize_country(country):
    c = _clean_text(country)
    if not c:
        return "Unknown"
    key = re.sub(r"[.,]", " ", c.lower())
    key = re.sub(r"\s+", " ", key).strip()
    return COUNTRY_NORMALIZE.get(key, c)

def split_country_region(raw):
    raw = _clean_text(raw)
    if not raw:
        return ("Unknown", None)

    raw = raw.strip('"').strip("'")

    if ":" in raw:
        c, r = raw.split(":", 1)
        return (normalize_country(c), r.strip())

    if "," in raw:
        c, r = raw.split(",", 1)
        return (normalize_country(c), r.strip())

    return (normalize_country(raw), None)

def try_extract_lat_lon(text_blob):
    if not text_blob:
        return (None, None)

    m = re.search(r"(-?\d+\.\d+)\s*,\s*(-?\d+\.\d+)", str(text_blob))
    if m:
        return float(m.group(1)), float(m.group(2))

    return (None, None)

def extract_country_from_isolate(isolate):
    isolate = _clean_text(isolate)
    if not isolate:
        return ("Unknown", None)

    patterns = [
        r"\bUSA\b", r"\bUnited States\b", r"\bChina\b",
        r"\bUK\b", r"\bNetherlands\b", r"\bSpain\b", r"\bBrazil\b"
    ]

    for p in patterns:
        m = re.search(p, isolate, flags=re.IGNORECASE)
        if m:
            return (normalize_country(m.group(0)), None)

    return ("Unknown", None)

# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

def main():
    print("🗺️  Extracting geo metadata (robust mode)...")
    print(f"🕒 Run timestamp: {RUN_TIMESTAMP}")

    with engine.connect() as conn:
        df = pd.read_sql(
            text("""
                SELECT id, accession, country, region,
                       latitude, longitude, isolate
                FROM geo_metadata
                ORDER BY id DESC
                LIMIT 200000
            """),
            conn
        )

    if df.empty:
        print("⚠️ geo_metadata is empty.")
        return

    fixed = []

    for _, row in df.iterrows():
        country = row["country"]
        region = row["region"]
        lat = row["latitude"]
        lon = row["longitude"]

        # Country / region normalization
        if country and str(country).lower() not in ["unknown", "none", "nan"]:
            c, r = split_country_region(country)
            country = c
            if not region and r:
                region = r
        else:
            c2, r2 = extract_country_from_isolate(row["isolate"])
            country = c2
            if not region:
                region = r2

        # Lat/Lon fallback
        if lat is None or lon is None:
            lat2, lon2 = try_extract_lat_lon(row["isolate"])
            lat = lat if lat is not None else lat2
            lon = lon if lon is not None else lon2

        fixed.append({
            "id": row["id"],
            "country": country,
            "region": region,
            "latitude": lat,
            "longitude": lon,
        })

    print("✅ Writing geo corrections back into MySQL...")

    with engine.begin() as conn:
        for r in fixed:
            conn.execute(
                text("""
                    UPDATE geo_metadata
                    SET country = :country,
                        region = :region,
                        latitude = :latitude,
                        longitude = :longitude,
                        run_timestamp = :run_ts
                    WHERE id = :id
                """),
                {
                    **r,
                    "run_ts": RUN_TIMESTAMP
                }
            )

    print("✅ Geo metadata updated with new run_timestamp")

    with engine.connect() as conn:
        ts = conn.execute(
            text("SELECT MAX(run_timestamp) FROM geo_metadata")
        ).scalar()

    print(f"📊 geo_metadata latest timestamp: {ts}")
    print("🎉 Geo metadata refresh complete.")

if __name__ == "__main__":
    main()

