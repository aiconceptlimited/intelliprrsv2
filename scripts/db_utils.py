#!/usr/bin/env python3
# ============================================================
# 🧬 IntelliPRRSV2 — Central MySQL Database Utility (Append Mode)
# ============================================================

import pandas as pd
from sqlalchemy import create_engine
from datetime import datetime
from .config import DB_CONFIG

# ============================================================
# DATABASE CONFIGURATION
# ============================================================
ENGINE_URL = (
    f"mysql+mysqlconnector://{DB_CONFIG['user']}:{DB_CONFIG['password']}"
    f"@{DB_CONFIG['host']}/{DB_CONFIG['database']}"
)

# ============================================================
# FUNCTION: insert_dataframe (append instead of overwrite)
# ============================================================
def insert_dataframe(df: pd.DataFrame, table_name: str):
    """
    Inserts a DataFrame into MySQL using append mode.
    Adds a 'run_timestamp' column for tracking pipeline runs.
    """
    if df.empty:
        print(f"⚠️ Skipped insert — DataFrame for `{table_name}` is empty.")
        return

    try:
        engine = create_engine(ENGINE_URL)
        df["run_timestamp"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        df.to_sql(table_name, con=engine, if_exists='append', index=False)
        print(f"✅ {len(df)} new rows appended to `{table_name}` at {datetime.now()}.")

    except Exception as e:
        print(f"❌ MySQL insert failed for `{table_name}`: {e}")

# ============================================================
# UTILITY: execute_query
# ============================================================
def execute_query(query):
    """Run a custom SQL query (use with caution)."""
    try:
        engine = create_engine(ENGINE_URL)
        with engine.connect() as conn:
            conn.execute(query)
        print("✅ Query executed successfully.")
    except Exception as e:
        print(f"❌ SQL query failed: {e}")

# ============================================================
# UTILITY: fetch_dataframe
# ============================================================
def fetch_dataframe(query):
    """Fetches a SQL query result into a pandas DataFrame."""
    try:
        engine = create_engine(ENGINE_URL)
        df = pd.read_sql(query, engine)
        return df
    except Exception as e:
        print(f"❌ Failed to fetch data: {e}")
        return pd.DataFrame()

# ============================================================
# SELF-TEST
# ============================================================
if __name__ == "__main__":
    try:
        engine = create_engine(ENGINE_URL)
        with engine.connect() as conn:
            conn.execute("SELECT 1")
        print("🟢 MySQL connection successful!")
    except Exception as e:
        print(f"🔴 MySQL connection failed: {e}")

