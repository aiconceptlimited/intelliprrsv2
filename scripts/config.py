"""
Runtime configuration for IntelliPRRSV2.

Database credentials and deployment-specific settings are supplied
through environment variables and are never stored in source code.
"""

import os


DB_HOST = os.getenv("INTELLIPRRSV2_DB_HOST", "localhost")
DB_PORT = int(os.getenv("INTELLIPRRSV2_DB_PORT", "3306"))
DB_USER = os.getenv("INTELLIPRRSV2_DB_USER", "")
DB_PASSWORD = os.getenv("INTELLIPRRSV2_DB_PASSWORD", "")
DB_NAME = os.getenv("INTELLIPRRSV2_DB_NAME", "")

NCBI_API_KEY = os.getenv("NCBI_API_KEY") or None


DB_CONFIG = {
    "host": DB_HOST,
    "port": DB_PORT,
    "user": DB_USER,
    "password": DB_PASSWORD,
    "database": DB_NAME,
}


def get_db_config():
    return DB_CONFIG.copy()


def validate_db_config():
    missing = []

    if not DB_USER:
        missing.append("INTELLIPRRSV2_DB_USER")
    if not DB_PASSWORD:
        missing.append("INTELLIPRRSV2_DB_PASSWORD")
    if not DB_NAME:
        missing.append("INTELLIPRRSV2_DB_NAME")

    if missing:
        raise RuntimeError(
            "Missing required database environment variables: "
            + ", ".join(missing)
        )
