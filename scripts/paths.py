"""
Project-relative paths for IntelliPRRSV2.

All filesystem paths are derived from the repository root so the
project is portable across development, testing, and deployment
environments.
"""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
CLEAN_DIR = DATA_DIR / "clean"
METADATA_DIR = DATA_DIR / "metadata"

REFERENCE_DIR = PROJECT_ROOT / "reference"
LINEAGE_REFERENCE_DIR = REFERENCE_DIR / "lineages"
MIRNA_REFERENCE_FILE = REFERENCE_DIR / "mirna" / "sus_mirna.fa"

RESULTS_DIR = PROJECT_ROOT / "results"
PHYLO_RESULTS_DIR = RESULTS_DIR / "phylo"
MIRNA_RESULTS_DIR = RESULTS_DIR / "mirna"
SUMMARY_RESULTS_DIR = RESULTS_DIR / "summary"

TEMP_DIR = PROJECT_ROOT / "temp"
LOG_DIR = PROJECT_ROOT / "logs"

# Specific files
MUTATION_RESULT_FILE = RESULTS_DIR / "mutation_hotspots_fast.csv"
MIRNA_RESULT_FILE = RESULTS_DIR / "mirna_interactions_real.csv"
CLEAN_METADATA_FILE = DATA_DIR / "clean_metadata.csv"
LINEAGE_RESULT_FILE = PHYLO_RESULTS_DIR / "real_lineage_table.csv"
MIRNA_SUMMARY_INPUT = MIRNA_RESULTS_DIR / "mirna_interactions.csv"
FETCHED_SEQUENCES_DIR = DATA_DIR / "fetched_sequences"
UID_MAP_LOG = LOG_DIR / "uid_map.log"


def ensure_runtime_directories():
    """Create runtime directories when they do not already exist."""
    for directory in (
        RAW_DIR,
        CLEAN_DIR,
        METADATA_DIR,
        RESULTS_DIR,
        PHYLO_RESULTS_DIR,
        MIRNA_RESULTS_DIR,
        SUMMARY_RESULTS_DIR,
        TEMP_DIR,
        LOG_DIR,
    ):
        directory.mkdir(parents=True, exist_ok=True)
