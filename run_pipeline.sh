#!/bin/bash
# ==========================================================
# 🧠 IntelliPRRSV2 — Complete Real PRRSV Analysis Pipeline
# Author: AI Concepts Limited
# Version: 2026.01
# ==========================================================

# Set up environment
source venv/bin/activate
LOG_DIR="logs"
mkdir -p $LOG_DIR

echo "🚀 Starting IntelliPRRSV2 pipeline at $(date)"
echo "Logs will be saved to $LOG_DIR"

# ----------------------------------------------------------
# STEP 1 — FETCH SEQUENCES + METADATA TO MYSQL
# ----------------------------------------------------------
echo "📡 Step 1: Fetching PRRSV sequences from NCBI..."
python3 scripts/fetch_1000_to_mysql.py > $LOG_DIR/fetch.log 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Step 1 complete: Sequences fetched successfully"
else
    echo "❌ Step 1 failed — check $LOG_DIR/fetch.log"
    exit 1
fi

# ----------------------------------------------------------
# STEP 2 — DATA CLEANING & STANDARDIZATION
# ----------------------------------------------------------
echo "🧹 Step 2: Cleaning and standardizing sequences..."
python3 scripts/clean_sequences_parallel.py > $LOG_DIR/clean.log 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Step 2 complete: Cleaned FASTA files ready"
else
    echo "❌ Step 2 failed — check $LOG_DIR/clean.log"
    exit 1
fi

# ----------------------------------------------------------
# STEP 3 — MUTATION HOTSPOT ANALYSIS
# ----------------------------------------------------------
echo "🔥 Step 3: Detecting mutation hotspots..."
python3 scripts/mutation_hotspot.py > $LOG_DIR/mutation.log 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Step 3 complete: Mutation hotspots identified"
else
    echo "❌ Step 3 failed — check $LOG_DIR/mutation.log"
    exit 1
fi

# ----------------------------------------------------------
# STEP 4 — PHYLOGENETIC TREE & LINEAGE ASSIGNMENT
# ----------------------------------------------------------
echo "🌳 Step 4: Building phylogenetic tree and assigning lineages..."
python3 scripts/build_phylogeny.py > $LOG_DIR/phylo.log 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Step 4 complete: Lineage data generated"
else
    echo "❌ Step 4 failed — check $LOG_DIR/phylo.log"
    exit 1
fi

# ----------------------------------------------------------
# STEP 5 — VACCINE ESCAPE PREDICTION
# ----------------------------------------------------------
echo "💉 Step 5: Computing vaccine escape similarity..."
python3 scripts/vaccine_escape_predictor_fast.py > $LOG_DIR/vaccine.log 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Step 5 complete: Vaccine similarity results ready"
else
    echo "❌ Step 5 failed — check $LOG_DIR/vaccine.log"
    exit 1
fi

# ----------------------------------------------------------
# STEP 6 — MIRNA INTERACTION MAPPING
# ----------------------------------------------------------
echo "🧠 Step 6: Analyzing miRNA–virus binding interactions..."
python3 scripts/mirna_interaction_fast.py > $LOG_DIR/mirna.log 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Step 6 complete: miRNA interaction file created"
else
    echo "❌ Step 6 failed — check $LOG_DIR/mirna.log"
    exit 1
fi

# ----------------------------------------------------------
# STEP 7 — DASHBOARD SETUP
# ----------------------------------------------------------
echo "📊 Step 7: Preparing Streamlit dashboard..."
echo "To launch dashboard, run: streamlit run dashboard/app.py"
echo "✅ Dashboard ready — data linked to results folder."

# ----------------------------------------------------------
# DONE
# ----------------------------------------------------------
echo "🎯 Pipeline completed successfully at $(date)"
deactivate

