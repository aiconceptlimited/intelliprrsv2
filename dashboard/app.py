#!/usr/bin/env python3
# ============================================================
# 🌌 IntelliPRRSV2 — Neon Bioinformatics Dashboard (Optimized v2.4)
# ============================================================

from pathlib import Path
import pandas as pd
import numpy as np
import streamlit as st
from scripts.config import DB_HOST, DB_USER, DB_PASSWORD, DB_NAME


import plotly.express as px
import plotly.graph_objects as go
from sqlalchemy import create_engine, text
from io import StringIO
from Bio import Phylo
import matplotlib.pyplot as plt
import time
from datetime import datetime
import json
import networkx as nx

# ============================================================
# ⚙️ DATABASE CONFIGURATION
# ============================================================
DB_CONFIG = {
    "host": DB_HOST,
    "user": DB_USER,
    "password": DB_PASSWORD,
    "database": DB_NAME,
}

ENGINE_URL = (
    f"mysql+mysqlconnector://{DB_CONFIG['user']}:{DB_CONFIG['password']}"
    f"@{DB_CONFIG['host']}/{DB_CONFIG['database']}"
)
engine = create_engine(ENGINE_URL)
# ============================================================
# 🧬 PIPELINE HEALTH — AUTHORITATIVE SOURCE
# ============================================================

def get_pipeline_health():
    with engine.connect() as conn:
        return conn.execute(
            text("""
                SELECT run_id, start_time, status
                FROM pipeline_runs
                WHERE status = 'completed'
                ORDER BY run_id DESC
                LIMIT 1
            """)
        ).fetchone()

PIPELINE_RUN = get_pipeline_health()
PIPELINE_RUN_ID = PIPELINE_RUN[0]
PIPELINE_RUN_TS = PIPELINE_RUN[1]


# ============================================================
# 🕒 PIPELINE RUN TIMESTAMP (AUTHORITATIVE)
# ============================================================

# (Removed PIPELINE_RUN_TS / PIPELINE_RUN_STR logic)


# ============================================================
# 🧭 STREAMLIT PAGE SETUP
# ============================================================

st.set_page_config(
    page_title="🌌 IntelliPRRSV2 Dashboard",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ✅ DEBUG BUILD (safe location)
st.write("✅ DEBUG BUILD LOADED:", datetime.now())

# ============================================================
# 🎨 CUSTOM NEON CSS
# ============================================================
st.markdown("""
<style>
body {
    background-color: #0d1117;
    color: #e6edf3;
    font-family: 'Inter', sans-serif;
}
.block-container {
    padding-top: 1.5rem;
    max-width: 1300px;
    margin: auto;
}
h1 { color: #00f0ff; text-shadow: 0 0 15px #00f0ff; text-align: center; }
h2, h3 { color: #00e6b8; }
.stTabs [data-baseweb="tab-list"] {
    justify-content: center; flex-wrap: wrap; background-color: #0d1117;
    padding: 6px 0; border-bottom: 2px solid #00f0ff30;
}
.stTabs [data-baseweb="tab"] {
    font-weight: 600; color: #fff; background-color: #161b22;
    border-radius: 8px 8px 0 0; margin: 0 5px; padding: 8px 14px;
    box-shadow: 0 0 8px #00f0ff40;
}
.stTabs [aria-selected="true"] {
    background-color: #00f0ff !important; color: #000 !important;
    box-shadow: 0 0 15px #00f0ff;
}
</style>
""", unsafe_allow_html=True)
# ============================================================
# 🧠 HELPER FUNCTIONS — PIPELINE-AWARE LOADER
# ============================================================

@st.cache_data(ttl=300, show_spinner=False)
def load_from_mysql(table_name: str):
    """
    Load ONLY data from the latest pipeline run (snapshot consistency).
    Handles:
      - geo_metadata: avoids NULL run_timestamp breaking latest-run selection
      - mirna_interactions: uses 'timestamp' instead of run_timestamp
      - phylogenetic_tree: no run_timestamp column, uses latest id
    """
    try:
        limit_map = {
            "sequences": 100000,
            "mutation_hotspots": 100000,
            "vaccine_escape_predictions": 50000,
            "lineage_assignments": 100000,
            "mirna_interactions": 50000,
            "geo_metadata": 20000,
            "phylogenetic_tree": 5,
        }
        limit = limit_map.get(table_name, 50000)

        # Choose timestamp column
        if table_name in ("mirna_interactions", "vaccine_escape_predictions"):
            ts_col = "timestamp"
        elif table_name == "phylogenetic_tree":
            ts_col = None
        else:
            ts_col = "run_timestamp"

        # Build query
        if ts_col is None:
            # phylogenetic_tree: no run_timestamp
            query = f"""
                SELECT *
                FROM `{table_name}`
                ORDER BY id DESC
                LIMIT {limit}
            """

        elif table_name == "geo_metadata":
            # geo_metadata: protect against NULL run_timestamp
            query = f"""
                SELECT *
                FROM `{table_name}`
                WHERE run_timestamp IS NOT NULL
                  AND run_timestamp = (
                      SELECT MAX(run_timestamp)
                      FROM `{table_name}`
                      WHERE run_timestamp IS NOT NULL
                  )
                ORDER BY run_timestamp DESC
                LIMIT {limit}
            """

        else:
            # Default latest-run snapshot load
            query = f"""
                SELECT *
                FROM `{table_name}`
                WHERE `{ts_col}` IS NOT NULL
                  AND `{ts_col}` = (
                      SELECT MAX(`{ts_col}`)
                      FROM `{table_name}`
                      WHERE `{ts_col}` IS NOT NULL
                  )
                ORDER BY `{ts_col}` DESC
                LIMIT {limit}
            """

        # Execute query
        with engine.connect() as conn:
            df = pd.read_sql(query, conn)

        return df

    except Exception:
        return pd.DataFrame()


def test_mysql_connection():
    """
    Test MySQL connectivity and return latency in milliseconds.
    """
    try:
        start = time.time()
        with engine.connect() as conn:
            conn.execute(text("SELECT 1"))
        return round((time.time() - start) * 1000, 2)
    except Exception:
        return None


def format_latency(ms):
    """
    Format latency for dashboard display.
    """
    if ms is None:
        return "❌ No Connection"
    if ms < 100:
        return f"🟢 {ms} ms"
    if ms < 500:
        return f"🟡 {ms} ms"
# ============================================================
# 🕒 TABLE FRESHNESS (READ-ONLY METADATA)
# ============================================================

def get_table_freshness(table_name):
    with engine.connect() as conn:
        ts_col = "timestamp" if table_name == "vaccine_escape_predictions" else "run_timestamp"
        return conn.execute(
            text(f"SELECT MAX({ts_col}) FROM {table_name}")
        ).scalar()


# ============================================================
# 🧩 SAFE PLOTTING WRAPPER (Fixes serialization + width warnings)
# ============================================================
def safe_plotly(fig):
    """Safely render Plotly charts with dtype + period fixes."""
    import pandas as pd

    for trace in fig.data:
        # Normalize Period and Timestamp objects for Plotly
        for attr in ["x", "y"]:
            if hasattr(trace, attr):
                data = getattr(trace, attr)
                if isinstance(data, pd.Series) or isinstance(data, list):
                    data = pd.Series(data)
                    if pd.api.types.is_period_dtype(data):
                        data = data.astype(str)
                    elif pd.api.types.is_datetime64_any_dtype(data):
                        data = data.dt.strftime("%Y-%m-%d")
                    setattr(trace, attr, data.tolist())
    st.plotly_chart(fig, width="stretch")


# ============================================================
# 🧩 SAFE PLOTLY WRAPPER — replaces deprecated use_container_width
# ============================================================
def safe_plotly(fig, width="stretch"):
    """
    Render Plotly charts safely and JSON-compatibly.
    - Converts Period/Timestamp objects automatically
    - Uses new width argument instead of use_container_width
    """
    import plotly.io as pio

    try:
        _ = json.loads(pio.to_json(fig, validate=False))
        st.plotly_chart(fig, width=width)

    except TypeError as e:
        if "Period" in str(e):
            st.warning("Auto-fixing non-serializable Period data → string conversion")
            for trace in fig.data:
                if hasattr(trace, "x"):
                    trace.x = [str(v) for v in trace.x]
            st.plotly_chart(fig, width=width)
        else:
            st.error(f"❌ Plotly rendering error: {e}")
# ============================================================
# ✅ DATABASE VALIDATION
# ============================================================
required_tables = [
    "sequences", "mutation_hotspots", "vaccine_escape_predictions",
    "lineage_assignments",
    "phylogenetic_tree", "geo_metadata"
]
with engine.connect() as conn:
    tables = set(pd.read_sql("SHOW TABLES", conn).iloc[:, 0].tolist())
missing = [t for t in required_tables if t not in tables]
if missing:
    st.warning(f"⚠️ Missing tables: {', '.join(missing)}")
else:
    pass
# ============================================================
# 🌌 HEADER METRICS
# ============================================================
st.markdown("<h1>🧬 IntelliPRRSV2 Intelligence Dashboard</h1>",
            unsafe_allow_html=True)
st.markdown("---")

col1, col2, col3, col4 = st.columns(4)
latency = test_mysql_connection()
col3.metric("📡 Status", "🟢 Online" if latency else "🔴 Offline")
col4.metric("⚡ Latency", format_latency(latency))

if latency:
    with engine.connect() as conn:
        seqs = pd.read_sql("SELECT COUNT(*) c FROM sequences", conn).iloc[0, 0]
        preds = pd.read_sql("SELECT COUNT(*) c FROM vaccine_escape_predictions", conn).iloc[0, 0]
    col1.metric("🧬 Sequences", f"{seqs:,}")
    col2.metric("💉 Escape Predictions", f"{preds:,}")

st.markdown("---")
# ============================================================
# 🧭 NAVIGATION
# ============================================================
tabs = st.tabs([
    "🏠 Overview", 
    "🌿 Lineage Trends", 
    "🧬 Mutation Hotspots",
    "💉 Vaccine Escape", 
    "🧠 miRNA Interactions",
    "🌳 Phylogenetic Tree", 
    "🗺️ Geographical Map"
])

# ============================================================
# 🏠 TAB 1 — OVERVIEW
# ============================================================

with tabs[0]:

    # --------------------------------------------------------
    # SECTION HEADER
    # --------------------------------------------------------
    st.markdown(
        """
        <div style="margin-bottom:0.6rem;">
            <h2 style="margin-bottom:0.15rem;">
                📊 PRRSV Genomic Overview
            </h2>
            <div style="color:#9aa4b2; font-size:0.9rem;">
                Current validated surveillance snapshot
                · Alignment-based mutation analysis
            </div>
        </div>
        """,
        unsafe_allow_html=True
    )

    # --------------------------------------------------------
    # LOAD LATEST ANALYTICAL SNAPSHOTS
    # --------------------------------------------------------
    df_mut = load_from_mysql("mutation_hotspots")
    df_lin = load_from_mysql("lineage_assignments")
    df_esc = load_from_mysql("vaccine_escape_predictions")

    # --------------------------------------------------------
    # NORMALIZE DATA
    # --------------------------------------------------------
    if not df_mut.empty:
        df_mut["position"] = pd.to_numeric(
            df_mut["position"], errors="coerce"
        )
        df_mut["mutation_frequency"] = pd.to_numeric(
            df_mut["mutation_frequency"], errors="coerce"
        )
        df_mut = df_mut.dropna(
            subset=["position", "mutation_frequency"]
        )

    if not df_lin.empty and "assigned_lineage" in df_lin.columns:
        df_lin["assigned_lineage"] = (
            df_lin["assigned_lineage"]
            .astype(str)
            .str.strip()
        )

    if (
        not df_esc.empty
        and "escape_probability" in df_esc.columns
    ):
        df_esc["escape_probability"] = pd.to_numeric(
            df_esc["escape_probability"],
            errors="coerce"
        )

    # --------------------------------------------------------
    # VALIDATED LINEAGE SNAPSHOT
    # --------------------------------------------------------
    if not df_lin.empty and "assigned_lineage" in df_lin.columns:

        lineage_counts = (
            df_lin["assigned_lineage"]
            .value_counts()
            .rename_axis("Lineage")
            .reset_index(name="Count")
        )

        lineage_counts["Percentage"] = (
            lineage_counts["Count"]
            / lineage_counts["Count"].sum()
            * 100
        )

        if "accession" in df_lin.columns:
            validated_isolates = int(
                df_lin["accession"].nunique()
            )
        else:
            validated_isolates = int(len(df_lin))

    else:
        lineage_counts = pd.DataFrame(
            columns=["Lineage", "Count", "Percentage"]
        )
        validated_isolates = 0

    # --------------------------------------------------------
    # MUTATION METRICS
    # --------------------------------------------------------
    mutation_records = int(len(df_mut))

    analyzed_positions = (
        int(df_mut["position"].nunique())
        if not df_mut.empty
        else 0
    )

    if not df_mut.empty:

        # One position-level value for the Overview.
        position_frequency = (
            df_mut.groupby("position", as_index=False)
            ["mutation_frequency"]
            .max()
            .sort_values("position")
        )

        mean_mutation_frequency = float(
            position_frequency["mutation_frequency"].mean()
        )

    else:

        position_frequency = pd.DataFrame(
            columns=["position", "mutation_frequency"]
        )

        mean_mutation_frequency = np.nan

    # --------------------------------------------------------
    # VACCINE ESCAPE METRIC
    # --------------------------------------------------------
    if (
        not df_esc.empty
        and "escape_probability" in df_esc.columns
    ):
        mean_escape = float(
            df_esc["escape_probability"].mean()
        )
    else:
        mean_escape = np.nan

    # ========================================================
    # PRIMARY ANALYTICAL KPIs
    # ========================================================
    st.markdown("### Surveillance Snapshot")

    k1, k2, k3, k4 = st.columns(4)

    k1.metric(
        "🧬 Validated Isolates",
        f"{validated_isolates:,}"
    )

    k2.metric(
        "📍 Alignment Positions",
        f"{analyzed_positions:,}"
    )

    k3.metric(
        "🔬 Mutation Records",
        f"{mutation_records:,}"
    )

    k4.metric(
        "💉 Mean Escape Score",
        (
            f"{mean_escape:.3f}"
            if not np.isnan(mean_escape)
            else "N/A"
        )
    )

    st.markdown("<div style='height:0.5rem'></div>", unsafe_allow_html=True)

    # ========================================================
    # MAIN ANALYTICAL VISUALS
    # ========================================================

    # --------------------------------------------------------
    # LINEAGE DISTRIBUTION — FULL WIDTH
    # --------------------------------------------------------
    st.markdown("### 🌿 Current Lineage Distribution")
    st.caption(
        "Validated PRRSV-2 lineage composition within the current surveillance snapshot"
    )

    if not lineage_counts.empty:

        # Keep the biologically relevant L1–L5 ordering while
        # retaining any additional categories if present.
        canonical_order = ["L1", "L2", "L3", "L4", "L5"]

        lineage_plot = lineage_counts.copy()

        lineage_plot["Lineage"] = lineage_plot["Lineage"].astype(str)

        known = [
            x for x in canonical_order
            if x in lineage_plot["Lineage"].values
        ]

        extra = [
            x for x in lineage_plot["Lineage"].tolist()
            if x not in canonical_order
        ]

        ordered = known + extra

        lineage_plot["Lineage"] = pd.Categorical(
            lineage_plot["Lineage"],
            categories=ordered,
            ordered=True
        )

        lineage_plot = lineage_plot.sort_values("Lineage")

        fig_lineage = px.bar(
            lineage_plot,
            x="Count",
            y="Lineage",
            orientation="h",
            text="Percentage",
            template="plotly_dark"
        )

        fig_lineage.update_traces(
            texttemplate="%{text:.1f}%",
            textposition="outside",
            cliponaxis=False,
            marker=dict(
                line=dict(width=0)
            ),
            hovertemplate=(
                "<b>%{y}</b><br>"
                "Validated isolates: %{x:,}<br>"
                "Proportion: %{text:.2f}%"
                "<extra></extra>"
            )
        )

        max_count = int(lineage_plot["Count"].max())
        x_upper = max_count * 1.18 if max_count > 0 else 1

        fig_lineage.update_layout(
            height=390,
            showlegend=False,
            margin=dict(
                l=65,
                r=95,
                t=25,
                b=60
            ),
            xaxis=dict(
                title="Validated Isolates",
                range=[0, x_upper],
                showgrid=True,
                zeroline=False,
                tickformat=","
            ),
            yaxis=dict(
                title=None,
                categoryorder="array",
                categoryarray=ordered
            ),
            hoverlabel=dict(
                bgcolor="#161b22"
            )
        )

        safe_plotly(fig_lineage)

        # Compact analytical interpretation
        dominant = lineage_plot.loc[
            lineage_plot["Count"].idxmax()
        ]

        st.info(
            f"**Dominant lineage:** {dominant['Lineage']} — "
            f"{int(dominant['Count']):,} validated isolates "
            f"({float(dominant['Percentage']):.1f}% of the current snapshot)."
        )

        st.caption(
            "Lineage assignment uses the current L1–L5 PRRSV-2 "
            "nomenclature and should not be interpreted as a "
            "genotype-independent classifier."
        )

    else:
        st.info("No validated lineage data available.")


    st.markdown("---")


    # --------------------------------------------------------
    # MUTATION-FREQUENCY LANDSCAPE — FULL WIDTH
    # --------------------------------------------------------
    st.markdown("### 🧬 Mutation-Frequency Landscape")
    st.caption(
        "Position-level variation summarized across the "
        "PRRSV multiple-sequence alignment"
    )

    if not position_frequency.empty:

        # ----------------------------------------------------
        # Aggregate into 100-alignment-position bins for the
        # Overview display. Underlying mutation records remain
        # unchanged.
        # ----------------------------------------------------
        overview_bins = 100

        plot_df = position_frequency.copy()

        plot_df["bin"] = (
            plot_df["position"] // overview_bins
        ) * overview_bins

        binned = (
            plot_df.groupby("bin", as_index=False)
            ["mutation_frequency"]
            .median()
        )

        binned["alignment_coordinate"] = (
            binned["bin"] + overview_bins / 2
        )

        fig_mutation = go.Figure()

        fig_mutation.add_trace(
            go.Scatter(
                x=binned["alignment_coordinate"],
                y=binned["mutation_frequency"],
                mode="lines",
                name="100-nt median",
                line=dict(width=3),
                hovertemplate=(
                    "Alignment coordinate: %{x:,.0f}<br>"
                    "Median mutation frequency: %{y:.3f}"
                    "<extra></extra>"
                )
            )
        )

        fig_mutation.add_hline(
            y=mean_mutation_frequency,
            line_dash="dash",
            annotation_text=(
                f"Overall mean = {mean_mutation_frequency:.3f}"
            ),
            annotation_position="top right"
        )

        fig_mutation.update_layout(
            height=460,
            template="plotly_dark",
            margin=dict(
                l=60,
                r=40,
                t=25,
                b=60
            ),
            hovermode="x",
            xaxis=dict(
                title="Alignment Coordinate (nt)",
                showgrid=True,
                zeroline=False
            ),
            yaxis=dict(
                title="Mutation Frequency",
                range=[0, 1],
                showgrid=True,
                zeroline=False
            ),
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.01,
                xanchor="left",
                x=0
            ),
            hoverlabel=dict(
                bgcolor="#161b22"
            )
        )

        safe_plotly(fig_mutation)

        st.caption(
            "The Overview uses 100-alignment-position median bins "
            "to improve readability. Mutation records in the "
            "underlying database are not modified."
        )

    else:
        st.info("No mutation data available.")


    # COMPACT ANALYTICAL SUMMARY
        # ========================================================
        st.markdown("---")
        st.markdown("### 📋 Analytical Summary")

        a1, a2, a3 = st.columns(3)

        a1.metric(
            "Mean Mutation Frequency",
            (
                f"{mean_mutation_frequency:.3f}"
                if not np.isnan(mean_mutation_frequency)
                else "N/A"
            )
        )

        a2.metric(
            "Alignment Coordinate Range",
            (
                f"{int(position_frequency['position'].min()):,}"
                f"–{int(position_frequency['position'].max()):,}"
                if not position_frequency.empty
                else "N/A"
            )
        )

        a3.metric(
            "Lineage Categories",
            f"{len(lineage_counts):,}"
        )

        st.caption(
            "Mutation frequency represents nucleotide discordance "
            "relative to the reference sequence at each analyzed "
            "alignment position; it is not an evolutionary "
            "substitution-rate estimate."
        )

# 🌿 TAB 2 — LINEAGE TRENDS
# ============================================================

with tabs[1]:
    lin_ts = get_table_freshness("lineage_assignments")
    st.caption(
        f"🕒 Last updated: {lin_ts.strftime('%b %d, %Y — %H:%M:%S')}"
    )

    st.subheader("🌿 Lineage Composition & Surveillance History")

    # --------------------------------------------------------
    # HISTORICAL LINEAGE COMPOSITION
    # --------------------------------------------------------
    # load_from_mysql() intentionally returns only the latest
    # snapshot. For historical lineage analysis, use a
    # read-only aggregated query over all snapshots.
    historical_counts = pd.DataFrame()

    try:
        conn = mysql.connector.connect(**DB_CONFIG)

        query = """
            SELECT
                run_timestamp,
                assigned_lineage,
                COUNT(*) AS Count
            FROM lineage_assignments
            WHERE run_timestamp IS NOT NULL
              AND assigned_lineage IN ('L1', 'L2', 'L3', 'L4', 'L5')
            GROUP BY run_timestamp, assigned_lineage
            ORDER BY run_timestamp
        """

        cursor = conn.cursor(dictionary=True)
        cursor.execute(query)
        historical_counts = pd.DataFrame(cursor.fetchall())

        cursor.close()
        conn.close()

    except Exception as e:
        st.error(
            "Unable to retrieve historical lineage snapshots."
        )
        historical_counts = pd.DataFrame()

    if not historical_counts.empty:

        historical_counts["run_timestamp"] = pd.to_datetime(
            historical_counts["run_timestamp"],
            errors="coerce"
        )

        historical_counts["assigned_lineage"] = (
            historical_counts["assigned_lineage"]
            .astype(str)
            .str.strip()
        )

        historical_counts["Count"] = pd.to_numeric(
            historical_counts["Count"],
            errors="coerce"
        )

        historical_counts = historical_counts.dropna(
            subset=[
                "run_timestamp",
                "assigned_lineage",
                "Count"
            ]
        )

        # ----------------------------------------------------
        # BUILD ONE COMPOSITION VECTOR PER SNAPSHOT
        # ----------------------------------------------------
        canonical = ["L1", "L2", "L3", "L4", "L5"]

        composition = (
            historical_counts
            .pivot(
                index="run_timestamp",
                columns="assigned_lineage",
                values="Count"
            )
            .fillna(0)
        )

        for lineage in canonical:
            if lineage not in composition.columns:
                composition[lineage] = 0

        extra = [
            col for col in composition.columns
            if col not in canonical
        ]

        ordered_columns = canonical + extra
        composition = composition[ordered_columns]

        # Collapse repeated identical lineage-composition
        # snapshots. These represent repeated pipeline states,
        # not independent biological observations.
        composition_unique = composition.loc[
            composition.ne(composition.shift()).any(axis=1)
        ].copy()

        composition_unique = composition_unique.sort_index()

        if not composition_unique.empty:

            history = composition_unique.reset_index()

            history["Snapshot"] = [
                ts.strftime("%d %b %Y")
                for ts in history["run_timestamp"]
            ]

            history["State"] = [
                f"State {i + 1}"
                for i in range(len(history))
            ]

            history["Total"] = history[ordered_columns].sum(axis=1)

            history_long = history.melt(
                id_vars=[
                    "run_timestamp",
                    "Snapshot",
                    "State",
                    "Total"
                ],
                value_vars=ordered_columns,
                var_name="Lineage",
                value_name="Count"
            )

            history_long["Percentage"] = (
                history_long["Count"]
                / history_long["Total"]
                * 100
            )

            # ------------------------------------------------
            # COMPOSITION HISTORY
            # ------------------------------------------------
            st.markdown(
                "### 📊 Lineage Composition Across Surveillance Snapshots"
            )

            st.caption(
                "Distinct lineage-composition states identified "
                "from historical pipeline snapshots; repeated "
                "identical snapshots are collapsed."
            )

            fig_history = px.bar(
                history_long,
                x="Snapshot",
                y="Percentage",
                color="Lineage",
                category_orders={
                    "Lineage": ordered_columns
                },
                text="Percentage",
                custom_data=[
                    "Count",
                    "Total",
                    "run_timestamp"
                ],
                template="plotly_dark"
            )

            fig_history.update_traces(
                texttemplate="%{text:.1f}%",
                textposition="inside",
                hovertemplate=(
                    "<b>%{x}</b><br>"
                    "Lineage: %{fullData.name}<br>"
                    "Share: %{y:.1f}%<br>"
                    "Validated isolates: %{customdata[0]:,}<br>"
                    "Snapshot total: %{customdata[1]:,}<br>"
                    "Timestamp: %{customdata[2]}"
                    "<extra></extra>"
                )
            )

            fig_history.update_layout(
                barmode="stack",
                barnorm="percent",
                height=430,
                margin=dict(
                    l=65,
                    r=45,
                    t=30,
                    b=75
                ),
                xaxis=dict(
                    title="Distinct Surveillance Snapshot",
                    type="category"
                ),
                yaxis=dict(
                    title="Lineage Composition (%)",
                    range=[0, 100],
                    ticksuffix="%"
                ),
                legend=dict(
                    title="Lineage",
                    orientation="h",
                    yanchor="bottom",
                    y=1.02,
                    xanchor="left",
                    x=0
                ),
                hoverlabel=dict(
                    bgcolor="#161b22"
                )
            )

            safe_plotly(fig_history)

            # ------------------------------------------------
            # SURVEILLANCE COHORT SUMMARY
            # ------------------------------------------------
            if len(history) > 1:

                first_total = int(history.iloc[0]["Total"])
                last_total = int(history.iloc[-1]["Total"])
                isolate_change = last_total - first_total

                change_text = (
                    f"+{isolate_change:,}"
                    if isolate_change >= 0
                    else f"{isolate_change:,}"
                )

                st.info(
                    f"**Surveillance cohort:** "
                    f"{first_total:,} → {last_total:,} validated "
                    f"isolates ({change_text} overall). "
                    f"**{len(history)} distinct lineage-composition "
                    f"states** were observed."
                )

            else:

                st.info(
                    f"One distinct lineage-composition state is "
                    f"available ({int(history.iloc[0]['Total']):,} "
                    f"validated isolates)."
                )

            st.caption(
                "Lineage composition describes the distribution of "
                "assigned PRRSV-2 lineages within each surveillance "
                "snapshot. Changes between snapshots should not be "
                "interpreted as direct measurements of evolutionary rate."
            )

            st.caption(
                f"📦 Historical aggregation: "
                f"{len(historical_counts):,} lineage-snapshot "
                f"records collapsed into {len(history):,} distinct "
                f"composition states."
            )

        else:
            st.info(
                "No distinct lineage-composition states could be "
                "derived from the available historical data."
            )

    else:
        st.info(
            "No historical lineage data available."
        )

# ============================================================
# 🧬 TAB 3 — MUTATION HOTSPOTS
# ============================================================

with tabs[2]:

    mut_ts = get_table_freshness("mutation_hotspots")
    st.caption(
        f"🕒 Last updated: {mut_ts.strftime('%b %d, %Y — %H:%M:%S')}"
    )

    st.markdown(
        "<h2 style='color:#00f0ff;'>"
        "🧬 Mutation Hotspot Intelligence"
        "</h2>",
        unsafe_allow_html=True
    )

    st.caption(
        "Latest validated mutation snapshot. Mutation frequency "
        "represents nucleotide discordance relative to the "
        "reference sequence at each analyzed alignment position."
    )

    df_mut = load_from_mysql("mutation_hotspots")

    if not df_mut.empty:

        # ----------------------------------------------------
        # DATA NORMALIZATION
        # ----------------------------------------------------
        df_mut["position"] = pd.to_numeric(
            df_mut["position"],
            errors="coerce"
        )

        df_mut["mutation_frequency"] = pd.to_numeric(
            df_mut["mutation_frequency"],
            errors="coerce"
        )

        df_mut = df_mut.dropna(
            subset=["position", "mutation_frequency"]
        ).copy()

        df_mut = df_mut.sort_values("position")

        # One value per alignment position.
        position_df = (
            df_mut.groupby("position", as_index=False)
            ["mutation_frequency"]
            .max()
            .sort_values("position")
        )

        mean_freq = float(
            position_df["mutation_frequency"].mean()
        )

        analyzed_positions_tab = int(
            position_df["position"].nunique()
        )

        mutation_records_tab = int(len(df_mut))

        # ----------------------------------------------------
        # SUMMARY METRICS
        # ----------------------------------------------------
        st.markdown("### 📌 Mutation Snapshot")

        m1, m2, m3 = st.columns(3)

        m1.metric(
            "Mutation Records",
            f"{mutation_records_tab:,}"
        )

        m2.metric(
            "Analyzed Positions",
            f"{analyzed_positions_tab:,}"
        )

        m3.metric(
            "Mean Mutation Frequency",
            f"{mean_freq:.3f}"
        )

        st.markdown("<div style='height:0.4rem'></div>",
                    unsafe_allow_html=True)

        # ----------------------------------------------------
        # 1. MUTATION-FREQUENCY LANDSCAPE
        # ----------------------------------------------------
        st.markdown(
            "### 📈 Mutation-Frequency Landscape"
        )

        st.caption(
            "Position-level nucleotide discordance across the "
            "PRRSV multiple-sequence alignment. Points represent "
            "analyzed alignment positions; the line shows the "
            "local rolling median."
        )

        landscape = position_df.copy()

        # Rolling window is defined in ordered analyzed
        # positions rather than physical genomic coordinates.
        landscape["rolling_median"] = (
            landscape["mutation_frequency"]
            .rolling(
                window=200,
                min_periods=20,
                center=True
            )
            .median()
        )

        fig_landscape = go.Figure()

        fig_landscape.add_trace(
            go.Scattergl(
                x=landscape["position"],
                y=landscape["mutation_frequency"],
                mode="markers",
                name="Analyzed positions",
                marker=dict(
                    size=4,
                    opacity=0.45
                ),
                hovertemplate=(
                    "Alignment coordinate: %{x:,.0f}<br>"
                    "Mutation frequency: %{y:.3f}"
                    "<extra></extra>"
                )
            )
        )

        fig_landscape.add_trace(
            go.Scatter(
                x=landscape["position"],
                y=landscape["rolling_median"],
                mode="lines",
                name="200-position rolling median",
                line=dict(width=3),
                hovertemplate=(
                    "Alignment coordinate: %{x:,.0f}<br>"
                    "Rolling median: %{y:.3f}"
                    "<extra></extra>"
                )
            )
        )

        fig_landscape.add_hline(
            y=mean_freq,
            line_dash="dash",
            annotation_text=(
                f"Mean = {mean_freq:.3f}"
            ),
            annotation_position="top right"
        )

        fig_landscape.update_layout(
            height=500,
            template="plotly_dark",
            margin=dict(
                l=65,
                r=45,
                t=35,
                b=65
            ),
            hovermode="closest",
            xaxis=dict(
                title="Alignment Coordinate (nt)",
                showgrid=True,
                zeroline=False
            ),
            yaxis=dict(
                title="Mutation Frequency",
                range=[0, 1],
                showgrid=True,
                zeroline=False
            ),
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.01,
                xanchor="left",
                x=0
            ),
            hoverlabel=dict(
                bgcolor="#161b22"
            )
        )

        safe_plotly(fig_landscape)

        st.caption(
            "Alignment coordinates are based on the multiple-"
            "sequence alignment and are not equivalent to ungapped "
            "reference-genome coordinates. Mutation frequency is "
            "not an evolutionary substitution-rate estimate."
        )

        st.markdown("---")

        # ----------------------------------------------------
        # 2. MUTATION RECORD BURDEN BY GENOMIC REGION
        # ----------------------------------------------------
        if "orf" in df_mut.columns:

            st.markdown(
                "### 🧬 Mutation Record Burden by Genomic Region"
            )

            orf_summary = (
                df_mut.groupby("orf", dropna=False)
                .agg(
                    mutation_records=("mutation_frequency", "size"),
                    analyzed_positions=("position", "nunique"),
                    mean_frequency=("mutation_frequency", "mean")
                )
                .reset_index()
            )

            orf_summary["orf"] = (
                orf_summary["orf"]
                .fillna("Unannotated")
                .astype(str)
            )

            orf_summary = orf_summary.sort_values(
                "mutation_records",
                ascending=True
            )

            fig_burden = px.bar(
                orf_summary,
                x="mutation_records",
                y="orf",
                orientation="h",
                text="mutation_records",
                custom_data=[
                    "analyzed_positions",
                    "mean_frequency"
                ],
                template="plotly_dark"
            )

            fig_burden.update_traces(
                texttemplate="%{text:,}",
                textposition="outside",
                cliponaxis=False,
                hovertemplate=(
                    "<b>%{y}</b><br>"
                    "Mutation records: %{x:,}<br>"
                    "Analyzed positions: %{customdata[0]:,}<br>"
                    "Mean frequency: %{customdata[1]:.3f}"
                    "<extra></extra>"
                )
            )

            fig_burden.update_layout(
                height=470,
                showlegend=False,
                margin=dict(
                    l=105,
                    r=95,
                    t=30,
                    b=60
                ),
                xaxis=dict(
                    title="Mutation Records",
                    showgrid=True,
                    zeroline=False,
                    tickformat=","
                ),
                yaxis=dict(
                    title=None
                )
            )

            safe_plotly(fig_burden)

            st.caption(
                "Mutation-record burden reflects the number of "
                "mutation records assigned to each annotated "
                "genomic region; it does not represent an "
                "evolutionary mutation rate."
            )

            st.markdown("---")

            # ------------------------------------------------
            # 3. MEAN MUTATION FREQUENCY BY REGION
            # ------------------------------------------------
            st.markdown(
                "### 📊 Mean Mutation Frequency by Genomic Region"
            )

            orf_frequency = orf_summary.sort_values(
                "mean_frequency",
                ascending=True
            )

            fig_orf_frequency = px.bar(
                orf_frequency,
                x="mean_frequency",
                y="orf",
                orientation="h",
                text="mean_frequency",
                custom_data=[
                    "mutation_records",
                    "analyzed_positions"
                ],
                template="plotly_dark"
            )

            fig_orf_frequency.update_traces(
                texttemplate="%{text:.3f}",
                textposition="outside",
                cliponaxis=False,
                hovertemplate=(
                    "<b>%{y}</b><br>"
                    "Mean mutation frequency: %{x:.3f}<br>"
                    "Mutation records: %{customdata[0]:,}<br>"
                    "Analyzed positions: %{customdata[1]:,}"
                    "<extra></extra>"
                )
            )

            fig_orf_frequency.update_layout(
                height=470,
                showlegend=False,
                margin=dict(
                    l=105,
                    r=80,
                    t=30,
                    b=60
                ),
                xaxis=dict(
                    title="Mean Mutation Frequency",
                    range=[
                        0,
                        min(
                            0.25,
                            max(
                                0.25,
                                float(
                                    orf_frequency[
                                        "mean_frequency"
                                    ].max()
                                ) * 1.15
                            )
                        )
                    ],
                    showgrid=True,
                    zeroline=False
                ),
                yaxis=dict(
                    title=None
                )
            )

            safe_plotly(fig_orf_frequency)

            st.caption(
                "Mean mutation frequency summarizes the average "
                "reference-discordance frequency among mutation "
                "records assigned to each region."
            )

    else:
        st.info("No mutation hotspot data available.")

# ============================================================
# 💉 TAB 4 — VACCINE ESCAPE
# ============================================================

with tabs[3]:

    esc_ts = get_table_freshness("vaccine_escape_predictions")

    if esc_ts is not None:
        st.caption(
            f"🕒 Last updated: {esc_ts.strftime('%b %d, %Y — %H:%M:%S')}"
        )
    else:
        st.caption("🕒 Last updated: unavailable")

    st.markdown(
        "<h2 style='color:#00f0ff;'>"
        "💉 Vaccine Escape Intelligence"
        "</h2>",
        unsafe_allow_html=True
    )

    st.caption(
        "Model-derived assessment of sequence match and predicted "
        "escape probability across the latest validated "
        "vaccine-reference snapshot."
    )

    df_esc = load_from_mysql("vaccine_escape_predictions")

    if not df_esc.empty:

        # ----------------------------------------------------
        # NORMALIZE DATA
        # ----------------------------------------------------
        df_esc["escape_probability"] = pd.to_numeric(
            df_esc["escape_probability"],
            errors="coerce"
        )

        df_esc["match_score_%"] = pd.to_numeric(
            df_esc["match_score_%"],
            errors="coerce"
        )

        df_esc = df_esc.dropna(
            subset=["escape_probability", "match_score_%"]
        ).copy()

        # ----------------------------------------------------
        # SNAPSHOT METRICS
        # ----------------------------------------------------
        n_records = len(df_esc)
        n_isolates = df_esc["isolate"].nunique()
        n_vaccines = df_esc["vaccine"].nunique()
        mean_escape = float(df_esc["escape_probability"].mean())

        st.markdown("### 📌 Vaccine Escape Snapshot")

        c1, c2, c3, c4 = st.columns(4)

        c1.metric(
            "Prediction Records",
            f"{n_records:,}"
        )

        c2.metric(
            "PRRSV Isolates",
            f"{n_isolates:,}"
        )

        c3.metric(
            "Vaccine References",
            f"{n_vaccines:,}"
        )

        c4.metric(
            "Mean Escape Probability",
            f"{mean_escape:.3f}"
        )

        st.markdown(
            "<div style='height:0.35rem'></div>",
            unsafe_allow_html=True
        )

        # ----------------------------------------------------
        # RISK DISTRIBUTION
        # ----------------------------------------------------
        st.markdown("### ⚠️ Predicted Escape-Risk Distribution")

        risk_summary = (
            df_esc.groupby("predicted_escape_risk")
            .agg(
                records=("escape_probability", "size"),
                isolates=("isolate", "nunique"),
                mean_probability=("escape_probability", "mean")
            )
            .reset_index()
        )

        risk_order = {
            "High": 0,
            "Moderate": 1,
            "Low": 2
        }

        risk_summary["sort_order"] = (
            risk_summary["predicted_escape_risk"]
            .map(risk_order)
            .fillna(99)
        )

        risk_summary = risk_summary.sort_values("sort_order")

        fig_risk = px.bar(
            risk_summary,
            x="predicted_escape_risk",
            y="records",
            text="records",
            custom_data=[
                "isolates",
                "mean_probability"
            ],
            template="plotly_dark"
        )

        fig_risk.update_traces(
            texttemplate="%{text:,}",
            textposition="outside",
            cliponaxis=False,
            hovertemplate=(
                "<b>%{x}</b><br>"
                "Prediction records: %{y:,}<br>"
                "Unique isolates: %{customdata[0]:,}<br>"
                "Mean escape probability: %{customdata[1]:.3f}"
                "<extra></extra>"
            )
        )

        fig_risk.update_layout(
            height=430,
            showlegend=False,
            margin=dict(
                l=65,
                r=65,
                t=35,
                b=65
            ),
            xaxis=dict(
                title="Predicted Escape Risk"
            ),
            yaxis=dict(
                title="Prediction Records",
                tickformat=",",
                showgrid=True,
                zeroline=False
            )
        )

        safe_plotly(fig_risk)

        st.caption(
            "Risk categories summarize model predictions across "
            "20,526 isolate–vaccine-reference comparisons. They "
            "should not be interpreted as experimentally confirmed "
            "vaccine failure."
        )

        st.markdown("---")

        # ----------------------------------------------------
        # VACCINE REFERENCE COMPARISON
        # ----------------------------------------------------
        st.markdown(
            "### 🧪 Mean Predicted Escape Probability by Vaccine Reference"
        )

        vaccine_summary = (
            df_esc.groupby("vaccine")
            .agg(
                records=("escape_probability", "size"),
                isolates=("isolate", "nunique"),
                mean_escape=("escape_probability", "mean"),
                mean_match=("match_score_%", "mean")
            )
            .reset_index()
            .sort_values(
                "mean_escape",
                ascending=True
            )
        )

        fig_vaccine = px.bar(
            vaccine_summary,
            x="mean_escape",
            y="vaccine",
            orientation="h",
            text="mean_escape",
            custom_data=[
                "records",
                "isolates",
                "mean_match"
            ],
            template="plotly_dark"
        )

        fig_vaccine.update_traces(
            texttemplate="%{text:.3f}",
            textposition="outside",
            cliponaxis=False,
            hovertemplate=(
                "<b>%{y}</b><br>"
                "Mean escape probability: %{x:.3f}<br>"
                "Prediction records: %{customdata[0]:,}<br>"
                "Unique isolates: %{customdata[1]:,}<br>"
                "Mean match score: %{customdata[2]:.2f}%"
                "<extra></extra>"
            )
        )

        fig_vaccine.update_layout(
            height=560,
            showlegend=False,
            margin=dict(
                l=230,
                r=95,
                t=35,
                b=65
            ),
            xaxis=dict(
                title="Mean Predicted Escape Probability",
                range=[0, 0.45],
                showgrid=True,
                zeroline=False
            ),
            yaxis=dict(
                title=None
            )
        )

        safe_plotly(fig_vaccine)

        st.caption(
            "Each vaccine reference is represented by 1,866 "
            "isolate-level comparisons. Higher values indicate "
            "higher model-derived predicted escape probability "
            "within this dataset."
        )

        st.markdown("---")

        # ----------------------------------------------------
        # MATCH SCORE / ESCAPE RELATIONSHIP
        # ----------------------------------------------------
        st.markdown(
            "### 🔬 Match Score vs Predicted Escape Probability"
        )

        st.caption(
            "Each point represents one vaccine reference summarized "
            "across the 1,866 PRRSV isolates."
        )

        fig_relationship = px.scatter(
            vaccine_summary,
            x="mean_match",
            y="mean_escape",
            text="vaccine",
            size="records",
            hover_name="vaccine",
            template="plotly_dark"
        )

        fig_relationship.update_traces(
            textposition="top center",
            marker=dict(
                size=14,
                opacity=0.85
            ),
            hovertemplate=(
                "<b>%{hovertext}</b><br>"
                "Mean match score: %{x:.2f}%<br>"
                "Mean escape probability: %{y:.3f}"
                "<extra></extra>"
            )
        )

        fig_relationship.update_layout(
            height=500,
            showlegend=False,
            margin=dict(
                l=75,
                r=50,
                t=35,
                b=65
            ),
            xaxis=dict(
                title="Mean Match Score (%)",
                range=[
                    max(
                        0,
                        float(vaccine_summary["mean_match"].min()) - 5
                    ),
                    min(
                        100,
                        float(vaccine_summary["mean_match"].max()) + 5
                    )
                ],
                showgrid=True,
                zeroline=False
            ),
            yaxis=dict(
                title="Mean Predicted Escape Probability",
                range=[0, 0.45],
                showgrid=True,
                zeroline=False
            )
        )

        safe_plotly(fig_relationship)

        st.markdown("---")

        # ----------------------------------------------------
        # DATA COMPLETENESS
        # ----------------------------------------------------
        missing_lineage = int(
            df_esc["lineage"].isna().sum()
        ) if "lineage" in df_esc.columns else len(df_esc)

        missing_accession = int(
            df_esc["accession"].isna().sum()
        ) if "accession" in df_esc.columns else len(df_esc)

        st.info(
            "Lineage-stratified vaccine-escape analysis is not "
            "displayed because the current validated vaccine "
            "snapshot contains no lineage values. Accession "
            "identifiers are also absent from this snapshot."
        )

        st.caption(
            f"Data completeness: {missing_lineage:,} of "
            f"{len(df_esc):,} records lack lineage values; "
            f"{missing_accession:,} lack accession identifiers."
        )

        st.markdown(
            "#### ⚠️ Interpretation"
        )

        st.caption(
            "Escape probability and predicted escape risk are "
            "computational model outputs. They do not constitute "
            "experimental evidence of immune escape, vaccine "
            "failure, or clinical protection failure. Experimental "
            "and epidemiological validation would be required to "
            "establish biological or field significance."
        )

    else:
        st.info("No vaccine escape prediction data available.")


# ============================================================
# 🧠 TAB 5 — miRNA–PRRSV INTERACTION ANALYTICS
# ============================================================

with tabs[4]:

    st.markdown(
        "<h2 style='color:#00f0ff;'>"
        "🧠 miRNA–PRRSV Interaction Intelligence"
        "</h2>",
        unsafe_allow_html=True
    )

    st.caption(
        "Validated computational miRNA–PRRSV interaction dataset "
        "generated by the real-miRNA analysis pipeline."
    )

    # --------------------------------------------------------
    # LOAD VALIDATED RESULT — NOT THE 89.9M-ROW MYSQL TABLE
    # --------------------------------------------------------
    from scripts.paths import MIRNA_RESULT_FILE

    MIRNA_FILE = MIRNA_RESULT_FILE

    @st.cache_data(ttl=600, show_spinner=False)
    def load_validated_mirna_results(path):
        df = pd.read_csv(path)

        df["binding_score"] = pd.to_numeric(
            df["binding_score"],
            errors="coerce"
        )

        df["energy_kcal"] = pd.to_numeric(
            df["energy_kcal"],
            errors="coerce"
        )

        df = df.dropna(
            subset=[
                "miRNA",
                "target_genome",
                "binding_score",
                "energy_kcal"
            ]
        ).copy()

        return df

    if not MIRNA_FILE.exists():

        st.error(
            "Validated miRNA result file was not found: "
            f"{MIRNA_FILE}"
        )

    else:

        df_mir = load_validated_mirna_results(
            str(MIRNA_FILE)
        )

        # ----------------------------------------------------
        # SNAPSHOT METRICS
        # ----------------------------------------------------
        n_interactions = len(df_mir)
        n_mirnas = df_mir["miRNA"].nunique()
        n_targets = df_mir["target_genome"].nunique()

        mean_score = float(
            df_mir["binding_score"].mean()
        )

        mean_energy = float(
            df_mir["energy_kcal"].mean()
        )

        st.markdown("### 📌 Interaction Snapshot")

        c1, c2, c3, c4 = st.columns(4)

        c1.metric(
            "Predicted Interactions",
            f"{n_interactions:,}"
        )

        c2.metric(
            "Unique miRNAs",
            f"{n_mirnas:,}"
        )

        c3.metric(
            "Target Genomes",
            f"{n_targets:,}"
        )

        c4.metric(
            "Mean Binding Score",
            f"{mean_score:.3f}"
        )

        st.caption(
            f"Mean predicted binding energy: "
            f"{mean_energy:.3f} kcal/mol"
        )

        st.markdown(
            "<div style='height:0.35rem'></div>",
            unsafe_allow_html=True
        )

        # ----------------------------------------------------
        # TOP miRNAs
        # ----------------------------------------------------
        st.markdown(
            "### 🎯 Top miRNAs by Predicted Interaction Coverage"
        )

        top_mir = (
            df_mir.groupby("miRNA")
            .agg(
                interactions=("miRNA", "size"),
                target_genomes=("target_genome", "nunique"),
                mean_score=("binding_score", "mean"),
                mean_energy=("energy_kcal", "mean")
            )
            .reset_index()
            .sort_values(
                "interactions",
                ascending=False
            )
            .head(15)
            .sort_values(
                "interactions",
                ascending=True
            )
        )

        fig_top_mir = px.bar(
            top_mir,
            x="interactions",
            y="miRNA",
            orientation="h",
            text="interactions",
            custom_data=[
                "target_genomes",
                "mean_score",
                "mean_energy"
            ],
            template="plotly_dark"
        )

        fig_top_mir.update_traces(
            texttemplate="%{text:,}",
            textposition="outside",
            cliponaxis=False,
            hovertemplate=(
                "<b>%{y}</b><br>"
                "Predicted interactions: %{x:,}<br>"
                "Target genomes: %{customdata[0]:,}<br>"
                "Mean binding score: %{customdata[1]:.4f}<br>"
                "Mean energy: %{customdata[2]:.3f} kcal/mol"
                "<extra></extra>"
            )
        )

        fig_top_mir.update_layout(
            height=560,
            showlegend=False,
            margin=dict(
                l=190,
                r=90,
                t=30,
                b=65
            ),
            xaxis=dict(
                title="Predicted Interaction Count",
                tickformat=",",
                showgrid=True,
                zeroline=False
            ),
            yaxis=dict(
                title=None
            )
        )

        safe_plotly(fig_top_mir)

        st.caption(
            "Interaction coverage is the number of predicted "
            "miRNA–target-genome pairs. Target-genome counts "
            "represent distinct PRRSV genomes."
        )

        st.markdown("---")

        # ----------------------------------------------------
        # BINDING ENERGY DISTRIBUTION
        # ----------------------------------------------------
        st.markdown(
            "### 🧪 Predicted Binding-Energy Distribution"
        )

        fig_energy = px.histogram(
            df_mir,
            x="energy_kcal",
            nbins=40,
            template="plotly_dark"
        )

        fig_energy.update_traces(
            hovertemplate=(
                "Energy: %{x:.2f} kcal/mol<br>"
                "Predicted interactions: %{y:,}"
                "<extra></extra>"
            )
        )

        fig_energy.update_layout(
            height=430,
            showlegend=False,
            margin=dict(
                l=70,
                r=45,
                t=30,
                b=65
            ),
            xaxis=dict(
                title="Predicted Binding Energy (kcal/mol)",
                showgrid=True,
                zeroline=False
            ),
            yaxis=dict(
                title="Predicted Interactions",
                tickformat=",",
                showgrid=True,
                zeroline=False
            )
        )

        safe_plotly(fig_energy)

        st.caption(
            "More negative predicted binding energies indicate "
            "stronger predicted interaction energetics within "
            "the computational scoring framework. These values "
            "are predictions, not experimental binding measurements."
        )

        st.markdown("---")

        # ----------------------------------------------------
        # SCORE VS ENERGY — DETERMINISTIC BINNED VIEW
        # ----------------------------------------------------
        st.markdown(
            "### 🔬 Binding Score vs Predicted Binding Energy"
        )

        st.caption(
            "The visualization uses deterministic score–energy "
            "bins to summarize the full validated dataset without "
            "plotting hundreds of thousands of individual points."
        )

        score_bins = pd.cut(
            df_mir["binding_score"],
            bins=12,
            include_lowest=True
        )

        energy_bins = pd.cut(
            df_mir["energy_kcal"],
            bins=12,
            include_lowest=True
        )

        heat_df = (
            df_mir.assign(
                score_bin=score_bins,
                energy_bin=energy_bins
            )
            .groupby(
                ["energy_bin", "score_bin"],
                observed=False
            )
            .size()
            .reset_index(name="interactions")
        )

        heat_df["score_label"] = (
            heat_df["score_bin"]
            .map(lambda x: f"{x.left:.3f}–{x.right:.3f}")
        )

        heat_df["energy_label"] = (
            heat_df["energy_bin"]
            .map(lambda x: f"{x.left:.2f}–{x.right:.2f}")
        )

        fig_heat = px.density_heatmap(
            heat_df,
            x="score_label",
            y="energy_label",
            z="interactions",
            text_auto=True,
            template="plotly_dark"
        )

        fig_heat.update_layout(
            height=500,
            margin=dict(
                l=105,
                r=55,
                t=30,
                b=95
            ),
            xaxis=dict(
                title="Binding Score Range"
            ),
            yaxis=dict(
                title="Binding Energy Range (kcal/mol)"
            ),
            coloraxis_colorbar=dict(
                title="Interactions"
            )
        )

        safe_plotly(fig_heat)

        st.markdown("---")

        # ----------------------------------------------------
        # TOP-miRNA SUMMARY TABLE
        # ----------------------------------------------------
        st.markdown(
            "### 📋 Top miRNA Interaction Summary"
        )

        table_df = top_mir.sort_values(
            "interactions",
            ascending=False
        ).copy()

        table_df["mean_score"] = table_df[
            "mean_score"
        ].round(4)

        table_df["mean_energy"] = table_df[
            "mean_energy"
        ].round(3)

        table_df = table_df.rename(
            columns={
                "miRNA": "miRNA",
                "interactions": "Predicted interactions",
                "target_genomes": "Target genomes",
                "mean_score": "Mean binding score",
                "mean_energy": "Mean energy (kcal/mol)"
            }
        )

        st.dataframe(
            table_df,
            use_container_width=True,
            hide_index=True
        )

        st.markdown("---")

        # ----------------------------------------------------
        # INTERPRETATION
        # ----------------------------------------------------
        st.markdown(
            "#### ⚠️ Interpretation"
        )

        st.caption(
            "All miRNA–PRRSV interactions shown here are "
            "computational predictions from the validated "
            "miRNA analysis dataset. They do not constitute "
            "experimental evidence of miRNA binding, target "
            "regulation, antiviral activity, or biological "
            "effect. Experimental validation is required to "
            "establish biological significance."
        )

        st.caption(
            "Dataset provenance: "
            "results/mirna_interactions_real.csv · "
            f"{n_interactions:,} records · "
            f"{n_mirnas:,} miRNAs · "
            f"{n_targets:,} target genomes."
        )

# ============================================================
# 🌳 TAB 6 — PHYLOGENETIC TREE VISUALIZATION
# ============================================================

with tabs[5]:

    st.markdown(
        "<h2 style='color:#00f0ff;'>"
        "🌳 PRRSV Phylogenetic Tree"
        "</h2>",
        unsafe_allow_html=True
    )

    st.caption(
        "Interactive phylogram reconstructed from the authoritative "
        "Newick tree stored in the IntelliPRRSV2 database."
    )

    # --------------------------------------------------------
    # LOAD AUTHORITATIVE TREE
    # --------------------------------------------------------
    df_tree = load_from_mysql("phylogenetic_tree")

    if df_tree.empty or "newick_data" not in df_tree.columns:

        st.info("No phylogenetic tree data available in MySQL.")

    else:

        try:

            # The table contains one authoritative tree.
            newick_data = df_tree.iloc[0]["newick_data"]

            if not isinstance(newick_data, str) or not newick_data.strip():

                raise ValueError(
                    "The phylogenetic tree contains empty Newick data."
                )

            # ------------------------------------------------
            # PARSE NEWICK
            # ------------------------------------------------
            handle = StringIO(newick_data)
            phylo_tree = Phylo.read(handle, "newick")

            terminals = phylo_tree.get_terminals()
            internal_nodes = phylo_tree.get_nonterminals()

            branch_lengths = [
                clade.branch_length
                for clade in phylo_tree.find_clades()
                if clade.branch_length is not None
            ]

            support_values = [
                clade.confidence
                for clade in phylo_tree.find_clades()
                if clade.confidence is not None
            ]

            # ------------------------------------------------
            # TREE METRICS
            # ------------------------------------------------
            n_tips = len(terminals)
            n_internal = len(internal_nodes)
            n_branches = len(branch_lengths)
            n_support = len(support_values)

            max_branch = (
                max(branch_lengths)
                if branch_lengths else None
            )

            # ------------------------------------------------
            # SNAPSHOT METRICS
            # ------------------------------------------------
            st.markdown("### 📌 Tree Snapshot")

            c1, c2, c3, c4 = st.columns(4)

            c1.metric(
                "Terminal Isolates",
                f"{n_tips:,}"
            )

            c2.metric(
                "Internal Nodes",
                f"{n_internal:,}"
            )

            c3.metric(
                "Branches",
                f"{n_branches:,}"
            )

            c4.metric(
                "Nodes with Support",
                f"{n_support:,}"
            )

            st.caption(
                "Database tree timestamp: "
                f"{df_tree.iloc[0].get('timestamp', 'unavailable')}"
            )

            # ------------------------------------------------
            # TREE COORDINATES
            #
            # A rectangular phylogram is constructed directly
            # from the parsed Newick topology.
            #
            # x = cumulative branch length
            # y = terminal ordering / topology
            # ------------------------------------------------
            y_counter = [0]

            def assign_coordinates(clade, x_parent=0.0):

                branch = clade.branch_length or 0.0
                x_here = x_parent + branch

                if clade.is_terminal():

                    y_here = y_counter[0]
                    y_counter[0] += 1

                else:

                    child_coords = [
                        assign_coordinates(
                            child,
                            x_here
                        )
                        for child in clade.clades
                    ]

                    child_y = [
                        coords[1]
                        for coords in child_coords
                    ]

                    y_here = (
                        sum(child_y) / len(child_y)
                        if child_y else y_counter[0]
                    )

                return x_here, y_here

            assign_coordinates(phylo_tree.root)

            # ------------------------------------------------
            # SECOND PASS — BUILD BRANCH SEGMENTS
            # ------------------------------------------------
            branch_x = []
            branch_y = []

            terminal_x = []
            terminal_y = []
            terminal_names = []

            support_x = []
            support_y = []
            support_text = []

            y_counter = [0]

            def build_segments(clade, x_parent=0.0):

                branch = clade.branch_length or 0.0
                x_here = x_parent + branch

                if clade.is_terminal():

                    y_here = y_counter[0]
                    y_counter[0] += 1

                    terminal_x.append(x_here)
                    terminal_y.append(y_here)
                    terminal_names.append(
                        clade.name if clade.name else "Unnamed tip"
                    )

                    return x_here, y_here

                child_coords = []

                for child in clade.clades:

                    child_x, child_y = build_segments(
                        child,
                        x_here
                    )

                    child_coords.append(
                        (child_x, child_y)
                    )

                    # Horizontal branch
                    branch_x.extend(
                        [x_here, child_x, None]
                    )
                    branch_y.extend(
                        [child_y, child_y, None]
                    )

                child_y_values = [
                    item[1]
                    for item in child_coords
                ]

                y_here = (
                    sum(child_y_values) /
                    len(child_y_values)
                )

                # Vertical connector
                if child_y_values:

                    branch_x.extend(
                        [x_here, x_here, None]
                    )

                    branch_y.extend(
                        [
                            min(child_y_values),
                            max(child_y_values),
                            None
                        ]
                    )

                # Store internal node support.
                if clade.confidence is not None:

                    support_x.append(x_here)
                    support_y.append(y_here)
                    support_text.append(
                        f"Node support: {clade.confidence:.3f}"
                    )

                return x_here, y_here

            build_segments(phylo_tree.root)

            # ------------------------------------------------
            # PHYLOGRAM
            # ------------------------------------------------
            st.markdown(
                "### 🧬 Interactive Phylogenetic Tree"
            )

            fig_tree = go.Figure()

            # Tree branches
            fig_tree.add_trace(
                go.Scatter(
                    x=branch_x,
                    y=branch_y,
                    mode="lines",
                    line=dict(
                        width=0.8
                    ),
                    hoverinfo="skip",
                    name="Tree branches"
                )
            )

            # Terminal tips
            fig_tree.add_trace(
                go.Scatter(
                    x=terminal_x,
                    y=terminal_y,
                    mode="markers",
                    marker=dict(
                        size=4
                    ),
                    customdata=terminal_names,
                    hovertemplate=(
                        "<b>%{customdata}</b>"
                        "<br>Terminal isolate"
                        "<br>Branch coordinate: %{x:.6f}"
                        "<extra></extra>"
                    ),
                    name="Isolates"
                )
            )

            fig_tree.update_layout(
                template="plotly_dark",
                height=max(
                    850,
                    min(1800, 450 + n_tips * 0.45)
                ),
                margin=dict(
                    l=55,
                    r=45,
                    t=45,
                    b=70
                ),
                showlegend=True,
                legend=dict(
                    orientation="h",
                    yanchor="bottom",
                    y=1.01,
                    xanchor="left",
                    x=0
                ),
                plot_bgcolor="#0d1117",
                paper_bgcolor="#0d1117",
                xaxis=dict(
                    title="Branch Length (tree units)",
                    showgrid=True,
                    zeroline=False,
                    showticklabels=True
                ),
                yaxis=dict(
                    title="Tree Tip Order",
                    showgrid=False,
                    zeroline=False,
                    showticklabels=False,
                    autorange="reversed"
                )
            )

            safe_plotly(fig_tree)

            # ------------------------------------------------
            # TREE INTERPRETATION
            # ------------------------------------------------
            st.caption(
                "The topology represents an unrooted phylogenetic tree. "
                "Horizontal branch position reflects cumulative branch "
                "length encoded in the Newick tree; it should not be "
                "interpreted as sampling time. Individual isolate labels "
                "are available by hovering over terminal tips."
            )

            # ------------------------------------------------
            # OPTIONAL SUPPORT-VALUE VIEW
            # ------------------------------------------------
            if n_support > 0:

                with st.expander(
                    "🔎 Inspect encoded node-support values"
                ):

                    support_df = pd.DataFrame(
                        {
                            "Node position": support_x,
                            "Tree order": support_y,
                            "Support": [
                                float(
                                    text.split(": ")[1]
                                )
                                for text in support_text
                            ]
                        }
                    )

                    support_df = support_df.sort_values(
                        "Support",
                        ascending=False
                    )

                    st.dataframe(
                        support_df.head(100),
                        use_container_width=True,
                        hide_index=True
                    )

                    st.caption(
                        f"The Newick tree contains {n_support:,} "
                        "internal-node support values. The table shows "
                        "the 100 highest encoded values; no support "
                        "values are inferred or recalculated by the dashboard."
                    )

            # ------------------------------------------------
            # TREE QC
            # ------------------------------------------------
            st.markdown("### 🔍 Tree Quality-Control Summary")

            q1, q2, q3 = st.columns(3)

            q1.metric(
                "Unique Tip Labels",
                f"{len(set(terminal_names)):,}"
            )

            q2.metric(
                "Duplicate Tip Labels",
                f"{len(terminal_names) - len(set(terminal_names)):,}"
            )

            q3.metric(
                "Maximum Branch Length",
                f"{max_branch:.6f}"
                if max_branch is not None
                else "N/A"
            )

            if len(set(terminal_names)) == n_tips:

                st.success(
                    "✅ All terminal labels are unique."
                )

            else:

                st.warning(
                    "⚠️ Duplicate terminal labels detected."
                )

            # ------------------------------------------------
            # RAW NEWICK
            # ------------------------------------------------
            with st.expander("🧾 View Raw Newick Data"):

                st.caption(
                    f"{len(newick_data):,} characters"
                )

                st.code(
                    newick_data[:5000] + (
                        "\n..."
                        if len(newick_data) > 5000
                        else ""
                    )
                )

            # ------------------------------------------------
            # SCIENTIFIC DISCLAIMER
            # ------------------------------------------------
            st.markdown("#### ⚠️ Interpretation")

            st.caption(
                "This visualization displays the topology and branch "
                "lengths encoded in the stored Newick tree. The tree is "
                "unrooted, and the dashboard does not infer direction "
                "of transmission, ancestry through time, or geographic "
                "spread. Node-support values are displayed as encoded "
                "in the source tree and are not independently recalculated "
                "by the dashboard."
            )

        except Exception as e:

            st.error(
                f"⚠️ Error rendering phylogenetic tree: {e}"
            )

# ============================================================
# 🌍 TAB 7 — GLOBAL DISTRIBUTION MAP
# ============================================================

with tabs[6]:

    st.markdown(
        "<h2 style='color:#00f0ff;'>"
        "🌍 Global PRRSV Geographic Distribution"
        "</h2>",
        unsafe_allow_html=True
    )

    st.caption(
        "Country-level geographic distribution of PRRSV records "
        "with available geographic metadata."
    )

    # --------------------------------------------------------
    # LOAD LATEST GEOGRAPHIC SNAPSHOT
    # --------------------------------------------------------
    df_geo = load_from_mysql("geo_metadata")

    if df_geo.empty:

        st.info(
            "No geographic metadata are available "
            "for the latest pipeline snapshot."
        )

    else:

        df_geo.columns = [
            c.lower() for c in df_geo.columns
        ]

        # ----------------------------------------------------
        # NORMALIZE GEOGRAPHIC FIELDS
        # ----------------------------------------------------
        df_geo["country"] = (
            df_geo["country"]
            .fillna("Unknown")
            .astype(str)
            .str.strip()
        )

        df_geo.loc[
            df_geo["country"].eq(""),
            "country"
        ] = "Unknown"

        df_geo["latitude"] = pd.to_numeric(
            df_geo.get("latitude"),
            errors="coerce"
        )

        df_geo["longitude"] = pd.to_numeric(
            df_geo.get("longitude"),
            errors="coerce"
        )

        # ----------------------------------------------------
        # VALID COORDINATE PAIRS
        # ----------------------------------------------------
        valid_coords = df_geo[
            df_geo["latitude"].between(-90, 90, inclusive="both")
            &
            df_geo["longitude"].between(-180, 180, inclusive="both")
        ].copy()

        n_records = len(df_geo)

        n_countries = df_geo[
            df_geo["country"].ne("Unknown")
        ]["country"].nunique()

        n_coord_records = len(valid_coords)

        n_mapped_countries = valid_coords[
            valid_coords["country"].ne("Unknown")
        ]["country"].nunique()

        n_unmapped_records = (
            n_records - n_coord_records
        )

        n_unmapped_countries = (
            n_countries - n_mapped_countries
        )

        # ----------------------------------------------------
        # SNAPSHOT
        # ----------------------------------------------------
        st.markdown("### 📌 Geographic Snapshot")

        c1, c2, c3, c4 = st.columns(4)

        c1.metric(
            "Geographic Records",
            f"{n_records:,}"
        )

        c2.metric(
            "Countries",
            f"{n_countries:,}"
        )

        c3.metric(
            "Records with Coordinates",
            f"{n_coord_records:,}"
        )

        c4.metric(
            "Mapped Countries",
            f"{n_mapped_countries:,}"
        )

        coord_pct = (
            100 * n_coord_records / n_records
            if n_records else 0
        )

        st.caption(
            f"Coordinate coverage: {coord_pct:.1f}% of geographic "
            "records. Coordinates are country-level representations "
            "rather than individual sampling locations."
        )

        st.markdown("---")

        # ----------------------------------------------------
        # COUNTRY SUMMARY
        # ----------------------------------------------------
        country_summary = (
            df_geo.groupby("country")
            .size()
            .reset_index(name="records")
            .sort_values(
                "records",
                ascending=False
            )
        )

        country_summary["percentage"] = (
            100 *
            country_summary["records"] /
            n_records
        )

        # ----------------------------------------------------
        # COUNTRY-LEVEL COORDINATE MAP
        # ----------------------------------------------------
        st.markdown(
            "### 🗺️ Country-Level Geographic Distribution"
        )

        if not valid_coords.empty:

            # One coordinate per country in the validated dataset.
            # The audit established that coordinates correspond to
            # country-level geographic representations.
            mapped_summary = (
                valid_coords[
                    valid_coords["country"].ne("Unknown")
                ]
                .groupby("country")
                .agg(
                    records=("country", "size"),
                    latitude=("latitude", "first"),
                    longitude=("longitude", "first")
                )
                .reset_index()
            )

            mapped_summary["percentage"] = (
                100 *
                mapped_summary["records"] /
                n_records
            )

            fig_geo = px.scatter_geo(
                mapped_summary,
                lat="latitude",
                lon="longitude",
                size="records",
                hover_name="country",
                custom_data=[
                    "records",
                    "percentage"
                ],
                projection="natural earth",
                template="plotly_dark"
            )

            fig_geo.update_traces(
                marker=dict(
                    sizemin=6,
                    opacity=0.78,
                    line=dict(
                        width=0.7
                    )
                ),
                hovertemplate=(
                    "<b>%{hovertext}</b><br>"
                    "Records: %{customdata[0]:,}<br>"
                    "Share of geographic records: "
                    "%{customdata[1]:.1f}%"
                    "<extra></extra>"
                )
            )

            fig_geo.update_geos(
                projection_type="natural earth",
                showland=True,
                showcountries=True,
                showcoastlines=True,
                showocean=True,
                lataxis_range=[-10, 70],
                lonaxis_range=[-120, 160]
            )

            fig_geo.update_layout(
                height=650,
                margin=dict(
                    l=20,
                    r=20,
                    t=20,
                    b=20
                ),
                showlegend=False,
                paper_bgcolor="#0d1117",
                plot_bgcolor="#0d1117",
                geo=dict(
                    showland=True,
                    showcountries=True,
                    showcoastlines=True
                )
            )

            safe_plotly(fig_geo)

            st.caption(
                f"{len(mapped_summary):,} countries are represented "
                "on the map. Marker size reflects the number of "
                "geographic records assigned to each country."
            )

        else:

            st.info(
                "No valid coordinate pairs are available "
                "for the current geographic snapshot."
            )

        st.markdown("---")

        # ----------------------------------------------------
        # COUNTRY RANKING
        # ----------------------------------------------------
        st.markdown(
            "### 📊 PRRSV Records by Country"
        )

        ranking = country_summary.sort_values(
            "records",
            ascending=True
        )

        fig_country = px.bar(
            ranking,
            x="records",
            y="country",
            orientation="h",
            text="records",
            custom_data=["percentage"],
            template="plotly_dark"
        )

        fig_country.update_traces(
            texttemplate="%{text:,}",
            textposition="outside",
            cliponaxis=False,
            hovertemplate=(
                "<b>%{y}</b><br>"
                "Records: %{x:,}<br>"
                "Share: %{customdata[0]:.1f}%"
                "<extra></extra>"
            )
        )

        fig_country.update_layout(
            height=max(
                600,
                300 + len(ranking) * 24
            ),
            margin=dict(
                l=150,
                r=85,
                t=25,
                b=60
            ),
            showlegend=False,
            xaxis=dict(
                title="Geographic Records",
                tickformat=",",
                showgrid=True,
                zeroline=False
            ),
            yaxis=dict(
                title=None
            )
        )

        safe_plotly(fig_country)

        st.markdown("---")

        # ----------------------------------------------------
        # METADATA COVERAGE
        # ----------------------------------------------------
        st.markdown(
            "### 🔎 Geographic Metadata Coverage"
        )

        q1, q2, q3 = st.columns(3)

        q1.metric(
            "Coordinate-bearing Records",
            f"{n_coord_records:,}"
        )

        q2.metric(
            "Records without Coordinates",
            f"{n_unmapped_records:,}"
        )

        q3.metric(
            "Countries without Coordinates",
            f"{n_unmapped_countries:,}"
        )

        # ----------------------------------------------------
        # COLLECTION DATE COMPLETENESS
        # ----------------------------------------------------
        if "collection_date" in df_geo.columns:

            date_values = (
                df_geo["collection_date"]
                .fillna("")
                .astype(str)
                .str.strip()
                .str.lower()
            )

            missing_dates = date_values.isin(
                ["", "none", "nan", "nat"]
            ).sum()

            date_pct = (
                100 *
                (n_records - missing_dates) /
                n_records
                if n_records else 0
            )

            st.caption(
                f"Collection-date completeness: "
                f"{date_pct:.1f}% "
                f"({n_records - missing_dates:,} of "
                f"{n_records:,} records)."
            )

        # ----------------------------------------------------
        # COUNTRIES WITHOUT COORDINATES
        # ----------------------------------------------------
        mapped_country_set = set(
            mapped_summary["country"]
            if not valid_coords.empty
            else []
        )

        unmapped_country_summary = (
            country_summary[
                ~country_summary["country"].isin(
                    mapped_country_set
                )
            ]
            .sort_values(
                "records",
                ascending=False
            )
        )

        if not unmapped_country_summary.empty:

            with st.expander(
                "🌐 Countries without mapped coordinates"
            ):

                display_unmapped = (
                    unmapped_country_summary
                    .rename(
                        columns={
                            "country": "Country",
                            "records": "Records",
                            "percentage": "Share (%)"
                        }
                    )
                    .copy()
                )

                display_unmapped["Share (%)"] = (
                    display_unmapped["Share (%)"]
                    .round(1)
                )

                st.dataframe(
                    display_unmapped[
                        [
                            "Country",
                            "Records",
                            "Share (%)"
                        ]
                    ],
                    use_container_width=True,
                    hide_index=True
                )

        # ----------------------------------------------------
        # INTERPRETATION
        # ----------------------------------------------------
        st.markdown(
            "#### ⚠️ Interpretation"
        )

        st.caption(
            "The geographic coordinates in this dataset represent "
            "country-level geographic assignments. They should not "
            "be interpreted as precise sampling locations. The map "
            "therefore summarizes the geographic distribution of "
            "available PRRSV records rather than individual collection "
            "sites."
        )

        st.caption(
            "The geographic metadata do not by themselves establish "
            "prevalence, transmission pathways, temporal spread, or "
            "causal geographic associations. Records without coordinate "
            "metadata remain included in the country-level summary."
        )

        # ----------------------------------------------------
        # DATA PROVENANCE
        # ----------------------------------------------------
        latest_geo_ts = (
            df_geo["run_timestamp"].max()
            if "run_timestamp" in df_geo.columns
            else "unavailable"
        )

        st.caption(
            "Data provenance: geo_metadata · "
            f"latest snapshot {latest_geo_ts} · "
            f"{n_records:,} geographic records."
        )

# ============================================================
# 🌌 FOOTER — LIVE PIPELINE METADATA

# ============================================================

st.markdown("<hr>", unsafe_allow_html=True)
st.markdown(f"""
<div style='text-align:center; color:#00f0ff;'>
    <b>🧬 IntelliPRRSV2 — Neon Dashboard</b><br>
    Last Pipeline Run: {last_sync}<br>
    Version: 2.4 · Stable Neon Release<br>
    © {datetime.now().year} IntelliPRRSV Research Group
</div>
""", unsafe_allow_html=True)

