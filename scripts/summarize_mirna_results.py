import pandas as pd
import matplotlib.pyplot as plt
import os

# ----------------------------------------------------------
# Paths
# ----------------------------------------------------------
INPUT_FILE = str(MIRNA_SUMMARY_INPUT)
OUT_DIR = str(MIRNA_RESULTS_DIR)

# ----------------------------------------------------------
# Load data
# ----------------------------------------------------------
if not os.path.exists(INPUT_FILE):
    print(f"❌ File not found: {INPUT_FILE}")
    exit()

df = pd.read_csv(INPUT_FILE)
print(f"✅ Loaded {len(df)} miRNA–virus interactions")

# ----------------------------------------------------------
# Summary 1: Top 10 miRNAs with most viral targets
# ----------------------------------------------------------
top_mirnas = df["miRNA"].value_counts().head(10)
print("\n📊 Top 10 most active miRNAs:")
print(top_mirnas)

plt.figure(figsize=(10, 6))
top_mirnas.plot(kind="bar")
plt.title("Top 10 Most Active Sus scrofa miRNAs Targeting PRRSV")
plt.xlabel("miRNA")
plt.ylabel("Number of Targeted PRRSV Sequences")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "top_mirnas.png"))
plt.close()

# ----------------------------------------------------------
# Summary 2: Binding score distribution
# ----------------------------------------------------------
plt.figure(figsize=(8, 5))
df["binding_score"].hist(bins=20)
plt.title("Distribution of miRNA–PRRSV Binding Scores")
plt.xlabel("Binding Score")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "binding_score_distribution.png"))
plt.close()

# ----------------------------------------------------------
# Summary 3: Strongest interactions
# ----------------------------------------------------------
top_hits = df.sort_values("binding_score", ascending=False).head(20)
top_hits.to_csv(os.path.join(OUT_DIR, "top_strong_interactions.csv"), index=False)

print(f"\n✅ Saved summary plots and tables to {OUT_DIR}")
print("   - top_mirnas.png")
print("   - binding_score_distribution.png")
print("   - top_strong_interactions.csv")

