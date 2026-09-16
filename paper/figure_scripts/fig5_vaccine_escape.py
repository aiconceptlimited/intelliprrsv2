#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

os.makedirs('../figures', exist_ok=True)

# Load frozen data
df = pd.read_csv('../figure_data/Fig5_escape.csv')

# Create figure
fig, ax = plt.subplots(figsize=(8, 6))

# Bar plot
colors = ['#C73E1D' if x == df['avg_escape'].max() else '#2E86AB' for x in df['avg_escape']]
bars = ax.bar(df['lineage'], df['avg_escape'], color=colors)

# Add value labels on bars
for bar, val in zip(bars, df['avg_escape']):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 0.02,
            f'{val:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

# Labels and title
ax.set_xlabel('Lineage', fontsize=12)
ax.set_ylabel('Average Escape Probability', fontsize=12)
ax.set_title('Vaccine Escape Probability by Lineage', fontsize=14, fontweight='bold')
ax.set_ylim(0, 0.5)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('../figures/Fig5_v1.0.png', dpi=300, bbox_inches='tight')
plt.savefig('../figures/Fig5_v1.0.pdf', bbox_inches='tight')
plt.close()

print("✅ Figure 5 generated: Fig5_v1.0.png and Fig5_v1.0.pdf")
