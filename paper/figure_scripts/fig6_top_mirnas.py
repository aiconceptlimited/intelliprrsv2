#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

os.makedirs('../figures', exist_ok=True)

# Load frozen data
df = pd.read_csv('../figure_data/Fig6_miRNA_REPRODUCIBLE.csv')

# Create figure
fig, ax = plt.subplots(figsize=(10, 6))

# Horizontal bar plot
colors = plt.cm.plasma(np.linspace(0.3, 0.9, len(df)))
bars = ax.barh(df['miRNA'], df['interaction_count'], color=colors)

# Add value labels
for bar, val in zip(bars, df['interaction_count']):
    ax.text(val + 1000, bar.get_y() + bar.get_height()/2.,
            f'{val:,}', va='center', fontsize=9)

# Labels and title
ax.set_xlabel('Interaction Count', fontsize=12)
ax.set_ylabel('miRNA', fontsize=12)
ax.set_title('Top 10 miRNAs by Interaction Count', fontsize=14, fontweight='bold')

# Format x-axis with commas
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('../figures/Fig6_v1.0.png', dpi=300, bbox_inches='tight')
plt.savefig('../figures/Fig6_v1.0.pdf', bbox_inches='tight')
plt.close()

print("✅ Figure 6 generated: Fig6_v1.0.png and Fig6_v1.0.pdf")
