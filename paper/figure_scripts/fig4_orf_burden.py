#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

os.makedirs('../figures', exist_ok=True)

# Load frozen data
df = pd.read_csv('../figure_data/Fig4_orf_burden.csv')

# Create figure with twin axes
fig, ax1 = plt.subplots(figsize=(10, 6))

# Bar plot on primary axis
bars = ax1.bar(df['orf'], df['mutation_count'], color='#2E86AB', alpha=0.7)
ax1.set_xlabel('ORF', fontsize=12)
ax1.set_ylabel('Mutation Count', fontsize=12, color='#2E86AB')
ax1.tick_params(axis='y', labelcolor='#2E86AB')

# Line plot on secondary axis
ax2 = ax1.twinx()
ax2.plot(df['orf'], df['avg_frequency'], color='#C73E1D', marker='o', 
         linewidth=2, markersize=8, label='Avg Frequency')
ax2.set_ylabel('Average Mutation Frequency', fontsize=12, color='#C73E1D')
ax2.tick_params(axis='y', labelcolor='#C73E1D')

# Add value labels on bars
for bar, count in zip(bars, df['mutation_count']):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height + 500,
             f'{count:,}', ha='center', va='bottom', fontsize=9)

# Title
ax1.set_title('ORF-Specific Mutation Burden', fontsize=14, fontweight='bold')

# Remove top and right spines
ax1.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)

plt.tight_layout()
plt.savefig('../figures/Fig4_v1.0.png', dpi=300, bbox_inches='tight')
plt.savefig('../figures/Fig4_v1.0.pdf', bbox_inches='tight')
plt.close()

print("✅ Figure 4 generated: Fig4_v1.0.png and Fig4_v1.0.pdf")
