"""Figure 3: Genome-Wide Mutation Landscape (Redesigned).

Uses Fig3_mutations.csv as the validated data source.
Coordinates are alignment-column positions, not ungapped genomic positions.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Load data from the validated source
DATA_FILE = Path(__file__).parent.parent / "figure_data" / "Fig3_mutations.csv"

def generate_fig3():
    """Generate redesigned Figure 3."""
    
    df = pd.read_csv(DATA_FILE)
    
    # Aggregate by position (alignment coordinates)
    df_agg = df.groupby('position', as_index=False)['mutation_frequency'].mean()
    
    print(f"Generating Figure 3: Genome-Wide Mutation Landscape")
    print(f"  Positions (alignment coordinates): {len(df_agg)}")
    print(f"  Range: {df_agg['position'].min()} – {df_agg['position'].max()}")
    print(f"  Mean frequency: {df_agg['mutation_frequency'].mean():.4f}")
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Scatter plot with semi-transparent points
    ax.scatter(df_agg['position'], df_agg['mutation_frequency'],
               s=2, color='#1f77b4', alpha=0.4, rasterized=True)
    
    # Rolling median
    df_sorted = df_agg.sort_values('position')
    window = 200
    rolling_median = df_sorted['mutation_frequency'].rolling(window, center=True, min_periods=1).median()
    ax.plot(df_sorted['position'], rolling_median, 
            color='#d7191c', linewidth=1.5, alpha=0.7,
            label='Rolling median (200 positions)')
    
    # Mean line
    mean_freq = df_agg['mutation_frequency'].mean()
    ax.axhline(y=mean_freq, color='#fdae61', linestyle='--', linewidth=1.5,
               label=f'Mean = {mean_freq:.3f}')
    
    ax.set_xlabel('Alignment Coordinate (nt)', fontweight='bold')
    ax.set_ylabel('Mutation Frequency', fontweight='bold')
    ax.set_title('Figure 3: Genome-Wide Mutation Landscape Across PRRSV Multiple-Sequence Alignment', fontweight='bold')
    ax.set_ylim(0, 1.05)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.1, linestyle='--')
    
    plt.tight_layout()
    
    output_dir = Path(__file__).parent.parent / "outputs" / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_base = output_dir / "fig3_mutation_landscape"
    plt.savefig(f"{output_base}.pdf", dpi=600, bbox_inches='tight')
    plt.savefig(f"{output_base}.png", dpi=600, bbox_inches='tight')
    
    print(f"✅ Figure 3 generated: {output_base}.pdf / .png")
    return fig

if __name__ == "__main__":
    generate_fig3()
