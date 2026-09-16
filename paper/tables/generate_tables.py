"""Generate all tables from IntelliPRRSV2 exports."""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent))
from statistics import (
    mutation_summary, orf_summary, lineage_summary,
    mirna_summary, vaccine_summary, geographic_summary
)
from config import TABLES_DIR


def generate_table_1():
    """Table 1: Dataset characteristics."""
    
    lineage = lineage_summary()
    geo = geographic_summary()
    
    table = f"""\\begin{{table}}[htbp]
\\centering
\\caption{{Characteristics of the PRRSV dataset analyzed in this study.}}
\\begin{{tabular}}{{lrr}}
\\hline
\\textbf{{Characteristic}} & \\textbf{{Category}} & \\textbf{{Count}} \\\\
\\hline
Total sequences & - & {lineage['total_sequences']} \\\\
Countries represented & - & {geo['total_countries']} \\\\
\\hline
\\multicolumn{{3}}{{l}}{{\\textbf{{Lineage Distribution}}}} \\\\
\\hline
"""
    
    for l in lineage['lineages']:
        table += f"    {l['lineage']} & - & {l['count']} \\\\\n"
    
    table += f"""
\\hline
\\multicolumn{{3}}{{l}}{{\\textbf{{Geographic Distribution}}}} \\\\
\\hline
"""
    
    for c in geo['top_5']:
        table += f"    {c['country']} & - & {c['count']} \\\\\n"
    
    table += """
\\hline
\\end{tabular}
\\end{table}
"""
    
    return table


def generate_table_2():
    """Table 2: Top mutation hotspots."""
    
    stats = mutation_summary()
    top = stats['top_hotspots'][:20]
    
    table = """
\\begin{table}[htbp]
\\centering
\\caption{Top 20 mutation hotspots identified in PRRSV genomes.}
\\begin{tabular}{rrrr}
\\hline
\\textbf{Rank} & \\textbf{Position} & \\textbf{Frequency} & \\textbf{Region} \\\\
\\hline
"""
    
    for i, row in enumerate(top, 1):
        table += f"    {i} & {int(row['position'])} & {row['mutation_frequency']:.4f} & {row['region']} \\\\\n"
    
    table += """
\\hline
\\end{tabular}
\\end{table}
"""
    
    return table


def generate_table_3():
    """Table 3: ORF mutation burden."""
    
    stats = orf_summary()
    
    table = """
\\begin{table}[htbp]
\\centering
\\caption{Mutation burden across PRRSV genomic regions.}
\\begin{tabular}{lr}
\\hline
\\textbf{ORF/Gene Region} & \\textbf{Mean Mutation Frequency} \\\\
\\hline
"""
    
    for row in stats['all_regions']:
        table += f"    {row['region']} & {row['mean_frequency']:.2f} \\\\\n"
    
    table += """
\\hline
\\end{tabular}
\\end{table}
"""
    
    return table


def generate_table_4():
    """Table 4: Top miRNA predictions."""
    
    stats = mirna_summary()
    
    table = """
\\begin{table}[htbp]
\\centering
\\caption{Top 10 predicted miRNA-PRRSV interactions.}
\\begin{tabular}{lrr}
\\hline
\\textbf{miRNA ID} & \\textbf{Binding Score} & \\textbf{Energy (kcal/mol)} \\\\
\\hline
"""
    
    for row in stats['top_10']:
        table += f"    {row['mirna']} & {row['binding_score']:.3f} & {row['energy_kcal']:.2f} \\\\\n"
    
    table += """
\\hline
\\end{tabular}
\\end{table}
"""
    
    return table


def generate_table_5():
    """Table 5: Vaccine escape summary."""
    
    stats = vaccine_summary()
    
    table = """
\\begin{table}[htbp]
\\centering
\\caption{Vaccine escape risk summary.}
\\begin{tabular}{lrr}
\\hline
\\textbf{Risk Category} & \\textbf{Count} & \\textbf{Percentage} \\\\
\\hline
"""
    
    for category in ['Low', 'Moderate', 'High']:
        count = stats['risk_distribution'].get(category, 0)
        pct = stats['risk_percentages'].get(category, 0)
        table += f"    {category} & {count} & {pct:.1f}\\% \\\\\n"
    
    table += f"""
\\hline
\\textbf{{Total}} & {stats['total_comparisons']} & 100.0\\% \\\\
\\hline
\\end{{tabular}}
\\end{{table}}
"""
    
    return table


def generate_all_tables():
    """Generate all tables."""
    
    tables = {
        'table1_dataset.tex': generate_table_1,
        'table2_hotspots.tex': generate_table_2,
        'table3_orf_burden.tex': generate_table_3,
        'table4_mirnas.tex': generate_table_4,
        'table5_vaccine.tex': generate_table_5
    }
    
    print("  Generating tables...")
    for filename, generator in tables.items():
        output_path = TABLES_DIR / filename
        with open(output_path, 'w') as f:
            f.write(generator())
        print(f"    ✓ {filename}")
    
    return tables


if __name__ == "__main__":
    generate_all_tables()
