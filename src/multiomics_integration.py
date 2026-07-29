import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import os

GWAS_PATH = r"D:\MULTIOMICS\data\top_loci_genes.csv"
RNASEQ_PATH = r"D:\MULTIOMICS\data\differential_expression.csv"
OUT_PATH = r"D:\MULTIOMICS\outputs\multiomics_overlap.csv"
FIG_PATH = r"D:\MULTIOMICS\figures\venn_diagram.png"

def load_gwas_genes():
    print("Loading GWAS annotated genes...")
    df = pd.read_csv(GWAS_PATH)
    all_genes = []
    for genes_str in df['Nearby_Genes']:
        if genes_str not in ["N/A", "No protein-coding gene nearby"]:
            all_genes.extend([g.strip() for g in str(genes_str).split(",")])
    unique_genes = set(all_genes)
    print(f"GWAS unique genes: {len(unique_genes)}")
    return unique_genes, df

def load_rnaseq_genes(top_n=200):
    print(f"\nLoading top {top_n} RNA-seq genes by raw p-value...")
    df = pd.read_csv(RNASEQ_PATH)
    df = df.dropna(subset=['Symbol'])
    df_sorted = df.sort_values('p_value').head(top_n)
    unique_genes = set(df_sorted['Symbol'].unique())
    print(f"RNA-seq candidate genes: {len(unique_genes)}")
    return unique_genes, df_sorted

def find_overlap(gwas_genes, rnaseq_genes, gwas_df, rnaseq_df):
    overlap = gwas_genes.intersection(rnaseq_genes)
    print(f"\n=== OVERLAP RESULT ===")
    print(f"Genes in BOTH GWAS and RNA-seq candidate lists: {len(overlap)}")
    print(f"Overlapping genes: {sorted(overlap)}")

    # Build detailed overlap table
    records = []
    for gene in overlap:
        gwas_rows = gwas_df[gwas_df['Nearby_Genes'].str.contains(gene, na=False)]
        rnaseq_row = rnaseq_df[rnaseq_df['Symbol'] == gene]

        gwas_pval = gwas_rows['PVAL'].min() if len(gwas_rows) > 0 else np.nan
        rnaseq_pval = rnaseq_row['p_value'].values[0] if len(rnaseq_row) > 0 else np.nan
        rnaseq_log2fc = rnaseq_row['log2FC'].values[0] if len(rnaseq_row) > 0 else np.nan

        records.append({
            'Gene': gene,
            'GWAS_min_pvalue': gwas_pval,
            'RNAseq_pvalue': rnaseq_pval,
            'RNAseq_log2FC': rnaseq_log2fc
        })

    overlap_df = pd.DataFrame(records).sort_values('GWAS_min_pvalue')
    os.makedirs(r"D:\MULTIOMICS\outputs", exist_ok=True)
    overlap_df.to_csv(OUT_PATH, index=False)
    print(f"\nSaved overlap table to {OUT_PATH}")
    return overlap_df

def plot_venn(gwas_genes, rnaseq_genes, overlap):
    fig, ax = plt.subplots(figsize=(8, 8))
    fig.patch.set_facecolor('#0D1117')

    v = venn2([gwas_genes, rnaseq_genes],
              set_labels=('GWAS Loci Genes', 'RNA-seq Candidate Genes'),
              ax=ax)

    for text in v.set_labels:
        if text:
            text.set_color('white')
            text.set_fontsize(12)
    for text in v.subset_labels:
        if text:
            text.set_color('white')
            text.set_fontsize(14)
            text.set_fontweight('bold')

    if v.get_patch_by_id('10'):
        v.get_patch_by_id('10').set_color('#2E86C1')
        v.get_patch_by_id('10').set_alpha(0.6)
    if v.get_patch_by_id('01'):
        v.get_patch_by_id('01').set_color('#E74C3C')
        v.get_patch_by_id('01').set_alpha(0.6)
    if v.get_patch_by_id('11'):
        v.get_patch_by_id('11').set_color('#9B59B6')
        v.get_patch_by_id('11').set_alpha(0.8)

    ax.set_title(f'Multi-Omics Gene Overlap\nGWAS ∩ RNA-seq = {len(overlap)} genes',
                 color='white', fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()
    os.makedirs(r"D:\MULTIOMICS\figures", exist_ok=True)
    plt.savefig(FIG_PATH, dpi=150, bbox_inches='tight', facecolor='#0D1117')
    print(f"Saved Venn diagram to {FIG_PATH}")
    plt.show()

if __name__ == "__main__":
    gwas_genes, gwas_df = load_gwas_genes()
    rnaseq_genes, rnaseq_df = load_rnaseq_genes(top_n=200)

    overlap_df = find_overlap(gwas_genes, rnaseq_genes, gwas_df, rnaseq_df)
    plot_venn(gwas_genes, rnaseq_genes, set(overlap_df['Gene']))

    print("\nDone.")