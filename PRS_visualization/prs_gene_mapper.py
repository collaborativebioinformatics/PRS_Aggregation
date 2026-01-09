"""
prs_gene_mapper.py

Helper module to map PRS SNPs to Genes and visualize gene tracks.
Features:
- Auto-downloads UCSC RefGene (hg19/GRCh37) annotations if missing.
- Maps genes to PRS SNPs with a configurable flanking window (default 25kb).
- Includes a progress bar for overlap calculation.
"""
import os
import sys
import requests
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from typing import List, Dict, Optional, Set, Any, Tuple

# -----------------------------------------------------------------------------
# Configuration & Data Loading
# -----------------------------------------------------------------------------

# We use hg19 (GRCh37) as it is the standard for most PRS.
GENE_DATA_URL = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"
GENE_CACHE_FILE = "refGene_hg19.txt.gz"

REFGENE_COLS = [
    "bin", "name", "chrom", "strand", "txStart", "txEnd",
    "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
    "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"
]

_GENE_DF_CACHE: Optional[pd.DataFrame] = None
_TOP_GENES_CACHE: List[Dict[str, str]] = []  # Stores the precomputed top 100 list

def _download_gene_data():
    """Download gene annotation if not present."""
    if not os.path.exists(GENE_CACHE_FILE):
        print(f"\n[gene_mapper] Downloading gene annotations from UCSC (hg19)...")
        try:
            r = requests.get(GENE_DATA_URL, stream=True)
            r.raise_for_status()
            with open(GENE_CACHE_FILE, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            print("[gene_mapper] Download complete.")
        except Exception as e:
            print(f"[gene_mapper] Error downloading gene data: {e}")

def _load_genes() -> pd.DataFrame:
    """Load and cache the gene dataframe."""
    global _GENE_DF_CACHE
    if _GENE_DF_CACHE is not None:
        return _GENE_DF_CACHE

    _download_gene_data()
    
    if not os.path.exists(GENE_CACHE_FILE):
        return pd.DataFrame(columns=REFGENE_COLS)

    try:
        # Load only necessary columns to save memory
        df = pd.read_csv(
            GENE_CACHE_FILE, 
            sep='\t', 
            header=None, 
            names=REFGENE_COLS, 
            usecols=["chrom", "txStart", "txEnd", "name2", "strand"],
            compression='gzip'
        )
        # Clean chromosome names
        df['chrom'] = df['chrom'].astype(str).str.replace('^chr', '', regex=True)
        # Filter to standard chromosomes
        valid_chrs = set([str(i) for i in range(1, 23)] + ['X', 'Y', 'M', 'MT'])
        df = df[df['chrom'].isin(valid_chrs)].copy()
        
        # Keep longest transcript per gene symbol to simplify visualization/counting
        df['length'] = df['txEnd'] - df['txStart']
        df = df.sort_values('length', ascending=False).drop_duplicates(subset=['name2'])
        df = df.sort_values(['chrom', 'txStart'])
        
        _GENE_DF_CACHE = df
        return df
    except Exception as e:
        print(f"[gene_mapper] Error parsing gene file: {e}")
        return pd.DataFrame(columns=REFGENE_COLS)

# -----------------------------------------------------------------------------
# Core Mapping & Progress Logic
# -----------------------------------------------------------------------------

def precompute_top_genes(tracks: List[Any], flank_bp: int = 25000):
    """
    Calculates overlap for all genes and stores the top 100 in a cache.
    Prints a text-based progress bar to stdout.
    """
    global _TOP_GENES_CACHE
    genes = _load_genes()
    if genes.empty:
        return

    # Extract all SNPs from tracks
    all_snps = []
    for t in tracks:
        if not t.df.empty:
            all_snps.append(t.df[['hm_chr', 'hm_pos']].copy())
    
    if not all_snps:
        return
    
    snps = pd.concat(all_snps).drop_duplicates()
    snps['hm_pos'] = pd.to_numeric(snps['hm_pos'], errors='coerce')
    snps = snps.dropna()

    unique_chrs = snps['hm_chr'].unique()
    total_steps = len(unique_chrs)
    
    gene_counts = {}
    
    print(f"\n[info] Calculating gene coverage (±{flank_bp//1000}kb) for {len(snps)} unique SNPs...")

    # Progress loop
    for i, chrom in enumerate(unique_chrs):
        # Progress Bar
        percent = (i + 1) / total_steps
        bar_len = 30
        filled_len = int(bar_len * percent)
        bar = '=' * filled_len + '-' * (bar_len - filled_len)
        sys.stdout.write(f'\r[{bar}] {int(percent * 100)}% (processing chr{chrom})')
        sys.stdout.flush()

        # Calculation
        chrom_str = str(chrom)
        genes_chr = genes[genes['chrom'] == chrom_str]
        snps_chr = snps[snps['hm_chr'] == chrom_str]
        
        if genes_chr.empty or snps_chr.empty:
            continue
            
        g_starts = genes_chr['txStart'].values - flank_bp
        g_ends = genes_chr['txEnd'].values + flank_bp
        g_names = genes_chr['name2'].values
        
        s_pos = np.sort(snps_chr['hm_pos'].values)
        
        # Vectorized search
        idx_start = np.searchsorted(s_pos, g_starts, side='left')
        idx_end = np.searchsorted(s_pos, g_ends, side='right')
        counts = idx_end - idx_start
        
        for name, count in zip(g_names, counts):
            if count > 0:
                gene_counts[name] = gene_counts.get(name, 0) + count

    sys.stdout.write('\n')  # End progress bar line
    print("[info] Gene coverage calculation finished.")

    sorted_genes = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)[:100]
    _TOP_GENES_CACHE = [
        {'label': f"{sym} ({count} SNPs)", 'value': sym}
        for sym, count in sorted_genes
    ]

def get_top_genes_options() -> List[Dict[str, str]]:
    """Returns the cached top genes list for the UI dropdown."""
    return _TOP_GENES_CACHE

def genes_with_snps_in_window(tracks: List[Any], chrom: str, xmin: int, xmax: int) -> List[Dict[str, str]]:
    """Returns genes located within the specific viewing window."""
    genes = _load_genes()
    if genes.empty:
        return []
        
    chrom = str(chrom).replace('chr', '')
    # Overlap logic: gene_start < window_end AND gene_end > window_start
    mask = (genes['chrom'] == chrom) & (genes['txStart'] <= xmax) & (genes['txEnd'] >= xmin)
    subset = genes[mask].sort_values('txStart')
    
    return [{'label': row['name2'], 'value': row['name2']} for _, row in subset.iterrows()]

def fetch_gene_by_id(gene_symbol: str) -> Optional[Dict[str, Any]]:
    genes = _load_genes()
    row = genes[genes['name2'] == gene_symbol]
    if row.empty:
        return None
    dat = row.iloc[0]
    return {
        "seq_region_name": dat['chrom'],
        "start": int(dat['txStart']),
        "end": int(dat['txEnd']),
        "strand": dat['strand'],
        "id": dat['name2']
    }

# -----------------------------------------------------------------------------
# Visualization
# -----------------------------------------------------------------------------

def make_gene_track_figure(chrom: str, xmin: int, xmax: int, only_gene_ids: Optional[Set[str]] = None) -> go.Figure:
    fig = go.Figure()
    genes = _load_genes()
    chrom = str(chrom).replace('chr', '')
    
    if only_gene_ids:
        subset = genes[(genes['chrom'] == chrom) & (genes['name2'].isin(only_gene_ids))]
    else:
        subset = genes[(genes['chrom'] == chrom) & (genes['txEnd'] >= xmin) & (genes['txStart'] <= xmax)]

    if subset.empty:
        fig.update_layout(
            title="No genes found in region",
            xaxis=dict(range=[xmin, xmax]),
            yaxis=dict(showticklabels=False, showgrid=False),
            height=100, template="plotly_white"
        )
        return fig

    # Stagger genes Y=0, 1, 2
    y_levels = 3
    
    for i, (_, row) in enumerate(subset.iterrows()):
        g_start, g_end, g_name, strand = row['txStart'], row['txEnd'], row['name2'], row['strand']
        y_pos = i % y_levels
        color = "#d62728" if only_gene_ids and g_name in only_gene_ids else "navy" # Red if selected, Navy otherwise
        
        # Gene Line
        fig.add_trace(go.Scatter(
            x=[g_start, g_end], y=[y_pos, y_pos], mode="lines",
            line=dict(color=color, width=3),
            hoverinfo="text", hovertext=f"{g_name} ({strand})", showlegend=False
        ))
        
        # Label
        fig.add_trace(go.Scatter(
            x=[(g_start+g_end)/2], y=[y_pos + 0.3], mode="text",
            text=[g_name], textposition="top center",
            textfont=dict(size=10, color=color), showlegend=False
        ))
        
        # TSS Marker
        tss = g_start if strand == '+' else g_end
        symbol = "triangle-right" if strand == '+' else "triangle-left"
        fig.add_trace(go.Scatter(
            x=[tss], y=[y_pos], mode="markers",
            marker=dict(symbol=symbol, size=8, color=color), hoverinfo="skip", showlegend=False
        ))

    fig.update_layout(
        title=dict(text=f"Genes on chr{chrom} (RefSeq hg19)", font=dict(size=12), y=0.95),
        xaxis=dict(range=[xmin, xmax], showgrid=True, tickformat=","),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-0.5, y_levels+1]),
        height=180, margin=dict(l=30, r=20, t=30, b=20), template="plotly_white"
    )
    return fig
