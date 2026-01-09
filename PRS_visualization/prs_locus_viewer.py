#!/usr/bin/env python3
"""
PRS Locus Viewer (Dash)
-----------------------
A small Dash app that reads 2–10 PRS/PGS scoring files (PGS Catalog-style .txt)
and lets you pick an rsID + window size to visualize SNP locations and effect
weights across multiple PRSs.

Key features:
- Robust reader that skips top annotation rows (lines starting with '#')
- Extracts PRS/PGS id/name from header metadata if present
- Union rsID search (type to search; doesn't render millions of options at once)
- Locus plot: one horizontal "track" per PRS (| ticks). Selected SNP is colored
  by effect_weight and highlighted with a red triangle. Nearby SNPs are gray.
- Optional: color ALL SNPs in window by effect_weight.
- Optional: restrict to SNPs shared across all displayed PRSs.

Run:
  python prs_locus_viewer.py PGS1.txt PGS2.txt [PGS3.txt ...]
or just:
  python prs_locus_viewer.py
and you'll get a file picker dialog.

Dependencies:
  pip install dash pandas plotly numpy
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import requests

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Dash, Input, Output, State, ctx, dcc, html

import prs_info_parser as info_parser
import prs_gene_mapper as gene_mapper

# ----------------------------
# IO helpers
# ----------------------------

_META_KV_RE = re.compile(r"^#\s*([^=\s]+)\s*=\s*(.*)\s*$")
_WINDOW_KB_OPTIONS = [10, 20, 50, 100, 250, 500, 1000, 2000, 10000]
_PGS_ID_RE = re.compile(r"PGS\d{6}", re.IGNORECASE)
_MEGA_PEARL = [
    "#F3C7D4",  # soft rose
    "#C3D7F0",  # pale blue
    "#B8A2D1",  # lavender
    "#FBE5C8",  # cream
    "#8FC9C5",  # teal
    "#F08A9A",  # vivid pink
]
_SET3 = [
    "#8DD3C7",
    "#FFFFB3",
    "#BEBADA",
    "#FB8072",
    "#80B1D3",
    "#FDB462",
    "#B3DE69",
    "#FCCDE5",
    "#D9D9D9",
    "#BC80BD",
    "#CCEBC5",
    "#FFED6F",
]
_ANCESTRY_PRIORITY = ["AFR", "AMR", "ASN", "EAS", "EUR", "GME", "SAS", "MID", "OCE", "OTH", "UNK"]


def _build_ancestry_palette() -> Dict[str, str]:
    palette_cycle = _MEGA_PEARL + _SET3
    colors: Dict[str, str] = {}
    for idx, code in enumerate(_ANCESTRY_PRIORITY):
        colors[code] = palette_cycle[idx % len(palette_cycle)]
    return colors


_ANCESTRY_COLORS = _build_ancestry_palette()
_MARKER_SIZE_COLOR = 31  # reduced ~30%
_MARKER_SIZE_GRAY = 22   # reduced ~30%
_MARKER_SIZE_SELECTED = int(round(_MARKER_SIZE_COLOR * 1.5))  # selected rsID stands out


def _normalize_chr(chrom: Optional[str]) -> Optional[str]:
    if not chrom:
        return None
    c = re.sub(r"\s+", "", str(chrom)).lower()
    if c.startswith("chr"):
        c = c[3:]
    if c in ("all", "*"):
        return None
    return c


def _open_text_maybe_gzip(path: str):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def _build_pgs_score_url(pgs_id: Optional[str]) -> Optional[str]:
    if not pgs_id:
        return None
    pid = pgs_id.upper().strip()
    if not _PGS_ID_RE.match(pid):
        return None
    return f"https://www.pgscatalog.org/score/{pid}/"


def _format_count_short(value: Optional[int], max_len: int = 6) -> str:
    """
    Format counts so displayed strings stay compact (<= max_len-ish characters).
    Example: 100000, ~10.6k, ~1.2M
    """
    if value is None:
        return "n/a"
    val = int(value)
    plain = str(val)
    if len(plain) <= max_len:
        return plain

    suffixes = [(1_000_000_000, "B"), (1_000_000, "M"), (1_000, "k")]
    for threshold, suffix in suffixes:
        if val >= threshold:
            scaled = val / threshold
            # Try decreasing precision until it fits the target length budget
            for decimals in [2, 1, 0]:
                num_part = f"{scaled:.{decimals}f}".rstrip("0").rstrip(".")
                approx = f"~{num_part}{suffix}"
                if len(approx) <= max_len:
                    return approx
            return f"~{int(round(scaled))}{suffix}"

    return plain[:max_len]


def read_pgs_scoring_file(path: str) -> Tuple[Dict[str, str], pd.DataFrame]:
    """
    Read a PRS/PGS scoring file and return: (metadata_dict, dataframe)

    Supports PGS Catalog-style files (comment/metadata lines start with '#') and
    more generic "annotation rows on top" tables by *detecting the header line*.

    Output dataframe columns (normalized):
      - rsID (string)
      - hm_chr (string, without 'chr' prefix)
      - hm_pos (int)
      - effect_weight (float)
    """
    meta: Dict[str, str] = {}

    header_idx: Optional[int] = None
    header_line: Optional[str] = None

    # 1) Scan the top of the file:
    #    - collect metadata key=value pairs from lines starting with '#'
    #    - detect the first plausible header line (contains required-ish columns)
    try:
        with _open_text_maybe_gzip(path) as f:
            for idx, line in enumerate(f):
                if idx > 5000:  # safety for weird files
                    break

                raw = line.rstrip("\n")
                if not raw.strip():
                    continue

                if raw.startswith("#"):
                    m = _META_KV_RE.match(raw)
                    if m:
                        k, v = m.group(1), m.group(2)
                        meta[k] = v
                    continue

                # Candidate header (tab or whitespace separated)
                tokens = re.split(r"\t|\s+", raw.strip())
                low = {t.lower() for t in tokens}

                has_rsid = any(x in low for x in ["rsid", "rs_id", "hm_rsid", "variant_id", "snp", "snpid"])
                has_chr = any(x in low for x in ["hm_chr", "chr", "chr_name", "chrom", "chromosome", "hm_chrom"])
                has_pos = any(x in low for x in ["hm_pos", "pos", "position", "chr_position", "bp", "hm_position"])
                has_w = any(x in low for x in ["effect_weight", "weight", "beta", "effect", "effect_size", "effectsize"])

                if has_rsid and has_chr and has_pos and has_w:
                    header_idx = idx
                    header_line = raw
                    break
    except Exception:
        # metadata/header detection is best-effort; we'll fall back to pandas parsing
        header_idx = None
        header_line = None

    # 2) Read the table
    if header_idx is not None:
        # Decide delimiter from the detected header line
        sep = "\t" if ("\t" in (header_line or "")) else r"\s+"
        try:
            df = pd.read_csv(
                path,
                sep=sep,
                skiprows=header_idx,
                dtype=str,
                engine="python" if sep != "\t" else None,
                low_memory=False,
            )
        except Exception:
            # last resort: try tab with comment filtering
            df = pd.read_csv(path, sep="\t", comment="#", dtype=str, low_memory=False)
    else:
        # Common case for PGS Catalog: comment rows start with '#'
        try:
            df = pd.read_csv(path, sep="\t", comment="#", dtype=str, low_memory=False)
        except Exception:
            df = pd.read_csv(path, sep=r"\s+", comment="#", dtype=str, engine="python", low_memory=False)

    if df.empty:
        raise ValueError(f"No data rows detected in file: {path}")

    # 3) Normalize column names (case-insensitive)
    colmap = {c.lower(): c for c in df.columns}

    def pick_col(*candidates: str) -> Optional[str]:
        for c in candidates:
            if c.lower() in colmap:
                return colmap[c.lower()]
        return None

    rs_col = pick_col("rsid", "rs_id", "rs", "hm_rsid", "variant_id", "snpid", "snp")
    chr_col = pick_col("hm_chr", "chrom", "chr", "chr_name", "chromosome", "hm_chrom")
    pos_col = pick_col("hm_pos", "pos", "position", "chr_position", "bp", "hm_position")
    w_col = pick_col("effect_weight", "weight", "beta", "effect", "effectsize", "effect_size")

    missing = [name for name, col in [("rsID", rs_col), ("hm_chr", chr_col), ("hm_pos", pos_col), ("effect_weight", w_col)] if col is None]
    if missing:
        raise ValueError(
            f"File {os.path.basename(path)} is missing required column(s): {missing}. "
            f"Columns found: {list(df.columns)}"
        )

    out = df[[rs_col, chr_col, pos_col, w_col]].copy()
    out.columns = ["rsID", "hm_chr", "hm_pos", "effect_weight"]

    # 4) Clean/convert
    out["rsID"] = out["rsID"].astype(str).str.strip()
    out["hm_chr"] = out["hm_chr"].astype(str).str.strip().str.replace("^chr", "", regex=True)
    out["hm_pos"] = pd.to_numeric(out["hm_pos"], errors="coerce").astype("Int64")
    out["effect_weight"] = pd.to_numeric(out["effect_weight"], errors="coerce")

    out = out.dropna(subset=["rsID", "hm_chr", "hm_pos", "effect_weight"]).copy()
    out["hm_chr"] = out["hm_chr"].astype(str)
    out["hm_pos"] = out["hm_pos"].astype(int)

    # Some files may contain duplicate rsIDs; keep the first occurrence
    out = out.sort_values(["hm_chr", "hm_pos"]).drop_duplicates(subset=["rsID"], keep="first").reset_index(drop=True)
    return meta, out


def _get_meta_value(meta: Dict[str, str], *keys: str) -> str:
    if not meta:
        return ""
    lowered = {k.lower(): v for k, v in meta.items()}
    for key in keys:
        if key.lower() in lowered:
            return str(lowered[key.lower()]).strip()
    return ""


def _fetch_publication_info(pgs_id: Optional[str]) -> Tuple[Optional[str], Optional[str]]:
    if not pgs_id:
        return None, None

    try:
        score_url = f"https://www.pgscatalog.org/rest/score/{pgs_id}"
        resp = requests.get(score_url, timeout=6)
        if not resp.ok:
            return None, None
        score = resp.json()
        pub = score.get("publication") or {}
        author = pub.get("firstauthor") or "Unknown"
        journal = pub.get("journal") or "Unknown Journal"
        year = str(pub.get("date_publication", "0000")).split("-")[0]
        citation = f"{author} et al. {journal} ({year})"
        doi = pub.get("doi")
        pub_url = None
        if doi:
            pub_url = doi if str(doi).startswith("http") else f"https://doi.org/{doi}"
        return citation, pub_url
    except Exception:
        return None, None


@dataclass(frozen=True)
class ScoreTrack:
    label: str          # e.g. "PGS000125"
    name: str           # e.g. "Qi_T2D_2017"
    path: str
    df: pd.DataFrame
    total_snps: int
    meta: Dict[str, str]
    trait_reported: str
    pgs_id: str
    gwas_ancestry: Optional[Dict[str, float]]
    gwas_count: Optional[int]
    dev_count: Optional[int]
    publication_citation: Optional[str]
    publication_url: Optional[str]
    pgs_score_url: Optional[str]


def _infer_track_label(meta: Dict[str, str], path: str) -> Tuple[str, str]:
    """
    Returns (label, name) for display.
    """
    pgs_id = _get_meta_value(meta, "pgs_id", "pgs", "PGS_ID")
    pgs_name = _get_meta_value(meta, "pgs_name", "PGS_NAME")
    if pgs_id:
        return pgs_id, pgs_name
    # fallback: filename stem
    stem = os.path.splitext(os.path.basename(path))[0]
    return stem, ""


def _symmetric_cmin_cmax(values: Iterable[float]) -> Tuple[float, float]:
    vals = np.asarray([v for v in values if v is not None], dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return (-1.0, 1.0)
    m = float(np.max(np.abs(vals)))
    if m == 0:
        m = 1.0
    return (-m, m)


def _resolve_locus(tracks: List[ScoreTrack], rsid: str) -> Tuple[Optional[str], Optional[int]]:
    rsid = (rsid or "").strip()
    if not rsid:
        return None, None
    for tr in tracks:
        hit = tr.df.loc[tr.df["rsID"] == rsid, ["hm_chr", "hm_pos"]]
        if not hit.empty:
            return str(hit.iloc[0]["hm_chr"]), int(hit.iloc[0]["hm_pos"])
    return None, None


def _infer_trait(meta: Dict[str, str]) -> str:
    return info_parser.infer_trait(meta)


def _summarize_traits(tracks: List[ScoreTrack]) -> str:
    traits = [t.trait_reported for t in tracks if t.trait_reported]
    if not traits:
        return "Trait: Unknown"
    counts: Dict[str, int] = {}
    for tr in traits:
        counts[tr] = counts.get(tr, 0) + 1
    if len(counts) == 1:
        return f"Trait: {next(iter(counts))}"
    total = sum(counts.values())
    parts = [f"{trait} ({round(100 * count / total)}%)" for trait, count in sorted(counts.items())]
    return "Traits: " + " / ".join(parts)


def _extract_ancestry_from_meta(
    meta: Dict[str, str],
    pgs_id: Optional[str] = None,
) -> Tuple[Optional[Dict[str, float]], Optional[int], Optional[int]]:
    return info_parser.extract_ancestry_from_meta(meta, pgs_id=pgs_id)


def _order_ancestry_labels(dist: Dict[str, float]) -> List[str]:
    priority = _ANCESTRY_PRIORITY
    labels_upper = {k.upper(): k for k in dist.keys()}
    ordered = [labels_upper[p] for p in priority if p in labels_upper]
    remaining = [labels_upper[k] for k in labels_upper if k not in priority]
    return ordered + sorted(remaining)


def _ancestry_color(label: str) -> str:
    key = label.upper()
    if key in _ANCESTRY_COLORS:
        return _ANCESTRY_COLORS[key]
    palette_cycle = _MEGA_PEARL + _SET3
    idx = abs(hash(key)) % len(palette_cycle)
    return palette_cycle[idx]


def _make_ancestry_bar(track: ScoreTrack) -> Optional[go.Figure]:
    dist = track.gwas_ancestry or {}
    if not dist:
        return None
    labels = _order_ancestry_labels(dist)
    bars = []
    for label in labels:
        val = float(dist[label])
        bars.append(
            go.Bar(
                x=[val],
                y=[""],
                name=label,
                orientation="h",
                marker=dict(color=_ancestry_color(label)),
                hovertemplate=f"{label}: {val:.1f}%<extra></extra>",
                width=0.4,
            )
        )
    total = sum(dist.values())
    xmax = max(100.0, total) * 1.02
    fig = go.Figure(data=bars)
    fig.update_layout(
        barmode="stack",
        height=120,
        margin=dict(l=10, r=30, t=30, b=10),
        title=dict(text="Ancestry distribution", font=dict(size=12)),
        xaxis=dict(automargin=True, range=[0, xmax], tickformat="~s", nticks=4, tickangle=0, tickfont=dict(size=10)),
        yaxis=dict(showticklabels=False, showgrid=False),
        template="plotly_white",
        showlegend=True,
        legend=dict(
            orientation="v",
            yanchor="middle",
            y=0.5,
            xanchor="left",
            x=1.02,
            font=dict(size=9),
        ),
    )
    return fig


def _make_single_bar(
    value: Optional[int],
    title: str,
    xmax_override: Optional[float] = None,
    text_formatter=None,
) -> go.Figure:
    val = int(value) if value is not None else 0
    xmax = max(val * 1.25, 1) if xmax_override is None else xmax_override
    display_text = text_formatter(val) if (text_formatter and value is not None) else ("n/a" if value is None else str(val))
    fig = go.Figure(
        data=[
            go.Bar(
                x=[val],
                y=[""],
                orientation="h",
                text=[display_text],
                textposition="outside",
                width=0.4,
                textfont=dict(size=10),
                marker=dict(color="rgba(70,130,180,0.7)"),
            )
        ]
    )
    fig.update_layout(
        height=120,
        margin=dict(l=10, r=30, t=30, b=10),
        title=dict(text=title, font=dict(size=12)),
        xaxis=dict(automargin=True, range=[0, xmax], tickformat="~s", nticks=4, tickangle=0, tickfont=dict(size=10)),
        yaxis=dict(showticklabels=False, showgrid=False),
        template="plotly_white",
    )
    return fig


def _add_vline(fig: go.Figure, x: int, n_tracks: int, gap: float = 0.25) -> None:
    if n_tracks <= 0:
        return
    
    # Define the full vertical range (from bottom of first track to top of last)
    y_min = -0.5
    y_max = n_tracks - 0.5

    # Draw ONE continuous dashed line
    fig.add_shape(
        type="line",
        x0=x,
        x1=x,
        y0=y_min,
        y1=y_max,
        # Keep dash="dot" for the dashed style
        line=dict(color="red", width=1, dash="dot"), 
        layer="below",
    )


# ----------------------------
# Plot builder
# ----------------------------

def make_locus_figure(
    tracks: List[ScoreTrack],
    rsid: str,
    window_kb: int,
    show_track_labels: bool = True,
    color_all_snps: bool = False,
    shared_only: bool = False,
    locus_override: Optional[Tuple[str, int]] = None,
) -> Tuple[go.Figure, str]:
    """
    Returns (figure, status_text).
    """
    rsid = (rsid or "").strip()
    if not rsid and locus_override is None:
        fig = go.Figure()
        fig.update_layout(
            height=220 + 35 * len(tracks),
            margin=dict(l=30, r=20, t=40, b=40),
            title="Select an rsID to begin",
        )
        return fig, "Waiting for rsID input."

    # Determine locus (chr, pos) using first track that contains rsid
    if locus_override is not None:
        locus_chr, locus_pos = locus_override
    else:
        locus_chr, locus_pos = _resolve_locus(tracks, rsid)

    if locus_chr is None or locus_pos is None:
        fig = go.Figure()
        fig.update_layout(
            height=220 + 35 * len(tracks),
            margin=dict(l=30, r=20, t=40, b=40),
            title=f"rsID not found: {rsid}",
        )
        return fig, f"rsID '{rsid}' was not found in any loaded PRS file."

    half_window = int(window_kb) * 1000
    xmin, xmax = locus_pos - half_window, locus_pos + half_window

    # Precompute shared-only rsID set (within the window & chromosome)
    window_sets = []
    for tr in tracks:
        d = tr.df
        w = d[(d["hm_chr"] == locus_chr) & (d["hm_pos"].between(xmin, xmax))]
        window_sets.append(set(w["rsID"].tolist()))
    shared_rsids = set.intersection(*window_sets) if (shared_only and window_sets) else None

    # Determine color scale range
    if color_all_snps or not rsid:
        all_vals = []
        for tr in tracks:
            d = tr.df
            w = d[(d["hm_chr"] == locus_chr) & (d["hm_pos"].between(xmin, xmax))]
            if shared_rsids is not None:
                w = w[w["rsID"].isin(shared_rsids)]
            all_vals.extend(w["effect_weight"].astype(float).tolist())
        cmin, cmax = _symmetric_cmin_cmax(all_vals)
    else:
        sel_vals = []
        for tr in tracks:
            hit = tr.df.loc[tr.df["rsID"] == rsid, "effect_weight"]
            if not hit.empty:
                sel_vals.append(float(hit.iloc[0]))
        cmin, cmax = _symmetric_cmin_cmax(sel_vals)

    fig = go.Figure()

    y_vals = list(range(len(tracks)))
    y_ticktext = [tr.label if tr.name == "" else tr.label for tr in tracks]

    if rsid:
        status_lines = [f"Selected locus: chr{locus_chr}:{locus_pos:,} (±{window_kb} kb)"]
    else:
        status_lines = [f"Gene window: chr{locus_chr}:{xmin:,}-{xmax:,} (±{window_kb} kb)"]

    selected_traces: List[go.Scatter] = []
    # Add one track at a time
    for i, tr in enumerate(tracks):
        d = tr.df
        w = d[(d["hm_chr"] == locus_chr) & (d["hm_pos"].between(xmin, xmax))].copy()

        if shared_rsids is not None:
            w = w[w["rsID"].isin(shared_rsids)]

        n_snps = w.shape[0]
        status_lines.append(f"{tr.label}: {n_snps} SNPs in window")

        if n_snps == 0:
            continue

        # Grey SNPs (context)
        w_other = w[w["rsID"] != rsid]
        if not w_other.empty:
            if color_all_snps:
                fig.add_trace(go.Scattergl(
                    x=w_other["hm_pos"],
                    y=[i] * len(w_other),
                    mode="markers",
                    marker=dict(
                        symbol="line-ns",
                        size=_MARKER_SIZE_COLOR,
                        color=w_other["effect_weight"].astype(float),
                        coloraxis="coloraxis",
                        line=dict(width=0),
                    ),
                    hovertemplate=(
                        f"{tr.label}<br>"
                        "rsID=%{customdata[0]}<br>"
                        "pos=%{x:,}<br>"
                        "effect_weight=%{customdata[1]:.4g}<extra></extra>"
                    ),
                    customdata=np.column_stack([w_other["rsID"], w_other["effect_weight"].astype(float)]),
                    showlegend=False,
                ))
            else:
                fig.add_trace(go.Scattergl(
                    x=w_other["hm_pos"],
                    y=[i] * len(w_other),
                    mode="markers",
                    marker=dict(
                        symbol="line-ns",
                        size=_MARKER_SIZE_GRAY,
                        color="rgba(120,120,120,0.55)",
                        line=dict(width=1),
                    ),
                    hovertemplate=(
                        f"{tr.label}<br>"
                        "rsID=%{customdata[0]}<br>"
                        "pos=%{x:,}<br>"
                        "effect_weight=%{customdata[1]:.4g}<extra></extra>"
                    ),
                    customdata=np.column_stack([w_other["rsID"], w_other["effect_weight"].astype(float)]),
                    showlegend=False,
                ))

        # Selected rsID in this track (if present)
        if rsid:
            w_sel = w[w["rsID"] == rsid]
            if not w_sel.empty:
                ew = float(w_sel.iloc[0]["effect_weight"])
                selected_traces.append(go.Scattergl(
                    x=[int(w_sel.iloc[0]["hm_pos"])],
                    y=[i],
                    mode="markers",
                    marker=dict(
                        symbol="line-ns",
                        size=_MARKER_SIZE_GRAY,    # Height: matches background (22)
                        color=[ew],                # Fill: maps to effect weight
                        coloraxis="coloraxis",
                        line=dict(
                            width=2,               # Width: 2px (thicker than background's 1px)
                            color=[ew],            # Border Color: MUST match effect weight
                            coloraxis="coloraxis"  # Ensure border uses the heatmap scale
                        ),
                    ),
                    hovertemplate=(
                        f"{tr.label}<br>"
                        f"**SELECTED** {rsid}<br>"
                        "pos=%{x:,}<br>"
                        f"effect_weight={ew:.6g}<extra></extra>"
                    ),
                    showlegend=False,
                ))
                # keep selected SNP as a thicker "|" marker only

    # Vertical reference line at selected position (with gaps for track markers)
    if rsid:
        _add_vline(fig, locus_pos, len(tracks))
        for tr_sel in selected_traces:
            fig.add_trace(tr_sel)
        fig.add_trace(go.Scatter(
            x=[locus_pos],
            y=[len(tracks) - 0.45],
            mode="markers",
            marker=dict(symbol="triangle-up", size=10, color="red"),
            hovertemplate=f"Selected locus chr{locus_chr}:{locus_pos:,}<extra></extra>",
            showlegend=False,
        ))

    # Layout polish
    fig.update_layout(
        height=220 + 35 * max(1, len(tracks)),
        margin=dict(l=120 if show_track_labels else 30, r=20, t=50, b=55),
        title="PRS locus view",
        xaxis=dict(title=f"Genomic position (bp) on chr{locus_chr}", range=[xmin, xmax], tickformat=","),
        yaxis=dict(
            title="",
            tickmode="array",
            tickvals=y_vals,
            ticktext=y_ticktext if show_track_labels else ["" for _ in y_vals],
            autorange="reversed",
        ),
        coloraxis=dict(
            colorscale="Spectral",
            reversescale=True,
            cmin=cmin,
            cmax=cmax,
            colorbar=dict(title="effect_weight"),
        ),
        template="plotly_white",
    )

    return fig, " | ".join(status_lines)


# ----------------------------
# Dash app factory
# ----------------------------

def create_app(tracks: List[ScoreTrack], chrom_filter: Optional[str] = None) -> Dash:
    # Global rsID universe for searching
    all_rsids = sorted(set().union(*[set(tr.df["rsID"].tolist()) for tr in tracks]))
    trait_summary = _summarize_traits(tracks)
    chrom_label = f"chr{chrom_filter}" if chrom_filter else "all chromosomes"

    app = Dash(__name__)
    app.title = "PRS Locus Viewer"

    app.layout = html.Div(
        style={"maxWidth": "1200px", "margin": "0 auto", "padding": "16px"},
        children=[
            html.H2("PRS locus view", style={"fontSize": "26px"}),
            html.Div(
                style={"marginBottom": "8px", "fontWeight": "600", "fontSize": "24px"},
                children=trait_summary,
            ),
            html.Div(
                style={"marginBottom": "8px", "color": "#555"},
                children=[
                    html.Div(f"Loaded {len(tracks)} PRS file(s): " + ", ".join([t.label for t in tracks])),
                    html.Div("Pick an rsID and a window size to view variant locations and effect weights."),
                    html.Div(
                        children=[
                            "Viewing chromosome: ",
                            html.B(chrom_label),
                        ]
                    ),
                ],
            ),

            html.Hr(),

            html.Div(
                style={"display": "flex", "gap": "16px", "flexWrap": "wrap", "alignItems": "flex-end"},
                children=[
                    html.Div(
                        style={"minWidth": "420px", "flex": "3"},
                        children=[
                            html.Label("rsID"),
                            dcc.Dropdown(
                                id="rsid",
                                options=[{"label": r, "value": r} for r in all_rsids[:50]],
                                placeholder="Type to search (e.g. rs7903146)",
                                value=all_rsids[0] if all_rsids else None,
                                searchable=True,
                                clearable=True,
                            ),
                            html.Div(
                                style={"fontSize": "12px", "color": "#666", "marginTop": "4px"},
                                children="Tip: start typing 'rs' + digits. The dropdown will filter."
                            ),
                            html.Div(
                                style={"marginTop": "10px"},
                                children=[
                                    html.Label("Window size (kb, flanking region)"),
                                    dcc.Slider(
                                        id="window-kb",
                                        min=0, max=len(_WINDOW_KB_OPTIONS) - 1, step=1, value=2,
                                        marks={i: str(v) for i, v in enumerate(_WINDOW_KB_OPTIONS)},
                                        tooltip={"placement": "bottom", "always_visible": False},
                                    ),
                                ],
                            ),
                        ],
                    ),

                    html.Div(
                        style={"minWidth": "260px", "flex": "1"},
                        children=[
                            html.Label("View mode"),
                            dcc.RadioItems(
                                id="view-mode",
                                options=[
                                    {"label": "SNP", "value": "snp"},
                                    {"label": "Gene", "value": "gene"},
                                ],
                                value="snp",
                                labelStyle={"display": "block"},
                            ),
                            html.Div(
                                style={"marginTop": "8px"},
                                children=[
                                    html.Label("Gene"),
                                    dcc.Dropdown(id="gene-select", placeholder="Select a gene"),
                                ],
                            ),
                        ],
                    ),

                    html.Div(
                        style={"minWidth": "260px", "flex": "1"},
                        children=[
                            html.Label("Display options"),
                            dcc.Checklist(
                                id="opts",
                                options=[
                                    {"label": "Color ALL SNPs by effect_weight", "value": "color_all"},
                                    {"label": "Only show SNPs shared across all loaded PRSs (within window)", "value": "shared_only"},
                                    {"label": "Show track labels on y-axis", "value": "show_labels"},
                                ],
                                value=["show_labels"],
                            ),
                        ],
                    ),
                ],
            ),


            html.Div(
                style={
                    "marginTop": "10px",
                    "padding": "10px",
                    "border": "1px solid #eee",
                    "borderRadius": "8px",
                    "background": "#fafafa",
                },
                children=[
                    html.Label("PRS tracks to display"),
                    dcc.Checklist(
                        id="track-select",
                        options=[
                            {"label": t.label, "value": t.path}
                            for t in tracks
                        ],
                        value=[t.path for t in tracks],
                        labelStyle={"display": "block"},
                    ),
                    html.Div(
                        style={"fontSize": "12px", "color": "#666", "marginTop": "4px"},
                        children="Tip: uncheck tracks to focus the view; rsID search still spans all loaded tracks.",
                    ),
                ],
            ),

            html.Div(style={"marginTop": "12px", "color": "#444"}, id="status"),

            dcc.Graph(id="locus-plot", style={"marginTop": "8px"}),
            dcc.Graph(id="gene-plot", style={"marginTop": "8px"}),

            html.Div(
                style={
                    "marginTop": "16px",
                    "padding": "12px",
                    "border": "1px solid #eee",
                    "borderRadius": "8px",
                    "background": "#fafafa",
                },
                children=[
                    html.Label("Track details (#SNPS, #individual, ancestry, publication)", style={"fontWeight": "600"}),
                    html.Div(id="track-details", style={"marginTop": "8px"}),
                ],
            ),

            html.Hr(),

            html.Details(
                open=False,
                children=[
                    html.Summary("Show loaded track details"),
                    html.Ul([
                        html.Li(f"{t.label}" + (f" — {t.name}" if t.name else "") + f" — {os.path.basename(t.path)}")
                        for t in tracks
                    ])
                ],
            ),
        ],
    )

    # Dynamic dropdown options (prevents huge option payloads)
    @app.callback(
        Output("rsid", "options"),
        Input("rsid", "search_value"),
        State("rsid", "value"),
        prevent_initial_call=True,
    )
    def _update_rsid_options(search_value: Optional[str], current_value: Optional[str]):
        if not all_rsids:
            return []
        if not search_value:
            # Small default list, but keep current value if set
            opts = all_rsids[:50]
        else:
            s = search_value.strip().lower()
            # contains-match; cap for performance
            opts = [r for r in all_rsids if s in r.lower()][:200]
            if not opts:
                # allow exact typed rsID that may not exist (user may be pasting)
                opts = [search_value.strip()]
        if current_value and current_value not in opts:
            opts = [current_value] + opts
        return [{"label": r, "value": r} for r in opts]

    @app.callback(
        Output("gene-select", "options"),
        Output("gene-select", "value"),
        Input("view-mode", "value"),
        Input("rsid", "value"),
        Input("window-kb", "value"),
        Input("track-select", "value"),
        State("gene-select", "value"),
    )
    def _update_gene_options(
        view_mode: str,
        rsid_value: Optional[str],
        window_kb_value: int,
        selected_paths: List[str],
        current_gene: Optional[str],
    ):
        selected_paths = selected_paths or []
        selected_set = set(selected_paths)
        view_tracks = [t for t in tracks if t.path in selected_set]
        window_kb = _WINDOW_KB_OPTIONS[int(window_kb_value)] if window_kb_value is not None else 50

        options = []

        if view_mode == "gene":
            # UPDATED: Use the precomputed cache from the new module
            options = gene_mapper.get_top_genes_options()
        else:
            locus_chr, locus_pos = _resolve_locus(view_tracks, str(rsid_value) if rsid_value else "")
            if locus_chr is None or locus_pos is None:
                return [], None

            half_window = int(window_kb) * 1000
            xmin, xmax = locus_pos - half_window, locus_pos + half_window

            # UPDATED: Call without 'flank_bp' and assign directly (module now formats the dicts)
            options = gene_mapper.genes_with_snps_in_window(view_tracks, locus_chr, xmin, xmax)

        # Logic to persist the selected gene if it exists in the new options
        if current_gene and any(o["value"] == current_gene for o in options):
            return options, current_gene

        # Default to the first gene if options exist, otherwise None
        return options, options[0]["value"] if options else None

    @app.callback(
        Output("locus-plot", "figure"),
        Output("status", "children"),
        Output("gene-plot", "figure"),
        Output("track-details", "children"),
        Input("rsid", "value"),
        Input("window-kb", "value"),
        Input("opts", "value"),
        Input("track-select", "value"),
        Input("view-mode", "value"),
        Input("gene-select", "value"),
    )
    def _update_plot(
        rsid_value: Optional[str],
        window_kb_value: int,
        opts_value: List[str],
        selected_paths: List[str],
        view_mode: str,
        gene_value: Optional[str],
    ):
        opts_value = opts_value or []
        selected_paths = selected_paths or []
        selected_set = set(selected_paths)
        view_tracks = [t for t in tracks if t.path in selected_set]
        if not view_tracks:
            fig = go.Figure()
            fig.update_layout(
                height=260,
                margin=dict(l=30, r=20, t=50, b=40),
                title="No PRS tracks selected",
                template="plotly_white",
            )
            empty_gene = go.Figure()
            empty_gene.update_layout(
                height=180,
                margin=dict(l=30, r=20, t=40, b=30),
                title="Gene window view (GRCh37)",
                template="plotly_white",
            )
            return fig, "No PRS tracks selected (use the checklist above to enable at least one).", empty_gene, []

        window_kb = _WINDOW_KB_OPTIONS[int(window_kb_value)] if window_kb_value is not None else 50
        locus_override = None

        if view_mode == "gene" and not gene_value:
            fig = go.Figure()
            fig.update_layout(
                height=260,
                margin=dict(l=30, r=20, t=50, b=40),
                title="Select a gene to view the window",
                template="plotly_white",
            )
            gene_fig = go.Figure()
            gene_fig.update_layout(
                height=180,
                margin=dict(l=30, r=20, t=40, b=30),
                title="Gene window view (GRCh37)",
                template="plotly_white",
            )
            return fig, "Select a gene to view SNPs in the flanking window.", gene_fig, []

        if view_mode == "gene" and gene_value:
            gene_info = gene_mapper.fetch_gene_by_id(gene_value)
            if gene_info:
                gene_chr = str(gene_info.get("seq_region_name"))
                exons = gene_info.get("Exon") or []
                if isinstance(exons, list) and exons:
                    starts = [int(e.get("start", 0)) for e in exons if e.get("start") is not None]
                    ends = [int(e.get("end", 0)) for e in exons if e.get("end") is not None]
                    if starts and ends:
                        coding_start = min(starts)
                        coding_end = max(ends)
                    else:
                        coding_start = int(gene_info.get("start", 0))
                        coding_end = int(gene_info.get("end", 0))
                else:
                    coding_start = int(gene_info.get("start", 0))
                    coding_end = int(gene_info.get("end", 0))
                flank_bp = window_kb * 1000
                xmin, xmax = coding_start - flank_bp, coding_end + flank_bp
                locus_override = (gene_chr, int((coding_start + coding_end) / 2))
                gene_fig = gene_mapper.make_gene_track_figure(gene_chr, xmin, xmax, only_gene_ids={gene_value})
            else:
                gene_fig = go.Figure()
                gene_fig.update_layout(
                    height=180,
                    margin=dict(l=30, r=20, t=40, b=30),
                    title="Gene window view (GRCh37)",
                    template="plotly_white",
                )
        else:
            locus_chr, locus_pos = _resolve_locus(view_tracks, str(rsid_value) if rsid_value else "")
            if locus_chr is not None and locus_pos is not None:
                half_window = window_kb * 1000
                xmin, xmax = locus_pos - half_window, locus_pos + half_window
                only_gene_ids = {gene_value} if gene_value else None
                gene_fig = gene_mapper.make_gene_track_figure(locus_chr, xmin, xmax, only_gene_ids=only_gene_ids)
            else:
                gene_fig = go.Figure()
                gene_fig.update_layout(
                    height=180,
                    margin=dict(l=30, r=20, t=40, b=30),
                    title="Gene window view (GRCh37)",
                    template="plotly_white",
                )

        fig, status = make_locus_figure(
            tracks=view_tracks,
            rsid=str(rsid_value) if (view_mode == "snp" and rsid_value) else "",
            window_kb=int(window_kb),
            show_track_labels=("show_labels" in opts_value),
            color_all_snps=("color_all" in opts_value),
            shared_only=("shared_only" in opts_value),
            locus_override=locus_override,
        )

        max_snp_count = max([int(t.total_snps) for t in view_tracks] or [1])
        max_ind_count = max(
            [
                v
                for t in view_tracks
                for v in [t.gwas_count if t.gwas_count is not None else t.dev_count]
                if v is not None
            ]
            or [1]
        )
        snp_xmax = max(max_snp_count * 1.25, 1)
        ind_xmax = max(max_ind_count * 1.25, 1)

        details_children = []
        for tr in view_tracks:
            snp_count = int(tr.total_snps)
            gwas_val = tr.gwas_count if tr.gwas_count is not None else tr.dev_count
            snp_fig = _make_single_bar(snp_count, "#SNPS", xmax_override=snp_xmax, text_formatter=_format_count_short)
            gwas_fig = _make_single_bar(gwas_val, "#individual", xmax_override=ind_xmax)
            ancestry_fig = _make_ancestry_bar(tr)
            ancestry_graph = (
                dcc.Graph(figure=ancestry_fig, config={"displayModeBar": False})
                if ancestry_fig is not None
                else html.Div(style={"color": "#666", "paddingTop": "8px"}, children="Ancestry not available")
            )
            label_link = (
                html.A(
                    tr.label,
                    href=tr.pgs_score_url,
                    target="_blank",
                    rel="noopener noreferrer",
                    style={"color": "#0b5394", "textDecoration": "none"},
                )
                if tr.pgs_score_url
                else html.Span(tr.label)
            )
            pub_child = None
            if tr.publication_citation:
                pub_child = html.A(
                    tr.publication_citation,
                    href=tr.publication_url or tr.pgs_score_url,
                    target="_blank",
                    rel="noopener noreferrer",
                    style={"color": "#0b5394", "textDecoration": "none"},
                )
            else:
                pub_child = html.Span("Publication not available", style={"color": "#666"})
            details_children.append(
                html.Div(
                    style={
                        "display": "flex",
                        "gap": "12px",
                        "flexWrap": "wrap",
                        "alignItems": "flex-start",
                        "marginTop": "8px",
                        "fontSize": "12px",
                    },
                    children=[
                        html.Div(
                            style={"minWidth": "160px", "fontWeight": "700", "fontSize": "14.4px"},
                            children=label_link,
                        ),
                        html.Div(style={"minWidth": "216px"}, children=dcc.Graph(figure=snp_fig, config={"displayModeBar": False})),
                        html.Div(style={"minWidth": "240px"}, children=dcc.Graph(figure=gwas_fig, config={"displayModeBar": False})),
                        html.Div(style={"minWidth": "264px"}, children=ancestry_graph),
                        html.Div(
                            style={"minWidth": "240px", "fontSize": "14.4px"},
                            children=pub_child,
                        ),
                    ],
                )
            )
        return fig, status, gene_fig, details_children
    
    return app


# ----------------------------
# File selection & launch
# ----------------------------

def choose_files_dialog() -> List[str]:
    """
    Native OS file picker (works when running locally with a GUI).
    """
    try:
        import tkinter as tk
        from tkinter import filedialog
    except Exception as e:
        raise RuntimeError(
            "tkinter is not available (common on minimal servers). "
            "Please pass file paths on the command line instead."
        ) from e

    root = tk.Tk()
    root.withdraw()
    paths = filedialog.askopenfilenames(
        title="Select PRS/PGS scoring files",
        filetypes=[
            ("PGS scoring files", "*.txt *.tsv *.tab *.gz"),
            ("All files", "*.*"),
        ],
    )
    return list(paths)


def load_tracks(paths: List[str], chrom_filter: Optional[str] = None) -> List[ScoreTrack]:
    target_chr = _normalize_chr(chrom_filter)
    tracks: List[ScoreTrack] = []
    for p in paths:
        meta, df = read_pgs_scoring_file(p)
        total_snps = df.shape[0]
        if target_chr:
            df = df[df["hm_chr"].str.lower() == target_chr]
            if df.empty:
                print(f"[warn] {os.path.basename(p)} has no variants on chr{target_chr}; keeping metadata only.")
        label, name = _infer_track_label(meta, p)
        trait_reported = _infer_trait(meta)
        pgs_id = info_parser.parse_pgs_id(meta, p) or label
        if _PGS_ID_RE.match(label) is None and pgs_id:
            label = pgs_id
        gwas_ancestry, gwas_count, dev_count = _extract_ancestry_from_meta(meta, pgs_id=pgs_id)
        publication_citation, publication_url = _fetch_publication_info(pgs_id)
        pgs_score_url = _build_pgs_score_url(pgs_id)
        tracks.append(
            ScoreTrack(
                label=label,
                name=name,
                path=p,
                df=df,
                total_snps=total_snps,
                meta=meta,
                trait_reported=trait_reported,
                pgs_id=pgs_id,
                gwas_ancestry=gwas_ancestry,
                gwas_count=gwas_count,
                dev_count=dev_count,
                publication_citation=publication_citation,
                publication_url=publication_url,
                pgs_score_url=pgs_score_url,
            )
        )
    if not tracks:
        raise ValueError("No input files provided.")
    # Order tracks by SNP count (descending) so denser scores appear first in the UI
    tracks = sorted(tracks, key=lambda t: t.total_snps, reverse=True)
    if len(tracks) > 10:
        print(f"[warn] You provided {len(tracks)} files; showing the first 10.")
        tracks = tracks[:10]
    return tracks


def launch_prs_viewer(
    file_paths: Optional[List[str]] = None,
    host: str = "127.0.0.1",
    port: int = 8050,
    debug: bool = True,
    chrom_filter: Optional[str] = "1",
):
    chrom_filter = _normalize_chr(chrom_filter)
    if file_paths is None or len(file_paths) == 0:
        file_paths = choose_files_dialog()
    tracks = load_tracks(file_paths, chrom_filter=chrom_filter)

    chrom_label = f"chr{chrom_filter}" if chrom_filter else "all chromosomes"
    flank_bp = 25000

    print("\n" + "=" * 60)
    print(f"  PRS LOCUS VIEWER | Chromosome: {chrom_label}")
    print("=" * 60 + "\n")

    print(f"[1/3] Loading data for {chrom_label}...")

    unique_snp_count = 0
    try:
        all_snps = pd.concat([t.df[["hm_chr", "hm_pos"]] for t in tracks], ignore_index=True)
        unique_snp_count = int(all_snps.drop_duplicates().shape[0])
    except Exception:
        unique_snp_count = sum(len(t.df) for t in tracks)

    print(f"[2/3] Calculating gene coverage (Window: ±{flank_bp // 1000}kb)...")
    if unique_snp_count:
        print(f"      Targeting {unique_snp_count:,} unique SNPs")

    gene_mapper.precompute_top_genes(tracks, flank_bp=flank_bp)

    print(f"[3/3] Calculation complete.\n")

    print("-" * 60)
    print(f" [READY] Viewer active at: http://{host}:{port}/")
    print("         (Press CTRL+C to stop)")
    print("-" * 60 + "\n")

    app = create_app(tracks, chrom_filter=chrom_filter)

    def _announce_chrom_filter():
        # Redundant print removed since we print above, but keeping hook if needed
        pass

    if hasattr(app.server, "before_serving"):
        app.server.before_serving(_announce_chrom_filter)
    
    app.run(host=host, port=port, debug=debug, use_reloader=False, dev_tools_hot_reload=False)

def main():
    ap = argparse.ArgumentParser(description="Dash PRS locus viewer")
    ap.add_argument("files", nargs="*", help="2–10 PRS/PGS scoring .txt files (PGS Catalog format supported)")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", default=8050, type=int)
    ap.add_argument("--no-debug", action="store_true")
    ap.add_argument("--chr", dest="chrom", default="1", help="Chromosome to load (default chr1; use 'all' to keep all). Example: --chr chr5")
    args = ap.parse_args()

    chrom_filter = _normalize_chr(args.chrom)
    file_paths = args.files if args.files else None
    launch_prs_viewer(file_paths=file_paths, host=args.host, port=args.port, debug=not args.no_debug, chrom_filter=chrom_filter)


if __name__ == "__main__":
    main()
