# 02_build_sparse_matrices.py
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import scipy.sparse as sp


@dataclass
class Config:
    processed_dir: Path = Path("./data/processed")
    tensor_dir: Path = Path("./data/tensors")

    master_annotated_parquet: Path = Path("./data/processed/pgs_variants_master_annotated.parquet")

    # Base numeric feature
    weight_col: str = "weight_scaled"

    # Annotation columns (we will canonicalize from *_y / *_x)
    use_n_genes: bool = True
    n_genes_col: str = "n_genes"
    n_genes_scale: float = 1.0

    # Canonical names we will create in-memory
    gene_region_out: str = "gene_region"
    mutation_type_out: str = "mutation_type"

    # Optional top-k by abs(weight_scaled) per PGS
    use_topk: bool = True
    topk: int = 300

    # Outputs
    out_csr_npz: Path = Path("./data/tensors/X_pgs_feat_csr.npz")
    out_pgs_ids_json: Path = Path("./data/tensors/pgs_ids.json")
    out_locus_index_json: Path = Path("./data/tensors/locus_index.json")
    out_feature_schema_json: Path = Path("./data/tensors/feature_schema.json")


def ensure_dirs(cfg: Config):
    cfg.tensor_dir.mkdir(parents=True, exist_ok=True)


def safe_str(x) -> str:
    if x is None:
        return ""
    s = str(x)
    if s.lower() in ["nan", "none"]:
        return ""
    return s


def pick_first_existing(df: pd.DataFrame, candidates: List[str], default):
    for c in candidates:
        if c in df.columns:
            return df[c]
    return default


def canonicalize_annotation_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Your file has gene_region_x/gene_region_y etc. We choose *_y (VEP) if present,
    else fallback to *_x, else default.
    """
    df = df.copy()

    gene_region = pick_first_existing(df, ["gene_region_y", "gene_region_x", "gene_region"], pd.Series(["other"] * len(df)))
    mutation_type = pick_first_existing(df, ["mutation_type_y", "mutation_type_x", "mutation_type"], pd.Series(["other"] * len(df)))

    df["gene_region"] = gene_region.map(safe_str).replace("", "other")
    df["mutation_type"] = mutation_type.map(safe_str).replace("", "other")

    # n_genes already exists in your columns list; normalize anyway
    if "n_genes" in df.columns:
        df["n_genes"] = pd.to_numeric(df["n_genes"], errors="coerce").fillna(0).astype(float)
    else:
        df["n_genes"] = 0.0

    return df


def main():
    cfg = Config()
    ensure_dirs(cfg)

    print("Config:", json.dumps({k: str(v) for k, v in asdict(cfg).items()}, indent=2))

    master = pd.read_parquet(cfg.master_annotated_parquet).copy()

    # Required columns
    for req in ["pgs_id", "locus_id", cfg.weight_col]:
        if req not in master.columns:
            raise ValueError(f"Missing required column '{req}' in {cfg.master_annotated_parquet}")

    master["pgs_id"] = master["pgs_id"].astype(str)
    master["locus_id"] = master["locus_id"].astype(str)

    # Canonicalize annotations (handles *_x/*_y)
    master = canonicalize_annotation_columns(master)

    # Determine categories globally (shared schema)
    gene_region_cats = sorted(master["gene_region"].unique().tolist())
    mutation_type_cats = sorted(master["mutation_type"].unique().tolist())

    # Build feature index within-locus
    feat_names: List[str] = []
    feat_index: Dict[str, int] = {}

    def add_feat(name: str):
        feat_index[name] = len(feat_names)
        feat_names.append(name)

    add_feat("weight_scaled")
    if cfg.use_n_genes:
        add_feat("n_genes")

    for gr in gene_region_cats:
        add_feat(f"gene_region={gr}")

    for mt in mutation_type_cats:
        add_feat(f"mutation_type={mt}")

    F = len(feat_names)
    print(f"Feature channels per locus: F={F}")
    print(f"gene_region categories: {gene_region_cats}")
    print(f"mutation_type categories: {mutation_type_cats}")

    # Index PGS and loci
    pgs_ids = sorted(master["pgs_id"].unique().tolist())
    pgs_to_row = {p: i for i, p in enumerate(pgs_ids)}

    loci = sorted(master["locus_id"].unique().tolist())
    locus_index = {l: i for i, l in enumerate(loci)}
    B, V = len(pgs_ids), len(loci)

    rows: List[int] = []
    cols: List[int] = []
    data: List[float] = []

    for pgs_id, df in master.groupby("pgs_id", sort=False):
        r = pgs_to_row[pgs_id]
        df = df.dropna(subset=["locus_id", cfg.weight_col]).copy()

        if cfg.use_topk and len(df) > cfg.topk:
            w = df[cfg.weight_col].astype(float).to_numpy()
            keep = np.argpartition(np.abs(w), -cfg.topk)[-cfg.topk:]
            df = df.iloc[keep]

        for _, row in df.iterrows():
            lid = str(row["locus_id"])
            j = locus_index.get(lid)
            if j is None:
                continue

            base = j * F

            # weight_scaled
            w = float(row[cfg.weight_col])
            if w != 0.0:
                rows.append(r)
                cols.append(base + feat_index["weight_scaled"])
                data.append(w)

            # n_genes (optional numeric)
            if cfg.use_n_genes:
                ng = float(row.get("n_genes", 0.0)) * float(cfg.n_genes_scale)
                if ng != 0.0:
                    rows.append(r)
                    cols.append(base + feat_index["n_genes"])
                    data.append(ng)

            # gene_region onehot
            gr = safe_str(row.get("gene_region", "other")) or "other"
            gr_key = f"gene_region={gr}"
            if gr_key in feat_index:
                rows.append(r)
                cols.append(base + feat_index[gr_key])
                data.append(1.0)

            # mutation_type onehot
            mt = safe_str(row.get("mutation_type", "other")) or "other"
            mt_key = f"mutation_type={mt}"
            if mt_key in feat_index:
                rows.append(r)
                cols.append(base + feat_index[mt_key])
                data.append(1.0)

    X_csr = sp.csr_matrix(
        (np.array(data, dtype=np.float32), (np.array(rows), np.array(cols))),
        shape=(B, V * F),
    )
    sp.save_npz(cfg.out_csr_npz, X_csr)

    cfg.out_pgs_ids_json.write_text(json.dumps(pgs_ids, indent=2))
    cfg.out_locus_index_json.write_text(json.dumps(locus_index, indent=2))

    feature_schema = {
        "F": F,
        "feat_names": feat_names,
        "feat_index": feat_index,
        "gene_region_categories": gene_region_cats,
        "mutation_type_categories": mutation_type_cats,
        "layout": "columns are (locus_index * F + feature_index)",
        "numeric_features": ["weight_scaled"] + (["n_genes"] if cfg.use_n_genes else []),
        "note": "gene_region/mutation_type canonicalized from *_y first, then *_x.",
    }
    cfg.out_feature_schema_json.write_text(json.dumps(feature_schema, indent=2))

    print("Saved multi-feature CSR:", cfg.out_csr_npz.resolve())
    print(" - shape:", X_csr.shape, "nnz:", X_csr.nnz)
    print("Saved feature_schema:", cfg.out_feature_schema_json.resolve())


if __name__ == "__main__":
    main()
