from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Tuple, List, Any, Optional

import numpy as np
import pandas as pd


@dataclass
class Config:
    raw_dir: Path = Path("./data/raw")
    processed_dir: Path = Path("./data/processed")
    cache_dir: Path = Path("./data/cache")

    master_parquet: Path = Path("./data/processed/pgs_variants_master.parquet")
    meta_parquet: Path = Path("./data/processed/pgs_metadata.parquet")
    report_json: Path = Path("./data/processed/validation_report.json")

    feature_col: str = "weight_scaled"

    # Validation
    drop_missing_locus: bool = True
    drop_missing_weight: bool = True

    # Optional gene annotation stage (kept as placeholder columns by default)
    annotate_genes: bool = False
    annotation_cache_json: Path = Path("./data/cache/variant_annotation_cache.json")


def ensure_dirs(cfg: Config):
    cfg.processed_dir.mkdir(parents=True, exist_ok=True)
    cfg.cache_dir.mkdir(parents=True, exist_ok=True)


def safe_float_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def try_parse_jsonish(value: str) -> Any:
    """
    Some header values look like JSON: {"True": 282, "False": 0}
    Others are plain strings. Try JSON parsing; fallback to string.
    """
    v = value.strip()
    if v.startswith("{") or v.startswith("["):
        try:
            return json.loads(v)
        except Exception:
            return v
    return v


def parse_pgs_file(file_path: Path) -> Tuple[Dict[str, Any], pd.DataFrame]:
    """
    Parse a PGS Catalog scoring file (format_version ~2.0).

    Header:
      - ignore "###" banner and "##" section markers
      - parse "#key=value" lines into metadata

    Data:
      - first non-# line is the tab header
      - remaining lines are rows (may be empty)
    """
    metadata: Dict[str, Any] = {"source_file": file_path.name}
    header_cols: Optional[List[str]] = None
    rows: List[List[str]] = []

    with open(file_path, "r") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line.strip():
                continue

            if line.startswith("###"):
                continue
            if line.startswith("##"):
                # section marker
                continue

            if line.startswith("#"):
                # metadata
                if "=" in line:
                    k, v = line[1:].split("=", 1)
                    metadata[k.strip()] = try_parse_jsonish(v)
                else:
                    # comment-only header line
                    pass
                continue

            # data section: first non-# line = header
            if header_cols is None:
                header_cols = [c.strip() for c in line.split("\t")]
            else:
                rows.append(line.split("\t"))

    if header_cols is None:
        # no table at all
        return metadata, pd.DataFrame()

    df = pd.DataFrame(rows, columns=header_cols)

    # Normalize common naming
    if "rsID" in df.columns and "variant_id" not in df.columns:
        df = df.rename(columns={"rsID": "variant_id"})
    if "rsid" in df.columns and "variant_id" not in df.columns:
        df = df.rename(columns={"rsid": "variant_id"})
    if "variant_id" not in df.columns:
        df["variant_id"] = None

    return metadata, df


def infer_ancestry(metadata: Dict[str, Any], file_path: Path) -> str:
    """
    Infer ancestry from pgs_name token patterns or filename words.
    Examples:
      - pgs_name includes "-SAS-" => South Asian
      - pgs_name includes "T2Dafr" => African
    """
    name = str(metadata.get("pgs_name", "")).lower()
    fname = file_path.stem.lower()

    # token map
    token_map = {
        "sas": "South Asian",
        "eas": "East Asian",
        "eur": "European",
        "afr": "African",
        "his": "Hispanic/Latino",
        "amr": "Hispanic/Latino",
        "lat": "Hispanic/Latino",
    }

    # look for -sas- / _sas_ / sas etc (carefully)
    for tok, label in token_map.items():
        if f"-{tok}-" in name or f"_{tok}_" in name or name.endswith(tok) or f"{tok}-" in name:
            return label
        if f"{tok}" in name and ("t2d" in name or "grs" in name or "prs" in name):
            # fallback heuristic for things like "GRS582_T2Dafr"
            if tok in ["afr", "his", "sas", "eas", "eur"]:
                return label

    # filename heuristics
    if "southasian" in fname or "south_asian" in fname:
        return "South Asian"
    if "eastasian" in fname or "east_asian" in fname:
        return "East Asian"
    if "european" in fname or "eur" in fname:
        return "European"
    if "african" in fname or "afr" in fname:
        return "African"
    if "hispanic" in fname or "latino" in fname or "his" in fname:
        return "Hispanic/Latino"

    return ""


def derive_weight(df: pd.DataFrame) -> pd.DataFrame:
    """
    Robust weight extraction irrespective of metadata weight_type.
    Preference:
      effect_weight > beta > OR(log)
    """
    df = df.copy()

    for c in ["effect_weight", "beta", "OR"]:
        if c in df.columns:
            df[c] = safe_float_series(df[c])

    if "effect_weight" in df.columns and df["effect_weight"].notna().any():
        df["weight_raw"] = df["effect_weight"]
        df["weight_source"] = "effect_weight"
    elif "beta" in df.columns and df["beta"].notna().any():
        df["weight_raw"] = df["beta"]
        df["weight_source"] = "beta"
    elif "OR" in df.columns and df["OR"].notna().any():
        orv = df["OR"].replace(0, np.nan)
        df["weight_raw"] = np.log(orv).replace([np.inf, -np.inf], np.nan).fillna(0.0)
        df["weight_source"] = "log(OR)"
    else:
        df["weight_raw"] = np.nan
        df["weight_source"] = "missing"

    df["weight_raw"] = safe_float_series(df["weight_raw"])
    return df


def robust_scale_within_pgs(weights: np.ndarray) -> np.ndarray:
    """
    Median/MAD robust z-score + tanh squash to [-1, 1]
    """
    w = weights.astype(np.float32)
    med = np.nanmedian(w)
    mad = np.nanmedian(np.abs(w - med))
    if mad > 1e-8:
        z = (w - med) / (1.4826 * mad)
    else:
        std = np.nanstd(w)
        z = (w - np.nanmean(w)) / (std + 1e-8)
    return np.tanh(z).astype(np.float32)


def build_locus_id(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build locus_id:
      1) hm_chr + hm_pos (+ alleles if available)
      2) chr_name + chr_position (+ other/effect_allele)
      3) variant_id (rsID) fallback
    """
    df = df.copy()

    def norm_chr(x) -> str:
        if pd.isna(x):
            return ""
        s = str(x).strip()
        s = s.replace("chr", "").replace("CHR", "")
        return s

    # Harmonized coordinates are your best cross-build key.
    if "hm_chr" in df.columns and "hm_pos" in df.columns:
        c = df["hm_chr"].map(norm_chr)
        p = df["hm_pos"].astype(str)
        # PGS scoring files usually do not include hm_ref/hm_alt, so we keep chr:pos.
        df["locus_id"] = c + ":" + p
        df["coord_source"] = "harmonized"
        return df

    # Raw coordinates
    if "chr_name" in df.columns and "chr_position" in df.columns:
        c = df["chr_name"].map(norm_chr)
        p = df["chr_position"].astype(str)
        if "other_allele" in df.columns and "effect_allele" in df.columns:
            ref = df["other_allele"].astype(str)
            alt = df["effect_allele"].astype(str)
            df["locus_id"] = c + ":" + p + ":" + ref + ":" + alt
        else:
            df["locus_id"] = c + ":" + p
        df["coord_source"] = "raw"
        return df

    # Fallback
    df["locus_id"] = df["variant_id"].astype(str)
    df["coord_source"] = "rsid"
    return df


def add_annotation_placeholders(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keeps pipeline stable even before you implement real annotation.
    """
    df = df.copy()
    for c in ["genes", "gene_region", "consequence", "mutation_type"]:
        if c not in df.columns:
            df[c] = ""
    return df


def build_metadata_record(meta: Dict[str, Any], df: pd.DataFrame, ancestry: str) -> Dict[str, Any]:
    # Normalize key names you’ll actually use later
    pgs_id = str(meta.get("pgs_id", ""))
    out = {
        "pgs_id": pgs_id,
        "pgs_name": str(meta.get("pgs_name", "")),
        "trait_reported": str(meta.get("trait_reported", "")),
        "trait_mapped": str(meta.get("trait_mapped", "")),
        "trait_efo": str(meta.get("trait_efo", "")),
        "genome_build": str(meta.get("genome_build", "")),
        "weight_type": str(meta.get("weight_type", "")),
        "variants_number_header": meta.get("variants_number", None),
        "HmPOS_build": str(meta.get("HmPOS_build", "")),
        "HmPOS_date": str(meta.get("HmPOS_date", "")),
        "source_file": str(meta.get("source_file", "")),
        "ancestry": ancestry,
        "n_rows_table": int(len(df)),
        "n_rows_nonnull_weight": int(df["weight_raw"].notna().sum()) if "weight_raw" in df.columns else 0,
        "n_rows_nonnull_hm": int(df["hm_chr"].notna().sum()) if "hm_chr" in df.columns else 0,
    }
    return out


def main():
    cfg = Config()
    ensure_dirs(cfg)

    raw_files = sorted(cfg.raw_dir.glob("*.txt"))
    if not raw_files:
        raise FileNotFoundError(f"No .txt files found in {cfg.raw_dir.resolve()}")

    all_meta: List[Dict[str, Any]] = []
    all_tables: List[pd.DataFrame] = []

    report = {
        "files_total": len(raw_files),
        "files_parsed": 0,
        "files_skipped_no_rows": 0,
        "files_errors": 0,
        "errors": [],
    }

    for fp in raw_files:
        try:
            meta, df = parse_pgs_file(fp)

            pgs_id = str(meta.get("pgs_id", fp.stem))
            meta["pgs_id"] = pgs_id

            # If table has 0 rows, record metadata and skip variants
            if df.empty:
                ancestry = infer_ancestry(meta, fp)
                all_meta.append(build_metadata_record(meta, df, ancestry))
                report["files_skipped_no_rows"] += 1
                continue

            ancestry = infer_ancestry(meta, fp)

            # Derive weight and locus
            df = derive_weight(df)
            df = build_locus_id(df)

            # Drop missing essentials
            if cfg.drop_missing_locus:
                df = df[df["locus_id"].notna() & (df["locus_id"].astype(str).str.len() > 0)]
            if cfg.drop_missing_weight:
                df = df[df["weight_raw"].notna()]

            # Scale within PGS
            if len(df) > 0:
                df[cfg.feature_col] = robust_scale_within_pgs(df["weight_raw"].to_numpy())
            else:
                df[cfg.feature_col] = np.array([], dtype=np.float32)

            # Keep enrichment-friendly master schema
            df_master = pd.DataFrame({
                "pgs_id": pgs_id,
                "pgs_name": str(meta.get("pgs_name", "")),
                "trait_mapped": str(meta.get("trait_mapped", "")),
                "genome_build": str(meta.get("genome_build", "")),
                "weight_type": str(meta.get("weight_type", "")),
                "ancestry": ancestry,

                "locus_id": df["locus_id"].astype(str),
                "coord_source": df["coord_source"].astype(str),

                "variant_id": df["variant_id"].astype(str) if "variant_id" in df.columns else "",
                "chr_name": df["chr_name"].astype(str) if "chr_name" in df.columns else "",
                "chr_position": df["chr_position"].astype(str) if "chr_position" in df.columns else "",
                "effect_allele": df["effect_allele"].astype(str) if "effect_allele" in df.columns else "",
                "other_allele": df["other_allele"].astype(str) if "other_allele" in df.columns else "",

                "hm_chr": df["hm_chr"].astype(str) if "hm_chr" in df.columns else "",
                "hm_pos": df["hm_pos"].astype(str) if "hm_pos" in df.columns else "",

                "weight_raw": df["weight_raw"].astype(float),
                "weight_source": df["weight_source"].astype(str),
                cfg.feature_col: df[cfg.feature_col].astype(float),
            })

            df_master = add_annotation_placeholders(df_master)

            all_tables.append(df_master)
            all_meta.append(build_metadata_record(meta, df, ancestry))

            report["files_parsed"] += 1

        except Exception as e:
            report["files_errors"] += 1
            report["errors"].append({"file": fp.name, "error": str(e)})

    # Concatenate all variant rows
    master_df = pd.concat(all_tables, ignore_index=True) if all_tables else pd.DataFrame()
    meta_df = pd.DataFrame(all_meta)

    # Final validation report stats
    if not master_df.empty:
        report["rows_master"] = int(len(master_df))
        report["n_pgs"] = int(master_df["pgs_id"].nunique())
        report["n_unique_loci"] = int(master_df["locus_id"].nunique())
        report["coord_source_counts"] = master_df["coord_source"].value_counts(dropna=False).to_dict()
        report["ancestry_counts"] = master_df["ancestry"].value_counts(dropna=False).to_dict()
    else:
        report["rows_master"] = 0
        report["n_pgs"] = int(meta_df["pgs_id"].nunique()) if not meta_df.empty else 0
        report["n_unique_loci"] = 0

    # Save
    master_df.to_parquet(cfg.master_parquet, index=False)
    meta_df.to_parquet(cfg.meta_parquet, index=False)
    cfg.report_json.write_text(json.dumps(report, indent=2))

    print("Saved:", cfg.master_parquet.resolve())
    print("Saved:", cfg.meta_parquet.resolve())
    print("Saved:", cfg.report_json.resolve())
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
