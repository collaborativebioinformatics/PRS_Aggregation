# 01b_annotate_variants_ensembl_vep.py
from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Any, List, Tuple

import pandas as pd
import requests


@dataclass
class Config:
    # Input from Script 1
    master_parquet: Path = Path("./data/processed/pgs_variants_master.parquet")

    # Outputs
    out_master_annotated: Path = Path("./data/processed/pgs_variants_master_annotated.parquet")
    out_locus_table: Path = Path("./data/processed/locus_annotations.parquet")

    # Cache so you don't re-hit VEP
    cache_json: Path = Path("./data/cache/vep_locus_cache.json")

    # VEP GRCh37 server + endpoint
    server: str = "https://grch37.rest.ensembl.org"
    endpoint: str = "/vep/homo_sapiens/region"

    # Request tuning
    batch_size: int = 200
    sleep_s: float = 0.2
    timeout_s: int = 60
    max_retries: int = 5

    # Add query params to improve output content
    # (These are optional flags supported by VEP REST.)
    query_params: str = "?canonical=1&numbers=1&variant_class=1&hgvs=1&vcf_string=1&minimal=1"


def ensure_dirs(cfg: Config):
    cfg.cache_json.parent.mkdir(parents=True, exist_ok=True)
    cfg.out_master_annotated.parent.mkdir(parents=True, exist_ok=True)


def load_cache(path: Path) -> Dict[str, Any]:
    if path.exists():
        return json.loads(path.read_text())
    return {}


def save_cache(path: Path, cache: Dict[str, Any]):
    path.write_text(json.dumps(cache, indent=2))


def norm_chr(chr_str: str) -> str:
    s = str(chr_str).strip()
    s = s.replace("chr", "").replace("CHR", "")
    return s


def build_vep_variant_string(chr_: str, pos: str, rsid: str, ref: str, alt: str) -> str:
    """
    Ensembl VEP region POST expects entries like:
      "21  26960070  rs116645811 G A . . ."
    (chrom, position, identifier, ref, alt, then VCF-style placeholder fields) :contentReference[oaicite:1]{index=1}
    """
    chr_ = norm_chr(chr_)
    pos = str(pos).strip()
    rsid = str(rsid).strip() if rsid else "."
    ref = str(ref).strip().upper() if ref else "."
    alt = str(alt).strip().upper() if alt else "."

    # Keep it robust: if we’re missing alleles, we can still attempt with "."
    return f"{chr_} {pos} {rsid} {ref} {alt} . . ."


def vep_post(cfg: Config, variant_strings: List[str]) -> List[Dict[str, Any]]:
    url = cfg.server + cfg.endpoint + cfg.query_params
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    payload = {"variants": variant_strings}

    last_err = None
    for attempt in range(1, cfg.max_retries + 1):
        try:
            r = requests.post(url, headers=headers, json=payload, timeout=cfg.timeout_s)
            if r.status_code == 429:
                # rate-limited
                time.sleep(cfg.sleep_s * attempt)
                continue
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last_err = e
            time.sleep(cfg.sleep_s * attempt)

    raise RuntimeError(f"VEP request failed after {cfg.max_retries} retries. Last error: {last_err}")


def map_gene_region(most_severe: str) -> str:
    s = (most_severe or "").lower()
    if "intergenic" in s:
        return "intergenic"
    if "upstream_gene_variant" in s or "downstream_gene_variant" in s:
        return "upstream_downstream"
    if "intron" in s:
        return "intronic"
    if "utr" in s:
        return "utr"
    if "regulatory" in s:
        return "regulatory"
    # exonic-ish consequences
    exonic_keys = [
        "missense", "synonymous", "stop_gained", "stop_lost",
        "start_lost", "frameshift", "inframe", "protein_altering",
        "splice", "coding_sequence_variant"
    ]
    if any(k in s for k in exonic_keys):
        return "exonic"
    return "other"


def map_mutation_type(most_severe: str) -> str:
    s = (most_severe or "").lower()
    if "missense" in s:
        return "missense"
    if "synonymous" in s:
        return "synonymous"
    if "stop_gained" in s or "stop_lost" in s or "start_lost" in s:
        return "stop_start"
    if "frameshift" in s:
        return "frameshift"
    if "inframe" in s:
        return "inframe"
    if "splice" in s:
        return "splice"
    if "utr" in s:
        return "utr"
    if "intron" in s:
        return "intronic"
    if "intergenic" in s:
        return "intergenic"
    if "regulatory" in s:
        return "regulatory"
    if "upstream_gene_variant" in s or "downstream_gene_variant" in s:
        return "upstream_downstream"
    return "other"


def extract_locus_annotation(vep_rec: Dict[str, Any]) -> Dict[str, Any]:
    """
    Pull a compact locus-level annotation from VEP JSON.
    We use:
      - most_severe_consequence (top-level)
      - transcript_consequences gene_symbol
    """
    most_severe = vep_rec.get("most_severe_consequence", "") or ""

    genes = set()
    tcs = vep_rec.get("transcript_consequences", []) or []
    for tc in tcs:
        gs = tc.get("gene_symbol")
        if gs:
            genes.add(str(gs))

    genes_sorted = sorted(genes)
    return {
        "genes": ",".join(genes_sorted),
        "n_genes": int(len(genes_sorted)),
        "most_severe_consequence": most_severe,
        "gene_region": map_gene_region(most_severe),
        "mutation_type": map_mutation_type(most_severe),
    }


def chunked(lst: List[Any], n: int):
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def main():
    cfg = Config()
    ensure_dirs(cfg)

    master = pd.read_parquet(cfg.master_parquet)

    # We expect harmonized coords available from Script 1
    # locus_id is like "1:20688352" and coord_source "harmonized"
    loci = master[["locus_id", "hm_chr", "hm_pos", "variant_id", "other_allele", "effect_allele"]].copy()
    loci = loci.drop_duplicates(subset=["locus_id"]).reset_index(drop=True)

    # Build variant strings for VEP
    loci["hm_chr"] = loci["hm_chr"].astype(str)
    loci["hm_pos"] = loci["hm_pos"].astype(str)
    loci["variant_id"] = loci["variant_id"].astype(str)
    loci["other_allele"] = loci["other_allele"].astype(str)
    loci["effect_allele"] = loci["effect_allele"].astype(str)

    cache = load_cache(cfg.cache_json)

    # Determine what needs annotation
    to_query = []
    for _, row in loci.iterrows():
        lid = row["locus_id"]
        if lid in cache:
            continue
        # If hm coords missing, skip (should be rare if coord_source is harmonized)
        if row["hm_chr"] in ["", "nan", "None"] or row["hm_pos"] in ["", "nan", "None"]:
            cache[lid] = {
                "genes": "",
                "n_genes": 0,
                "most_severe_consequence": "",
                "gene_region": "other",
                "mutation_type": "other",
                "status": "missing_harmonized_coords",
            }
            continue
        to_query.append(lid)

    print(f"Unique loci: {len(loci)}")
    print(f"Cached loci: {len(cache)}")
    print(f"To query: {len(to_query)}")

    # Index rows by locus_id for quick access
    loci_by_id = loci.set_index("locus_id").to_dict(orient="index")

    # Query in batches
    for batch_ids in chunked(to_query, cfg.batch_size):
        variant_strings = []
        for lid in batch_ids:
            r = loci_by_id[lid]
            variant_strings.append(
                build_vep_variant_string(
                    chr_=r["hm_chr"],
                    pos=r["hm_pos"],
                    rsid=r.get("variant_id", "."),
                    ref=r.get("other_allele", "."),
                    alt=r.get("effect_allele", "."),
                )
            )

        results = vep_post(cfg, variant_strings)

        # Results come back in the same order as the input variants (typical for VEP REST)
        for lid, rec in zip(batch_ids, results):
            ann = extract_locus_annotation(rec)
            ann["status"] = "ok"
            cache[lid] = ann

        save_cache(cfg.cache_json, cache)
        time.sleep(cfg.sleep_s)

        print(f"Annotated {len(cache)} / {len(loci)} loci...")

    # Build locus annotation table
    ann_rows = []
    for lid, ann in cache.items():
        ann_rows.append({"locus_id": lid, **ann})
    ann_df = pd.DataFrame(ann_rows).drop_duplicates(subset=["locus_id"])

    # Merge back into master
    master2 = master.merge(ann_df, on="locus_id", how="left")

    # Fill any missing annotations
    for c, default in [
        ("genes", ""),
        ("n_genes", 0),
        ("most_severe_consequence", ""),
        ("gene_region", "other"),
        ("mutation_type", "other"),
        ("status", "missing"),
    ]:
        if c in master2.columns:
            master2[c] = master2[c].fillna(default)

    master2.to_parquet(cfg.out_master_annotated, index=False)
    ann_df.to_parquet(cfg.out_locus_table, index=False)

    print("Saved annotated master:", cfg.out_master_annotated.resolve())
    print("Saved locus table:", cfg.out_locus_table.resolve())
    print("Saved cache:", cfg.cache_json.resolve())


if __name__ == "__main__":
    main()