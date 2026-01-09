"""
Local metadata parser helpers for PRS/PGS scoring files.
"""
from __future__ import annotations

import json
import os
import re
from typing import Dict, Optional, Tuple

import requests

_PGS_CACHE: Dict[str, Dict[str, object]] = {}
_PGS_ID_RE = re.compile(r"(PGS\d{6})", re.IGNORECASE)

def get_meta_value(meta: Dict[str, str], *keys: str) -> str:
    if not meta:
        return ""
    lowered = {k.lower(): v for k, v in meta.items()}
    for key in keys:
        if key.lower() in lowered:
            return str(lowered[key.lower()]).strip()
    return ""


def infer_trait(meta: Dict[str, str]) -> str:
    return get_meta_value(meta, "trait_reported", "trait")


def _coerce_float(value: object) -> Optional[float]:
    if isinstance(value, str):
        value = value.replace(",", "").strip()
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _coerce_int(value: object) -> Optional[int]:
    val = _coerce_float(value)
    return int(val) if val is not None else None


def _parse_ancestry_blob(raw: str) -> Optional[object]:
    if not raw:
        return None
    text = raw.strip()
    if (text.startswith("{") and text.endswith("}")) or (text.startswith("[") and text.endswith("]")):
        try:
            return json.loads(text)
        except json.JSONDecodeError:
            try:
                return json.loads(text.replace("'", "\""))
            except json.JSONDecodeError:
                return None
    return None


def _classify_stage(entry: Dict[str, object]) -> str:
    stage_raw = str(
        entry.get("sample_type")
        or entry.get("study_stage")
        or entry.get("ancestry_type")
        or entry.get("stage")
        or ""
    ).lower()
    if any(token in stage_raw for token in ["eval", "validation", "replication", "test", "follow"]):
        return "eval"
    return "gwas"


def _extract_label(entry: Dict[str, object]) -> str:
    label = (
        entry.get("ancestry_broad")
        or entry.get("ancestry_category")
        or entry.get("ancestry")
        or entry.get("ancestry_country")
        or entry.get("ancestry_label")
        or "Unknown"
    )
    label_str = str(label).strip()
    return label_str or "Unknown"


def _extract_count(entry: Dict[str, object]) -> Optional[float]:
    raw = (
        entry.get("number")
        or entry.get("sample_number")
        or entry.get("sample_size")
        or entry.get("samples_number")
        or entry.get("ancestry_number")
        or entry.get("count")
    )
    return _coerce_float(raw)


def _sanitize_stage_blob(stage_blob: object) -> Dict[str, object]:
    if not isinstance(stage_blob, dict):
        return {}
    dist_raw = stage_blob.get("dist")
    dist: Optional[Dict[str, float]] = None
    if isinstance(dist_raw, dict):
        cleaned: Dict[str, float] = {}
        for k, v in dist_raw.items():
            coerced = _coerce_float(v)
            if coerced is not None:
                cleaned[str(k)] = coerced
        if cleaned:
            dist = cleaned
    count_val = _coerce_int(stage_blob.get("count"))
    stage: Dict[str, object] = {}
    if dist:
        stage["dist"] = dist
    if count_val is not None:
        stage["count"] = count_val
    return stage


def _aggregate_ancestry_entries(entries: object) -> Dict[str, Dict[str, object]]:
    if not isinstance(entries, list):
        return {}
    stage_counts: Dict[str, Dict[str, float]] = {"gwas": {}, "eval": {}}
    stage_pct: Dict[str, Dict[str, float]] = {"gwas": {}, "eval": {}}
    stage_totals: Dict[str, float] = {"gwas": 0.0, "eval": 0.0}

    for entry in entries:
        if not isinstance(entry, dict):
            continue
        stage = _classify_stage(entry)
        label = _extract_label(entry)
        count = _extract_count(entry)
        pct = _coerce_float(entry.get("percentage") or entry.get("percent"))
        if count is not None:
            stage_counts[stage][label] = stage_counts[stage].get(label, 0.0) + float(count)
            stage_totals[stage] += float(count)
        elif pct is not None:
            stage_pct[stage][label] = stage_pct[stage].get(label, 0.0) + float(pct)

    results: Dict[str, Dict[str, object]] = {}
    for stage in ("gwas", "eval"):
        dist = None
        if stage_counts[stage]:
            total = stage_totals[stage] if stage_totals[stage] > 0 else sum(stage_counts[stage].values())
            if total > 0:
                dist = {k: 100 * v / total for k, v in stage_counts[stage].items()}
        elif stage_pct[stage]:
            total_pct = sum(stage_pct[stage].values())
            if total_pct > 0:
                dist = {k: 100 * v / total_pct for k, v in stage_pct[stage].items()}

        stage_out: Dict[str, object] = {}
        if dist:
            stage_out["dist"] = dist
        count_total = int(stage_totals[stage]) if stage_totals[stage] > 0 else None
        if count_total is not None:
            stage_out["count"] = count_total
        if stage_out:
            results["gwas" if stage == "gwas" else "eval"] = stage_out
    return results


def _normalize_ancestry_payload(data: object) -> Dict[str, Dict[str, object]]:
    if isinstance(data, dict):
        out: Dict[str, Dict[str, object]] = {}
        gwas_stage = _sanitize_stage_blob(data.get("gwas"))
        eval_stage = _sanitize_stage_blob(data.get("eval"))
        dev_stage = _sanitize_stage_blob(data.get("dev"))
        if gwas_stage:
            out["gwas"] = gwas_stage
        if eval_stage:
            out["eval"] = eval_stage
        if dev_stage:
            out["dev"] = dev_stage
        if out:
            return out
    return _aggregate_ancestry_entries(data)


def _unpack_ancestry(ancestry: Dict[str, Dict[str, object]]) -> Tuple[Optional[Dict[str, float]], Optional[int], Optional[int]]:
    gwas_stage = ancestry.get("gwas") or {}
    dev_stage = ancestry.get("dev") or {}
    eval_stage = ancestry.get("eval") or {}

    primary_stage = gwas_stage if gwas_stage else (dev_stage if dev_stage else eval_stage)
    gwas_dist = primary_stage.get("dist") if isinstance(primary_stage.get("dist"), dict) else None
    gwas_count = _coerce_int(primary_stage.get("count"))

    dev_count = _coerce_int(dev_stage.get("count")) or _coerce_int(eval_stage.get("count"))
    return gwas_dist, gwas_count, dev_count


def get_pgs_ancestry(pgs_id: str) -> Dict[str, Dict[str, object]]:
    if not pgs_id:
        return {}
    if pgs_id in _PGS_CACHE:
        return _PGS_CACHE[pgs_id]

    url = f"https://www.pgscatalog.org/rest/score/{pgs_id}"
    try:
        resp = requests.get(url, headers={"Accept": "application/json"}, timeout=60)
        resp.raise_for_status()
        payload = resp.json()
    except (requests.RequestException, ValueError):
        _PGS_CACHE[pgs_id] = {}
        return {}

    ancestry = _normalize_ancestry_payload(payload.get("ancestry_distribution"))
    _PGS_CACHE[pgs_id] = ancestry
    return ancestry


def extract_ancestry_from_meta(
    meta: Dict[str, str],
    pgs_id: Optional[str] = None,
) -> Tuple[Optional[Dict[str, float]], Optional[int], Optional[int]]:
    candidates = [
        get_meta_value(meta, "ancestry_distribution"),
        get_meta_value(meta, "ancestry"),
        get_meta_value(meta, "pgs_ancestry_distribution"),
    ]
    for raw in candidates:
        parsed = _parse_ancestry_blob(raw)
        ancestry_from_meta = _normalize_ancestry_payload(parsed) if parsed is not None else {}
        gwas_dist, gwas_count, dev_count = _unpack_ancestry(ancestry_from_meta)
        if gwas_dist or gwas_count or dev_count:
            return gwas_dist, gwas_count, dev_count

    ancestry = get_pgs_ancestry(pgs_id) if pgs_id else {}
    return _unpack_ancestry(ancestry)


def parse_pgs_id(meta: Dict[str, str], path: str) -> str:
    in_meta = get_meta_value(meta, "pgs_id", "pgs", "pgs accession", "pgsid")
    if in_meta:
        m = _PGS_ID_RE.search(in_meta)
        if m:
            return m.group(1).upper()
    base = os.path.basename(path)
    m = _PGS_ID_RE.search(base)
    if m:
        return m.group(1).upper()
    return ""
