#!/usr/bin/env python3
"""Local RNA-seq analysis web app backed by R/Bioconductor DESeq2."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from email.parser import BytesParser
from email.policy import default as email_policy
import io
import json
import math
import os
import random
import socket
import subprocess
import tempfile
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


APP_DIR = Path(__file__).resolve().parent

DEFAULT_RECIPE: dict[str, Any] = {
    "deg": {
        "padj_cutoff": 0.05,
        "log2fc_cutoff": 1.0,
        "gene_types": ["protein_coding"],
        "protein_coding_only_for_labels": True,
        "protein_coding_only_for_heatmap": True,
        "protein_coding_only_for_enrichment": True,
    },
    "volcano": {
        "xlim": [-7, 7],
        "max_labels": 20,
        "colors": {"UP": "#D62728", "DOWN": "#1F77B4", "NS": "#B3B3B3"},
    },
    "heatmap": {
        "top_genes": 50,
        "clip": 2.0,
        "colors": {"low": "#1F77B4", "mid": "#F8FAFC", "high": "#D62728"},
        "group_palette": ["#4D4D4D", "#D62728", "#1F77B4", "#047C7B", "#B7791F", "#28724F", "#8B4F9F"],
        "group_colors": {"control": "#4D4D4D", "KO1": "#D62728", "KO2": "#1F77B4"},
    },
    "enrichment": {
        "min_overlap": 2,
        "max_terms": 80,
    },
}


@dataclass
class UploadedFile:
    filename: str
    file: io.BytesIO


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged: dict[str, Any] = {}
    for key, value in base.items():
        if isinstance(value, dict):
            merged[key] = deep_merge(value, {})
        else:
            merged[key] = value
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def parse_recipe(params: dict[str, Any]) -> dict[str, Any]:
    raw = params.get("analysis_recipe")
    if not raw:
        recipe = deep_merge(DEFAULT_RECIPE, {})
    elif isinstance(raw, str):
        recipe = deep_merge(DEFAULT_RECIPE, json.loads(raw))
    elif isinstance(raw, dict):
        recipe = deep_merge(DEFAULT_RECIPE, raw)
    else:
        raise ValueError("分析 recipe 格式不正确。")

    recipe["deg"]["padj_cutoff"] = float(recipe["deg"].get("padj_cutoff", 0.05))
    recipe["deg"]["log2fc_cutoff"] = float(recipe["deg"].get("log2fc_cutoff", 1.0))
    gene_types = recipe["deg"].get("gene_types", ["protein_coding"])
    if not isinstance(gene_types, list):
        gene_types = []
    recipe["deg"]["gene_types"] = [str(item).strip() for item in gene_types if str(item).strip()]
    recipe["volcano"]["max_labels"] = int(recipe["volcano"].get("max_labels", 20))
    recipe["heatmap"]["top_genes"] = int(recipe["heatmap"].get("top_genes", 50))
    recipe["heatmap"]["clip"] = float(recipe["heatmap"].get("clip", 2.0))
    recipe["enrichment"]["min_overlap"] = int(recipe["enrichment"].get("min_overlap", 2))
    recipe["enrichment"]["max_terms"] = int(recipe["enrichment"].get("max_terms", 80))
    normalize_colors(recipe)
    return recipe


def is_hex_color(value: Any) -> bool:
    text = str(value or "")
    return len(text) == 7 and text.startswith("#") and all(ch in "0123456789abcdefABCDEF" for ch in text[1:])


def normalize_colors(recipe: dict[str, Any]) -> None:
    defaults = DEFAULT_RECIPE
    for key in ("UP", "DOWN", "NS"):
        value = recipe["volcano"].setdefault("colors", {}).get(key)
        if not is_hex_color(value):
            recipe["volcano"]["colors"][key] = defaults["volcano"]["colors"][key]
    for key in ("low", "mid", "high"):
        value = recipe["heatmap"].setdefault("colors", {}).get(key)
        if not is_hex_color(value):
            recipe["heatmap"]["colors"][key] = defaults["heatmap"]["colors"][key]
    palette = recipe["heatmap"].get("group_palette")
    if not isinstance(palette, list):
        palette = []
    clean_palette = [str(color) for color in palette if is_hex_color(color)]
    if len(clean_palette) < 3:
        clean_palette = defaults["heatmap"]["group_palette"].copy()
    recipe["heatmap"]["group_palette"] = clean_palette
    group_colors = recipe["heatmap"].get("group_colors")
    if not isinstance(group_colors, dict):
        group_colors = {}
    recipe["heatmap"]["group_colors"] = {str(group): str(color) for group, color in group_colors.items() if is_hex_color(color)}


def json_response(handler: SimpleHTTPRequestHandler, payload: dict[str, Any], status: int = 200) -> None:
    body = json.dumps(payload, ensure_ascii=False, allow_nan=False).encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "application/json; charset=utf-8")
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def binary_response(
    handler: SimpleHTTPRequestHandler,
    content: bytes,
    filename: str,
    content_type: str,
    status: int = 200,
) -> None:
    handler.send_response(status)
    handler.send_header("Content-Type", content_type)
    handler.send_header("Content-Disposition", f'attachment; filename="{filename}"')
    handler.send_header("Content-Length", str(len(content)))
    handler.end_headers()
    handler.wfile.write(content)


def error_response(handler: SimpleHTTPRequestHandler, message: str, status: int = 400) -> None:
    json_response(handler, {"ok": False, "error": message}, status)


def local_network_addresses(port: int) -> list[str]:
    addresses: set[str] = set()
    try:
        hostname = socket.gethostname()
        for item in socket.getaddrinfo(hostname, port, family=socket.AF_INET, type=socket.SOCK_STREAM):
            ip = item[4][0]
            if ip and not ip.startswith("127."):
                addresses.add(f"http://{ip}:{port}")
    except OSError:
        pass
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as sock:
            sock.connect(("8.8.8.8", 80))
            ip = sock.getsockname()[0]
            if ip and not ip.startswith("127."):
                addresses.add(f"http://{ip}:{port}")
    except OSError:
        pass
    return sorted(addresses)


def safe_sheet_name(name: str) -> str:
    cleaned = "".join(ch if ch not in '[]:*?/\\' else "_" for ch in str(name))
    return (cleaned[:31] or "sheet").strip()


def export_deseq_workbook(payload: dict[str, Any]) -> bytes:
    output = io.BytesIO()
    summary = payload.get("summary") or {}
    params = payload.get("params") or {}
    recipe = payload.get("recipe") or {}
    selected_gene_types = ((recipe.get("deg") or {}).get("gene_types") or [])
    contrasts = payload.get("contrasts") or []
    enrichment = payload.get("enrichment") or {}

    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        pd.DataFrame(
            [
                {"key": key, "value": json.dumps(value, ensure_ascii=False) if isinstance(value, (dict, list)) else value}
                for key, value in {**summary, **params}.items()
            ]
        ).to_excel(writer, index=False, sheet_name="summary")

        if recipe:
            pd.DataFrame({"analysis_recipe_json": [json.dumps(recipe, ensure_ascii=False, indent=2)]}).to_excel(
                writer, index=False, sheet_name="recipe"
            )

        flat_frames = []
        used_names: set[str] = set()
        for contrast in contrasts:
            label = contrast.get("label", "contrast")
            rows = contrast.get("results") or []
            if not rows:
                continue
            frame = pd.DataFrame(rows)
            frame.insert(0, "contrast", label)
            flat_frames.append(frame)
            group = str(contrast.get("group") or label.split()[0] or "contrast")
            sheet = safe_sheet_name(f"{group}_all")
            base = sheet
            suffix = 2
            while sheet in used_names:
                sheet = safe_sheet_name(f"{base}_{suffix}")
                suffix += 1
            used_names.add(sheet)
            sheet_frame = frame
            if selected_gene_types and "gene_type" in sheet_frame.columns:
                sheet_frame = sheet_frame[sheet_frame["gene_type"].isin(selected_gene_types)]
            sheet_frame.to_excel(writer, index=False, sheet_name=sheet)

        if flat_frames:
            pd.concat(flat_frames, ignore_index=True).to_excel(writer, index=False, sheet_name="all_DESeq_results")

        for key, sheet in (("go", "GO_enrichment"), ("kegg", "KEGG_enrichment")):
            rows = enrichment.get(key) or []
            if rows:
                pd.DataFrame(rows).to_excel(writer, index=False, sheet_name=sheet)

    return output.getvalue()


def clean_header(value: Any) -> str:
    text = str(value).strip()
    return text if text and text.lower() != "nan" else "Unnamed"


def read_uploaded_table(file_item: UploadedFile) -> tuple[pd.DataFrame, pd.DataFrame | None, pd.DataFrame | None, dict[str, Any]]:
    filename = (file_item.filename or "").lower()
    content = file_item.file.read()
    if not content:
        raise ValueError("上传文件为空。")

    metadata: pd.DataFrame | None = None
    annotations: pd.DataFrame | None = None
    info: dict[str, Any] = {"filename": file_item.filename or "uploaded file", "sheets": []}

    if filename.endswith((".xlsx", ".xlsm", ".xltx", ".xltm", ".xls")):
        excel = pd.ExcelFile(io.BytesIO(content))
        info["sheets"] = excel.sheet_names
        lowered = {name.lower(): name for name in excel.sheet_names}
        count_sheet = excel.sheet_names[0]
        for candidate in ("counts", "raw_counts", "raw counts", "count", "matrix"):
            if candidate in lowered:
                count_sheet = lowered[candidate]
                break

        counts = pd.read_excel(excel, sheet_name=count_sheet)
        for name in excel.sheet_names:
            lower = name.lower()
            if lower in {"metadata", "meta", "sample", "samples", "design"}:
                metadata = pd.read_excel(excel, sheet_name=name)
            if any(token in lower for token in ("annotation", "geneset", "gene_set", "go", "kegg", "enrichment")):
                annotations = pd.read_excel(excel, sheet_name=name)
        info["count_sheet"] = count_sheet
    elif filename.endswith((".csv", ".txt", ".tsv")):
        sep = "\t" if filename.endswith((".tsv", ".txt")) else ","
        counts = pd.read_csv(io.BytesIO(content), sep=sep)
        info["count_sheet"] = "file"
    else:
        raise ValueError("请上传 .xlsx、.xls、.csv 或 .tsv 文件。")

    counts.columns = [clean_header(c) for c in counts.columns]
    if metadata is not None:
        metadata.columns = [clean_header(c) for c in metadata.columns]
    if annotations is not None:
        annotations.columns = [clean_header(c) for c in annotations.columns]
    return counts, metadata, annotations, info


def numeric_sample_columns(df: pd.DataFrame, gene_col: str | None = None) -> list[str]:
    columns: list[str] = []
    for col in df.columns:
        if col == gene_col:
            continue
        values = pd.to_numeric(df[col], errors="coerce")
        if values.notna().sum() >= max(3, len(values) * 0.5):
            columns.append(col)
    return columns


def infer_metadata(metadata: pd.DataFrame | None, sample_cols: list[str]) -> dict[str, Any] | None:
    if metadata is None or metadata.empty:
        return None
    lower = {c.lower(): c for c in metadata.columns}
    sample_col = lower.get("sample") or lower.get("sample_id") or lower.get("samples")
    condition_col = lower.get("condition") or lower.get("group") or lower.get("treatment")
    if sample_col is None or condition_col is None:
        return None

    rows = []
    for _, row in metadata.iterrows():
        sample = str(row[sample_col]).strip()
        condition = str(row[condition_col]).strip()
        if sample in sample_cols and condition:
            rows.append({"sample": sample, "condition": condition})
    conditions = sorted({row["condition"] for row in rows})
    return {"rows": rows, "conditions": conditions}


def parse_multipart_form(headers: Any, body: bytes) -> tuple[list[UploadedFile], dict[str, str]]:
    content_type = headers.get("content-type", "")
    raw_message = (
        f"Content-Type: {content_type}\r\n"
        "MIME-Version: 1.0\r\n"
        "\r\n"
    ).encode("utf-8") + body
    message = BytesParser(policy=email_policy).parsebytes(raw_message)
    if not message.is_multipart():
        raise ValueError("请使用表单上传文件。")

    files: list[UploadedFile] = []
    fields: dict[str, str] = {}
    for part in message.iter_parts():
        if part.get_content_disposition() != "form-data":
            continue
        name = part.get_param("name", header="content-disposition")
        if not name:
            continue
        payload = part.get_payload(decode=True) or b""
        filename = part.get_filename()
        if filename:
            files.append(UploadedFile(filename=filename, file=io.BytesIO(payload)))
        else:
            charset = part.get_content_charset() or "utf-8"
            fields[name] = payload.decode(charset, errors="replace")
    return files, fields


def infer_group_from_sample(sample: str) -> str:
    cleaned = sample.strip()
    for sep in ("__", "_", "-", "."):
        if sep in cleaned:
            head = cleaned.split(sep)[0].strip()
            if head:
                return head
    while cleaned and cleaned[-1].isdigit():
        cleaned = cleaned[:-1]
    return cleaned.strip("_-. ") or sample


def metadata_group_map(metadata: pd.DataFrame | None, sample_cols: list[str]) -> dict[str, str]:
    inferred = infer_metadata(metadata, sample_cols)
    if inferred is None:
        return {}
    return {row["sample"]: row["condition"] for row in inferred["rows"]}


def unique_sample_name(name: str, batch: str, used: set[str]) -> str:
    candidate = name
    if candidate in used:
        candidate = f"{batch}__{name}"
    suffix = 2
    base = candidate
    while candidate in used:
        candidate = f"{base}_{suffix}"
        suffix += 1
    used.add(candidate)
    return candidate


def prepare_uploaded_dataset(
    file_items: list[UploadedFile],
    batch_overrides: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    if not file_items:
        raise ValueError("请选择一个或多个数据文件。")

    tables: list[pd.DataFrame] = []
    all_meta: list[dict[str, str]] = []
    annotations: list[pd.DataFrame] = []
    file_info: list[dict[str, Any]] = []
    used_samples: set[str] = set()

    batch_overrides = batch_overrides or []
    for index, item in enumerate(file_items, start=1):
        counts, metadata, sheet_annotation, info = read_uploaded_table(item)
        override = batch_overrides[index - 1] if index - 1 < len(batch_overrides) else {}
        batch = str(override.get("batch") or Path(info["filename"]).stem or f"batch{index}").strip()
        gene_col = "gene" if "gene" in counts.columns else counts.columns[0]
        sample_cols = numeric_sample_columns(counts, gene_col)
        if not sample_cols:
            raise ValueError(f"{info['filename']} 没有检测到 raw count 样本列。")

        group_map = metadata_group_map(metadata, sample_cols)
        rename: dict[str, str] = {}
        for sample in sample_cols:
            renamed = unique_sample_name(sample, batch, used_samples)
            rename[sample] = renamed
            all_meta.append(
                {
                    "sample": renamed,
                    "original_sample": sample,
                    "batch": batch,
                    "group": group_map.get(sample) or infer_group_from_sample(sample),
                    "file": info["filename"],
                }
            )

        annotation_cols = [
            col
            for col in counts.columns
            if col not in sample_cols and col != gene_col
        ]
        keep_cols = [gene_col] + annotation_cols + sample_cols
        table = counts[keep_cols].copy().rename(columns={gene_col: "gene", **rename})
        tables.append(table)

        local_annotation_cols = [col for col in ("gene", "gene_name", "gene_type", "entrezgene_id") if col in table.columns]
        if len(local_annotation_cols) > 1:
            annotations.append(table[local_annotation_cols].drop_duplicates("gene"))
        if sheet_annotation is not None:
            standardized = standardize_annotation_columns(sheet_annotation)
            if standardized is not None:
                annotations.append(standardized)

        file_info.append(
            {
                "filename": info["filename"],
                "batch": batch,
                "count_sheet": info.get("count_sheet", "file"),
                "samples": [rename[sample] for sample in sample_cols],
                "original_samples": sample_cols,
            }
        )

    common_genes = set(tables[0]["gene"].astype(str))
    for table in tables[1:]:
        common_genes &= set(table["gene"].astype(str))
    if not common_genes:
        raise ValueError("多个文件之间没有共同 gene，无法合并分析。")

    merged: pd.DataFrame | None = None
    annotation_frames = []
    for table in tables:
        table = table.copy()
        table["gene"] = table["gene"].astype(str)
        table = table[table["gene"].isin(common_genes)].drop_duplicates("gene", keep="first")
        sample_cols = [row["sample"] for row in all_meta if row["sample"] in table.columns]
        annotation_cols = [col for col in table.columns if col not in sample_cols and col != "gene"]
        if annotation_cols:
            annotation_frames.append(table[["gene"] + annotation_cols])
        compact = table[["gene"] + sample_cols]
        merged = compact if merged is None else merged.merge(compact, on="gene", how="inner")

    if merged is None:
        raise ValueError("没有可合并的数据。")

    gene_annotation = None
    if annotation_frames:
        gene_annotation = pd.concat(annotation_frames, ignore_index=True)
        gene_annotation = gene_annotation.replace("", np.nan).groupby("gene", as_index=False).first()
        merged = merged.merge(gene_annotation, on="gene", how="left")

    term_annotation = None
    term_frames = [frame for frame in annotations if {"gene", "term_id", "term_name", "category"}.issubset(frame.columns)]
    if term_frames:
        term_annotation = pd.concat(term_frames, ignore_index=True).drop_duplicates()

    sample_cols = [row["sample"] for row in all_meta]
    for col in sample_cols:
        merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0).clip(lower=0)

    groups = sorted({row["group"] for row in all_meta}, key=lambda value: (value.lower() != "control", value.lower()))
    return {
        "counts": merged,
        "sample_cols": sample_cols,
        "sample_meta": all_meta,
        "groups": groups,
        "gene_annotation": gene_annotation,
        "term_annotation": term_annotation,
        "file_info": file_info,
        "common_gene_count": len(common_genes),
    }


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=float)
    adjusted = np.ones_like(pvalues)
    finite = np.isfinite(pvalues)
    if not finite.any():
        return adjusted
    idx = np.where(finite)[0]
    order = idx[np.argsort(pvalues[idx])]
    ranked = pvalues[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted[order] = np.clip(ranked, 0, 1)
    return adjusted


def normal_two_sided_p(z: np.ndarray) -> np.ndarray:
    return np.array([math.erfc(abs(float(value)) / math.sqrt(2)) for value in z], dtype=float)


def median_ratio_size_factors(counts: pd.DataFrame) -> pd.Series:
    matrix = counts.to_numpy(dtype=float)
    positive_all = (matrix > 0).all(axis=1)
    usable = matrix[positive_all]
    if usable.shape[0] >= 2:
        geo_means = np.exp(np.mean(np.log(usable), axis=1))
        ratios = usable / geo_means[:, None]
        factors = np.median(ratios, axis=0)
    else:
        totals = matrix.sum(axis=0)
        median_total = np.median(totals[totals > 0]) if (totals > 0).any() else 1.0
        factors = totals / median_total
    factors = np.where(np.isfinite(factors) & (factors > 0), factors, 1.0)
    return pd.Series(factors, index=counts.columns)


def run_differential_expression(
    counts: pd.DataFrame,
    gene_col: str,
    group_a: list[str],
    group_b: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected = group_a + group_b
    raw = counts[[gene_col] + selected].copy()
    raw[gene_col] = raw[gene_col].astype(str).str.strip()
    raw = raw[raw[gene_col].ne("")].drop_duplicates(subset=[gene_col], keep="first")
    for col in selected:
        raw[col] = pd.to_numeric(raw[col], errors="coerce").fillna(0).clip(lower=0)

    count_matrix = raw[selected]
    size_factors = median_ratio_size_factors(count_matrix)
    normalized = count_matrix.div(size_factors, axis=1)

    a = normalized[group_a]
    b = normalized[group_b]
    mean_a = a.mean(axis=1)
    mean_b = b.mean(axis=1)
    base_mean = normalized.mean(axis=1)
    log2_fc = np.log2((mean_b + 1.0) / (mean_a + 1.0))

    log_norm = np.log2(normalized + 1.0)
    log_a = log_norm[group_a]
    log_b = log_norm[group_b]
    var_a = log_a.var(axis=1, ddof=1).fillna(0)
    var_b = log_b.var(axis=1, ddof=1).fillna(0)
    poisson_floor = (1 / (mean_a + 1) + 1 / (mean_b + 1)) / (math.log(2) ** 2)
    se = np.sqrt(var_a / max(len(group_a), 1) + var_b / max(len(group_b), 1) + poisson_floor)
    se = np.where(np.isfinite(se) & (se > 1e-9), se, 1.0)
    stat = log2_fc / se
    pvalue = normal_two_sided_p(stat)
    padj = bh_adjust(pvalue)

    dispersion = ((normalized.var(axis=1, ddof=1) - base_mean) / np.maximum(base_mean ** 2, 1e-9)).clip(lower=0)

    results = pd.DataFrame(
        {
            "gene": raw[gene_col].to_numpy(),
            "baseMean": base_mean.to_numpy(),
            "meanA": mean_a.to_numpy(),
            "meanB": mean_b.to_numpy(),
            "log2FoldChange": log2_fc.to_numpy(),
            "lfcSE": se,
            "stat": stat,
            "pvalue": pvalue,
            "padj": padj,
            "dispersion": dispersion.to_numpy(),
        }
    ).sort_values(["padj", "pvalue", "baseMean"], ascending=[True, True, False])

    normalized.insert(0, "gene", raw[gene_col].to_numpy())
    return results, normalized


def normalize_all_counts(counts: pd.DataFrame, sample_cols: list[str]) -> pd.DataFrame:
    count_matrix = counts[sample_cols].copy()
    for col in sample_cols:
        count_matrix[col] = pd.to_numeric(count_matrix[col], errors="coerce").fillna(0).clip(lower=0)
    size_factors = median_ratio_size_factors(count_matrix)
    normalized = count_matrix.div(size_factors, axis=1)
    normalized.insert(0, "gene", counts["gene"].astype(str).to_numpy())
    return normalized


def annotate_de_results(results: pd.DataFrame, counts: pd.DataFrame, recipe: dict[str, Any]) -> pd.DataFrame:
    annotation_cols = [col for col in ("gene_name", "gene_type", "entrezgene_id") if col in counts.columns]
    if annotation_cols:
        annotation = counts[["gene"] + annotation_cols].drop_duplicates("gene")
        results = results.merge(annotation, on="gene", how="left")
    if "gene_name" not in results.columns:
        results["gene_name"] = results["gene"]
    results["gene_name"] = results["gene_name"].fillna("").astype(str)
    results["display_gene"] = np.where(results["gene_name"].str.strip().ne(""), results["gene_name"], results["gene"])
    if "gene_type" not in results.columns:
        results["gene_type"] = ""
    results["gene_type"] = results["gene_type"].fillna("").astype(str)
    results["regulation"] = "NS"
    fc_cutoff = float(recipe["deg"]["log2fc_cutoff"])
    padj_cutoff = float(recipe["deg"]["padj_cutoff"])
    results.loc[(results["padj"] < padj_cutoff) & (results["log2FoldChange"] > fc_cutoff), "regulation"] = "UP"
    results.loc[(results["padj"] < padj_cutoff) & (results["log2FoldChange"] < -fc_cutoff), "regulation"] = "DOWN"
    return results


def batch_correct_log_matrix(normalized: pd.DataFrame, sample_meta: list[dict[str, str]], sample_cols: list[str]) -> pd.DataFrame:
    log_matrix = np.log2(normalized.set_index("gene")[sample_cols].astype(float) + 1.0)
    batches = pd.Series({row["sample"]: row["batch"] for row in sample_meta})
    overall = log_matrix.mean(axis=1)
    corrected = log_matrix.copy()
    for batch in batches.loc[sample_cols].unique():
        batch_samples = [sample for sample in sample_cols if batches[sample] == batch]
        if not batch_samples:
            continue
        batch_mean = log_matrix[batch_samples].mean(axis=1)
        corrected[batch_samples] = log_matrix[batch_samples].sub(batch_mean, axis=0).add(overall, axis=0)
    return corrected


def build_multi_heatmap(
    normalized: pd.DataFrame,
    contrast_results: list[dict[str, Any]],
    sample_meta: list[dict[str, str]],
    control_group: str,
    ko_groups: list[str],
    sample_cols: list[str],
    recipe: dict[str, Any],
) -> dict[str, Any]:
    frames = []
    selected_gene_types = recipe["deg"].get("gene_types") or []
    for contrast in contrast_results:
        df = pd.DataFrame(contrast["results"])
        if df.empty:
            continue
        keep = (df["regulation"] != "NS")
        if selected_gene_types and "gene_type" in df.columns and df["gene_type"].astype(str).str.strip().any():
            keep &= df["gene_type"].isin(selected_gene_types)
        frames.append(df.loc[keep, ["gene", "display_gene", "padj"]])
    if not frames:
        return {"genes": [], "gene_ids": [], "samples": sample_cols, "groups": [], "matrix": []}

    top = (
        pd.concat(frames, ignore_index=True)
        .sort_values("padj")
        .drop_duplicates("gene")
        .head(int(recipe["heatmap"]["top_genes"]))
    )
    if top.empty:
        return {"genes": [], "gene_ids": [], "samples": sample_cols, "groups": [], "matrix": []}

    group_rank = {control_group: 0}
    group_rank.update({group: index + 1 for index, group in enumerate(ko_groups)})
    meta_map = {row["sample"]: row for row in sample_meta}
    ordered_samples = sorted(
        [sample for sample in sample_cols if meta_map[sample]["group"] in group_rank],
        key=lambda sample: (group_rank.get(meta_map[sample]["group"], 99), meta_map[sample]["batch"], sample),
    )
    corrected = batch_correct_log_matrix(normalized, sample_meta, ordered_samples)
    expr = corrected.loc[top["gene"], ordered_samples]
    row_mean = expr.mean(axis=1)
    row_std = expr.std(axis=1).replace(0, 1)
    clip = float(recipe["heatmap"]["clip"])
    z = expr.sub(row_mean, axis=0).div(row_std, axis=0).clip(-clip, clip)
    return {
        "genes": top["display_gene"].astype(str).tolist(),
        "gene_ids": top["gene"].astype(str).tolist(),
        "samples": ordered_samples,
        "groups": [meta_map[sample]["group"] for sample in ordered_samples],
        "batches": [meta_map[sample]["batch"] for sample in ordered_samples],
        "matrix": z.round(4).values.tolist(),
    }


def build_heatmap(normalized: pd.DataFrame, results: pd.DataFrame, selected_samples: list[str], limit: int = 50) -> dict[str, Any]:
    top_genes = results[results["baseMean"] > 0].head(limit)["gene"].tolist()
    if not top_genes:
        return {"genes": [], "samples": selected_samples, "matrix": []}
    matrix = normalized.set_index("gene").loc[top_genes, selected_samples].astype(float)
    log_matrix = np.log2(matrix + 1.0)
    row_mean = log_matrix.mean(axis=1)
    row_std = log_matrix.std(axis=1).replace(0, 1)
    z = log_matrix.sub(row_mean, axis=0).div(row_std, axis=0).clip(-2.5, 2.5)
    return {
        "genes": top_genes,
        "samples": selected_samples,
        "matrix": z.round(4).values.tolist(),
    }


def standardize_annotation_columns(annotations: pd.DataFrame | None) -> pd.DataFrame | None:
    if annotations is None or annotations.empty:
        return None
    lower = {c.lower(): c for c in annotations.columns}
    gene_col = lower.get("gene") or lower.get("gene_id") or lower.get("symbol")
    term_col = lower.get("term_id") or lower.get("term") or lower.get("go") or lower.get("kegg") or lower.get("pathway_id")
    name_col = lower.get("term_name") or lower.get("name") or lower.get("pathway") or lower.get("description")
    category_col = lower.get("category") or lower.get("ontology") or lower.get("source") or lower.get("database")
    if gene_col is None or term_col is None:
        return None
    out = pd.DataFrame(
        {
            "gene": annotations[gene_col].astype(str).str.strip(),
            "term_id": annotations[term_col].astype(str).str.strip(),
            "term_name": annotations[name_col].astype(str).str.strip() if name_col else annotations[term_col].astype(str).str.strip(),
            "category": annotations[category_col].astype(str).str.strip() if category_col else "Annotation",
        }
    )
    out = out[out["gene"].ne("") & out["term_id"].ne("")]
    return out.drop_duplicates()


def log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k: int, population: int, successes: int, draws: int) -> float:
    upper = min(successes, draws)
    logs = []
    for x in range(k, upper + 1):
        logs.append(log_choose(successes, x) + log_choose(population - successes, draws - x) - log_choose(population, draws))
    if not logs:
        return 1.0
    max_log = max(logs)
    return min(1.0, math.exp(max_log) * sum(math.exp(value - max_log) for value in logs))


def run_enrichment(
    results: pd.DataFrame,
    annotations: pd.DataFrame | None,
    fc_cutoff: float,
    padj_cutoff: float,
) -> dict[str, Any]:
    standardized = standardize_annotation_columns(annotations)
    if standardized is None:
        return {
            "available": False,
            "message": "未找到注释表。可在 Excel 中加入 annotation/go/kegg sheet，列名包含 gene、term_id、term_name、category。",
            "terms": [],
        }

    universe = set(results["gene"].astype(str))
    sig = set(
        results[
            (results["padj"] <= padj_cutoff)
            & (results["log2FoldChange"].abs() >= fc_cutoff)
            & (results["baseMean"] > 0)
        ]["gene"].astype(str)
    )
    ann = standardized[standardized["gene"].isin(universe)]
    population = len(universe)
    draws = len(sig)
    if population == 0 or draws == 0:
        return {"available": True, "message": "当前阈值下没有显著基因。", "terms": []}

    rows: list[dict[str, Any]] = []
    grouped = ann.groupby(["category", "term_id", "term_name"])
    for (category, term_id, term_name), frame in grouped:
        term_genes = set(frame["gene"])
        successes = len(term_genes)
        overlap = sorted(term_genes & sig)
        k = len(overlap)
        if k < 2:
            continue
        pvalue = hypergeom_sf(k, population, successes, draws)
        rows.append(
            {
                "category": category,
                "term_id": term_id,
                "term_name": term_name,
                "overlap": k,
                "term_size": successes,
                "significant_genes": draws,
                "pvalue": pvalue,
                "genes": ", ".join(overlap[:30]),
            }
        )

    if not rows:
        return {"available": True, "message": "没有达到最小重叠数的富集条目。", "terms": []}
    enriched = pd.DataFrame(rows)
    enriched["padj"] = bh_adjust(enriched["pvalue"].to_numpy())
    enriched["gene_ratio"] = enriched["overlap"] / enriched["significant_genes"]
    enriched = enriched.sort_values(["padj", "pvalue", "overlap"], ascending=[True, True, False]).head(80)
    return {"available": True, "message": "", "terms": clean_records(enriched)}


def run_multi_enrichment(
    contrast_results: list[dict[str, Any]],
    annotations: pd.DataFrame | None,
    recipe: dict[str, Any],
) -> dict[str, Any]:
    standardized = standardize_annotation_columns(annotations)
    if standardized is None:
        return {
            "available": False,
            "message": "未找到 GO/KEGG 注释表。可在 Excel 中加入 annotation/go/kegg sheet，列名包含 gene、term_id、term_name、category。",
            "terms": [],
        }

    rows: list[dict[str, Any]] = []
    fc_cutoff = float(recipe["deg"]["log2fc_cutoff"])
    padj_cutoff = float(recipe["deg"]["padj_cutoff"])
    min_overlap = int(recipe["enrichment"]["min_overlap"])
    max_terms = int(recipe["enrichment"]["max_terms"])
    selected_gene_types = recipe["deg"].get("gene_types") or []
    for contrast in contrast_results:
        results = pd.DataFrame(contrast["results"])
        if results.empty:
            continue
        universe = set(results["gene"].astype(str))
        ann = standardized[standardized["gene"].isin(universe)]
        population = len(universe)
        if population == 0:
            continue
        for direction, sign in (("UP", 1), ("DOWN", -1)):
            if sign > 0:
                sig_frame = results[
                    (results["padj"] <= padj_cutoff)
                    & (results["log2FoldChange"] >= fc_cutoff)
                ]
            else:
                sig_frame = results[
                    (results["padj"] <= padj_cutoff)
                    & (results["log2FoldChange"] <= -fc_cutoff)
                ]
            if selected_gene_types and "gene_type" in sig_frame.columns and sig_frame["gene_type"].astype(str).str.strip().any():
                sig_frame = sig_frame[sig_frame["gene_type"].isin(selected_gene_types)]
            sig = set(sig_frame["gene"].astype(str))
            draws = len(sig)
            if draws == 0:
                continue
            for (category, term_id, term_name), frame in ann.groupby(["category", "term_id", "term_name"]):
                term_genes = set(frame["gene"])
                successes = len(term_genes)
                overlap = sorted(term_genes & sig)
                k = len(overlap)
                if k < min_overlap:
                    continue
                pvalue = hypergeom_sf(k, population, successes, draws)
                rows.append(
                    {
                        "contrast": contrast["label"],
                        "direction": direction,
                        "category": category,
                        "term_id": term_id,
                        "term_name": term_name,
                        "overlap": k,
                        "term_size": successes,
                        "significant_genes": draws,
                        "pvalue": pvalue,
                        "genes": ", ".join(overlap[:30]),
                    }
                )

    if not rows:
        return {"available": True, "message": "当前阈值下没有达到最小重叠数的富集条目。", "terms": [], "go": [], "kegg": []}
    enriched = pd.DataFrame(rows)
    enriched["padj"] = bh_adjust(enriched["pvalue"].to_numpy())
    enriched["gene_ratio"] = enriched["overlap"] / enriched["significant_genes"]
    enriched = enriched.sort_values(["padj", "pvalue", "overlap"], ascending=[True, True, False]).head(max_terms * 2)
    is_kegg = enriched["category"].astype(str).str.lower().str.contains("kegg")
    go = enriched[~is_kegg].head(max_terms)
    kegg = enriched[is_kegg].head(max_terms)
    return {
        "available": True,
        "message": "",
        "terms": clean_records(enriched),
        "go": clean_records(go),
        "kegg": clean_records(kegg),
    }


def clean_records(df: pd.DataFrame, limit: int | None = None) -> list[dict[str, Any]]:
    if limit is not None:
        df = df.head(limit)
    clean = df.replace([np.inf, -np.inf], np.nan)
    records = clean.where(pd.notnull(clean), None).to_dict(orient="records")
    for row in records:
        for key, value in list(row.items()):
            if isinstance(value, np.generic):
                row[key] = value.item()
    return records


def run_r_deseq_pipeline(
    counts: pd.DataFrame,
    sample_cols: list[str],
    sample_meta: list[dict[str, Any]],
    control_group: str,
    ko_groups: list[str],
    recipe: dict[str, Any],
    species: str = "human",
) -> dict[str, Any]:
    script = APP_DIR / "deseq_pipeline.R"
    if not script.exists():
        raise ValueError("找不到 R 分析脚本 deseq_pipeline.R。")

    annotation_cols = [col for col in ("gene_name", "gene_type", "entrezgene_id") if col in counts.columns]
    export_cols = ["gene"] + sample_cols + annotation_cols
    export_counts = counts[export_cols].copy()
    export_counts["gene"] = export_counts["gene"].astype(str)
    for col in sample_cols:
        export_counts[col] = pd.to_numeric(export_counts[col], errors="coerce").fillna(0).round().astype(int)

    meta_df = pd.DataFrame(sample_meta)
    meta_df = meta_df[["sample", "batch", "group"]]
    params = {
        "control_group": control_group,
        "ko_groups": ko_groups,
        "analysis_recipe": recipe,
        "species": species,
    }

    with tempfile.TemporaryDirectory(prefix="rna_seq_deseq_") as tmp:
        tmp_path = Path(tmp)
        counts_path = tmp_path / "counts.csv"
        metadata_path = tmp_path / "metadata.csv"
        params_path = tmp_path / "params.json"
        output_path = tmp_path / "result.json"
        export_counts.to_csv(counts_path, index=False)
        meta_df.to_csv(metadata_path, index=False)
        params_path.write_text(json.dumps(params, ensure_ascii=False), encoding="utf-8")

        cmd = ["Rscript", str(script), str(counts_path), str(metadata_path), str(params_path), str(output_path)]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if proc.returncode != 0:
            message = (proc.stderr or proc.stdout or "R DESeq2 分析失败").strip()
            raise ValueError(f"R DESeq2 分析失败：{message[-2000:]}")
        if not output_path.exists():
            raise ValueError("R DESeq2 分析没有生成结果。")
        return json.loads(output_path.read_text(encoding="utf-8"))


def analyze_payload(file_items: list[UploadedFile], params: dict[str, Any]) -> dict[str, Any]:
    recipe = parse_recipe(params)
    dataset = prepare_uploaded_dataset(file_items, params.get("batch_overrides") or [])
    counts = dataset["counts"]
    sample_meta = params.get("sample_meta") or dataset["sample_meta"]
    if not isinstance(sample_meta, list):
        raise ValueError("样本分组格式不正确。")
    meta_map = {row["sample"]: row for row in sample_meta if row.get("sample") in dataset["sample_cols"]}
    if len(meta_map) != len(dataset["sample_cols"]):
        raise ValueError("样本分组信息不完整。")
    sample_meta = [meta_map[sample] for sample in dataset["sample_cols"]]
    for row in sample_meta:
        row["group"] = str(row.get("group") or "").strip()
        if not row["group"]:
            raise ValueError(f"{row['sample']} 没有 group。")

    control_group = str(params.get("control_group") or "control").strip()
    species = str(params.get("species") or "human").strip().lower()
    if species not in {"human", "mouse"}:
        raise ValueError("Species 只能选择 human 或 mouse。")
    ko_groups = params.get("ko_groups") or []
    if not isinstance(ko_groups, list):
        raise ValueError("KO 分组格式不正确。")
    ko_groups = [str(group).strip() for group in ko_groups if str(group).strip() and str(group).strip() != control_group]
    if not ko_groups:
        raise ValueError("请至少选择一个 KO / treatment group。")

    control_samples = [row["sample"] for row in sample_meta if row["group"] == control_group]
    if not control_samples:
        raise ValueError(f"找不到 control group：{control_group}")
    for group in ko_groups:
        if not [row for row in sample_meta if row["group"] == group]:
            raise ValueError(f"找不到 KO / treatment group：{group}")

    payload = run_r_deseq_pipeline(
        counts=counts,
        sample_cols=dataset["sample_cols"],
        sample_meta=sample_meta,
        control_group=control_group,
        ko_groups=ko_groups,
        recipe=recipe,
        species=species,
    )
    payload["info"] = {
        "files": dataset["file_info"],
        "common_gene_count": dataset["common_gene_count"],
    }
    payload["params"] = {
        "control_group": control_group,
        "ko_groups": ko_groups,
        "species": species,
        "gene_types": recipe["deg"]["gene_types"],
        "fc_cutoff": recipe["deg"]["log2fc_cutoff"],
        "padj_cutoff": recipe["deg"]["padj_cutoff"],
    }
    if "gene_type" in counts.columns:
        payload["gene_types"] = sorted([str(item) for item in counts["gene_type"].dropna().unique() if str(item).strip()])
    return payload


def inspect_payload(file_items: list[UploadedFile], params: dict[str, Any] | None = None) -> dict[str, Any]:
    params = params or {}
    dataset = prepare_uploaded_dataset(file_items, params.get("batch_overrides") or [])
    counts = dataset["counts"]
    preview_cols = ["gene"]
    for col in ("gene_name", "gene_type"):
        if col in counts.columns:
            preview_cols.append(col)
    preview_cols.extend(dataset["sample_cols"][:8])
    preview = counts[preview_cols].head(6).replace([np.inf, -np.inf], np.nan).where(pd.notnull(counts[preview_cols].head(6)), None)
    gene_types = []
    if "gene_type" in counts.columns:
        gene_types = sorted([str(item) for item in counts["gene_type"].dropna().unique() if str(item).strip()])
    return {
        "ok": True,
        "info": {
            "files": dataset["file_info"],
            "common_gene_count": dataset["common_gene_count"],
        },
        "columns": list(counts.columns),
        "gene_col": "gene",
        "sample_cols": dataset["sample_cols"],
        "sample_meta": dataset["sample_meta"],
        "groups": dataset["groups"],
        "gene_types": gene_types,
        "has_gene_name": "gene_name" in counts.columns,
        "has_gene_type": "gene_type" in counts.columns,
        "has_annotation": dataset["term_annotation"] is not None,
        "recipe": DEFAULT_RECIPE,
        "preview": preview.to_dict(orient="records"),
    }


def make_demo_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = random.Random(11)
    genes = [f"Gene{str(i).zfill(4)}" for i in range(1, 701)]
    samples = ["B1_control", "B1_KO1", "B1_KO2", "B2_control", "B2_KO1", "B2_KO2"]
    rows = []
    for index, gene in enumerate(genes):
        base = rng.randint(30, 1600)
        ko1_fold = 1.0
        ko2_fold = 1.0
        if index < 42:
            ko1_fold = rng.uniform(2.2, 5.0)
        elif 42 <= index < 82:
            ko1_fold = rng.uniform(0.18, 0.5)
        if 70 <= index < 118:
            ko2_fold = rng.uniform(2.1, 5.6)
        elif 118 <= index < 160:
            ko2_fold = rng.uniform(0.16, 0.48)

        batch_shift = [1.0, 1.0, 1.0, 1.18, 1.18, 1.18]
        means = [base, base * ko1_fold, base * ko2_fold, base, base * ko1_fold, base * ko2_fold]
        counts = [max(0, int(rng.gauss(mean * shift, mean * shift * 0.18))) for mean, shift in zip(means, batch_shift)]
        gene_type = "protein_coding" if index % 9 != 0 else "lncRNA"
        rows.append(
            {
                "gene": gene,
                "gene_name": f"GN{index + 1}",
                "gene_type": gene_type,
                **dict(zip(samples, counts)),
            }
        )
    counts = pd.DataFrame(rows)

    terms = []
    term_defs = [
        ("GO", "GO:0006955", "immune response", genes[:70]),
        ("GO", "GO:0008283", "cell proliferation", genes[35:130]),
        ("GO", "GO:0006915", "apoptotic process", genes[80:160]),
        ("KEGG", "hsa04630", "JAK-STAT signaling pathway", genes[:55] + genes[300:340]),
        ("KEGG", "hsa04110", "Cell cycle", genes[50:120] + genes[360:390]),
        ("KEGG", "hsa04010", "MAPK signaling pathway", genes[120:210]),
    ]
    for category, term_id, name, gene_list in term_defs:
        for gene in gene_list:
            terms.append({"gene": gene, "term_id": term_id, "term_name": name, "category": category})
    return counts, pd.DataFrame(terms)


def demo_payload() -> dict[str, Any]:
    counts, _annotations = make_demo_data()
    recipe = parse_recipe({})
    sample_meta = [
        {"sample": "B1_control", "original_sample": "control", "batch": "batch1", "group": "control", "file": "demo_batch1.xlsx"},
        {"sample": "B1_KO1", "original_sample": "KO1", "batch": "batch1", "group": "KO1", "file": "demo_batch1.xlsx"},
        {"sample": "B1_KO2", "original_sample": "KO2", "batch": "batch1", "group": "KO2", "file": "demo_batch1.xlsx"},
        {"sample": "B2_control", "original_sample": "control", "batch": "batch2", "group": "control", "file": "demo_batch2.xlsx"},
        {"sample": "B2_KO1", "original_sample": "KO1", "batch": "batch2", "group": "KO1", "file": "demo_batch2.xlsx"},
        {"sample": "B2_KO2", "original_sample": "KO2", "batch": "batch2", "group": "KO2", "file": "demo_batch2.xlsx"},
    ]
    sample_cols = [row["sample"] for row in sample_meta]
    control_group = "control"
    ko_groups = ["KO1", "KO2"]
    payload = run_r_deseq_pipeline(
        counts=counts,
        sample_cols=sample_cols,
        sample_meta=sample_meta,
        control_group=control_group,
        ko_groups=ko_groups,
        recipe=recipe,
        species="human",
    )
    payload["info"] = {
        "files": [
            {"filename": "demo_batch1.xlsx", "batch": "batch1", "samples": sample_cols[:3]},
            {"filename": "demo_batch2.xlsx", "batch": "batch2", "samples": sample_cols[3:]},
        ],
        "common_gene_count": int(len(counts)),
    }
    payload["params"] = {
        "control_group": control_group,
        "ko_groups": ko_groups,
        "species": "human",
        "fc_cutoff": recipe["deg"]["log2fc_cutoff"],
        "padj_cutoff": recipe["deg"]["padj_cutoff"],
    }
    return payload


class AppHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, directory=str(APP_DIR), **kwargs)

    def log_message(self, format: str, *args: Any) -> None:
        print(format % args)

    def do_GET(self) -> None:
        if self.path == "/api/demo":
            try:
                json_response(self, demo_payload())
            except Exception as exc:  # pragma: no cover - user-facing guard
                error_response(self, str(exc), 500)
            return
        if self.path == "/":
            self.path = "/index.html"
        super().do_GET()

    def do_POST(self) -> None:
        if self.path == "/api/export_deseq":
            try:
                length = int(self.headers.get("content-length", "0"))
                payload = json.loads(self.rfile.read(length).decode("utf-8") or "{}")
                workbook = export_deseq_workbook(payload)
                binary_response(
                    self,
                    workbook,
                    "DESeq2_all_results.xlsx",
                    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                )
            except Exception as exc:
                error_response(self, str(exc), 400)
            return

        if self.path not in {"/api/inspect", "/api/analyze"}:
            error_response(self, "Unknown endpoint.", 404)
            return
        content_type = self.headers.get("content-type", "")
        if "multipart/form-data" not in content_type:
            error_response(self, "请使用表单上传文件。")
            return
        length = int(self.headers.get("content-length", "0") or "0")
        max_upload_mb = int(os.environ.get("MAX_UPLOAD_MB", "200"))
        max_upload_bytes = max_upload_mb * 1024 * 1024
        if length <= 0:
            error_response(self, "没有收到上传内容。")
            return
        if length > max_upload_bytes:
            error_response(self, f"上传内容超过 {max_upload_mb} MB，请减少文件数量或调高服务器 MAX_UPLOAD_MB。", 413)
            return
        body = self.rfile.read(length)
        try:
            file_items, fields = parse_multipart_form(self.headers, body)
        except Exception as exc:
            error_response(self, f"上传表单解析失败：{exc}", 400)
            return
        if not file_items:
            error_response(self, "请选择一个或多个数据文件。")
            return
        try:
            params = json.loads(fields.get("params") or "{}")
            if self.path == "/api/inspect":
                payload = inspect_payload(file_items, params)
            else:
                payload = analyze_payload(file_items, params)
            json_response(self, payload)
        except Exception as exc:
            error_response(self, str(exc), 400)


def choose_port(preferred: int) -> int:
    for port in range(preferred, preferred + 50):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            if sock.connect_ex(("127.0.0.1", port)) != 0:
                return port
    return preferred


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default=os.environ.get("HOST", "127.0.0.1"))
    parser.add_argument("--port", type=int, default=int(os.environ.get("PORT", "8787")))
    args = parser.parse_args()
    port = choose_port(args.port)
    server = ThreadingHTTPServer((args.host, port), AppHandler)
    if args.host in {"0.0.0.0", "::"}:
        print(f"RNA-seq web app running locally at http://127.0.0.1:{port}")
        for address in local_network_addresses(port):
            print(f"LAN access: {address}")
    else:
        print(f"RNA-seq web app running at http://{args.host}:{port}")
        print("For LAN access, start with: HOST=0.0.0.0 PORT=%s python3 server.py" % port)
    server.serve_forever()


if __name__ == "__main__":
    main()
