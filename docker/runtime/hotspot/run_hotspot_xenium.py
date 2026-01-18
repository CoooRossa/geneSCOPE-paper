#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Uses the official Hotspot implementation to compute local-correlation p-values and FDR.

Tutorial reference (Hotspot Spatial Tutorial / Slide-seq demo):
https://hotspot.readthedocs.io/en/latest/Spatial_Tutorial.html
"""
import argparse
import hashlib
import json
import os
import re
import sys
import time
from datetime import datetime
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def parse_bool(val: Optional[str], default: bool = False) -> bool:
    if val is None:
        return default
    if isinstance(val, bool):
        return val
    s = str(val).strip().lower()
    if s in {"1", "true", "t", "yes", "y"}:
        return True
    if s in {"0", "false", "f", "no", "n"}:
        return False
    return default


HOTSPOT_LOCALCORR_PVALUE_METHOD = "hotspot_localcorr_z_official_p_fdr"

DEFAULT_EXCLUDE_FEATURE_LABELS = (
    "Unassigned,NegControl,Background,DeprecatedCodeword,SystemControl,Negative"
)

FEATURE_LABEL_CANDIDATE_COLS: Tuple[str, ...] = (
    "feature_type",
    "feature_types",
    "feature_category",
    "feature_class",
    "probe_type",
    "probe_class",
    "target_type",
    "target_category",
    "category",
    "type",
)

FEATURE_LABEL_FALLBACK_DELIMS: Tuple[str, ...] = ("_", "-", ":", ".", "|", " ")

DEFAULT_EXCLUDE_PREFIX: Tuple[str, ...] = (
    "negcontrol",
    "unassigned",
    "blank",
    "control",
    "negatives",
    "ambiguous",
    "background",
    "deprecatedcodeword",
    "systemcontrol",
    "negative",
)


def _hotspot_official_localcorr_p_fdr(
    z_df: pd.DataFrame,
    genes: Sequence[str],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Re-use Hotspot's module logic (norm.sf + multipletests) to turn pair-wise Z
    scores into p-values and BH-FDR values so we keep the official implementation.
    """
    gene_list = [str(x) for x in genes]
    if len(gene_list) == 0:
        return (
            np.array([], dtype=float),
            np.array([], dtype=float),
            np.empty((0, 0), dtype=float),
        )

    try:
        from hotspot import modules as hotspot_modules
    except Exception as exc:
        raise RuntimeError(
            "Cannot import hotspot.modules for official local-correlation p/FDR calculation"
        ) from exc

    try:
        norm = hotspot_modules.norm
        multipletests = hotspot_modules.multipletests
    except AttributeError as exc:
        raise RuntimeError(
            "Hotspot modules missing norm/multipletests; upgrade the hotspot package"
        ) from exc

    subset = z_df.reindex(index=gene_list, columns=gene_list)
    mat = subset.to_numpy(dtype=float)
    sym = (mat + mat.T) / 2.0
    np.fill_diagonal(sym, 0.0)

    try:
        from scipy.spatial.distance import squareform
    except Exception as exc:
        raise RuntimeError(
            "Missing scipy.spatial.distance.squareform needed for official p/FDR computation"
        ) from exc

    try:
        vec = squareform(sym)
    except Exception as exc:
        raise RuntimeError(
            "Failed to vectorize Hotspot local correlation Z matrix"
        ) from exc

    p_vec = np.full(vec.shape, np.nan, dtype=float)
    fdr_vec = np.full(vec.shape, np.nan, dtype=float)
    finite = np.isfinite(vec)
    if finite.any():
        p_vals = norm.sf(vec[finite])
        fdr_vals = multipletests(p_vals, method="fdr_bh")[1]
        p_vec[finite] = p_vals
        fdr_vec[finite] = fdr_vals

    n = sym.shape[0]
    tri = np.triu_indices(n, k=1)
    fdr_mat = np.full(sym.shape, np.nan, dtype=float)
    fdr_mat[tri] = fdr_vec
    fdr_mat[(tri[1], tri[0])] = fdr_vec
    np.fill_diagonal(fdr_mat, 0.0)

    return p_vec, fdr_vec, fdr_mat


def _write_empty_gene_graph_edges(outdir: str) -> None:
    cols = ["gene_a", "gene_b", "weight", "p_value", "fdr", "n_obs", "edge_type"]
    pd.DataFrame(columns=cols).to_csv(
        os.path.join(outdir, "gene_graph_edges.tsv"), sep="\t", index=False
    )


def emit_hotspot_gene_graph_edges(
    local_corr: Optional[object],
    local_corr_z: Optional[object],
    local_corr_c: Optional[object],
    outdir: str,
    n_obs: int,
    fdr_threshold: float,
) -> Dict[str, object]:
    edge_type = "hotspot_local_corr"
    edges_path = os.path.join(outdir, "gene_graph_edges.tsv")

    def _as_df(obj: Optional[object]) -> Optional[pd.DataFrame]:
        if obj is None:
            return None
        if isinstance(obj, pd.DataFrame):
            return obj.copy()
        try:
            return pd.DataFrame(obj)
        except Exception:
            return None

    z_df = _as_df(local_corr_z)
    c_df = _as_df(local_corr_c)
    w_df = c_df if c_df is not None else _as_df(local_corr)

    if w_df is None and z_df is None:
        _write_empty_gene_graph_edges(outdir)
        return {
            "edges_total_n": 0,
            "edges_sig_n": 0,
            "edges_weight_type": "local_corr",
            "edges_pvalue_method": HOTSPOT_LOCALCORR_PVALUE_METHOD,
            "edges_fdr_threshold": float(fdr_threshold),
        }

    if w_df is None and z_df is not None:
        w_df = z_df

    if w_df is None:
        _write_empty_gene_graph_edges(outdir)
        return {
            "edges_total_n": 0,
            "edges_sig_n": 0,
            "edges_weight_type": "local_corr",
            "edges_pvalue_method": HOTSPOT_LOCALCORR_PVALUE_METHOD,
            "edges_fdr_threshold": float(fdr_threshold),
        }

    # Prefer the square matrix representation (Hotspot v0.9.1) to recover the
    # module-building gene–gene graph from local correlations.
    if w_df.shape[0] == w_df.shape[1]:
        genes_idx = w_df.index.astype(str)
        genes_col = w_df.columns.astype(str)
        if set(genes_idx) != set(genes_col):
            raise ValueError("Hotspot local correlation matrix missing gene names.")

        w_df = w_df.copy()
        w_df.index = genes_idx
        w_df.columns = genes_col
        if list(w_df.columns) != list(w_df.index):
            w_df = w_df.reindex(index=genes_idx, columns=genes_idx)

        if z_df is None:
            raise RuntimeError(
                "Hotspot local_correlation_z missing; cannot compute official p/FDR."
            )
        z_df = z_df.copy()
        z_df.index = z_df.index.astype(str)
        z_df.columns = z_df.columns.astype(str)
        z_df = z_df.reindex(index=genes_idx, columns=genes_idx)

        genes = genes_idx.to_numpy()
        w_mat = w_df.apply(pd.to_numeric, errors="coerce").to_numpy()
        p_vals, fdr_vals, _ = _hotspot_official_localcorr_p_fdr(z_df, genes)
        pvalue_method = HOTSPOT_LOCALCORR_PVALUE_METHOD

        tri = np.triu_indices(len(genes), k=1)
        w_vals = w_mat[tri]
        mask = np.isfinite(w_vals) & np.isfinite(p_vals) & np.isfinite(fdr_vals)

        edges_total_n = int(np.sum(mask))
        if edges_total_n == 0:
            _write_empty_gene_graph_edges(outdir)
            return {
                "edges_total_n": 0,
                "edges_sig_n": 0,
                "edges_weight_type": "local_corr",
                "edges_pvalue_method": pvalue_method,
                "edges_fdr_threshold": float(fdr_threshold),
            }

        w_vals = w_vals[mask]
        p_vals = p_vals[mask]
        fdr_vals = fdr_vals[mask]
        gene_a = genes[tri[0]][mask]
        gene_b = genes[tri[1]][mask]

        edges_df = pd.DataFrame(
            {
                "gene_a": gene_a.astype(str),
                "gene_b": gene_b.astype(str),
                "weight": w_vals.astype(float),
                "p_value": p_vals.astype(float),
                "fdr": fdr_vals.astype(float),
                "n_obs": int(n_obs),
                "edge_type": edge_type,
            }
        )

        # The Hotspot module graph is based on positively correlated gene pairs.
        edges_sig = edges_df[
            (edges_df["weight"] > 0)
            & np.isfinite(edges_df["fdr"])
            & (edges_df["fdr"] <= float(fdr_threshold))
        ].copy()

        edges_sig = edges_sig.sort_values(["fdr", "weight"], ascending=[True, False])
        edges_sig.to_csv(edges_path, sep="\t", index=False)

        return {
            "edges_total_n": int(edges_total_n),
            "edges_sig_n": int(edges_sig.shape[0]),
            "edges_weight_type": "local_corr" if c_df is not None else "local_corr_or_z",
            "edges_pvalue_method": pvalue_method,
            "edges_fdr_threshold": float(fdr_threshold),
        }

    # Fallback: unsupported representation (avoid emitting a dense g×g dump).
    _write_empty_gene_graph_edges(outdir)
    return {
        "edges_total_n": 0,
        "edges_sig_n": 0,
        "edges_weight_type": "local_corr",
        "edges_pvalue_method": HOTSPOT_LOCALCORR_PVALUE_METHOD,
        "edges_fdr_threshold": float(fdr_threshold),
    }


def file_md5(path: str) -> str:
    if not path or not os.path.exists(path):
        return ""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def file_sha256(path: str) -> str:
    if not path or not os.path.exists(path):
        return ""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_obs_ids(obs_ids: Iterable[str]) -> str:
    ids = sorted(str(x) for x in obs_ids)
    if not ids:
        return ""
    payload = "\n".join(ids)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def compute_roi_hash(roi_source_sha256: str, obs_id_hash: str) -> str:
    if not obs_id_hash:
        return ""
    payload = f"{roi_source_sha256}\n{obs_id_hash}"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def normalize_dataset_id(dataset_id: str, data_dir: str) -> str:
    if dataset_id:
        return dataset_id
    return os.path.basename(os.path.abspath(data_dir))


def normalize_roi_id(roi_id: str, roi_csv: str) -> str:
    if roi_id:
        return roi_id
    if roi_csv:
        return os.path.splitext(os.path.basename(roi_csv))[0]
    return "full"


def read_resource_stats(path: str) -> Dict[str, object]:
    resources: Dict[str, object] = {"stats_path": path}
    if not path or not os.path.exists(path):
        return resources
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return resources
    mem_cur_col = None
    for cand in ("mem_current_bytes", "memory_current_bytes", "mem_bytes"):
        if cand in df.columns:
            mem_cur_col = cand
            break
    mem_peak_col = None
    for cand in ("mem_peak_bytes", "memory_peak_bytes", "mem_max_bytes", "memory_max_bytes"):
        if cand in df.columns:
            mem_peak_col = cand
            break
    if mem_cur_col is not None:
        resources["max_mem_current_bytes"] = float(pd.to_numeric(df[mem_cur_col], errors="coerce").max())
    if mem_peak_col is not None:
        resources["max_mem_peak_bytes"] = float(pd.to_numeric(df[mem_peak_col], errors="coerce").max())
    elif mem_cur_col is not None and "max_mem_current_bytes" in resources:
        resources["max_mem_peak_bytes"] = resources["max_mem_current_bytes"]
    if "cpu_usage_usec" in df.columns:
        cpu_vals = pd.to_numeric(df["cpu_usage_usec"], errors="coerce").dropna()
        if not cpu_vals.empty:
            resources["cpu_usage_usec"] = float(cpu_vals.iloc[-1])
    return resources


def read_roi_polygon(csv_path: str) -> np.ndarray:
    df = pd.read_csv(csv_path, comment="#")
    if df.empty:
        raise ValueError(f"ROI CSV has no data rows after skipping comments. File={csv_path}")

    cols = {c.strip().lower(): c for c in df.columns}
    xcol = cols.get("x")
    ycol = cols.get("y")
    if xcol is None or ycol is None:
        if df.shape[1] < 2:
            raise ValueError(
                "ROI CSV must have at least 2 columns after skipping comments. "
                f"File={csv_path} Columns={list(df.columns)}"
            )
        xcol, ycol = df.columns[:2]

    xy = df[[xcol, ycol]].copy()
    xy[xcol] = pd.to_numeric(xy[xcol], errors="coerce")
    xy[ycol] = pd.to_numeric(xy[ycol], errors="coerce")
    xy = xy.dropna()

    pts = xy.to_numpy(dtype=float)

    if pts.shape[0] >= 2:
        keep = [True]
        for i in range(1, pts.shape[0]):
            keep.append(not np.allclose(pts[i], pts[i - 1]))
        pts = pts[np.array(keep, dtype=bool)]

    if pts.shape[0] >= 2 and np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]

    if pts.shape[0] < 3:
        raise ValueError(
            "ROI polygon needs >=3 vertices after cleaning. "
            f"File={csv_path} Columns={list(df.columns)}"
        )

    if not np.allclose(pts[0], pts[-1]):
        pts = np.vstack([pts, pts[0]])

    return pts


def clip_cells_to_roi(cells_df: pd.DataFrame, roi_vertices: np.ndarray) -> pd.DataFrame:
    from shapely.geometry import Polygon

    poly = Polygon(roi_vertices)
    xs = cells_df["x"].to_numpy()
    ys = cells_df["y"].to_numpy()

    mask = None
    try:
        import shapely
        if hasattr(shapely, "points") and hasattr(shapely, "intersects"):
            pts = shapely.points(xs, ys)
            mask = shapely.intersects(poly, pts)
        elif hasattr(shapely, "vectorized"):
            from shapely import vectorized
            mask = vectorized.contains(poly, xs, ys) | vectorized.touches(poly, xs, ys)
    except Exception:
        mask = None

    if mask is None:
        from shapely.geometry import Point
        mask = np.fromiter(
            (poly.intersects(Point(float(x), float(y))) for x, y in zip(xs, ys)),
            dtype=bool,
            count=len(xs),
        )

    return cells_df.loc[mask].copy()


def _infer_cell_id_col(df: pd.DataFrame) -> str:
    lc = {c.lower(): c for c in df.columns}
    for cand in ["cell_id", "barcode", "id"]:
        if cand in lc:
            return lc[cand]
    raise ValueError(
        "cells.parquet: cannot find a cell id column. "
        f"cols(head)={list(df.columns)[:80]}"
    )


def _infer_xy_cols(df: pd.DataFrame) -> Tuple[str, str]:
    lc = {c.lower(): c for c in df.columns}
    x_candidates = ["x_centroid", "x", "center_x", "x_um", "x_location"]
    y_candidates = ["y_centroid", "y", "center_y", "y_um", "y_location"]
    xcol = next((lc.get(c) for c in x_candidates if c in lc), None)
    ycol = next((lc.get(c) for c in y_candidates if c in lc), None)
    if xcol is None or ycol is None:
        raise ValueError(
            "cells.parquet: cannot infer x/y columns. "
            f"Need one of {x_candidates} and {y_candidates}. cols(head)={list(df.columns)[:80]}"
        )
    return xcol, ycol


def _iter_arrow_batches(dataset, tx_pq_path: str, columns: List[str], batch_size: int) -> Iterable:
    if hasattr(dataset, "scanner"):
        scanner = dataset.scanner(columns=columns, batch_size=batch_size)
        for batch in scanner.to_batches():
            yield batch
        return

    if hasattr(dataset, "to_batches"):
        for batch in dataset.to_batches(columns=columns, batch_size=batch_size):
            yield batch
        return

    import pyarrow.parquet as pq
    pf = pq.ParquetFile(tx_pq_path)
    for batch in pf.iter_batches(batch_size=batch_size, columns=columns):
        yield batch


def load_cells_df(cells_pq: str) -> pd.DataFrame:
    import pyarrow.parquet as pq

    cdf = pq.read_table(cells_pq).to_pandas()
    idcol = _infer_cell_id_col(cdf)
    xcol, ycol = _infer_xy_cols(cdf)

    cells_df = cdf[[idcol, xcol, ycol]].copy()
    cells_df.columns = ["cell_id", "x", "y"]
    cells_df["cell_id"] = cells_df["cell_id"].astype(str)
    return cells_df


def apply_roi_and_downsample(
    cells_df: pd.DataFrame,
    roi_csv: str,
    max_cells: int,
    seed: int,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    stats = {
        "n_cells_input": int(len(cells_df)),
        "n_cells_after_roi": int(len(cells_df)),
        "n_cells_used": int(len(cells_df)),
        "roi_applied": 0,
        "downsample_applied": 0,
    }
    if roi_csv:
        roi = read_roi_polygon(roi_csv)
        before = len(cells_df)
        cells_df = clip_cells_to_roi(cells_df, roi)
        stats["roi_applied"] = 1
        stats["n_cells_after_roi"] = int(len(cells_df))
        eprint(f"[INFO] ROI clip on cells: {before} -> {len(cells_df)}")

    if max_cells > 0 and len(cells_df) > max_cells:
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(cells_df), size=max_cells, replace=False)
        cells_df = cells_df.iloc[idx].copy()
        stats["downsample_applied"] = 1
        eprint(f"[INFO] Downsample cells to max_cells={max_cells}: n={len(cells_df)}")

    stats["n_cells_used"] = int(len(cells_df))
    stats["obs_id_hash"] = hash_obs_ids(cells_df["cell_id"].astype(str).tolist())
    return cells_df, stats


def build_adata_from_transcripts(
    cells_df: pd.DataFrame,
    tx_pq: str,
    exclude_prefix: Tuple[str, ...],
) -> "anndata.AnnData":
    import pyarrow.dataset as ds
    from scipy import sparse
    import anndata as ad

    selected_cell_ids = set(cells_df["cell_id"].astype(str).tolist())
    if not selected_cell_ids:
        raise ValueError("No cells left after ROI/downsample. Check ROI and parameters.")

    dataset = ds.dataset(tx_pq, format="parquet")
    schema_names = [f.name for f in dataset.schema]
    lmap = {c.lower(): c for c in schema_names}

    cell_id_col = lmap.get("cell_id") or lmap.get("barcode")
    gene_col = lmap.get("feature_name") or lmap.get("gene") or lmap.get("name")

    if cell_id_col is None or gene_col is None:
        raise ValueError(
            "transcripts.parquet: cannot infer cell_id/gene columns. "
            f"schema={schema_names}"
        )

    scanned = 0
    kept = 0
    counts: Dict[Tuple[str, str], int] = {}

    batch_size = 1_000_000
    columns = [cell_id_col, gene_col]

    for batch in _iter_arrow_batches(dataset, tx_pq, columns=columns, batch_size=batch_size):
        df = batch.to_pandas()
        scanned += len(df)
        df[cell_id_col] = df[cell_id_col].astype(str)
        df = df[df[cell_id_col].isin(selected_cell_ids)]
        if df.empty:
            continue

        g = df[gene_col].astype(str)
        g_lower = g.str.lower()
        mask = np.ones(len(df), dtype=bool)
        for p in exclude_prefix:
            mask &= ~g_lower.str.startswith(p)
        df = df.loc[mask]
        if df.empty:
            continue

        kept += len(df)
        grp = df.groupby([cell_id_col, gene_col]).size()
        for (cid, gene), val in grp.items():
            key = (str(cid), str(gene))
            counts[key] = counts.get(key, 0) + int(val)

    eprint(
        f"[INFO] transcripts aggregation: scanned_rows={scanned} kept_rows={kept} cells={len(selected_cell_ids)}"
    )

    if not counts:
        raise ValueError("No transcript counts after filtering. Check ROI/downsample.")

    cell_ids = cells_df["cell_id"].astype(str).tolist()
    genes = sorted({gene for _, gene in counts.keys()})

    cell_index = {c: i for i, c in enumerate(cell_ids)}
    gene_index = {g: j for j, g in enumerate(genes)}

    rows: List[int] = []
    cols: List[int] = []
    data: List[int] = []
    for (cid, gene), val in counts.items():
        rows.append(cell_index[cid])
        cols.append(gene_index[gene])
        data.append(val)

    X = sparse.csr_matrix(
        (data, (rows, cols)), shape=(len(cell_ids), len(genes)), dtype=np.float32
    )

    adata = ad.AnnData(X=X)
    adata.obs_names = pd.Index(cell_ids, dtype=str)
    adata.var_names = pd.Index(genes, dtype=str)

    coords = cells_df.set_index("cell_id").loc[adata.obs_names, ["x", "y"]]
    adata.obs["x"] = coords["x"].to_numpy()
    adata.obs["y"] = coords["y"].to_numpy()
    adata.obsm["spatial"] = coords[["x", "y"]].to_numpy()

    return adata


def load_h5ad(path: str) -> "anndata.AnnData":
    import anndata as ad

    return ad.read_h5ad(path)


def ensure_xy(adata: "anndata.AnnData") -> "anndata.AnnData":
    if "x" not in adata.obs.columns or "y" not in adata.obs.columns:
        if "spatial" in adata.obsm:
            adata.obs["x"] = adata.obsm["spatial"][:, 0]
            adata.obs["y"] = adata.obsm["spatial"][:, 1]
        else:
            raise ValueError("No spatial coordinates found. Need obs['x','y'] or obsm['spatial'].")
    adata.obsm["spatial"] = adata.obs[["x", "y"]].to_numpy()
    return adata


def _parse_exclude_feature_labels(raw: Optional[str]) -> List[str]:
    if raw is None:
        return []
    s = str(raw)
    if s.strip() == "":
        return []

    out: List[str] = []
    seen = set()
    for part in s.split(","):
        lab = str(part).strip()
        if not lab:
            continue
        key = lab.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append(lab)
    return out


def _series_contains_any(series: pd.Series, needles_lower: Sequence[str]) -> pd.Series:
    if series.empty:
        return pd.Series([], dtype=bool, index=series.index)
    needles = [x for x in needles_lower if x]
    if not needles:
        return pd.Series(False, index=series.index)
    pattern = "|".join(re.escape(x) for x in needles)
    return series.str.contains(pattern, regex=True)


def filter_features_by_label(
    adata: "anndata.AnnData",
    exclude_feature_labels_raw: Optional[str],
    candidate_cols: Sequence[str] = FEATURE_LABEL_CANDIDATE_COLS,
    fallback_delims: Sequence[str] = FEATURE_LABEL_FALLBACK_DELIMS,
) -> Tuple["anndata.AnnData", Dict[str, object]]:
    n_vars_input = int(adata.n_vars)
    labels = _parse_exclude_feature_labels(exclude_feature_labels_raw)
    labels_lower = [x.strip().lower() for x in labels]
    labels_set = {x for x in labels_lower if x}

    meta: Dict[str, object] = {
        "exclude_feature_labels": labels,
        "exclude_label_cols": [],
        "fallback_used": False,
        "n_vars_dropped": 0,
    }

    if not labels_set or n_vars_input == 0:
        eprint(
            "[INFO] Feature label filtering (substring match): "
            f"n_vars {n_vars_input} -> {n_vars_input} (dropped 0); "
            f"exclude_label_cols={meta['exclude_label_cols']}; "
            f"exclude_feature_labels={meta['exclude_feature_labels']}"
        )
        return adata, meta

    drop_mask = np.zeros(n_vars_input, dtype=bool)
    used_cols: List[str] = []

    for col in candidate_cols:
        if col not in adata.var.columns:
            continue
        s = adata.var[col]
        norm = s.astype(str).str.strip().str.lower()
        hit = _series_contains_any(norm, sorted(labels_set))
        if bool(hit.any()):
            used_cols.append(col)
            drop_mask |= hit.to_numpy(dtype=bool)

    fallback_used = False
    if not bool(drop_mask.any()):
        fallback_used = True
        names_lower = pd.Series(adata.var_names.astype(str)).astype(str).str.strip().str.lower()
        hit = _series_contains_any(names_lower, sorted(labels_set))
        drop_mask = hit.to_numpy(dtype=bool)
        used_cols = ["__var_names__"]

    n_vars_dropped = int(np.sum(drop_mask))
    if n_vars_dropped > 0:
        keep = ~drop_mask
        adata = adata[:, keep].copy()

    n_vars_used = int(adata.n_vars)
    meta = {
        "exclude_feature_labels": labels,
        "exclude_label_cols": used_cols,
        "fallback_used": bool(fallback_used),
        "n_vars_dropped": int(n_vars_input - n_vars_used),
    }

    eprint(
        "[INFO] Feature label filtering (substring match): "
        f"n_vars {n_vars_input} -> {n_vars_used} (dropped {meta['n_vars_dropped']}); "
        f"exclude_label_cols={meta['exclude_label_cols']}; "
        f"exclude_feature_labels={meta['exclude_feature_labels']}"
    )
    return adata, meta


def apply_roi_and_downsample_adata(
    adata: "anndata.AnnData",
    roi_csv: str,
    max_cells: int,
    seed: int,
) -> Tuple["anndata.AnnData", Dict[str, int]]:
    stats = {
        "n_cells_input": int(adata.n_obs),
        "n_cells_after_roi": int(adata.n_obs),
        "n_cells_used": int(adata.n_obs),
        "roi_applied": 0,
        "downsample_applied": 0,
    }
    if roi_csv:
        roi = read_roi_polygon(roi_csv)
        coords = adata.obsm.get("spatial")
        if coords is None:
            raise ValueError("No spatial coordinates for ROI filtering.")

        cells_df = pd.DataFrame({
            "cell_id": adata.obs_names.astype(str),
            "x": coords[:, 0],
            "y": coords[:, 1],
        })
        before = adata.n_obs
        cells_df = clip_cells_to_roi(cells_df, roi)
        keep = adata.obs_names.astype(str).isin(cells_df["cell_id"].astype(str))
        adata = adata[keep].copy()
        stats["roi_applied"] = 1
        stats["n_cells_after_roi"] = int(adata.n_obs)
        eprint(f"[INFO] ROI clip on cells: {before} -> {adata.n_obs}")

    if max_cells > 0 and adata.n_obs > max_cells:
        rng = np.random.default_rng(seed)
        idx = rng.choice(adata.n_obs, size=max_cells, replace=False)
        adata = adata[idx].copy()
        stats["downsample_applied"] = 1

    stats["n_cells_used"] = int(adata.n_obs)
    stats["obs_id_hash"] = hash_obs_ids(adata.obs_names.astype(str).tolist())
    return adata, stats


def prepare_adata_for_hotspot(adata: "anndata.AnnData") -> "anndata.AnnData":
    from scipy import sparse

    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    elif not sparse.isspmatrix_csr(adata.X):
        adata.X = adata.X.tocsr()

    adata.var_names_make_unique()
    adata.var_names = adata.var_names.astype(str)

    adata.layers["csc_counts"] = adata.X.tocsc()
    total_counts = np.asarray(adata.X.sum(axis=1)).ravel()
    # Hotspot's bernoulli model expects integer UMI counts for binning.
    if np.issubdtype(total_counts.dtype, np.floating):
        rounded = np.rint(total_counts)
        if not np.allclose(total_counts, rounded):
            eprint("[WARN] total_counts are not integers; rounding for Hotspot umi_counts.")
        total_counts = rounded.astype(np.int64, copy=False)
    else:
        total_counts = total_counts.astype(np.int64, copy=False)
    zero_mask = total_counts <= 0
    if np.any(zero_mask):
        eprint(
            f"[WARN] total_counts has {int(np.sum(zero_mask))} zeros; setting to 1 to avoid log10(0)."
        )
        total_counts = total_counts.copy()
        total_counts[zero_mask] = 1
    adata.obs["total_counts"] = total_counts
    return adata


def _call_jobs(func, ncores: int, *args, **kwargs):
    try:
        return func(*args, jobs=ncores, **kwargs)
    except TypeError:
        return func(*args, n_jobs=ncores, **kwargs)


def _call_create_modules(hs, min_gene_threshold: int, core_only: bool, fdr_threshold: float):
    try:
        return hs.create_modules(
            min_gene_threshold=min_gene_threshold,
            core_only=core_only,
            fdr_threshold=fdr_threshold,
        )
    except TypeError:
        return hs.create_modules(min_gene_threshold=min_gene_threshold, core_only=core_only)


def normalize_modules(modules_obj) -> Optional[pd.DataFrame]:
    if modules_obj is None:
        return None
    if isinstance(modules_obj, pd.Series):
        df = modules_obj.to_frame(name="module_id")
    elif isinstance(modules_obj, pd.DataFrame):
        df = modules_obj.copy()
    else:
        try:
            df = pd.DataFrame(modules_obj)
        except Exception:
            return None

    if df.empty:
        return None

    if "gene" not in df.columns:
        df = df.reset_index().rename(columns={"index": "gene"})

    module_col = None
    for cand in ["module_id", "module", "Module", "ModuleID", "moduleID"]:
        if cand in df.columns:
            module_col = cand
            break
    if module_col is None:
        for col in df.columns:
            if col != "gene":
                module_col = col
                break
    if module_col is None:
        return None

    out = df[["gene", module_col]].rename(columns={module_col: "module_id"})
    return out


def build_modules_table(modules_obj, all_genes: List[str]) -> pd.DataFrame:
    mod_df = normalize_modules(modules_obj)
    if mod_df is None or mod_df.empty:
        return pd.DataFrame({"gene": all_genes, "module_id": -1})

    mod_df = mod_df.copy()
    mod_df["gene"] = mod_df["gene"].astype(str)
    mod_df = mod_df.drop_duplicates(subset=["gene"], keep="first")

    all_df = pd.DataFrame({"gene": all_genes})
    merged = all_df.merge(mod_df, on="gene", how="left")

    mod_vals = pd.to_numeric(merged["module_id"], errors="coerce")
    if mod_vals.notna().any():
        merged["module_id"] = mod_vals.fillna(-1).astype(int)
    else:
        merged["module_id"] = merged["module_id"].fillna(-1)

    return merged


def force_k_modules_from_linkage(
    linkage_matrix: object,
    leaf_labels: Iterable[str],
    k_modules: int,
) -> Optional[pd.Series]:
    try:
        from scipy.cluster.hierarchy import fcluster
    except Exception:
        return None

    try:
        Z = np.asarray(linkage_matrix)
    except Exception:
        return None

    if Z.ndim != 2 or Z.shape[1] < 4:
        return None

    labels = pd.Index([str(x) for x in leaf_labels], dtype=str)
    n = int(labels.size)
    if n <= 0:
        return pd.Series([], dtype=int)

    k = int(k_modules)
    if k < 1:
        return None
    if k > n:
        k = n

    try:
        cl = fcluster(Z, t=k, criterion="maxclust")
    except Exception:
        return None

    cl = np.asarray(cl).astype(int)
    if cl.shape[0] != n:
        return None

    # Re-label clusters in a stable, size-desc order (ties by min gene name).
    df = pd.DataFrame({"gene": labels.to_numpy(), "cl": cl})
    sizes = df.groupby("cl")["gene"].agg(["size", "min"]).reset_index()
    sizes = sizes.sort_values(["size", "min", "cl"], ascending=[False, True, True])
    remap: Dict[int, int] = {}
    for i, row in enumerate(sizes.itertuples(index=False)):
        remap[int(row.cl)] = int(i + 1)
    cl2 = np.array([remap[int(x)] for x in cl], dtype=int)
    return pd.Series(cl2, index=labels, name="module_id")


def _write_matrixmarket(matrix: np.ndarray, out_path: str) -> None:
    from scipy.io import mmwrite

    mmwrite(out_path, matrix)


def _emit_gg_matrices_from_hotspot(
    hs: object,
    genes: List[str],
    outdir: str,
) -> Dict[str, object]:
    # Weight matrix: prefer correlation coefficients ("c"), fallback to z-scores.
    weight_df = getattr(hs, "local_correlation_c", None)
    weight_type = "local_correlation_c"
    if weight_df is None:
        weight_df = getattr(hs, "local_correlation_z", None)
        weight_type = "local_correlation_z"

    if weight_df is None:
        raise RuntimeError("Hotspot did not expose local correlation matrices on the object")

    gg_genes = [str(x) for x in list(getattr(weight_df, "index", genes))]
    if len(gg_genes) == 0:
        gg_genes = [str(x) for x in genes]

    gg_genes_path = os.path.join(outdir, "gg_genes.tsv")
    pd.DataFrame({"idx": np.arange(1, len(gg_genes) + 1), "gene": gg_genes}).to_csv(
        gg_genes_path, sep="\t", index=False
    )

    weight_mat = np.asarray(weight_df.loc[gg_genes, gg_genes].to_numpy())
    _write_matrixmarket(weight_mat, os.path.join(outdir, "gg_weight_matrix.mtx"))

    z_df = getattr(hs, "local_correlation_z", None)
    if z_df is None:
        raise RuntimeError(
            "Hotspot did not expose local_correlation_z for official p/FDR conversion"
        )

    _, _, fdr_mat = _hotspot_official_localcorr_p_fdr(
        pd.DataFrame(z_df), gg_genes
    )

    _write_matrixmarket(fdr_mat, os.path.join(outdir, "gg_fdr_matrix.mtx"))

    return {
        "gg_weight_type": weight_type,
        "gg_fdr_source": HOTSPOT_LOCALCORR_PVALUE_METHOD,
        "gg_genes_n": int(len(gg_genes)),
    }


def _init_hotspot(
    adata: "anndata.AnnData",
    model: str,
) -> "hotspot.Hotspot":
    import hotspot
    import inspect
    from scipy import sparse

    sig = inspect.signature(hotspot.Hotspot)
    if "layer_key" in sig.parameters:
        return hotspot.Hotspot(
            adata,
            layer_key="csc_counts",
            model=model,
            latent_obsm_key="spatial",
            umi_counts_obs_key="total_counts",
        )

    counts = adata.layers.get("csc_counts", adata.X)
    if sparse.issparse(counts):
        counts_t = counts.T.tocsr()
        counts_df = pd.DataFrame.sparse.from_spmatrix(
            counts_t,
            index=adata.var_names.astype(str),
            columns=adata.obs_names.astype(str),
        )
    else:
        counts_df = pd.DataFrame(
            counts.T,
            index=adata.var_names.astype(str),
            columns=adata.obs_names.astype(str),
        )

    latent = adata.obsm.get("spatial")
    if latent is None:
        raise ValueError("Hotspot requires spatial coordinates in adata.obsm['spatial']")
    latent_df = pd.DataFrame(
        latent,
        index=adata.obs_names.astype(str),
        columns=[f"latent_{i}" for i in range(latent.shape[1])],
    )
    umi = pd.Series(
        np.asarray(adata.obs["total_counts"]).ravel(),
        index=adata.obs_names.astype(str),
    )
    return hotspot.Hotspot(counts_df, model=model, latent=latent_df, umi_counts=umi)


def run_hotspot(
    adata: "anndata.AnnData",
    n_neighbors: int,
    model: str,
    fdr_autocorr: float,
    min_gene_threshold: int,
    k_modules: int,
    emit_edges: bool,
    edge_fdr: float,
    edge_top_n: int,
    emit_gg_matrices: bool,
    core_only: bool,
    ncores: int,
    outdir: str,
) -> Tuple[int, Dict[str, object]]:
    hs = _init_hotspot(adata, model)
    hs.create_knn_graph(weighted_graph=False, n_neighbors=n_neighbors)

    hs_results = _call_jobs(hs.compute_autocorrelations, ncores)
    hs_results = hs_results.copy()
    fdr_col = None
    for col in hs_results.columns:
        if col.lower() == "fdr":
            fdr_col = col
            break

    if fdr_col is None:
        eprint("[WARN] FDR column not found in autocorrelations; using all genes")
        sig_mask = np.ones(hs_results.shape[0], dtype=bool)
    else:
        sig_mask = hs_results[fdr_col] < fdr_autocorr

    hs_genes = hs_results.index[sig_mask].tolist()
    n_sig = len(hs_genes)

    modules_obj = None
    local_corr = None
    local_corr_z = None
    local_corr_c = None
    edges_meta: Dict[str, object] = {
        "edges_outfile": "gene_graph_edges.tsv",
        "edges_schema": ["gene_a", "gene_b", "weight", "p_value", "fdr", "n_obs", "edge_type"],
        "edges_total_n": 0,
        "edges_sig_n": 0,
        "edges_weight_type": "local_corr",
        "edges_pvalue_method": HOTSPOT_LOCALCORR_PVALUE_METHOD,
        "edges_fdr_threshold": float(edge_fdr),
    }
    gg_meta: Optional[Dict[str, object]] = None
    if n_sig > 0:
        local_corr = _call_jobs(hs.compute_local_correlations, ncores, hs_genes)
        local_corr_z = getattr(hs, "local_correlation_z", None)
        local_corr_c = getattr(hs, "local_correlation_c", None)
        try:
            if emit_edges:
                edges_meta_update = emit_hotspot_gene_graph_edges(
                    local_corr=local_corr,
                    local_corr_z=local_corr_z,
                    local_corr_c=local_corr_c,
                    outdir=outdir,
                    n_obs=int(adata.n_obs),
                    fdr_threshold=float(edge_fdr),
                )
                edges_meta = {**edges_meta, **edges_meta_update}
            else:
                _write_empty_gene_graph_edges(outdir)
        except Exception as exc:
            eprint(f"[WARN] gene_graph_edges.tsv write failed: {exc}")
            _write_empty_gene_graph_edges(outdir)
        modules_obj = _call_create_modules(hs, min_gene_threshold, core_only, fdr_autocorr)
        if emit_gg_matrices:
            gg_meta = _emit_gg_matrices_from_hotspot(hs=hs, genes=hs_genes, outdir=outdir)
        if int(k_modules) > 0:
            forced = None
            try:
                forced = force_k_modules_from_linkage(
                    linkage_matrix=getattr(hs, "linkage", None),
                    leaf_labels=getattr(modules_obj, "index", []),
                    k_modules=int(k_modules),
                )
            except Exception:
                forced = None

            if forced is None or forced.empty:
                eprint("[WARN] k_modules requested but linkage cut failed; using Hotspot default modules")
            else:
                modules_obj = forced
                eprint(
                    "[INFO] k_modules override: requested=%d actual=%d genes=%d",
                    int(k_modules),
                    int(pd.unique(forced.values).size),
                    int(forced.shape[0]),
                )
    else:
        eprint("[WARN] No genes pass FDR threshold; skipping modules")
        _write_empty_gene_graph_edges(outdir)

    all_genes = hs_results.index.astype(str).tolist()
    modules_df = build_modules_table(modules_obj, all_genes)

    modules_df.to_csv(os.path.join(outdir, "modules.tsv"), sep="\t", index=False)
    if gg_meta is not None:
        edges_meta = {**edges_meta, **gg_meta}
    return n_sig, edges_meta


def _hotspot_version() -> str:
    try:
        from importlib.metadata import version

        return version("hotspot")
    except Exception:
        return "unknown"


def _write_meta_json(
    args: argparse.Namespace,
    cell_stats: Dict[str, int],
    n_genes_input: int,
    n_genes_used: int,
    edges_meta: Optional[Dict[str, object]],
    wall_time_sec: float,
    gene_filtering_meta: Optional[Dict[str, object]] = None,
) -> None:
    roi_csv = args.coord_file or args.roi_csv or ""
    dataset_id = normalize_dataset_id(args.dataset_id or "", args.data_dir)
    roi_id = normalize_roi_id(args.roi_id or "", roi_csv)
    roi_source_sha256 = file_sha256(roi_csv) if roi_csv else ""
    obs_id_hash = str(cell_stats.get("obs_id_hash", ""))
    roi_hash = compute_roi_hash(roi_source_sha256, obs_id_hash)

    gene_filtering_meta = gene_filtering_meta or {}
    gene_filtering = {
        "exclude_prefix": list(DEFAULT_EXCLUDE_PREFIX),
        "exclude_feature_labels": list(gene_filtering_meta.get("exclude_feature_labels", [])),
        "exclude_label_cols": list(gene_filtering_meta.get("exclude_label_cols", [])),
        "fallback_used": bool(gene_filtering_meta.get("fallback_used", False)),
        "n_vars_dropped": int(gene_filtering_meta.get("n_vars_dropped", 0)),
        "fdr_autocorr": float(args.fdr_autocorr),
        "min_gene_threshold": int(args.min_gene_threshold),
        "core_only": parse_bool(args.core_only, False),
    }

    roi_applied = bool(int(cell_stats.get("roi_applied", 0)))
    downsample_applied = bool(int(cell_stats.get("downsample_applied", 0)))
    k_modules_enabled = int(args.k_modules) > 0

    meta = {
        "method": "hotspot",
        "method_version": _hotspot_version(),
        "run_id": f"hotspot_{int(time.time())}",
        "timestamp": datetime.now().isoformat(),
        "workflow_vendor": "hotspot",
        "workflow_tutorial": "Demo: Spatial data from Slide-seq",
        "workflow_tutorial_url": "https://hotspot.readthedocs.io/en/latest/Spatial_Tutorial.html",
        "tutorial_defaults": {
            "model": "bernoulli",
            "n_neighbors": 300,
            "core_only": False,
            "min_gene_threshold": 20,
            "fdr_autocorr": 0.05,
            "weighted_graph": False,
            "create_modules": {
                "core_only": False,
                "min_gene_threshold": 20,
                "fdr_threshold": 0.05,
            },
        },
        "tutorial_actual": {
            "model": str(args.model),
            "n_neighbors": int(args.n_neighbors),
            "core_only": parse_bool(args.core_only, False),
            "min_gene_threshold": int(args.min_gene_threshold),
            "fdr_autocorr": float(args.fdr_autocorr),
            "roi_applied": roi_applied,
            "downsample_applied": downsample_applied,
            "max_cells": int(args.max_cells),
            "k_modules_enabled": bool(k_modules_enabled),
            "k_modules": int(args.k_modules),
        },
        "gg_weight_type": str(edges_meta.get("gg_weight_type", "")) if edges_meta else "",
        "gg_fdr_source": str(edges_meta.get("gg_fdr_source", "")) if edges_meta else "",
        "gg_dim": {"n": int(edges_meta.get("gg_genes_n", 0))} if edges_meta else {"n": 0},
        "gg_genes_n": int(edges_meta.get("gg_genes_n", 0)) if edges_meta else 0,
        "input_dataset_id": dataset_id,
        "roi_id": roi_id,
        "roi_hash": roi_hash,
        "roi_source_path": roi_csv or "",
        "roi_source_sha256": roi_source_sha256,
        "obs_id_hash": obs_id_hash,
        "n_obs_raw": int(cell_stats.get("n_cells_input", 0)),
        "n_obs_roi": int(cell_stats.get("n_cells_after_roi", cell_stats.get("n_cells_input", 0))),
        "n_obs_input": int(
            cell_stats.get("n_cells_after_roi", cell_stats.get("n_cells_input", 0))
        ),
        "n_obs_used": int(cell_stats.get("n_cells_used", 0)),
        "n_vars_input": int(n_genes_input),
        "n_vars_used": int(n_genes_used),
        "edges_outfile": str(edges_meta.get("edges_outfile", "gene_graph_edges.tsv"))
        if edges_meta
        else "gene_graph_edges.tsv",
        "edges_schema": list(
            edges_meta.get(
                "edges_schema",
                ["gene_a", "gene_b", "weight", "p_value", "fdr", "n_obs", "edge_type"],
            )
        )
        if edges_meta
        else ["gene_a", "gene_b", "weight", "p_value", "fdr", "n_obs", "edge_type"],
        "edges_all_fullgene_disabled": True,
        "edges_total_n": int(edges_meta.get("edges_total_n", 0)) if edges_meta else 0,
        "edges_sig_n": int(edges_meta.get("edges_sig_n", 0)) if edges_meta else 0,
        "edges_weight_type": str(edges_meta.get("edges_weight_type", "")) if edges_meta else "",
        "edges_pvalue_method": str(edges_meta.get("edges_pvalue_method", "")) if edges_meta else "",
        "edges_fdr_threshold": float(edges_meta.get("edges_fdr_threshold", args.fdr_autocorr))
        if edges_meta
        else float(args.fdr_autocorr),
        "gene_filtering": gene_filtering,
        "random_seed": int(args.seed),
        "stochastic": True,
        "params": {
            "data_dir": args.data_dir,
            "coord_file": roi_csv,
            "max_cells": int(args.max_cells),
            "seed": int(args.seed),
            "n_neighbors": int(args.n_neighbors),
            "model": str(args.model),
            "fdr_autocorr": float(args.fdr_autocorr),
            "min_gene_threshold": int(args.min_gene_threshold),
            "k_modules": int(args.k_modules),
            "emit_gg_matrices": bool(int(getattr(args, "emit_gg_matrices", 0))),
            "emit_edges": parse_bool(getattr(args, "emit_edges", "1"), True),
            "edge_fdr": float(getattr(args, "edge_fdr", 0.05)),
            "edge_top_n": int(getattr(args, "edge_top_n", 1000)),
            "core_only": parse_bool(args.core_only, False),
            "write_module_scores": parse_bool(args.write_module_scores, True),
            "exclude_feature_labels": str(getattr(args, "exclude_feature_labels", "")),
            "ncores": int(args.ncores),
        },
        "runtime": {
            "wall_time_sec": float(wall_time_sec),
            "n_threads": int(args.ncores),
        },
        "resources": read_resource_stats(args.stats_file) if args.stats_file else {},
        "membership_type": "partition",
    }

    with open(os.path.join(args.outdir, "meta.json"), "w", encoding="utf-8") as f:
        import json

        json.dump(meta, f, indent=2)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Hotspot runner for Xenium outs with ROI clip + downsample"
    )
    parser.add_argument("--data_dir", required=True, help="Xenium outs directory or .h5ad")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--coord_file", default="", help="ROI CSV (QuPath format)")
    parser.add_argument("--roi_csv", default="", help="Alias for --coord_file")
    parser.add_argument("--max_cells", type=int, default=0)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--n_neighbors", type=int, default=300)
    parser.add_argument("--model", type=str, default="bernoulli")
    parser.add_argument("--fdr_autocorr", type=float, default=0.05)
    parser.add_argument("--min_gene_threshold", type=int, default=20)
    parser.add_argument("--k_modules", type=int, default=0, help="Force module count (0 disables)")
    parser.add_argument(
        "--emit_gg_matrices",
        type=int,
        default=0,
        help="Emit full gg matrices (0/1, default 0)",
    )
    parser.add_argument("--emit_edges", type=str, default="1")
    parser.add_argument("--edge_fdr", type=float, default=0.05)
    parser.add_argument("--edge_top_n", type=int, default=1000)
    parser.add_argument("--core_only", type=str, default="false")
    parser.add_argument("--write_module_scores", type=str, default="1")
    parser.add_argument("--ncores", type=int, default=8)
    parser.add_argument("--dataset_id", default="", help="Dataset identifier")
    parser.add_argument("--roi_id", default="", help="ROI identifier")
    parser.add_argument("--stats_file", default="", help="Cgroup stats TSV for resources")
    parser.add_argument(
        "--exclude_feature_labels",
        type=str,
        default=DEFAULT_EXCLUDE_FEATURE_LABELS,
        help=(
            "Comma-separated feature labels to exclude (case-insensitive substring match). "
            "Empty string disables label filtering."
        ),
    )

    args = parser.parse_args()

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    np.random.seed(args.seed)
    roi_csv = args.coord_file or args.roi_csv or ""

    exclude_prefix = DEFAULT_EXCLUDE_PREFIX

    start_time = time.time()
    cell_stats: Dict[str, int] = {}
    gene_filtering_meta: Dict[str, object] = {}
    n_genes_input = 0
    n_genes_used = 0

    if os.path.isfile(args.data_dir) and args.data_dir.endswith(".h5ad"):
        adata = load_h5ad(args.data_dir)
        adata = ensure_xy(adata)
        adata, cell_stats = apply_roi_and_downsample_adata(
            adata, roi_csv, args.max_cells, args.seed
        )
    else:
        cells_pq = os.path.join(args.data_dir, "cells.parquet")
        tx_pq = os.path.join(args.data_dir, "transcripts.parquet")
        if os.path.exists(cells_pq) and os.path.exists(tx_pq):
            cells_df = load_cells_df(cells_pq)
            cells_df, cell_stats = apply_roi_and_downsample(
                cells_df, roi_csv, args.max_cells, args.seed
            )
            adata = build_adata_from_transcripts(cells_df, tx_pq, exclude_prefix)
        else:
            h5ad_path = None
            for cand in ["adata_raw.h5ad", "adata.h5ad", "raw.h5ad"]:
                p = os.path.join(args.data_dir, cand)
                if os.path.exists(p):
                    h5ad_path = p
                    break
            if h5ad_path is None:
                raise ValueError(
                    "data_dir missing cells.parquet/transcripts.parquet and no .h5ad fallback found"
                )
            adata = load_h5ad(h5ad_path)
            adata = ensure_xy(adata)
            adata, cell_stats = apply_roi_and_downsample_adata(
                adata, roi_csv, args.max_cells, args.seed
            )

    adata = ensure_xy(adata)
    n_genes_input = int(adata.n_vars)
    adata, gene_filtering_meta = filter_features_by_label(
        adata, getattr(args, "exclude_feature_labels", DEFAULT_EXCLUDE_FEATURE_LABELS)
    )
    n_genes_used = int(adata.n_vars)
    adata = prepare_adata_for_hotspot(adata)

    _, edges_meta = run_hotspot(
        adata=adata,
        n_neighbors=args.n_neighbors,
        model=args.model,
        fdr_autocorr=args.fdr_autocorr,
        min_gene_threshold=args.min_gene_threshold,
        k_modules=args.k_modules,
        emit_edges=parse_bool(getattr(args, "emit_edges", "1"), True),
        edge_fdr=float(getattr(args, "edge_fdr", 0.05)),
        edge_top_n=int(getattr(args, "edge_top_n", 1000)),
        emit_gg_matrices=bool(int(args.emit_gg_matrices)),
        core_only=parse_bool(args.core_only, False),
        ncores=args.ncores,
        outdir=outdir,
    )

    wall_time_sec = time.time() - start_time
    try:
        _write_meta_json(
            args,
            cell_stats,
            n_genes_input,
            n_genes_used,
            edges_meta,
            wall_time_sec,
            gene_filtering_meta=gene_filtering_meta,
        )
    except Exception as exc:
        eprint(f"[WARN] meta.json write failed: {exc}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
