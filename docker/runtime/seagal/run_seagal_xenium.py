#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SEAGAL runner for Xenium outs with ROI clip + downsample + bench outputs.
"""
import argparse
import json
import os
import sys
import time
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd


def eprint(*args):
    print(*args, file=sys.stderr)


def _try_importlib_version(pkg: str) -> str:
    try:
        from importlib.metadata import version
        return version(pkg)
    except Exception:
        return "unknown"


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

    # drop consecutive duplicates
    if pts.shape[0] >= 2:
        keep = [True]
        for i in range(1, pts.shape[0]):
            keep.append(not np.allclose(pts[i], pts[i - 1]))
        pts = pts[np.array(keep, dtype=bool)]

    # drop closing point if equals first
    if pts.shape[0] >= 2 and np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]

    if pts.shape[0] < 3:
        raise ValueError(
            "ROI polygon needs >=3 vertices after cleaning. "
            f"File={csv_path} Columns={list(df.columns)}"
        )

    # ensure closed polygon
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


def _iter_arrow_batches(dataset, tx_pq_path: str, columns: List[str], batch_size: int):
    if hasattr(dataset, "scanner"):
        scanner = dataset.scanner(columns=columns, batch_size=batch_size)
        for b in scanner.to_batches():
            yield b
        return

    if hasattr(dataset, "to_batches"):
        for b in dataset.to_batches(columns=columns, batch_size=batch_size):
            yield b
        return

    import pyarrow.parquet as pq
    pf = pq.ParquetFile(tx_pq_path)
    for b in pf.iter_batches(batch_size=batch_size, columns=columns):
        yield b


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
) -> Tuple[pd.DataFrame, int]:
    roi_applied = 0
    if roi_csv:
        roi = read_roi_polygon(roi_csv)
        before = len(cells_df)
        cells_df = clip_cells_to_roi(cells_df, roi)
        roi_applied = 1
        eprint(f"[INFO] ROI clip on cells: {before} -> {len(cells_df)}")

    if max_cells > 0 and len(cells_df) > max_cells:
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(cells_df), size=max_cells, replace=False)
        cells_df = cells_df.iloc[idx].copy()
        eprint(f"[INFO] Downsample cells to max_cells={max_cells}: n={len(cells_df)}")

    return cells_df, roi_applied


def build_adata_from_transcripts(
    cells_df: pd.DataFrame,
    tx_pq: str,
    exclude_prefix: Tuple[str, ...],
) -> "anndata.AnnData":
    import pyarrow.dataset as ds
    from scipy import sparse
    import anndata as ad

    cell_ids = cells_df["cell_id"].astype(str).tolist()
    if not cell_ids:
        raise ValueError("No cells left after ROI/downsample. Check ROI and parameters.")

    cell_indexer = pd.Index(cell_ids, dtype=str)

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
    gene_index: Dict[str, int] = {}
    genes: List[str] = []
    row_chunks: List[np.ndarray] = []
    col_chunks: List[np.ndarray] = []
    data_chunks: List[np.ndarray] = []

    batch_size = 1_000_000
    columns = [cell_id_col, gene_col]

    for batch in _iter_arrow_batches(dataset, tx_pq, columns=columns, batch_size=batch_size):
        df = batch.to_pandas()
        scanned += len(df)
        cids = df[cell_id_col].astype(str)
        row_idx = cell_indexer.get_indexer(cids)
        mask_cell = row_idx >= 0
        if not mask_cell.any():
            continue

        row_idx = row_idx[mask_cell].astype(np.int32, copy=False)
        gser = df.loc[mask_cell, gene_col].astype(str)

        if exclude_prefix:
            mask_gene = ~gser.str.startswith(exclude_prefix)
            if not mask_gene.any():
                continue
            row_idx = row_idx[mask_gene.to_numpy(dtype=bool, copy=False)]
            gser = gser.loc[mask_gene]

        kept += len(gser)

        g_arr = gser.to_numpy(dtype=object, copy=False)
        uniq_genes, inv = np.unique(g_arr, return_inverse=True)
        uniq_genes = [str(g) for g in uniq_genes.tolist()]

        for g in uniq_genes:
            if g not in gene_index:
                gene_index[g] = len(genes)
                genes.append(g)

        cols_for_uniq = np.fromiter(
            (gene_index[g] for g in uniq_genes),
            dtype=np.int32,
            count=len(uniq_genes),
        )
        col_idx = cols_for_uniq[inv].astype(np.int32, copy=False)

        keys = (row_idx.astype(np.uint64) << np.uint64(32)) | col_idx.astype(np.uint64)
        uniq_keys, cnt = np.unique(keys, return_counts=True)
        row_u = (uniq_keys >> np.uint64(32)).astype(np.int32, copy=False)
        col_u = (uniq_keys & np.uint64(0xFFFFFFFF)).astype(np.int32, copy=False)

        row_chunks.append(row_u)
        col_chunks.append(col_u)
        data_chunks.append(cnt.astype(np.int32, copy=False))

    eprint(
        f"[INFO] transcripts aggregation: scanned_rows={scanned} kept_rows={kept} cells={len(cell_ids)}"
    )

    if not row_chunks:
        raise ValueError("No transcript counts after filtering. Check ROI/downsample.")

    rows = np.concatenate(row_chunks)
    cols = np.concatenate(col_chunks)
    data = np.concatenate(data_chunks)

    X = sparse.csr_matrix(
        sparse.coo_matrix(
            (data, (rows, cols)),
            shape=(len(cell_ids), len(genes)),
            dtype=np.float32,
        )
    )
    X.sum_duplicates()

    gene_arr = np.asarray(genes, dtype=str)
    if gene_arr.size > 1:
        order = np.argsort(gene_arr)
        if not np.array_equal(order, np.arange(gene_arr.size)):
            X = X[:, order]
            gene_arr = gene_arr[order]

    adata = ad.AnnData(X=X)
    adata.obs_names = pd.Index(cell_ids, dtype=str)
    adata.var_names = pd.Index(gene_arr.tolist(), dtype=str)

    coords = cells_df.set_index("cell_id").loc[adata.obs_names, ["x", "y"]]
    adata.obs["x"] = coords["x"].to_numpy()
    adata.obs["y"] = coords["y"].to_numpy()
    adata.obsm["spatial"] = coords[["x", "y"]].to_numpy()

    return adata


def load_from_h5(h5_path: str) -> "anndata.AnnData":
    import scanpy as sc

    adata = sc.read_10x_h5(h5_path, gex_only=False)
    if "feature_types" in adata.var.columns:
        mask = adata.var["feature_types"] == "Gene Expression"
        if mask.any():
            adata = adata[:, mask].copy()
    return adata


def load_from_h5ad(path: str) -> "anndata.AnnData":
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


def filter_top_genes(adata: "anndata.AnnData", top_genes: int) -> "anndata.AnnData":
    if top_genes <= 0 or adata.n_vars <= top_genes:
        return adata

    from scipy import sparse

    X = adata.X
    if sparse.issparse(X):
        means = np.asarray(X.mean(axis=0)).ravel()
    else:
        means = X.mean(axis=0)
    top_idx = np.argsort(means)[::-1][:top_genes]
    mask = np.zeros(adata.n_vars, dtype=bool)
    mask[top_idx] = True
    return adata[:, mask].copy()


def _estimate_dense_connectivities_gib(n_obs: int, *, bytes_per_entry: int = 8) -> float:
    if n_obs <= 0:
        return 0.0
    return float(n_obs) * float(n_obs) * float(bytes_per_entry) / float(1024**3)


def bin_adata_to_grid(adata: "anndata.AnnData", grid_um: float) -> "anndata.AnnData":
    if grid_um <= 0:
        return adata

    from scipy import sparse
    import anndata as ad

    adata = ensure_xy(adata)
    x = pd.to_numeric(adata.obs["x"], errors="coerce").to_numpy(dtype=float, copy=False)
    y = pd.to_numeric(adata.obs["y"], errors="coerce").to_numpy(dtype=float, copy=False)

    valid = np.isfinite(x) & np.isfinite(y)
    if not valid.all():
        adata = adata[valid].copy()
        x = x[valid]
        y = y[valid]
        if adata.n_obs == 0:
            raise ValueError("All observations have invalid x/y after grid binning filter.")

    x0 = float(np.nanmin(x))
    y0 = float(np.nanmin(y))
    gx = np.floor((x - x0) / float(grid_um)).astype(np.int64, copy=False)
    gy = np.floor((y - y0) / float(grid_um)).astype(np.int64, copy=False)

    pairs = np.stack([gx, gy], axis=1)
    uniq_pairs, inv = np.unique(pairs, axis=0, return_inverse=True)
    inv = inv.astype(np.int32, copy=False)

    n_bins = int(uniq_pairs.shape[0])
    n_cells = int(adata.n_obs)
    if n_bins <= 0:
        raise ValueError("Grid binning produced zero bins.")

    X = adata.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    else:
        X = X.tocsr()

    rows = inv
    cols = np.arange(n_cells, dtype=np.int32)
    data = np.ones(n_cells, dtype=np.int8)
    G = sparse.csr_matrix((data, (rows, cols)), shape=(n_bins, n_cells))
    Xb = (G @ X).tocsr()
    Xb.sum_duplicates()

    counts = np.bincount(rows, minlength=n_bins).astype(np.int64, copy=False)
    x_sum = np.bincount(rows, weights=x, minlength=n_bins)
    y_sum = np.bincount(rows, weights=y, minlength=n_bins)
    x_mean = x_sum / np.maximum(counts, 1)
    y_mean = y_sum / np.maximum(counts, 1)

    bin_ids = [f"bin_{int(a)}_{int(b)}" for a, b in uniq_pairs.tolist()]
    obs = pd.DataFrame(
        {"x": x_mean, "y": y_mean, "n_cells_in_bin": counts.astype(int)},
        index=pd.Index(bin_ids, dtype=str),
    )

    bdata = ad.AnnData(X=Xb, obs=obs)
    bdata.var_names = adata.var_names.astype(str)
    bdata.obsm["spatial"] = obs[["x", "y"]].to_numpy(dtype=float, copy=False)
    return bdata


def parse_extra_args(extra_args_json: str) -> Dict[str, object]:
    if not extra_args_json:
        return {}
    if os.path.exists(extra_args_json):
        with open(extra_args_json, "r", encoding="utf-8") as f:
            return json.load(f)
    return json.loads(extra_args_json)


def _safe_bool(v: object, default: bool) -> bool:
    if v is None:
        return bool(default)
    if isinstance(v, bool):
        return v
    if isinstance(v, (int, np.integer)):
        return bool(int(v))
    if isinstance(v, float):
        if not np.isfinite(v):
            return bool(default)
        return bool(int(v))
    s = str(v).strip().lower()
    if s in {"1", "true", "t", "yes", "y"}:
        return True
    if s in {"0", "false", "f", "no", "n"}:
        return False
    return bool(default)


def _safe_int(v: object, default: int) -> int:
    try:
        if v is None:
            return int(default)
        if isinstance(v, bool):
            return int(v)
        return int(v)
    except Exception:
        return int(default)


def _safe_float(v: object, default: float) -> float:
    try:
        if v is None:
            return float(default)
        return float(v)
    except Exception:
        return float(default)


_SEAGAL_REQUIRED_FUNCS = ("load_raw", "process_st", "spatial_pattern_genes")
_SEAGAL_STATE_CLASS_CANDIDATES = ("Seagal", "SEAGAL", "SeagalObj", "SeagalObject")


def _resolve_seagal_library_api() -> Tuple[object, str]:
    """
    Resolve SEAGAL *library* API without importing the CLI-style entrypoint
    (e.g. `seagal.seagal`) that may call `sys.exit(...)`.
    """
    import importlib
    import types

    # Prefer internal modules that define the core pipeline functions without CLI parsing.
    # Intentionally avoid importing `seagal.seagal` (CLI wrapper) and the package root (which may
    # re-export CLI symbols in some builds).
    candidates = (
        "seagal._utils",
        "seagal._gmod",
        "seagal.modules.gmod",
    )

    modules: Dict[str, object] = {}
    import_errors: Dict[str, str] = {}
    for mod in candidates:
        try:
            modules[mod] = importlib.import_module(mod)
        except Exception as exc:
            import_errors[mod] = f"{type(exc).__name__}: {exc}"

    composite = types.SimpleNamespace()
    for m in modules.values():
        for name in (
            "load_raw",
            "process_st",
            "spatial_process",
            "spatial_pattern_genes",
            "spatial_association",
            "group_adata_by_genes",
            "genemodules",
        ):
            fn = getattr(m, name, None)
            if callable(fn) and not hasattr(composite, name):
                setattr(composite, name, fn)

    required = ("process_st", "spatial_pattern_genes")
    missing = [fn for fn in required if not callable(getattr(composite, fn, None))]
    if missing:
        raise ImportError(
            "Cannot resolve SEAGAL library API. "
            f"Missing={missing} candidates={list(candidates)} import_errors={import_errors}"
        )

    entrypoint = "library(seagal._utils+seagal._gmod)"
    eprint(f"[INFO] SEAGAL library api resolved: {entrypoint}")
    return composite, entrypoint


def _load_count_meta_csv_as_adata(count_csv: str, meta_csv: str) -> "anndata.AnnData":
    import anndata as ad
    from scipy import sparse

    count = pd.read_csv(count_csv, index_col=0)
    meta = pd.read_csv(meta_csv, index_col=0)

    if count.empty:
        raise ValueError(f"count.csv is empty: {count_csv}")
    if meta.empty:
        raise ValueError(f"meta.csv is empty: {meta_csv}")
    if "x" not in meta.columns or "y" not in meta.columns:
        raise ValueError(f"meta.csv must contain x/y columns: cols={list(meta.columns)} file={meta_csv}")

    meta = meta.copy()
    meta["x"] = pd.to_numeric(meta["x"], errors="coerce")
    meta["y"] = pd.to_numeric(meta["y"], errors="coerce")
    if not np.isfinite(meta["x"]).all() or not np.isfinite(meta["y"]).all():
        raise ValueError("meta.csv contains non-finite x/y values after numeric coercion.")

    if not meta.index.equals(count.index):
        if set(meta.index) == set(count.index):
            meta = meta.reindex(count.index)
        else:
            inter = meta.index.intersection(count.index)
            if inter.empty:
                raise ValueError("count.csv and meta.csv have no overlapping cell ids.")
            eprint(
                f"[WARN] count/meta index mismatch; using intersection n={len(inter)} "
                f"(count={len(count)} meta={len(meta)})"
            )
            count = count.loc[inter]
            meta = meta.loc[inter]

    genes = [str(g) for g in count.columns.astype(str).tolist()]
    var_df = pd.DataFrame({"gene": genes}, index=pd.Index(genes, dtype=str))
    X = sparse.csr_matrix(count.to_numpy(dtype=float, copy=False))
    adata = ad.AnnData(X=X, obs=meta, var=var_df)
    adata.obsm["spatial"] = meta[["x", "y"]].to_numpy(dtype=float, copy=True)
    return adata


def resolve_seagal_api() -> Tuple[object, str, Optional[type]]:
    api, entrypoint = _resolve_seagal_library_api()
    return api, entrypoint, None


def _call_api(func, *args, **kwargs):
    import inspect

    try:
        sig = inspect.signature(func)
    except Exception:
        return func(*args, **kwargs)

    params = sig.parameters
    accepts_varkw = any(p.kind == inspect.Parameter.VAR_KEYWORD for p in params.values())

    alias: Dict[str, List[str]] = {
        "fdr_cutoff": ["FDR_cutoff", "fdr_cutoff", "fdrCutoff"],
        "l_cutoff": ["L_cutoff", "l_cutoff", "lCutoff"],
        "topK": ["topK", "top_k", "topk"],
        "topk": ["topK", "top_k", "topk"],
        "min_counts": ["min_counts", "minCounts"],
        "min_cells": ["min_cells", "minCells"],
        "use_pattern_genes": ["use_pattern_genes", "usePatternGenes"],
        "n_permutation": ["n_permutation", "n_permutations", "nPermutation", "n_perm"],
        "permute_ratio": ["permute_ratio", "permuteRatio", "perm_ratio"],
        "indep": ["indep", "independent", "independent_test"],
        "grouped_only": ["grouped_only", "groupedOnly"],
        "grouped": ["grouped", "grouped_only"],
        "use_grouped": ["use_grouped", "useGrouped"],
        "n_modules": ["n_modules", "nModules", "n_module"],
        "nmax": ["nmax", "n_max", "nMax", "max_k"],
    }

    filtered: Dict[str, object] = {}
    for k, v in kwargs.items():
        if accepts_varkw:
            filtered[k] = v
            continue
        if k in params:
            filtered[k] = v
            continue
        for alt in alias.get(k, []):
            if alt in params:
                filtered[alt] = v
                break

    return func(*args, **filtered)


def write_seagal_csv_inputs(adata: "anndata.AnnData", outdir: str) -> Tuple[str, str]:
    from scipy import sparse

    seagal_in_dir = os.path.join(outdir, "_seagal_input")
    os.makedirs(seagal_in_dir, exist_ok=True)

    adata = ensure_xy(adata)
    obs_names = adata.obs_names.astype(str)
    var_names = adata.var_names.astype(str)

    X = adata.X
    if sparse.issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X)

    count_df = pd.DataFrame(X, index=obs_names, columns=var_names)
    meta_df = pd.DataFrame(
        {
            "x": pd.to_numeric(adata.obs["x"], errors="coerce"),
            "y": pd.to_numeric(adata.obs["y"], errors="coerce"),
        },
        index=obs_names,
    )

    count_path = os.path.join(seagal_in_dir, "count.csv")
    meta_path = os.path.join(seagal_in_dir, "meta.csv")
    count_df.to_csv(count_path)
    meta_df.to_csv(meta_path)

    return count_path, meta_path


def run_seagal_official_api(
    count_csv: str,
    meta_csv: str,
    *,
    min_counts: int,
    min_cells: int,
    svg_i: Optional[float],
    svg_topk: Optional[int],
    use_pattern_genes: bool,
    n_permutation: int,
    permute_ratio: float,
    fdr_cutoff: float,
    l_cutoff: float,
    indep: bool,
    dense_max_obs: int,
    modules_enabled: bool,
    modules_nmax: int,
    modules_n: Optional[int],
) -> Tuple[
    object,
    pd.DataFrame,
    pd.DataFrame,
    Optional[Dict[str, object]],
    List[str],
    Dict[str, object],
]:
    import types
    import inspect

    def _is_adata(obj: object) -> bool:
        return all(hasattr(obj, a) for a in ("obs", "var", "X", "n_obs", "n_vars"))

    def _coerce_adata(ret: object, state: object) -> Optional[object]:
        if _is_adata(ret):
            return ret
        maybe = getattr(state, "adata", None)
        if _is_adata(maybe):
            return maybe
        return None

    def _find_df_with_cols(container: object, required_cols: List[str]) -> Optional[pd.DataFrame]:
        if container is None:
            return None
        if isinstance(container, pd.DataFrame):
            return container if all(c in container.columns for c in required_cols) else None
        if isinstance(container, dict):
            for v in container.values():
                hit = _find_df_with_cols(v, required_cols)
                if hit is not None:
                    return hit
        if hasattr(container, "uns") and isinstance(getattr(container, "uns"), dict):
            for v in container.uns.values():
                hit = _find_df_with_cols(v, required_cols)
                if hit is not None:
                    return hit
        return None

    def _find_module_dict(container: object) -> Optional[Dict[str, object]]:
        key_candidates = ("gene2mod", "gene_to_module", "gene2module", "gene_module", "gene_module_dict")

        def _is_mod_dict(d: object) -> bool:
            if not isinstance(d, dict):
                return False
            return any(isinstance(d.get(k), dict) and d.get(k) for k in key_candidates)

        if _is_mod_dict(container):
            return container  # type: ignore[return-value]
        if isinstance(container, dict):
            for v in container.values():
                hit = _find_module_dict(v)
                if hit is not None:
                    return hit
        if hasattr(container, "uns") and isinstance(getattr(container, "uns"), dict):
            return _find_module_dict(container.uns)
        return None

    def _func_accepts_kw(func, name: str) -> bool:
        try:
            sig = inspect.signature(func)
        except Exception:
            return False
        if name in sig.parameters:
            return True
        return any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values())

    def _is_sg_like(obj: object) -> bool:
        return obj is not None and hasattr(obj, "adata")

    def _first_positional_param_name(func) -> Optional[str]:
        try:
            sig = inspect.signature(func)
        except Exception:
            return None
        for p in sig.parameters.values():
            if p.kind in (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD):
                return str(p.name)
        return None

    def _spatial_association_call_modes(
        func,
        *,
        sg_obj: object,
        grouped_obj: object,
        adata_obj: object,
    ) -> List[Tuple[str, Tuple[object, ...]]]:
        """
        SEAGAL's spatial_association is expensive; avoid calling it multiple times with incompatible
        inputs. Prefer sg-like objects (with `.adata`) unless the signature clearly indicates AnnData.
        """
        first = (_first_positional_param_name(func) or "").strip().lower()
        adata_like_first = {"adata", "anndata", "data", "st"}

        modes: List[Tuple[str, Tuple[object, ...]]] = []
        if first in adata_like_first and _is_adata(adata_obj):
            modes.append(("adata", (adata_obj,)))
            if grouped_obj is not None and grouped_obj is not adata_obj and _is_adata(grouped_obj):
                modes.append(("grouped_result", (grouped_obj,)))
            if sg_obj is not None and sg_obj is not grouped_obj and _is_sg_like(sg_obj):
                modes.append(("sg", (sg_obj,)))
        else:
            if sg_obj is not None and _is_sg_like(sg_obj):
                modes.append(("sg", (sg_obj,)))
            if grouped_obj is not None and grouped_obj is not sg_obj and _is_sg_like(grouped_obj):
                modes.append(("grouped_result", (grouped_obj,)))

        # de-dup by identity
        seen: set[int] = set()
        out: List[Tuple[str, Tuple[object, ...]]] = []
        for mode, args in modes:
            if not args or args[0] is None:
                continue
            key = id(args[0])
            if key in seen:
                continue
            seen.add(key)
            out.append((mode, args))
        return out

    api, entrypoint = _resolve_seagal_library_api()

    meta: Dict[str, object] = {
        "seagal_entrypoint": str(entrypoint),
        "state_object": "namespace",
        "official_api_used": True,
        "official_api_error": None,
        "official_api_chain_effective": [],
        "official_api_chain_missing": [],
        "api_chain_effective": [],
        "api_chain_missing": [],
        "min_counts": int(min_counts),
        "min_cells": int(min_cells),
        "svg_I_used": "NA",
        "svg_topk_used": "NA",
        "n_obs_post_process": 0,
        "n_vars_post_process": 0,
        "use_pattern_genes": int(bool(use_pattern_genes)),
        "n_permutation": int(n_permutation),
        "permute_ratio": float(permute_ratio),
        "fdr_cutoff": float(fdr_cutoff),
        "l_cutoff": float(l_cutoff),
        "indep": int(bool(indep)),
        "dense_connectivities_max_obs": int(dense_max_obs),
        "modules_enabled": int(bool(modules_enabled)),
        "modules_nmax": int(modules_nmax),
        "modules_n": int(modules_n) if modules_n is not None else None,
    }
    def chain_ok(name: str) -> None:
        if name not in meta["api_chain_effective"]:
            meta["api_chain_effective"].append(name)
        if name not in meta["official_api_chain_effective"]:
            meta["official_api_chain_effective"].append(name)

    def chain_miss(name: str) -> None:
        if name not in meta["api_chain_missing"]:
            meta["api_chain_missing"].append(name)
        if name not in meta["official_api_chain_missing"]:
            meta["official_api_chain_missing"].append(name)

    sg = types.SimpleNamespace()

    # ---- count/meta -> AnnData (no SEAGAL CLI wrappers; no SystemExit) ----
    try:
        adata = _load_count_meta_csv_as_adata(count_csv, meta_csv)
        meta["load_raw_mode"] = "pandas_count_meta"
        chain_ok("load_raw")
    except BaseException as exc:
        meta["official_api_error"] = f"{type(exc).__name__}: {exc}"
        chain_miss("load_raw")
        raise

    sg.adata = adata

    # ---- process_st ----
    process_st = getattr(api, "process_st", None)
    if not callable(process_st):
        chain_miss("process_st")
        raise AttributeError(f"SEAGAL api missing process_st at entrypoint={entrypoint}")

    proc_exc: Optional[BaseException] = None
    for mode, args in [("sg_adata", (sg, adata)), ("adata", (adata,))]:
        try:
            ret = _call_api(process_st, *args, min_counts=min_counts, min_cells=min_cells)
            adata2 = _coerce_adata(ret, sg)
            if adata2 is not None:
                adata = adata2
            meta["process_st_mode"] = mode
            sg.adata = adata
            break
        except BaseException as exc:
            proc_exc = exc
            continue
    else:
        chain_miss("process_st")
        raise RuntimeError(f"SEAGAL process_st failed. entrypoint={entrypoint} err={proc_exc}")
    chain_ok("process_st")
    meta["n_obs_post_process"] = int(getattr(adata, "n_obs", 0))
    meta["n_vars_post_process"] = int(getattr(adata, "n_vars", 0))

    if adata.n_obs == 0 or adata.n_vars == 0:
        sg.adata = adata
        sg.co_expression = pd.DataFrame()
        sg.module_dict = None
        return sg, pd.DataFrame(), pd.DataFrame(), None, [], meta

    # ---- spatial_process (optional) ----
    spatial_process = getattr(api, "spatial_process", None)
    if callable(spatial_process):
        sp_exc: Optional[BaseException] = None
        for args in [(adata,), (sg, adata), (sg,)]:
            try:
                ret = _call_api(spatial_process, *args)
                adata2 = _coerce_adata(ret, sg)
                if adata2 is not None:
                    adata = adata2
                sg.adata = adata
                chain_ok("spatial_process")
                break
            except BaseException as exc:
                sp_exc = exc
                continue
        else:
            chain_miss("spatial_process")
            meta["spatial_process_error"] = str(sp_exc)
    else:
        chain_miss("spatial_process")

    # ---- spatial_pattern_genes (required) ----
    spatial_pattern_genes = getattr(api, "spatial_pattern_genes", None)
    if not callable(spatial_pattern_genes):
        chain_miss("spatial_pattern_genes")
        raise AttributeError(f"SEAGAL api missing spatial_pattern_genes at entrypoint={entrypoint}")

    n_vars = int(getattr(adata, "n_vars", 0))
    svg_topk_used: Optional[int] = None
    svg_i_used: Optional[float] = None

    if svg_topk is None and svg_i is None:
        svg_topk = min(1000, n_vars)

    pat_exc: Optional[BaseException] = None
    if svg_topk is not None:
        svg_topk_used = max(1, min(int(svg_topk), n_vars))
        try:
            var_sorted = adata.var.sort_values("moranI", ascending=False)
            svg_i_used = float(var_sorted["moranI"].iloc[svg_topk_used - 1])
        except Exception:
            svg_i_used = None
        pat_kwargs = {"I": None, "topK": svg_topk_used}
    else:
        svg_i_used = float(svg_i) if svg_i is not None else None
        pat_kwargs = {"I": svg_i_used, "topK": None}

    for mode, args in [("sg", (sg,)), ("adata", (adata,))]:
        try:
            ret = _call_api(spatial_pattern_genes, *args, **pat_kwargs)
            adata2 = _coerce_adata(ret, sg)
            if adata2 is not None:
                adata = adata2
            sg.adata = adata
            meta["spatial_pattern_genes_mode"] = mode
            chain_ok("spatial_pattern_genes")
            break
        except BaseException as exc:
            pat_exc = exc
            continue
    else:
        chain_miss("spatial_pattern_genes")
        raise RuntimeError(
            f"SEAGAL spatial_pattern_genes failed. entrypoint={entrypoint} err={pat_exc}"
        )

    meta["svg_I_used"] = float(svg_i_used) if svg_i_used is not None else "NA"
    meta["svg_topk_used"] = int(svg_topk_used) if svg_topk_used is not None else "NA"

    mask = None
    try:
        mask = getattr(sg, "adata", adata).var.get("high_pattern_genes")
    except Exception:
        mask = None
    n_svg = int(mask.sum()) if mask is not None else 0
    meta["n_svg"] = int(n_svg)

    meta["n_pattern_genes"] = int(n_svg)

    meta["edges_all_fullgene_disabled"] = True
    meta["n_assoc_genes_full"] = 0

    # ---- group_adata_by_genes(grouped=False) (optional/tutorial-faithful) ----
    group_adata_by_genes = getattr(api, "group_adata_by_genes", None)
    grouped_result = None
    if callable(group_adata_by_genes) and (
        _func_accepts_kw(group_adata_by_genes, "grouped")
        or _func_accepts_kw(group_adata_by_genes, "grouped_only")
    ):
        grp_exc: Optional[BaseException] = None
        for mode, args in [("sg", (sg,)), ("adata", (adata,))]:
            try:
                grouped_result = _call_api(group_adata_by_genes, *args, grouped=False)
                meta["group_adata_by_genes_mode"] = mode
                chain_ok("group_adata_by_genes")
                break
            except BaseException as exc:
                grp_exc = exc
                continue
        else:
            chain_miss("group_adata_by_genes")
            meta["group_adata_by_genes_skipped_reason"] = f"call_failed:{grp_exc}"
    else:
        chain_miss("group_adata_by_genes")
        meta["group_adata_by_genes_skipped_reason"] = (
            "not_found" if not callable(group_adata_by_genes) else "signature_no_grouped"
        )

    coexpr_full = pd.DataFrame(
        columns=[
            "gene_1",
            "gene_2",
            "L",
            "L.p_value",
            "L.FDR",
            "pair",
            "-log10(FDR)",
            "Association",
        ]
    )
    coexpr_pat = pd.DataFrame()
    module_dict_pat: Optional[Dict[str, object]] = None
    meta["full_edges_rows"] = 0
    try:
        sg.co_expression = coexpr_full
    except Exception:
        pass

    # ---- pattern-gene-only association + genemodules (modules.tsv) ----
    sg.module_dict = None
    meta["pattern_modules_enabled"] = int(bool(modules_enabled))
    pattern_genes_used: List[str] = []
    if not modules_enabled:
        meta["pattern_modules_skipped_reason"] = "modules_disabled"
    elif mask is None:
        meta["pattern_modules_skipped_reason"] = "high_pattern_genes_missing"
    elif int(n_svg) < 2:
        meta["pattern_modules_skipped_reason"] = "n_pattern_genes_lt2"
    else:
        try:
            adata_pat = getattr(sg, "adata", adata)[:, mask].copy()
            meta["n_pattern_genes_used"] = int(getattr(adata_pat, "n_vars", 0))
            try:
                pattern_genes_used = [str(g) for g in adata_pat.var_names.astype(str).tolist()]
            except Exception:
                pattern_genes_used = [str(g) for g in list(getattr(adata_pat, "var_names", []))]

            sg_pat = types.SimpleNamespace()
            sg_pat.adata = adata_pat

            grouped_result_pat = None
            if callable(group_adata_by_genes) and (
                _func_accepts_kw(group_adata_by_genes, "grouped")
                or _func_accepts_kw(group_adata_by_genes, "grouped_only")
            ):
                grp_exc_pat: Optional[BaseException] = None
                for mode, args in [("sg", (sg_pat,)), ("adata", (adata_pat,))]:
                    try:
                        grouped_result_pat = _call_api(group_adata_by_genes, *args, grouped=False)
                        meta["pattern_group_adata_by_genes_mode"] = mode
                        chain_ok("group_adata_by_genes")
                        break
                    except BaseException as exc:
                        grp_exc_pat = exc
                        continue
                else:
                    meta["pattern_group_adata_by_genes_skipped_reason"] = f"call_failed:{grp_exc_pat}"
            else:
                meta["pattern_group_adata_by_genes_skipped_reason"] = (
                    "not_found" if not callable(group_adata_by_genes) else "signature_no_grouped"
                )

            n_obs_pat = int(getattr(adata_pat, "n_obs", 0))
            if int(dense_max_obs) > 0 and n_obs_pat > int(dense_max_obs):
                est_gib = _estimate_dense_connectivities_gib(n_obs_pat)
                meta["pattern_modules_skipped_reason"] = "dense_connectivities_guard"
                meta["pattern_modules_skipped_detail"] = (
                    f"n_obs={n_obs_pat} dense_max_obs={int(dense_max_obs)} est_dense_connectivities_gib={est_gib:.1f}"
                )
                chain_miss("spatial_association")
                chain_miss("genemodules")
                if meta.get("official_api_error") in (None, "NA"):
                    meta["official_api_error"] = (
                        f"dense_connectivities_guard: n_obs={n_obs_pat} > dense_max_obs={int(dense_max_obs)} "
                        f"(est ~{est_gib:.1f} GiB for dense connectivities)"
                    )
                eprint(
                    "[WARN] official_api: skipping spatial_association/genemodules due to dense_connectivities_guard. "
                    f"n_obs={n_obs_pat} dense_max_obs={int(dense_max_obs)} est_dense_connectivities_gib={est_gib:.1f}"
                )
            else:
                spatial_association = getattr(api, "spatial_association", None)
                if int(getattr(adata_pat, "n_vars", 0)) >= 2 and callable(spatial_association):
                    assoc_exc_root_pat: Optional[BaseException] = None
                    assoc_errors_pat: Dict[str, str] = {}
                    for mode, args in _spatial_association_call_modes(
                        spatial_association,
                        sg_obj=sg_pat,
                        grouped_obj=grouped_result_pat,
                        adata_obj=adata_pat,
                    ):
                        try:
                            ret = _call_api(
                                spatial_association,
                                *args,
                                grouped_only=False,
                                use_pattern_genes=False,
                                genes=None,
                                n_permutation=n_permutation,
                                permute_ratio=permute_ratio,
                                fdr_cutoff=fdr_cutoff,
                                l_cutoff=l_cutoff,
                                FDR_cutoff=fdr_cutoff,
                                L_cutoff=l_cutoff,
                                indep=indep,
                            )
                            coexpr = _find_df_with_cols(ret, ["gene_1", "gene_2", "L"])
                            if coexpr is None:
                                coexpr = _find_df_with_cols(
                                    getattr(sg_pat, "co_expression", None), ["gene_1", "gene_2", "L"]
                                )
                            if coexpr is None:
                                coexpr = _find_df_with_cols(
                                    getattr(sg_pat, "adata", None), ["gene_1", "gene_2", "L"]
                                )
                            if coexpr is None:
                                raise RuntimeError(
                                    "Cannot locate SEAGAL SAG/co_expression table after spatial_association (pattern)"
                                )
                            sg_pat.co_expression = coexpr
                            coexpr_pat = coexpr.copy()
                            meta["spatial_association_pattern_mode"] = mode
                            chain_ok("spatial_association")
                            break
                        except BaseException as exc:
                            if assoc_exc_root_pat is None:
                                assoc_exc_root_pat = exc
                            assoc_errors_pat[mode] = f"{type(exc).__name__}: {exc}"
                            continue
                    else:
                        chain_miss("spatial_association")
                        meta["pattern_association_errors"] = assoc_errors_pat
                        if meta.get("official_api_error") in (None, "NA") and assoc_exc_root_pat is not None:
                            meta["official_api_error"] = (
                                f"{type(assoc_exc_root_pat).__name__}: {assoc_exc_root_pat}"
                            )
                else:
                    if not callable(spatial_association):
                        meta["pattern_association_skipped_reason"] = "spatial_association_not_found"
                    else:
                        meta["pattern_association_skipped_reason"] = "n_pattern_genes_lt2"

            meta["pattern_edges_rows"] = int(coexpr_pat.shape[0]) if isinstance(coexpr_pat, pd.DataFrame) else 0

            # ---- genemodules on pattern-gene universe ----
            if isinstance(coexpr_pat, pd.DataFrame) and not coexpr_pat.empty:
                n_genes_for_modules = int(getattr(adata_pat, "n_vars", 0))
                nmax_eff = min(int(modules_nmax), max(0, n_genes_for_modules - 1))
                if nmax_eff >= 2:
                    n_modules_eff = modules_n
                    if n_modules_eff is not None and (
                        n_modules_eff < 2 or n_modules_eff > (n_genes_for_modules - 1)
                    ):
                        n_modules_eff = None

                    genemodules = getattr(api, "genemodules", None)
                    if not callable(genemodules):
                        chain_miss("genemodules")
                        meta["pattern_modules_error"] = f"genemodules_not_found(entrypoint={entrypoint})"
                    else:
                        try:
                            ret = _call_api(
                                genemodules,
                                sg_pat,
                                nmax=int(nmax_eff),
                                use_grouped=False,
                                n_modules=n_modules_eff,
                            )
                            md = getattr(sg_pat, "module_dict", None)
                            if md is None:
                                md = _find_module_dict(ret)
                                if md is None:
                                    md = _find_module_dict(getattr(sg_pat, "adata", None))
                            module_dict_pat = md
                            chain_ok("genemodules")
                            meta["pattern_modules_nmax_used"] = int(nmax_eff)
                            meta["pattern_modules_n_used"] = (
                                int(n_modules_eff) if n_modules_eff is not None else None
                            )
                            try:
                                n_assigned = 0
                                n_mod = 0
                                if isinstance(md, dict):
                                    g2m = None
                                    for k in ("gene2mod", "gene_to_module", "gene2module", "gene_module", "gene_module_dict"):
                                        v = md.get(k)
                                        if isinstance(v, dict) and v:
                                            g2m = v
                                            break
                                    if g2m is None and md and all(not isinstance(v, dict) for v in md.values()):
                                        g2m = md
                                    if isinstance(g2m, dict):
                                        n_assigned = int(len(g2m))
                                        n_mod = int(len({str(x) for x in g2m.values()}))
                                label_nums: List[int] = []
                                if isinstance(g2m, dict):
                                    for lab in g2m.values():
                                        s = str(lab)
                                        if s.startswith("m") and s[1:].isdigit():
                                            label_nums.append(int(s[1:]))
                                        elif s.isdigit():
                                            label_nums.append(int(s))
                                label_min = min(label_nums) if label_nums else None
                                label_max = max(label_nums) if label_nums else None

                                meta["pattern_modules_genes_assigned_n"] = int(n_assigned)
                                meta["pattern_modules_n"] = int(n_mod)
                                meta["pattern_modules_label_min"] = label_min if label_min is not None else "NA"
                                meta["pattern_modules_label_max"] = label_max if label_max is not None else "NA"

                                eprint(
                                    "[INFO] official_api: genemodules succeeded. "
                                    f"n_modules={n_mod} module_genes_covered={n_assigned} "
                                    f"module_id_min={label_min if label_min is not None else 'NA'} "
                                    f"module_id_max={label_max if label_max is not None else 'NA'}"
                                )
                            except Exception:
                                pass
                        except BaseException as exc:
                            module_dict_pat = None
                            err = f"{type(exc).__name__}: {exc}"
                            meta["pattern_modules_error"] = err
                            if meta.get("official_api_error") in (None, "NA"):
                                meta["official_api_error"] = err
                            eprint(f"[WARN] official_api: genemodules failed; falling back to custom modules. error={err}")
                else:
                    meta["pattern_modules_skipped_reason"] = "n_pattern_genes_lt2"
            else:
                if not meta.get("pattern_modules_skipped_reason"):
                    meta["pattern_modules_skipped_reason"] = "pattern_association_empty"
        except BaseException as exc:
            module_dict_pat = None
            meta["pattern_modules_error"] = f"{type(exc).__name__}: {exc}"

    return sg, coexpr_full, coexpr_pat, module_dict_pat, pattern_genes_used, meta


def _extract_gene2label_from_module_dict(module_dict: Optional[object]) -> Dict[str, str]:
    if not module_dict:
        return {}

    if isinstance(module_dict, pd.DataFrame):
        gene_col = next((c for c in ["gene", "Gene"] if c in module_dict.columns), None)
        mod_col = next((c for c in ["module_id", "module", "Module"] if c in module_dict.columns), None)
        if gene_col and mod_col:
            return {
                str(g): str(m)
                for g, m in zip(
                    module_dict[gene_col].astype(str).tolist(),
                    module_dict[mod_col].tolist(),
                )
            }
        return {}

    key_candidates = (
        "gene2mod",
        "gene_to_module",
        "gene2module",
        "gene_module",
        "gene_module_dict",
    )

    def _find_gene2mod(d: object) -> Optional[Dict[str, object]]:
        if not isinstance(d, dict) or not d:
            return None
        for k in key_candidates:
            raw = d.get(k)
            if isinstance(raw, dict) and raw:
                return raw
        if d and all(not isinstance(v, dict) for v in d.values()):
            return d  # type: ignore[return-value]
        for v in d.values():
            hit = _find_gene2mod(v)
            if hit is not None:
                return hit
        return None

    gene2mod_raw = _find_gene2mod(module_dict)
    if not isinstance(gene2mod_raw, dict) or not gene2mod_raw:
        return {}
    return {str(k): str(v) for k, v in gene2mod_raw.items()}


def _label_to_module_id_map(gene2label: Dict[str, str]) -> Dict[str, int]:
    if not gene2label:
        return {}

    def _label_key(lab: str) -> Tuple[int, object]:
        s = str(lab)
        if s.startswith("m") and s[1:].isdigit():
            return (0, int(s[1:]))
        if s.isdigit():
            return (0, int(s))
        return (1, s)

    uniq_labels = sorted(set(gene2label.values()), key=_label_key)
    return {lab: i for i, lab in enumerate(uniq_labels, start=1)}


def build_modules_tsv_assigned(
    gene_universe: List[str],
    module_dict: Optional[object],
    *,
    include_unassigned: bool,
) -> pd.DataFrame:
    gene2label = _extract_gene2label_from_module_dict(module_dict)
    label_to_id = _label_to_module_id_map(gene2label)

    rows: List[Tuple[str, int]] = []
    for g in [str(x) for x in gene_universe]:
        lab = gene2label.get(g)
        if lab is None:
            if include_unassigned:
                rows.append((g, -1))
            continue
        mid = label_to_id.get(str(lab))
        if mid is None:
            if include_unassigned:
                rows.append((g, -1))
            continue
        rows.append((g, int(mid)))

    df = pd.DataFrame(rows, columns=["gene", "module_id"])
    if df.empty:
        return df
    return df.sort_values(["module_id", "gene"], ascending=[True, True], kind="mergesort").reset_index(drop=True)


def build_edges_all_tsv(coexpr: pd.DataFrame, *, n_obs: int) -> pd.DataFrame:
    cols = ["gene_a", "gene_b", "weight", "p_value", "fdr", "n_obs", "edge_type"]
    if coexpr is None or coexpr.empty:
        return pd.DataFrame(columns=cols)

    def _pick_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
        for c in candidates:
            if c in df.columns:
                return c
        lmap = {c.lower(): c for c in df.columns}
        for c in candidates:
            hit = lmap.get(c.lower())
            if hit is not None:
                return hit
        return None

    gene1_col = _pick_col(coexpr, ["gene_1", "gene1", "gene_a", "geneA", "gene"])
    gene2_col = _pick_col(coexpr, ["gene_2", "gene2", "gene_b", "geneB"])
    l_col = _pick_col(coexpr, ["L", "leeL", "global_L", "global_L_index"])
    p_col = _pick_col(coexpr, ["L.p_value", "L_p_value", "p_value", "pvalue", "P_value"])
    fdr_col = _pick_col(coexpr, ["L.FDR", "L_FDR", "fdr", "FDR", "q_value", "qvalue"])

    missing_required = [n for n, c in [("gene_1", gene1_col), ("gene_2", gene2_col), ("L", l_col)] if c is None]
    if missing_required:
        raise ValueError(f"SEAGAL co_expression missing required columns: {missing_required}")

    edges = pd.DataFrame(
        {
            "gene_a": coexpr[gene1_col].astype(str),
            "gene_b": coexpr[gene2_col].astype(str),
            "weight": pd.to_numeric(coexpr[l_col], errors="coerce"),
            "p_value": pd.to_numeric(coexpr[p_col], errors="coerce") if p_col else np.nan,
            "fdr": pd.to_numeric(coexpr[fdr_col], errors="coerce") if fdr_col else np.nan,
            "n_obs": int(n_obs),
            "edge_type": "seagal_global_L",
        }
    )

    edges = edges[edges["gene_a"] != edges["gene_b"]].copy()
    gene_a = edges[["gene_a", "gene_b"]].min(axis=1)
    gene_b = edges[["gene_a", "gene_b"]].max(axis=1)
    edges["gene_a"] = gene_a
    edges["gene_b"] = gene_b
    edges = edges.drop_duplicates(subset=["gene_a", "gene_b"])
    edges = edges[np.isfinite(edges["weight"])].copy()

    edges["_abs_weight"] = edges["weight"].abs()
    edges = (
        edges.sort_values(
            by=["_abs_weight", "gene_a", "gene_b"],
            ascending=[False, True, True],
            kind="mergesort",
        )
        .drop(columns=["_abs_weight"])
        .reset_index(drop=True)
    )

    return edges[cols]


def build_modules_tsv(gene_universe: List[str], module_dict: Optional[Dict[str, object]]) -> pd.DataFrame:
    modules = pd.DataFrame({"gene": [str(g) for g in gene_universe], "module_id": -1})
    if not module_dict:
        return modules

    gene2mod_raw: Optional[Dict[str, object]] = None
    if isinstance(module_dict, pd.DataFrame):
        gene_col = next((c for c in ["gene", "Gene"] if c in module_dict.columns), None)
        mod_col = next((c for c in ["module_id", "module", "Module"] if c in module_dict.columns), None)
        if gene_col and mod_col:
            gene2mod_raw = dict(
                zip(
                    module_dict[gene_col].astype(str).tolist(),
                    module_dict[mod_col].tolist(),
                )
            )
    elif isinstance(module_dict, dict):
        for key in [
            "gene2mod",
            "gene_to_module",
            "gene2module",
            "gene_module",
            "gene_module_dict",
        ]:
            raw = module_dict.get(key)
            if isinstance(raw, dict) and raw:
                gene2mod_raw = raw
                break
        if gene2mod_raw is None and module_dict and all(not isinstance(v, dict) for v in module_dict.values()):
            gene2mod_raw = module_dict  # type: ignore[assignment]

    if not isinstance(gene2mod_raw, dict) or not gene2mod_raw:
        return modules

    gene2label = {str(k): str(v) for k, v in gene2mod_raw.items()}

    def _label_key(lab: str) -> Tuple[int, object]:
        s = str(lab)
        if s.startswith("m") and s[1:].isdigit():
            return (0, int(s[1:]))
        if s.isdigit():
            return (0, int(s))
        return (1, s)

    uniq_labels = sorted(set(gene2label.values()), key=_label_key)
    label_to_id = {lab: i for i, lab in enumerate(uniq_labels, start=1)}

    modules["module_id"] = [
        int(label_to_id.get(gene2label.get(str(g)), -1)) for g in modules["gene"].tolist()
    ]
    return modules


def cluster_pattern_modules_from_weight_matrix(
    L_pat: np.ndarray,
    pattern_genes: List[str],
    *,
    nmax: int,
    n_modules: Optional[int],
) -> Optional[Dict[str, object]]:
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.metrics import silhouette_score

    genes = [str(g) for g in pattern_genes]
    n = int(len(genes))
    if n < 2:
        return None

    X = np.asarray(L_pat, dtype=np.float64)
    if X.shape != (n, n):
        raise ValueError(f"L_pat shape {X.shape} does not match pattern genes {n}.")
    X = X.copy()
    np.fill_diagonal(X, 0.0)
    X[~np.isfinite(X)] = 0.0

    nmax_eff = min(int(nmax), max(0, n - 1))
    if nmax_eff < 2:
        return None

    if n_modules is not None:
        nclust_opt = int(n_modules)
        if nclust_opt < 2 or nclust_opt > (n - 1):
            nclust_opt = 2
    else:
        range_n_clusters = list(range(2, nmax_eff + 1))
        sil_scores: List[float] = []
        for n_clusters in range_n_clusters:
            hc = AgglomerativeClustering(n_clusters=n_clusters)
            labels = hc.fit_predict(X)
            sil = float(silhouette_score(X, labels))
            sil_scores.append(sil)
            eprint(f"[INFO] modules silhouette n_clusters={n_clusters} score={sil:.6g}")
        sil_scores_r = sil_scores[::-1]
        i = len(sil_scores_r) - int(np.argmax(sil_scores_r)) - 1
        nclust_opt = int(range_n_clusters[i])

    eprint(f"[INFO] modules selected n_clusters={nclust_opt}")
    hc_opt = AgglomerativeClustering(n_clusters=nclust_opt)
    labels_opt = hc_opt.fit_predict(X)
    sil_opt = float(silhouette_score(X, labels_opt))

    labels_str = [f"m{int(c)}" for c in labels_opt.tolist()]
    gene2mod = dict(zip(genes, labels_str))
    mod2gene: Dict[str, List[str]] = {}
    for g, lab in gene2mod.items():
        mod2gene.setdefault(str(lab), []).append(g)

    return {
        "gene2mod": gene2mod,
        "mod2gene": mod2gene,
        "silhouette": sil_opt,
        "n_modules": int(nclust_opt),
    }

def _append_edges_all_from_seagal_coexpr(
    coexpr: pd.DataFrame,
    out_path: str,
    *,
    n_obs: int,
    edge_type: str,
) -> int:
    if not isinstance(coexpr, pd.DataFrame) or coexpr.empty:
        return 0

    if "gene_1" not in coexpr.columns or "gene_2" not in coexpr.columns or "L" not in coexpr.columns:
        raise ValueError(f"SEAGAL co_expression missing required columns: {list(coexpr.columns)}")

    p_col = "L.p_value" if "L.p_value" in coexpr.columns else None
    fdr_col = "L.FDR" if "L.FDR" in coexpr.columns else None

    out = pd.DataFrame(
        {
            "gene_a": coexpr["gene_1"].astype(str),
            "gene_b": coexpr["gene_2"].astype(str),
            "weight": pd.to_numeric(coexpr["L"], errors="coerce"),
            "p_value": pd.to_numeric(coexpr[p_col], errors="coerce") if p_col else 1.0,
            "fdr": pd.to_numeric(coexpr[fdr_col], errors="coerce") if fdr_col else 1.0,
            "n_obs": int(n_obs),
            "edge_type": str(edge_type),
        }
    )
    out = out[np.isfinite(out["weight"])].copy()
    out.to_csv(out_path, sep="\t", index=False, mode="a", header=False)
    return int(out.shape[0])


def _set_thread_env(n: int) -> None:
    n = int(max(1, n))
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = str(n)
    try:
        from threadpoolctl import threadpool_limits

        threadpool_limits(n)
    except Exception:
        pass


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data_dir", required=True, help="Xenium outs directory OR .h5ad file path")
    ap.add_argument("--outdir", required=True, help="Output directory (mounted from host)")
    ap.add_argument("--coord_file", default="", help="ROI polygon CSV (vertices).")
    ap.add_argument("--ncores", type=int, default=1, help="Parallel cores")
    ap.add_argument("--seed", type=int, default=1, help="Random seed")
    ap.add_argument("--max_cells", type=int, default=0, help="Downsample to at most this many cells (after ROI)")
    ap.add_argument("--top_genes", type=int, default=0, help="Keep top genes by mean (0 disables)")
    ap.add_argument("--grid_um", default="0", help="Grid binning step size (um; 0 disables)")
    ap.add_argument(
        "--modules_nmax",
        type=int,
        default=None,
        help="Max candidate module count (pattern-gene modules); overrides extra_args_json.modules_nmax when set",
    )
    ap.add_argument("--extra_args_json", default="", help="Extra JSON string or file path")
    args = ap.parse_args()

    _set_thread_env(int(args.ncores))

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    t0 = time.time()
    np.random.seed(args.seed)

    roi_csv = args.coord_file.strip()
    extra = parse_extra_args(args.extra_args_json)

    min_counts = _safe_int(extra.get("min_counts", 150), 150)
    min_cells = _safe_int(extra.get("min_cells", 10), 10)
    n_permutation = _safe_int(extra.get("n_permutation", 99), 99)
    permute_ratio = _safe_float(extra.get("permute_ratio", 0.2), 0.2)
    fdr_cutoff = _safe_float(extra.get("fdr_cutoff", 0.05), 0.05)
    l_cutoff = _safe_float(extra.get("l_cutoff", 0.1), 0.1)
    indep = _safe_bool(extra.get("indep", True), True)

    svg_i = extra.get("svg_I", extra.get("pattern_i"))
    svg_topk = extra.get("svg_topk", extra.get("pattern_topk"))
    if svg_i is not None:
        svg_i = _safe_float(svg_i, float("nan"))
        if not np.isfinite(svg_i):
            svg_i = None
    if svg_topk is not None:
        svg_topk = _safe_int(svg_topk, 0)
        if svg_topk < 1:
            svg_topk = None

    use_pattern_genes = _safe_bool(extra.get("use_pattern_genes", True), True)

    dense_max_obs = _safe_int(extra.get("seagal_dense_connectivities_max_obs", 20000), 20000)

    modules_enabled = _safe_bool(extra.get("modules_enabled", 1), True)
    modules_nmax = _safe_int(extra.get("modules_nmax", 6), 6)
    if getattr(args, "modules_nmax", None) is not None:
        modules_nmax = int(args.modules_nmax)
    modules_n = extra.get("modules_n", extra.get("modules_n_modules"))
    if modules_n is not None:
        modules_n = _safe_int(modules_n, 0)
        if modules_n < 2:
            modules_n = None
    modules_tsv_include_unassigned = _safe_bool(extra.get("modules_tsv_include_unassigned", 0), False)

    cells_df: Optional[pd.DataFrame] = None
    roi_applied = 0

    if os.path.isfile(args.data_dir) and args.data_dir.endswith(".h5ad"):
        eprint(f"[INFO] Loading from h5ad: {args.data_dir}")
        adata = load_from_h5ad(args.data_dir)
        adata = ensure_xy(adata)
        cells_df = pd.DataFrame(
            {
                "cell_id": adata.obs_names.astype(str),
                "x": adata.obs["x"].to_numpy(),
                "y": adata.obs["y"].to_numpy(),
            }
        )
        cells_df, roi_applied = apply_roi_and_downsample(cells_df, roi_csv, args.max_cells, args.seed)
        adata = adata[cells_df["cell_id"].astype(str).tolist()].copy()
    else:
        if not os.path.isdir(args.data_dir):
            raise FileNotFoundError(f"data_dir not found: {args.data_dir}")
        cells_pq = os.path.join(args.data_dir, "cells.parquet")
        tx_pq = os.path.join(args.data_dir, "transcripts.parquet")
        h5_path = os.path.join(args.data_dir, "cell_feature_matrix.h5")
        h5ad_path = os.path.join(args.data_dir, "adata_raw.h5ad")

        if os.path.exists(cells_pq) and os.path.exists(tx_pq):
            eprint("[INFO] Loading from Xenium parquet")
            cells_df = load_cells_df(cells_pq)
            cells_df, roi_applied = apply_roi_and_downsample(
                cells_df, roi_csv, args.max_cells, args.seed
            )
            adata = build_adata_from_transcripts(
                cells_df,
                tx_pq,
                exclude_prefix=(
                    "Unassigned",
                    "NegControl",
                    "Background",
                    "DeprecatedCodeword",
                    "SystemControl",
                    "Negative",
                ),
            )
        elif os.path.exists(h5_path):
            eprint("[INFO] Loading from cell_feature_matrix.h5")
            adata = load_from_h5(h5_path)
            if os.path.exists(cells_pq):
                cells_df = load_cells_df(cells_pq)
                cells_df, roi_applied = apply_roi_and_downsample(
                    cells_df, roi_csv, args.max_cells, args.seed
                )
                adata.obs_names = adata.obs_names.astype(str)
                cells_df = cells_df[cells_df["cell_id"].isin(adata.obs_names)].copy()
                if cells_df.empty:
                    raise ValueError("No overlapping cells between cells.parquet and h5 matrix")
                adata = adata[cells_df["cell_id"].tolist()].copy()
                adata.obs["x"] = cells_df["x"].to_numpy()
                adata.obs["y"] = cells_df["y"].to_numpy()
                adata.obsm["spatial"] = adata.obs[["x", "y"]].to_numpy()
            else:
                adata = ensure_xy(adata)
                cells_df = pd.DataFrame(
                    {
                        "cell_id": adata.obs_names.astype(str),
                        "x": adata.obs["x"].to_numpy(),
                        "y": adata.obs["y"].to_numpy(),
                    }
                )
                cells_df, roi_applied = apply_roi_and_downsample(
                    cells_df, roi_csv, args.max_cells, args.seed
                )
                adata = adata[cells_df["cell_id"].astype(str).tolist()].copy()
        elif os.path.exists(h5ad_path):
            eprint("[INFO] Loading from adata_raw.h5ad")
            adata = load_from_h5ad(h5ad_path)
            adata = ensure_xy(adata)
            cells_df = pd.DataFrame(
                {
                    "cell_id": adata.obs_names.astype(str),
                    "x": adata.obs["x"].to_numpy(),
                    "y": adata.obs["y"].to_numpy(),
                }
            )
            cells_df, roi_applied = apply_roi_and_downsample(
                cells_df, roi_csv, args.max_cells, args.seed
            )
            adata = adata[cells_df["cell_id"].astype(str).tolist()].copy()
        else:
            raise FileNotFoundError(
                "data_dir missing required files. Need one of: "
                "cells.parquet + transcripts.parquet, cell_feature_matrix.h5, or adata_raw.h5ad"
            )

    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("Loaded data has zero cells or genes after ROI/downsample.")

    if "spatial" not in adata.obsm:
        adata = ensure_xy(adata)

    adata = filter_top_genes(adata, args.top_genes)

    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("Loaded data has zero cells or genes after ROI/downsample/top_genes.")

    grid_um_raw = str(args.grid_um).strip()
    grid_um = _safe_float(grid_um_raw, 0.0)
    if not np.isfinite(grid_um) or grid_um <= 0:
        grid_um = 0.0

    n_obs_pre_grid = int(adata.n_obs)
    if grid_um > 0:
        eprint(f"[INFO] Grid binning enabled: grid_um={grid_um}")
        adata = bin_adata_to_grid(adata, grid_um)
        eprint(
            f"[INFO] Grid binning result: n_bins={adata.n_obs} (from n_obs_pre_grid={n_obs_pre_grid})"
        )

    eprint(
        f"[INFO] Prepared raw AnnData for SEAGAL CSV export: cells={adata.n_obs} genes={adata.n_vars} "
        f"roi_applied={roi_applied}"
    )

    count_csv, meta_csv = write_seagal_csv_inputs(adata, outdir)

    (
        sg,
        coexpr_full,
        coexpr_pat,
        module_dict_pat_official,
        pattern_genes_official,
        seagal_meta,
    ) = run_seagal_official_api(
        count_csv,
        meta_csv,
        min_counts=min_counts,
        min_cells=min_cells,
        svg_i=svg_i,
        svg_topk=svg_topk,
        use_pattern_genes=use_pattern_genes,
        n_permutation=n_permutation,
        permute_ratio=permute_ratio,
        fdr_cutoff=fdr_cutoff,
        l_cutoff=l_cutoff,
        indep=indep,
        dense_max_obs=dense_max_obs,
        modules_enabled=modules_enabled,
        modules_nmax=modules_nmax,
        modules_n=modules_n,
    )

    modules_df_official: Optional[pd.DataFrame] = None
    if modules_enabled and module_dict_pat_official and not pattern_genes_official:
        eprint("[WARN] genemodules returned module_dict but no pattern_genes; will fall back to custom_hclust")
        module_dict_pat_official = None
    if modules_enabled and module_dict_pat_official and pattern_genes_official:
        modules_df_official = build_modules_tsv_assigned(
            pattern_genes_official,
            module_dict_pat_official,
            include_unassigned=bool(modules_tsv_include_unassigned),
        )
        n_assigned_official = int((modules_df_official["module_id"] > 0).sum()) if not modules_df_official.empty else 0
        if n_assigned_official <= 0:
            err = seagal_meta.get("pattern_modules_error")
            if err:
                eprint(f"[WARN] genemodules returned no assignments; will fall back: {err}")
            else:
                eprint("[WARN] genemodules returned no assignments; will fall back to custom_hclust")
            module_dict_pat_official = None
            modules_df_official = None

    need_custom_modules = bool(modules_enabled and not module_dict_pat_official)

    adata_post = getattr(sg, "adata", None)
    n_obs_post = int(adata_post.n_obs) if adata_post is not None else 0
    n_vars_post = int(adata_post.n_vars) if adata_post is not None else 0

    # ---- edges_all.tsv (pattern genes only; full-gene all-to-all disabled) ----
    if adata_post is None or n_obs_post == 0 or n_vars_post == 0:
        raise ValueError("SEAGAL preprocessing returned empty AnnData (n_obs=0 or n_vars=0).")

    full_genes = sorted([str(g) for g in adata_post.var_names.astype(str).tolist()])

    edges_all_path = os.path.join(outdir, "edges_all.tsv")
    edges_all_fullgene_disabled = True
    edges_full_n = 0

    # ---- Pattern-gene modules (cluster on pattern submatrix) ----
    pattern_mask = None
    try:
        pattern_mask = adata_post.var.get("high_pattern_genes")
    except Exception:
        pattern_mask = None

    pattern_set: set[str] = set()
    if pattern_mask is not None:
        try:
            mask_arr = np.asarray(pattern_mask, dtype=bool)
            genes_orig = [str(g) for g in adata_post.var_names.astype(str).tolist()]
            pattern_set = {g for g, keep in zip(genes_orig, mask_arr.tolist()) if bool(keep)}
        except Exception:
            pattern_set = set()

    pattern_genes = sorted([g for g in full_genes if g in pattern_set])

    edges_pattern_n = 0
    module_dict_pat_custom: Optional[Dict[str, object]] = None
    modules_n_permutation = _safe_int(extra.get("seagal_modules_n_permutation", 0), 0)
    if pattern_genes and len(pattern_genes) >= 2:
        import importlib

        utils = importlib.import_module("seagal._utils")
        spatial_association = getattr(utils, "spatial_association", None)
        if not callable(spatial_association):
            raise ImportError("Cannot import seagal._utils.spatial_association")

        spatial_association(
            sg,
            grouped_only=False,
            use_pattern_genes=True,
            genes=None,
            n_permutation=int(modules_n_permutation),
            permute_ratio=float(permute_ratio),
            FDR_cutoff=float(fdr_cutoff),
            L_cutoff=float(l_cutoff),
            indep=bool(indep),
        )
        df_pat = getattr(sg, "co_expression", pd.DataFrame())

        with open(edges_all_path, "w", encoding="utf-8") as f:
            f.write("gene_a\tgene_b\tweight\tp_value\tfdr\tn_obs\tedge_type\n")
        if isinstance(df_pat, pd.DataFrame) and not df_pat.empty:
            edges_pattern_n = _append_edges_all_from_seagal_coexpr(
                df_pat, edges_all_path, n_obs=int(n_obs_post), edge_type="seagal_global_L_pattern"
            )

        if need_custom_modules and isinstance(df_pat, pd.DataFrame) and not df_pat.empty:
            df1 = df_pat[["gene_1", "gene_2", "L"]].copy()
            df2 = df1.rename(columns={"gene_1": "gene_2", "gene_2": "gene_1"})
            coexpr = pd.concat([df1, df2], ignore_index=True)
            mat = coexpr.pivot(index="gene_1", columns="gene_2", values="L")
            mat.index.name = None
            mat.columns.name = None
            mat = mat.reindex(index=pattern_genes, columns=pattern_genes).fillna(0.0)
            L_pat = mat.to_numpy(dtype=np.float64, copy=True)
            np.fill_diagonal(L_pat, 0.0)
            module_dict_pat_custom = cluster_pattern_modules_from_weight_matrix(
                L_pat,
                pattern_genes,
                nmax=int(modules_nmax),
                n_modules=modules_n,
            )
        try:
            setattr(sg, "co_expression", None)
        except Exception:
            pass
    else:
        with open(edges_all_path, "w", encoding="utf-8") as f:
            f.write("gene_a\tgene_b\tweight\tp_value\tfdr\tn_obs\tedge_type\n")

    modules_path = os.path.join(outdir, "modules.tsv")
    modules_source = "modules_disabled"
    modules_df = pd.DataFrame(columns=["gene", "module_id"])
    n_modules = 0
    n_module_genes = 0
    if modules_enabled:
        if modules_df_official is not None and not modules_df_official.empty:
            modules_source = "seagal_genemodules"
            modules_df = modules_df_official
        elif not module_dict_pat_official:
            if module_dict_pat_custom:
                modules_source = "custom_hclust_fallback"
                modules_df = build_modules_tsv_assigned(
                    pattern_genes,
                    module_dict_pat_custom,
                    include_unassigned=bool(modules_tsv_include_unassigned),
                )
            else:
                err = seagal_meta.get("pattern_modules_error")
                if err:
                    eprint(f"[WARN] genemodules unavailable; no fallback modules: {err}")
                modules_source = "custom_hclust_fallback" if pattern_genes else "pattern_genes_empty"
                modules_df = pd.DataFrame(columns=["gene", "module_id"])

        if not modules_df.empty:
            n_modules = int(modules_df.loc[modules_df["module_id"] > 0, "module_id"].nunique())
            n_module_genes = int((modules_df["module_id"] > 0).sum())
    else:
        modules_source = "modules_disabled"

    eprint(f"[INFO] modules_source={modules_source} n_modules={n_modules} n_module_genes={n_module_genes}")
    modules_df.to_csv(modules_path, sep="\t", index=False)

    module_k_best = int(n_modules)
    module_genes_covered = int(n_module_genes)

    # Optional: we no longer emit the full SEAGAL co_expression table by default (can be huge).
    out_tsv = os.path.join(outdir, "seagal_results.tsv")
    pd.DataFrame().to_csv(out_tsv, sep="\t", index=False)
    out_pat_tsv = os.path.join(outdir, "seagal_results_pattern.tsv")
    pd.DataFrame().to_csv(out_pat_tsv, sep="\t", index=False)

    n_sig = 0
    n_sig_full = 0
    n_sig_pattern = 0

    stats = {
        "n_cells": int(n_obs_post),
        "n_genes": int(n_vars_post),
        "n_sig": int(n_sig),
        "n_sig_full": int(n_sig_full),
        "n_sig_pattern": int(n_sig_pattern),
        "n_pattern_genes": int(len(pattern_genes)),
        "edges_all_fullgene_disabled": int(bool(edges_all_fullgene_disabled)),
        "edges_outfile": os.path.basename(edges_all_path),
        "edges_gene_set": "pattern_genes",
        "edges_full_n": 0,
        "edges_pattern_n": int(edges_pattern_n),
        "seagal_modules_n_permutation": int(modules_n_permutation),
        "roi_applied": int(roi_applied),
        "max_cells": int(args.max_cells),
        "top_genes": int(args.top_genes),
        "grid_um": float(grid_um),
        "n_obs_pre_grid": int(n_obs_pre_grid),
        "seed": int(args.seed),
        "modules_enabled": int(bool(modules_enabled)),
        "modules_source": str(modules_source),
        "n_modules": int(n_modules),
        "module_k_best": int(module_k_best),
        "n_module_genes": int(n_module_genes),
    }
    for k, v in seagal_meta.items():
        stats[f"seagal_{k}"] = v
    with open(os.path.join(outdir, "run_stats.tsv"), "w", encoding="utf-8") as f:
        for k, v in stats.items():
            f.write(f"{k}\t{v}\n")

    meta_json_path = os.path.join(outdir, "meta.json")
    try:
        seagal_entrypoint = str(seagal_meta.get("seagal_entrypoint", "seagal"))
        api_chain_effective = list(seagal_meta.get("api_chain_effective", []))
        api_chain_missing = list(seagal_meta.get("api_chain_missing", []))
        meta_payload = {
            "method": "seagal",
            "modules_source": str(modules_source),
            "n_modules": int(n_modules),
            "n_module_genes": int(n_module_genes),
            "edges_all_fullgene_disabled": True,
            "edges": {
                "outfile": "edges_all.tsv",
                "gene_set": "pattern_genes",
                "edge_type": "seagal_global_L_pattern",
                "edges_rows": int(edges_pattern_n),
            },
            "seagal_entrypoint": seagal_entrypoint,
            "state_object": seagal_meta.get("state_object", "unknown"),
            "api_chain_effective": api_chain_effective,
            "api_chain_missing": api_chain_missing,
            "api_chain": [f"{seagal_entrypoint}.{fn}" for fn in api_chain_effective],
            "seagal_version": _try_importlib_version("seagal"),
            "seagal_meta": seagal_meta,
            "full_edges": {
                "n_genes": int(n_vars_post),
                "use_pattern_genes": False,
                "disabled": True,
                "edges_rows": 0,
            },
            "pattern_modules": {
                "n_pattern_genes": int(len(pattern_genes)),
                "modules_enabled": bool(modules_enabled),
                "modules_source": str(modules_source),
                "n_module_genes": int(n_module_genes),
                "module_genes_covered": int(module_genes_covered),
                "n_modules": int(n_modules),
            },
            "params": {
                "min_counts": int(min_counts),
                "min_cells": int(min_cells),
                "grid_um": float(grid_um),
                "n_obs_pre_grid": int(n_obs_pre_grid),
                "edges_all_fullgene_disabled": True,
                "edges_outfile": "edges_all.tsv",
                "seagal_modules_n_permutation": int(modules_n_permutation),
                "seagal_dense_connectivities_max_obs": int(dense_max_obs),
                "svg_I": svg_i,
                "svg_topk": svg_topk,
                "use_pattern_genes": bool(use_pattern_genes),
                "n_permutation": int(n_permutation),
                "permute_ratio": float(permute_ratio),
                "fdr_cutoff": float(fdr_cutoff),
                "l_cutoff": float(l_cutoff),
                "indep": bool(indep),
                "modules_enabled": bool(modules_enabled),
                "modules_nmax": int(modules_nmax),
                "modules_n": int(modules_n) if modules_n is not None else None,
            },
            "counts": {
                "n_obs_post_process": int(n_obs_post),
                "n_vars_post_process": int(n_vars_post),
                "sag_pairs_n": int(edges_pattern_n),
                "sag_pairs_full_n": 0,
                "sag_pairs_pattern_n": int(edges_pattern_n),
                "n_pattern_genes": int(len(pattern_genes)),
                "modules_source": str(modules_source),
                "n_module_genes": int(n_module_genes),
                "module_genes_covered": int(module_genes_covered),
                "modules_n": int(n_modules),
            },
            "inputs": {
                "count_csv": os.path.relpath(count_csv, outdir),
                "meta_csv": os.path.relpath(meta_csv, outdir),
            },
            "outputs": {
                "seagal_results_tsv": "seagal_results.tsv",
                "seagal_results_pattern_tsv": "seagal_results_pattern.tsv",
                "edges_all_tsv": "edges_all.tsv",
                "modules_tsv": "modules.tsv",
            },
        }
        with open(meta_json_path, "w", encoding="utf-8") as f:
            json.dump(meta_payload, f, indent=2, ensure_ascii=False)
    except Exception as exc:
        eprint(f"[WARN] meta.json write failed: {exc}")

    wall = time.time() - t0
    eprint(f"[DONE] SEAGAL finished. wall_sec={wall:.2f}")
    eprint(f"[DONE] wrote: {out_tsv}")
    eprint(f"[DONE] wrote: {out_pat_tsv}")
    eprint(f"[DONE] wrote: {edges_all_path}")
    eprint(f"[DONE] wrote: {modules_path}")
    eprint(f"[DONE] wrote: {meta_json_path}")


if __name__ == "__main__":
    main()
