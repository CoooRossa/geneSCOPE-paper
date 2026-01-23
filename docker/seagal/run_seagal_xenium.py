#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SEAGAL runner (official CSV-format tutorial style).

Hard constraints for this runner:
- This file must not define any custom functions or classes (no `def`, no `class`).
- Only persistent outputs under --outdir are: all_edges.tsv and modules.tsv.
- Only call SEAGAL via public entrypoint: seagal.seagal (SEAGAL, spatial_pattern_genes, spatial_association, genemodules).
- Xenium support is via converting to SEAGAL CSV (count.csv/meta.csv) in a temporary directory.
"""

import argparse
import os
import sys
import tempfile

import numpy as np
import pandas as pd

os.environ.setdefault("MPLBACKEND", "Agg")

try:
    import multiprocessing as mp

    mp.set_start_method("fork", force=True)
except Exception:
    pass

ap = argparse.ArgumentParser()
ap.add_argument("--outdir", required=True, help="Output directory (writes only all_edges.tsv and modules.tsv).")

# Either provide count/meta directly, or provide Xenium/.h5ad via data_dir
ap.add_argument("--count_csv", default="", help="SEAGAL count.csv (rows=spots, cols=genes).")
ap.add_argument("--meta_csv", default="", help="SEAGAL meta.csv (index aligned with count; must include x,y).")
ap.add_argument("--data_dir", default="", help="Xenium outs directory OR .h5ad file path (used if count/meta not provided).")

# Xenium / h5ad preprocessing knobs
ap.add_argument("--coord_file", default="", help="Optional ROI polygon vertices CSV (x,y columns).")
ap.add_argument("--seed", type=int, default=1, help="Random seed.")
ap.add_argument("--max_cells", type=int, default=0, help="Downsample to at most this many cells (after ROI). 0 disables.")
ap.add_argument("--top_genes", type=int, default=0, help="Keep top genes by mean (0 disables). Recommended for Xenium.")
ap.add_argument(
    "--grid_um",
    type=float,
    default=100.0,
    help="Pseudo-spot grid size (microns). Default 100 (Visium-like pitch).",
)
ap.add_argument("--ncores", type=int, default=1, help="Set thread env (OMP/BLAS/MKL).")

# SEAGAL official API knobs
ap.add_argument("--svg_I", type=float, default=np.nan, help="Moran's I threshold (optional). NaN disables.")
ap.add_argument("--svg_topk", type=int, default=0, help="TopK Moran's I genes (optional). 0 disables.")
ap.add_argument("--n_permutation", type=int, default=99, help="Permutations for Global_L p-value.")
ap.add_argument("--permute_ratio", type=float, default=0.2, help="permute_ratio/percent used in permutation.")
ap.add_argument("--fdr_cutoff", type=float, default=0.05, help="FDR cutoff.")
ap.add_argument("--l_cutoff", type=float, default=0.1, help="Effect-size cutoff on L.")
ap.add_argument("--indep", type=int, default=1, choices=[0, 1], help="Use BH indep=True/False (1/0).")

ap.add_argument("--modules_nmax", type=int, default=6, help="Max candidate module count for silhouette search.")
ap.add_argument("--modules_n", type=int, default=0, help="Force number of modules if >=2; else auto.")

# Back-compat no-ops (accepted but intentionally ignored)
ap.add_argument("--extra_args_json", default="", help="(ignored) Back-compat with older wrappers.")
ap.add_argument("--repeat", type=int, default=1, help="(ignored) Back-compat with older wrappers.")
ap.add_argument("--sample_sec", type=float, default=0.0, help="(ignored) Back-compat with older wrappers.")

args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)

n = int(max(1, args.ncores))
for k in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[k] = str(n)

rng = np.random.default_rng(args.seed)

if args.extra_args_json.strip():
    print(
        "[WARN] --extra_args_json is accepted for back-compat but ignored; use explicit flags (e.g. --svg_topk).",
        file=sys.stderr,
    )

count_csv = args.count_csv.strip()
meta_csv = args.meta_csv.strip()
tmp_ctx = None

# --------------------------
# Path A: user-provided count/meta (pure SEAGAL CSV tutorial route)
# --------------------------
if count_csv and meta_csv:
    if not os.path.exists(count_csv):
        print(f"[ERROR] --count_csv not found: {count_csv}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(meta_csv):
        print(f"[ERROR] --meta_csv not found: {meta_csv}", file=sys.stderr)
        sys.exit(1)

    # SEAGAL load_raw enforces identical index AND order; pre-check here for clearer errors.
    try:
        _c = pd.read_csv(count_csv, index_col=0, nrows=5)
        _m = pd.read_csv(meta_csv, index_col=0, nrows=5)
        _ = (_c, _m)
    except Exception as e:
        print(f"[ERROR] failed to read count/meta csv: {type(e).__name__}: {e}", file=sys.stderr)
        sys.exit(1)

else:
    # --------------------------
    # Path B: Xenium outs dir or .h5ad -> build temporary count/meta
    # --------------------------
    if not args.data_dir:
        print("[ERROR] Provide either (--count_csv AND --meta_csv) OR --data_dir.", file=sys.stderr)
        sys.exit(1)

    data_dir = args.data_dir.strip()

    adata = None
    cells = None

    if os.path.isfile(data_dir) and data_dir.endswith(".h5ad"):
        try:
            import anndata as ad
        except Exception as e:
            print(f"[ERROR] import anndata failed: {type(e).__name__}: {e}", file=sys.stderr)
            sys.exit(1)

        adata = ad.read_h5ad(data_dir)
        if "x" not in adata.obs.columns or "y" not in adata.obs.columns:
            if "spatial" in adata.obsm:
                adata.obs["x"] = adata.obsm["spatial"][:, 0]
                adata.obs["y"] = adata.obsm["spatial"][:, 1]
            else:
                print("[ERROR] h5ad missing obs['x','y'] and obsm['spatial'].", file=sys.stderr)
                sys.exit(1)

        cells = pd.DataFrame(
            {"cell_id": adata.obs_names.astype(str), "x": adata.obs["x"].to_numpy(), "y": adata.obs["y"].to_numpy()}
        ).set_index("cell_id")

    else:
        if not os.path.isdir(data_dir):
            print(f"[ERROR] --data_dir not found or not a directory: {data_dir}", file=sys.stderr)
            sys.exit(1)

        # Load expression matrix
        try:
            import scanpy as sc
        except Exception as e:
            print(f"[ERROR] import scanpy failed: {type(e).__name__}: {e}", file=sys.stderr)
            sys.exit(1)

        mtx_dir = os.path.join(data_dir, "cell_feature_matrix")
        h5_path = os.path.join(data_dir, "cell_feature_matrix.h5")

        if os.path.isdir(mtx_dir):
            adata = sc.read_10x_mtx(mtx_dir, var_names="gene_symbols", make_unique=True)
        elif os.path.exists(h5_path):
            adata = sc.read_10x_h5(h5_path, gex_only=False)
            if "feature_types" in adata.var.columns:
                mask = adata.var["feature_types"] == "Gene Expression"
                if mask.any():
                    adata = adata[:, mask].copy()
        else:
            print("[ERROR] Xenium data_dir must contain cell_feature_matrix/ or cell_feature_matrix.h5.", file=sys.stderr)
            sys.exit(1)

        adata.obs_names = adata.obs_names.astype(str)
        adata.var_names = adata.var_names.astype(str)

        # Load cells coordinates
        cells_csv_gz = os.path.join(data_dir, "cells.csv.gz")
        cells_csv = os.path.join(data_dir, "cells.csv")
        cells_pq = os.path.join(data_dir, "cells.parquet")

        if os.path.exists(cells_csv_gz):
            cells_df = pd.read_csv(cells_csv_gz)
        elif os.path.exists(cells_csv):
            cells_df = pd.read_csv(cells_csv)
        elif os.path.exists(cells_pq):
            try:
                import pyarrow.parquet as pq
            except Exception as e:
                print(f"[ERROR] import pyarrow failed (needed for cells.parquet): {type(e).__name__}: {e}", file=sys.stderr)
                sys.exit(1)
            cells_df = pq.read_table(cells_pq).to_pandas()
        else:
            print("[ERROR] Xenium data_dir must contain cells.csv(.gz) or cells.parquet.", file=sys.stderr)
            sys.exit(1)

        cols_lc = {c.lower(): c for c in cells_df.columns}
        cid_col = cols_lc.get("cell_id") or cols_lc.get("barcode") or cols_lc.get("id")
        x_col = cols_lc.get("x_centroid") or cols_lc.get("x") or cols_lc.get("cell_centroid_x")
        y_col = cols_lc.get("y_centroid") or cols_lc.get("y") or cols_lc.get("cell_centroid_y")
        if cid_col is None or x_col is None or y_col is None:
            print(f"[ERROR] cells table missing required columns. cols={list(cells_df.columns)[:80]}", file=sys.stderr)
            sys.exit(1)

        cells_df = cells_df[[cid_col, x_col, y_col]].copy()
        cells_df.columns = ["cell_id", "x", "y"]
        cells_df["cell_id"] = cells_df["cell_id"].astype(str)
        cells_df["x"] = pd.to_numeric(cells_df["x"], errors="coerce")
        cells_df["y"] = pd.to_numeric(cells_df["y"], errors="coerce")
        cells_df = cells_df.dropna()
        cells = cells_df.set_index("cell_id")

    # Intersect cells with expression matrix
    common = pd.Index(adata.obs_names).intersection(pd.Index(cells.index.astype(str)))
    if len(common) == 0:
        print("[ERROR] No overlap between expression matrix barcodes and cells coordinate table.", file=sys.stderr)
        sys.exit(1)

    adata = adata[common].copy()
    cells = cells.loc[adata.obs_names].copy()

    # Optional ROI clip
    if args.coord_file.strip():
        roi_path = args.coord_file.strip()
        if not os.path.exists(roi_path):
            print(f"[ERROR] --coord_file not found: {roi_path}", file=sys.stderr)
            sys.exit(1)

        roi = pd.read_csv(roi_path, comment="#")
        roi_cols = {c.lower(): c for c in roi.columns}
        rx = roi_cols.get("x") or roi.columns[0]
        ry = roi_cols.get("y") or roi.columns[1]
        poly = roi[[rx, ry]].to_numpy(dtype=float, copy=False)
        if poly.shape[0] < 3:
            print("[ERROR] ROI polygon needs >=3 vertices.", file=sys.stderr)
            sys.exit(1)

        try:
            from matplotlib.path import Path
        except Exception as e:
            print(f"[ERROR] import matplotlib.path.Path failed: {type(e).__name__}: {e}", file=sys.stderr)
            sys.exit(1)

        pts = cells[["x", "y"]].to_numpy(dtype=float, copy=False)
        mask = Path(poly).contains_points(pts)
        if int(mask.sum()) == 0:
            print("[ERROR] ROI clip removed all cells.", file=sys.stderr)
            sys.exit(1)

        adata = adata[mask].copy()
        cells = cells.iloc[np.where(mask)[0]].copy()

    # Optional downsample
    if args.max_cells > 0 and adata.n_obs > args.max_cells:
        idx = rng.choice(adata.n_obs, size=args.max_cells, replace=False)
        adata = adata[idx].copy()
        cells = cells.iloc[idx].copy()

    # Attach coords
    adata.obs["x"] = cells["x"].to_numpy()
    adata.obs["y"] = cells["y"].to_numpy()
    adata.obsm["spatial"] = adata.obs[["x", "y"]].to_numpy()

    # Optional top genes
    if args.top_genes > 0 and adata.n_vars > args.top_genes:
        try:
            from scipy import sparse
        except Exception as e:
            print(f"[ERROR] import scipy failed (needed for top_genes): {type(e).__name__}: {e}", file=sys.stderr)
            sys.exit(1)

        X = adata.X
        if sparse.issparse(X):
            means = np.asarray(X.mean(axis=0)).ravel()
        else:
            means = np.asarray(X).mean(axis=0)
        top_idx = np.argsort(means)[::-1][: int(args.top_genes)]
        mask = np.zeros(adata.n_vars, dtype=bool)
        mask[top_idx] = True
        adata = adata[:, mask].copy()

    # Grid binning to pseudo-spots
    if args.grid_um > 0:
        try:
            from scipy import sparse
        except Exception as e:
            print(f"[ERROR] import scipy failed (needed for grid binning): {type(e).__name__}: {e}", file=sys.stderr)
            sys.exit(1)

        x = pd.to_numeric(adata.obs["x"], errors="coerce").to_numpy(dtype=float, copy=False)
        y = pd.to_numeric(adata.obs["y"], errors="coerce").to_numpy(dtype=float, copy=False)
        ok = np.isfinite(x) & np.isfinite(y)
        if not ok.all():
            adata = adata[ok].copy()
            x = x[ok]
            y = y[ok]
            if adata.n_obs == 0:
                print("[ERROR] All observations have invalid x/y.", file=sys.stderr)
                sys.exit(1)

        x0 = float(np.min(x))
        y0 = float(np.min(y))
        gx = np.floor((x - x0) / float(args.grid_um)).astype(np.int64, copy=False)
        gy = np.floor((y - y0) / float(args.grid_um)).astype(np.int64, copy=False)

        pairs = np.stack([gx, gy], axis=1)
        uniq_pairs, inv = np.unique(pairs, axis=0, return_inverse=True)
        inv = inv.astype(np.int32, copy=False)
        n_bins = int(uniq_pairs.shape[0])
        n_cells = int(adata.n_obs)

        X = adata.X
        if not sparse.issparse(X):
            X = sparse.csr_matrix(X)
        else:
            X = X.tocsr()

        rows = inv
        cols = np.arange(n_cells, dtype=np.int32)
        G = sparse.csr_matrix((np.ones(n_cells, dtype=np.int8), (rows, cols)), shape=(n_bins, n_cells))
        Xb = (G @ X).tocsr()
        Xb.sum_duplicates()

        counts = np.bincount(rows, minlength=n_bins).astype(np.int64, copy=False)
        x_mean = np.bincount(rows, weights=x, minlength=n_bins) / np.maximum(counts, 1)
        y_mean = np.bincount(rows, weights=y, minlength=n_bins) / np.maximum(counts, 1)

        bin_ids = [f"bin_{int(a)}_{int(b)}" for a, b in uniq_pairs.tolist()]

        # fail-fast for huge CSV
        n_genes = int(Xb.shape[1])
        if (n_bins * n_genes > 30_000_000) and (args.top_genes <= 0):
            print(
                f"[ERROR] count.csv would be too large (bins={n_bins}, genes={n_genes}). "
                "Set --top_genes (e.g. 3000-6000) and/or increase --grid_um.",
                file=sys.stderr,
            )
            sys.exit(1)

        tmp_ctx = tempfile.TemporaryDirectory(prefix="seagal_csv_")
        count_csv = os.path.join(tmp_ctx.name, "count.csv")
        meta_csv = os.path.join(tmp_ctx.name, "meta.csv")

        pd.DataFrame(Xb.toarray(), index=pd.Index(bin_ids, dtype=str), columns=adata.var_names.astype(str)).to_csv(
            count_csv
        )
        pd.DataFrame({"x": x_mean, "y": y_mean}, index=pd.Index(bin_ids, dtype=str)).to_csv(meta_csv)

    else:
        # No binning: use cells as spots (still needs fail-fast in Xenium scale)
        tmp_ctx = tempfile.TemporaryDirectory(prefix="seagal_csv_")
        count_csv = os.path.join(tmp_ctx.name, "count.csv")
        meta_csv = os.path.join(tmp_ctx.name, "meta.csv")

        try:
            from scipy import sparse
        except Exception:
            sparse = None

        X = adata.X
        if sparse is not None and sparse.issparse(X):
            if (adata.n_obs * adata.n_vars > 30_000_000) and (args.top_genes <= 0):
                print(
                    f"[ERROR] count.csv would be too large (obs={adata.n_obs}, genes={adata.n_vars}). "
                    "Set --top_genes and/or enable --grid_um.",
                    file=sys.stderr,
                )
                sys.exit(1)
            X_arr = X.toarray()
        else:
            X_arr = np.asarray(X)

        pd.DataFrame(X_arr, index=adata.obs_names.astype(str), columns=adata.var_names.astype(str)).to_csv(count_csv)
        pd.DataFrame(
            {"x": adata.obs["x"].to_numpy(), "y": adata.obs["y"].to_numpy()}, index=adata.obs_names.astype(str)
        ).to_csv(meta_csv)

# --------------------------
# SEAGAL official CSV route
# --------------------------
try:
    import scipy

    if not hasattr(scipy, "inf"):
        scipy.inf = np.inf
except Exception:
    pass

try:
    from seagal.seagal import SEAGAL, spatial_pattern_genes, spatial_association, genemodules
except Exception as e:
    print(f"[ERROR] import seagal.seagal failed: {type(e).__name__}: {e}", file=sys.stderr)
    sys.exit(1)

sg = SEAGAL(count_path=count_csv, meta_path=meta_csv, visium_path="__no_visium__")

svg_I = None if not np.isfinite(float(args.svg_I)) else float(args.svg_I)
svg_topk = None if int(args.svg_topk) < 1 else int(args.svg_topk)

if svg_I is None and svg_topk is None:
    print("[ERROR] Provide --svg_I or --svg_topk (SEAGAL spatial_pattern_genes requires one).", file=sys.stderr)
    sys.exit(2)

if svg_topk is not None:
    try:
        n_avail = int(getattr(getattr(sg, "adata", None), "n_vars", 0))
    except Exception:
        n_avail = 0
    if n_avail <= 0:
        print("[ERROR] SEAGAL produced empty adata (0 genes); cannot run spatial_pattern_genes.", file=sys.stderr)
        sys.exit(1)
    if svg_topk > n_avail:
        print(
            f"[WARN] --svg_topk={svg_topk} exceeds available genes after filtering (n_genes={n_avail}); "
            f"clamping to {n_avail}.",
            file=sys.stderr,
        )
        svg_topk = n_avail

spatial_pattern_genes(sg, I=svg_I, topK=svg_topk)

spatial_association(
    sg,
    grouped_only=False,
    use_pattern_genes=True,
    genes=None,
    n_permutation=int(args.n_permutation),
    permute_ratio=float(args.permute_ratio),
    FDR_cutoff=float(args.fdr_cutoff),
    L_cutoff=float(args.l_cutoff),
    indep=bool(int(args.indep)),
)

if getattr(sg, "co_expression", None) is None or not isinstance(sg.co_expression, pd.DataFrame) or sg.co_expression.empty:
    print("[ERROR] SEAGAL produced empty sg.co_expression; cannot write all_edges.tsv", file=sys.stderr)
    sys.exit(1)

coexpr = sg.co_expression.copy()
cols_map = {str(c).strip().lower(): c for c in coexpr.columns}

gene_a_col = None
for cand in ("gene_a", "gene_1", "gene1", "genea"):
    if cand in cols_map:
        gene_a_col = cols_map[cand]
        break
gene_b_col = None
for cand in ("gene_b", "gene_2", "gene2", "geneb"):
    if cand in cols_map:
        gene_b_col = cols_map[cand]
        break
weight_col = None
for cand in ("weight", "l"):
    if cand in cols_map:
        weight_col = cols_map[cand]
        break

if gene_a_col is None or gene_b_col is None or weight_col is None:
    print(f"[ERROR] SEAGAL co_expression missing required columns. cols={list(coexpr.columns)}", file=sys.stderr)
    sys.exit(1)

rename_map = {}
if gene_a_col != "gene_a":
    rename_map[gene_a_col] = "gene_a"
if gene_b_col != "gene_b":
    rename_map[gene_b_col] = "gene_b"
if weight_col != "weight":
    rename_map[weight_col] = "weight"

edges_df = coexpr.rename(columns=rename_map)
edges_df["gene_a"] = edges_df["gene_a"].astype(str)
edges_df["gene_b"] = edges_df["gene_b"].astype(str)
edges_df["weight"] = pd.to_numeric(edges_df["weight"], errors="coerce")
first_cols = ["gene_a", "gene_b", "weight"]
rest_cols = [c for c in edges_df.columns if c not in first_cols]
edges_df = edges_df[first_cols + rest_cols]

edges_path = os.path.join(args.outdir, "all_edges.tsv")
edges_df.to_csv(edges_path, sep="\t", index=False)

modules_n = None if int(args.modules_n) < 2 else int(args.modules_n)
genemodules(sg, nmax=int(args.modules_nmax), use_grouped=False, n_modules=modules_n)

if not hasattr(sg, "module_dict") or "gene2mod" not in sg.module_dict:
    print("[ERROR] SEAGAL did not produce module_dict['gene2mod']; cannot write modules.tsv", file=sys.stderr)
    sys.exit(1)

gene2mod = sg.module_dict["gene2mod"]
module_df = pd.DataFrame({"gene": list(gene2mod.keys()), "module_id": list(gene2mod.values())})
module_df["gene"] = module_df["gene"].astype(str)
module_id_raw = module_df["module_id"].astype(str)
module_id_is_int = module_id_raw.str.fullmatch(r"-?[0-9]+")
module_id_digits = module_id_raw.str.extract(r"^[^0-9]*([0-9]+)$", expand=False)
module_df["module_id"] = module_id_raw.where(module_id_is_int | module_id_digits.isna(), module_id_digits)
module_df = module_df[["gene", "module_id"]]

module_path = os.path.join(args.outdir, "modules.tsv")
module_df.to_csv(module_path, sep="\t", index=False)

print(f"[OK] Wrote: {edges_path}")
print(f"[OK] Wrote: {module_path}")
