#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os

import hotspot
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Hotspot v1+ minimal spatial workflow (AnnData route). "
            "Writes only modules.tsv and all_edges.tsv."
        )
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--xenium-outs", dest="xenium_outs", type=str)
    input_group.add_argument("--xenium_outs", dest="xenium_outs", type=str)
    input_group.add_argument("--h5ad", dest="h5ad", type=str)
    input_group.add_argument("--data_dir", dest="data_dir", type=str)

    parser.add_argument("--outdir", type=str, required=True, help="Output directory")
    parser.add_argument("--coord_file", type=str, default="")
    parser.add_argument("--roi_csv", type=str, default="")
    parser.add_argument("--model", type=str, default="bernoulli")
    parser.add_argument("--n-neighbors", dest="n_neighbors", type=int, default=300)
    parser.add_argument("--n_neighbors", dest="n_neighbors", type=int)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--threads", dest="jobs", type=int)
    parser.add_argument("--ncores", dest="jobs", type=int)
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--dataset_id", type=str, default="")
    parser.add_argument("--roi_id", type=str, default="")
    parser.add_argument("--max_cells", type=int, default=0)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--fdr_autocorr", type=float, default=0.05)
    parser.add_argument("--min_gene_threshold", type=int, default=20)
    parser.add_argument("--core_only", nargs="?", const="1", default="1")
    parser.add_argument("--stats_file", type=str, default="")

    args = parser.parse_args()

    if int(args.repeat) != 1:
        raise ValueError("--repeat must be 1 (this container only writes modules.tsv and all_edges.tsv).")

    outdir_root = args.outdir
    if os.path.basename(outdir_root).startswith("repeat_"):
        repdir = outdir_root
    else:
        repdir = os.path.join(outdir_root, "repeat_001")
    os.makedirs(repdir, exist_ok=True)
    module_path = os.path.join(repdir, "modules.tsv")
    edges_path = os.path.join(repdir, "all_edges.tsv")

    src_h5ad = args.h5ad
    src_outs = args.xenium_outs
    if src_h5ad is None and src_outs is None and args.data_dir is not None:
        if str(args.data_dir).lower().endswith(".h5ad"):
            src_h5ad = args.data_dir
        else:
            src_outs = args.data_dir

    if src_h5ad is not None:
        adata = sc.read_h5ad(src_h5ad)
        if "spatial" not in adata.obsm:
            raise ValueError("Input .h5ad missing obsm['spatial'].")
        if "counts" not in adata.layers:
            adata.layers["counts"] = adata.X
    else:
        outs_dir = src_outs
        matrix_path = os.path.join(outs_dir, "cell_feature_matrix.h5")
        cells_path = os.path.join(outs_dir, "cells.csv.gz")
        if not os.path.exists(matrix_path):
            raise FileNotFoundError(matrix_path)
        if not os.path.exists(cells_path):
            raise FileNotFoundError(cells_path)

        adata = sc.read_10x_h5(matrix_path)
        adata.var_names_make_unique()

        try:
            cells_df = pd.read_csv(
                cells_path,
                usecols=["cell_id", "x_centroid", "y_centroid"],
            )
        except ValueError as exc:
            raise ValueError(
                "cells.csv.gz must contain columns: cell_id, x_centroid, y_centroid"
            ) from exc

        cells_df["cell_id"] = cells_df["cell_id"].astype(str)
        cells_df = cells_df.set_index("cell_id")

        adata.obs_names = adata.obs_names.astype(str)
        common = adata.obs_names[adata.obs_names.isin(cells_df.index)]
        if common.size == 0:
            raise ValueError("No shared cell_id between matrix and cells.csv.gz.")

        adata = adata[common].copy()
        adata.obsm["spatial"] = cells_df.loc[
            common, ["x_centroid", "y_centroid"]
        ].to_numpy(dtype=float)
        adata.layers["counts"] = adata.X

    roi_csv = (args.coord_file or args.roi_csv or "").strip()
    if roi_csv:
        if not os.path.exists(roi_csv):
            raise FileNotFoundError(roi_csv)

        roi_df = pd.read_csv(roi_csv, comment="#")
        if roi_df.empty:
            raise ValueError(f"ROI CSV has no data rows after skipping comments: {roi_csv}")

        roi_cols = {c.strip().lower(): c for c in roi_df.columns}
        xcol = roi_cols.get("x")
        ycol = roi_cols.get("y")
        if xcol is None or ycol is None:
            if roi_df.shape[1] < 2:
                raise ValueError(
                    f"ROI CSV must have >=2 columns (or x/y columns): {roi_csv} cols={list(roi_df.columns)}"
                )
            xcol, ycol = roi_df.columns[:2]

        xy = roi_df[[xcol, ycol]].copy()
        xy[xcol] = pd.to_numeric(xy[xcol], errors="coerce")
        xy[ycol] = pd.to_numeric(xy[ycol], errors="coerce")
        xy = xy.dropna()
        roi_pts = xy.to_numpy(dtype=float)

        if roi_pts.shape[0] >= 2:
            keep = [True]
            for i in range(1, roi_pts.shape[0]):
                keep.append(not np.allclose(roi_pts[i], roi_pts[i - 1]))
            roi_pts = roi_pts[np.array(keep, dtype=bool)]

        if roi_pts.shape[0] >= 2 and np.allclose(roi_pts[0], roi_pts[-1]):
            roi_pts = roi_pts[:-1]

        if roi_pts.shape[0] < 3:
            raise ValueError(f"ROI polygon needs >=3 vertices after cleaning: {roi_csv}")

        if not np.allclose(roi_pts[0], roi_pts[-1]):
            roi_pts = np.vstack([roi_pts, roi_pts[0]])

        spatial = np.asarray(adata.obsm["spatial"], dtype=float)
        px = spatial[:, 0]
        py = spatial[:, 1]

        xv = roi_pts[:, 0]
        yv = roi_pts[:, 1]
        inside = np.zeros(px.shape[0], dtype=bool)

        for i in range(xv.shape[0] - 1):
            x0 = float(xv[i])
            y0 = float(yv[i])
            x1 = float(xv[i + 1])
            y1 = float(yv[i + 1])

            cond_y = (y0 > py) != (y1 > py)
            if not np.any(cond_y):
                continue
            x_int = (x1 - x0) * (py[cond_y] - y0) / (y1 - y0) + x0
            inside[cond_y] = np.logical_xor(inside[cond_y], px[cond_y] < x_int)

        adata = adata[inside].copy()
        if adata.n_obs == 0:
            pd.DataFrame({"gene": [], "module_id": []}).to_csv(
                module_path, sep="\t", index=False
            )
            pd.DataFrame({"gene_a": [], "gene_b": [], "weight": []}).to_csv(
                edges_path, sep="\t", index=False
            )
            raise SystemExit(0)

    if int(args.max_cells) > 0 and adata.n_obs > int(args.max_cells):
        rng = np.random.default_rng(int(args.seed))
        idx = rng.choice(adata.n_obs, size=int(args.max_cells), replace=False)
        adata = adata[idx].copy()

    counts = adata.layers["counts"]
    if not sp.issparse(counts):
        counts = sp.csr_matrix(counts)
    elif not sp.isspmatrix_csr(counts):
        counts = counts.tocsr()
    adata.layers["counts"] = counts

    if adata.n_vars == 0:
        pd.DataFrame({"gene": [], "module_id": []}).to_csv(module_path, sep="\t", index=False)
        pd.DataFrame({"gene_a": [], "gene_b": [], "weight": []}).to_csv(
            edges_path, sep="\t", index=False
        )
        raise SystemExit(0)

    gene_sum = np.asarray(adata.layers["counts"].sum(axis=0)).ravel().astype(float)
    gene_sum_sq = (
        np.asarray(adata.layers["counts"].multiply(adata.layers["counts"]).sum(axis=0))
        .ravel()
        .astype(float)
    )
    mean = gene_sum / float(adata.n_obs)
    mean_sq = gene_sum_sq / float(adata.n_obs)
    var = mean_sq - mean**2
    var_mask = var > 0.0

    if int(np.count_nonzero(var_mask)) == 0:
        pd.DataFrame({"gene": [], "module_id": []}).to_csv(module_path, sep="\t", index=False)
        pd.DataFrame({"gene_a": [], "gene_b": [], "weight": []}).to_csv(
            edges_path, sep="\t", index=False
        )
        raise SystemExit(0)

    if int(np.count_nonzero(var_mask)) != int(adata.n_vars):
        adata = adata[:, var_mask].copy()
        counts = adata.layers["counts"]
        if not sp.isspmatrix_csr(counts):
            adata.layers["counts"] = counts.tocsr()

    total_counts = np.asarray(adata.layers["counts"].sum(axis=1)).ravel()
    total_counts = np.maximum(total_counts, 1.0)
    adata.obs["total_counts"] = total_counts

    hs = hotspot.Hotspot(
        adata,
        layer_key="counts",
        model=args.model,
        latent_obsm_key="spatial",
        umi_counts_obs_key="total_counts",
    )
    hs.create_knn_graph(n_neighbors=int(args.n_neighbors), weighted_graph=False)

    hs_results = hs.compute_autocorrelations(jobs=int(args.jobs))
    if "FDR" not in hs_results.columns:
        raise ValueError("Hotspot autocorrelation results missing 'FDR' column.")
    hs_genes = hs_results.index[hs_results["FDR"] < float(args.fdr_autocorr)]
    hs_genes = hs_genes.astype(str)

    if hs_genes.size < 2:
        pd.DataFrame({"gene": [], "module_id": []}).to_csv(
            module_path, sep="\t", index=False
        )
        pd.DataFrame({"gene_a": [], "gene_b": [], "weight": []}).to_csv(
            edges_path, sep="\t", index=False
        )
        raise SystemExit(0)

    lcz = hs.compute_local_correlations(hs_genes, jobs=int(args.jobs))
    core_only_val = True
    core_only_str = str(args.core_only).strip().lower()
    if core_only_str in {"0", "false", "f", "no", "n"}:
        core_only_val = False
    elif core_only_str in {"1", "true", "t", "yes", "y"}:
        core_only_val = True
    else:
        raise ValueError("--core_only must be one of 0/1/true/false")

    modules = hs.create_modules(
        min_gene_threshold=int(args.min_gene_threshold),
        core_only=core_only_val,
        fdr_threshold=float(args.fdr_autocorr),
    )

    if isinstance(modules, pd.Series):
        module_series = modules.copy()
    elif isinstance(modules, pd.DataFrame) and modules.shape[1] == 1:
        module_series = modules.iloc[:, 0].copy()
    else:
        module_series = pd.Series(modules)

    module_series = module_series.dropna()
    module_series.index = module_series.index.astype(str)
    module_out = pd.DataFrame(
        {"gene": module_series.index.to_numpy(), "module_id": module_series.to_numpy()}
    )
    module_out.to_csv(module_path, sep="\t", index=False)

    if not isinstance(lcz, pd.DataFrame):
        lcz = pd.DataFrame(lcz)

    lcz = lcz.copy()
    lcz.index = lcz.index.astype(str)
    lcz.columns = lcz.columns.astype(str)
    if list(lcz.columns) != list(lcz.index):
        lcz = lcz.reindex(index=lcz.index, columns=lcz.index)

    genes = lcz.index.to_numpy()
    mat = lcz.to_numpy(dtype=float)
    tri = np.triu_indices(mat.shape[0], k=1)
    edges = pd.DataFrame(
        {
            "gene_a": genes[tri[0]],
            "gene_b": genes[tri[1]],
            "weight": mat[tri],
        }
    )
    edges = edges[np.isfinite(edges["weight"].to_numpy())]
    edges = edges[["gene_a", "gene_b", "weight"]]
    edges.to_csv(edges_path, sep="\t", index=False)
