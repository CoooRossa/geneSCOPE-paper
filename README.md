# geneSCOPE-paper

Companion materials for the **geneSCOPE** manuscript:

- ROI coordinate CSVs used in the paper
- Reproducible Docker runners (geneSCOPE + baseline methods)
- R scripts to generate benchmarking analyses (mapping to STRINGdb, edge-level, module-level, runtime panels)
- Docker execution with the corrected geneSCOPE v1.2.0 release

This repository does **not** include raw Xenium datasets; you must provide your own 10x Genomics Xenium `outs/` folder(s).

## Repository layout

- `benchmark.sh`: one-shot pipeline (run methods → stage outputs → run R benchmark scripts → runtime panels)
- `ROI-coordinate-files/`: ROI polygon CSVs (X/Y vertices) used in the paper
- `docker/`: Docker images and runners
  - `docker/genescope/`
  - `docker/giotto_grid/`
  - `docker/hotspot/`
  - `docker/seagal/`
- `benchmark-Rscripts/`: post-processing / benchmarking scripts
  - `mapping.R` (maps predicted edges to STRINGdb)
  - `edge-level.R`
  - `module-level.R`
  - `runtime-panels.R`
- `main-text-scripts/`: scripts used for main-text workflows/figures (edit paths before running)

## Inputs

### Xenium `outs/`

All runners expect a 10x Genomics Xenium `outs/` directory (host path mounted into containers), typically containing
`cell_feature_matrix.h5` and `cells.csv.gz`.

### ROI CSV format

The ROI CSVs in `ROI-coordinate-files/` represent a polygon in Xenium coordinates:

- Lines starting with `#` are treated as comments.
- Header must contain `X,Y` (case-insensitive). If missing, the first 2 columns are used.
- Remaining rows are polygon vertices (numeric `X` and `Y`).
- The polygon is auto-closed if the first/last vertex differ.

## Quick start (one command)

1. Select a sample and provide its Xenium `outs/` path. The matching ROI file
   under `ROI-coordinate-files/` is used by default.
2. Run:

```bash
bash benchmark.sh
```

Tip: for quick tests without editing the file, you can override parameters via environment variables:

```bash
OUTS="/path/to/xenium_outs" \
SAMPLE_ID="P5" \
THREADS=16 \
SEED=1 \
N_RUNS=1 \
bash benchmark.sh
```

## Outputs

`benchmark.sh` writes each sample under `final_out/<SAMPLE_ID>/`:

- Raw method outputs: `final_out/<SAMPLE_ID>/<method>/`
- Staged benchmark inputs: `final_out/<SAMPLE_ID>/for_compare/`
- STRINGdb and downstream panels: `final_out/<SAMPLE_ID>/string1step/` and related directories
- Runtime panels: `final_out/<SAMPLE_ID>/runtime_panels_input/`

### File-format guarantees (used by benchmarking scripts)

The pipeline standardizes outputs so the analysis scripts can read them:

- `all_edges.tsv`: first 3 columns are `gene_a`, `gene_b`, `weight` (tab-separated).
- `modules.tsv`: columns are `gene`, `module_id` (tab-separated).
- `mapping.R` reads `edges_all.tsv` (a staged copy of `all_edges.tsv`) under `final_out/<SAMPLE_ID>/for_compare/...`.

## R dependencies (benchmark scripts)

The analysis scripts in `benchmark-Rscripts/` require (at minimum):

- CRAN: `optparse`, `data.table`, `jsonlite`, `ggplot2`
- Bioconductor: `STRINGdb`

`mapping.R` downloads STRINGdb data on first run (internet required).
