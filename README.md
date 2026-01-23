# geneSCOPE-paper

Companion materials for the **geneSCOPE** manuscript:

- ROI coordinate CSVs used in the paper
- Reproducible Docker runners (geneSCOPE + baseline methods)
- R scripts to generate benchmarking analyses (mapping to STRINGdb, edge-level, module-level, runtime panels)

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

1. Edit the paths at the top of `benchmark.sh`:
   - `OUTS=/path/to/<xenium_outs_dir>`
   - `ROI_CSV=/path/to/<roi.csv>` (you can point to a file under `ROI-coordinate-files/`)
2. Run:

```bash
bash benchmark.sh
```

Tip: for quick tests without editing the file, you can override parameters via environment variables:

```bash
OUTS="/path/to/xenium_outs" \
ROI_CSV="/path/to/roi.csv" \
THREADS=16 \
SEED=1 \
DATASET_ID="GSE280314_P5" \
ROI_ID="P5_tumor_region_final" \
N_RUNS=1 \
bash benchmark.sh
```

## Outputs

`benchmark.sh` writes everything under `final_out/`:

- Raw method outputs:
  - `final_out/genescope/repeat_###/`
  - `final_out/giotto/repeat_###/`
  - `final_out/hotspot/repeat_###/`
  - `final_out/seagal/` (container-level `stats.tsv` + `run.log`) and `final_out/seagal/repeat_###/`
  - “All-edges” runs:
    - `final_out/giotto-alledges/repeat_001/`
    - `final_out/hotspot-alledges/repeat_001/`
    - `final_out/seagal-alledges/repeat_001/`
- Staged benchmark root for `mapping.R`:
  - `final_out/for_compare/<method>/repeat_001/edges_all.tsv`
  - `final_out/for_compare/<method>/repeat_001/modules.tsv`
- STRINGdb mapping + downstream panels:
  - `final_out/string1step/`
  - `final_out/precision_and_recall/`
  - `final_out/anchor-null2/`
- Runtime panels input + plots:
  - `final_out/runtime_panels_input/`
  - `final_out/runtime_panels_input/plots_runtime/`

### File-format guarantees (used by benchmarking scripts)

The pipeline standardizes outputs so the analysis scripts can read them:

- `all_edges.tsv`: first 3 columns are `gene_a`, `gene_b`, `weight` (tab-separated).
- `modules.tsv`: columns are `gene`, `module_id` (tab-separated).
- `mapping.R` reads `edges_all.tsv` (a staged copy of `all_edges.tsv`) under `final_out/for_compare/...`.

## R dependencies (benchmark scripts)

The analysis scripts in `benchmark-Rscripts/` require (at minimum):

- CRAN: `optparse`, `data.table`, `jsonlite`, `ggplot2`
- Bioconductor: `STRINGdb`

`mapping.R` downloads STRINGdb data on first run (internet required).
