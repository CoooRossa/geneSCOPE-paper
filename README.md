# geneSCOPE-paper

Companion materials for the **geneSCOPE** manuscript, including:

- ROI coordinate CSVs used in the paper
- Docker recipes/runners for reproducible execution (geneSCOPE + baseline methods)
- R scripts to generate the benchmarking analyses (edge-level, module-level, runtime)

This repository does **not** include raw Xenium datasets; you must provide your own 10x Genomics Xenium `outs/` folders.

## Repository layout

- `ROI-coordinate-files/`: ROI polygon CSVs (X/Y vertices) used in the paper
- `docker/`: Docker images and runners
  - `docker/runtime/`: repeat-run benchmarking containers (CPU/memory sampling + per-repeat outputs)
  - `docker/Full-edges/`: full-gene edge export containers (`edges_all.tsv`)
- `benchmark-Rscripts/`: benchmarking post-processing scripts (mapping, edge-level, module-level, runtime panels)
- `main-text-scripts/`: example scripts used for main-text workflows/figures (edit paths before running)
- `runtime-test.sh`: convenience wrapper to build+run runtime benchmark containers
- `full-edges-test.sh`: convenience wrapper to build+run full-edge export containers

## Inputs

### Xenium `outs/`

All runners expect a 10x Genomics Xenium `outs/` directory (container path), typically containing files like
`cells.parquet`, `transcripts.parquet`, etc.

### ROI CSV format

The ROI CSVs in `ROI-coordinate-files/` represent a polygon in Xenium coordinates:

- Lines starting with `#` are comments and can be ignored by runners.
- Header must contain `X,Y` (case-insensitive).
- Remaining rows are polygon vertices (numeric `X` and `Y`).
- The polygon is auto-closed if the first/last vertex differ.

## Quick start: runtime benchmark (repeat runs)

1. Edit `OUTS` and `ROI_CSV` in `runtime-test.sh` (the script will error if left as `/path/to/...`).
2. Run:

```bash
./runtime-test.sh
```

Outputs are written under:

`<this-repo>/<DATASET_ID>__<ROI_ID>/{genescope,giotto,hotspot,seagal}/repeat_*/`

Plot CPU/memory panels from the run output:

```bash
Rscript benchmark-Rscripts/runtime-panels.R \
  --run_root "<this-repo>/<DATASET_ID>__<ROI_ID>"
```

## Quick start: full-gene edge export (full edges)

1. Edit dataset paths at the bottom of `full-edges-test.sh` (the `*_OUTS` and `*_ROI_CSV` variables).
2. Optionally tune resources:

```bash
THREADS=64 MEM=256g ./full-edges-test.sh
```

Outputs are staged to:

`<BENCH_ROOT>/{P1,P2,P5,Lymph}/{genescope,giotto,hotspot,seagal}/repeat_001/`

See the method-specific Docker READMEs for details:

- `docker/Full-edges/genescope/`
- `docker/Full-edges/giotto_grid/`
- `docker/Full-edges/hotspot/`
- `docker/Full-edges/seagal/`

## Benchmark analysis (edge-level + module-level)

The analysis scripts in `benchmark-Rscripts/` expect a benchmark root with per-method outputs. For each method, the first
repeat directory should contain at least:

- `edges_all.tsv` (standardized edge list: `gene_a`, `gene_b`, `weight`, `fdr`, …)
- `modules.tsv` (gene → module assignment)
- `meta.json` (optional run metadata; used for provenance and some default thresholds)

### 1) Map predicted edges to STRINGdb

```bash
Rscript benchmark-Rscripts/mapping.R \
  --bench_root "<BENCH_ROOT>/P5" \
  --outdir "<BENCH_ROOT>/P5/stringdb-annotation" \
  --methods "genescope,giotto,hotspot,seagal" \
  --string_version "12.0" \
  --keep_subscores 1
```

Notes:

- `mapping.R` needs **geneSCOPE** R code to access shared STRING utilities: pass `--genescope_root /path/to/geneSCOPE`
  or install the `geneSCOPE` R package.
- `STRINGdb` downloads files on first run (internet required).

### 2) Edge-level benchmark panels

```bash
Rscript benchmark-Rscripts/edge-level.R \
  --map_dir "<BENCH_ROOT>/P5/stringdb-annotation" \
  --outdir "<BENCH_ROOT>/P5/edge-level" \
  --methods "genescope,giotto,hotspot,seagal"
```

### 3) Module-level benchmark panels

```bash
Rscript benchmark-Rscripts/module-level.R \
  --map_dir "<BENCH_ROOT>/P5/stringdb-annotation" \
  --outdir "<BENCH_ROOT>/P5/module-level" \
  --methods "genescope,giotto,hotspot,seagal" \
  --modules_tsv_by_method "<BENCH_ROOT>/P5/genescope/repeat_001/modules.tsv,<BENCH_ROOT>/P5/giotto/repeat_001/modules.tsv,<BENCH_ROOT>/P5/hotspot/repeat_001/modules.tsv,<BENCH_ROOT>/P5/seagal/repeat_001/modules.tsv"
```

## R dependencies (analysis scripts)

Minimum expected packages:

- `optparse`, `data.table`, `jsonlite`, `ggplot2`
- `STRINGdb` (Bioconductor)
