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
- `correction-analysis/`: frozen four-sample configuration, rerun entry point, and result gates

## Frozen geneSCOPE source

The audited geneSCOPE `main` line is version **1.0.2**. The workflows below pin
an exact commit from that line rather than relying on a moving branch name.

The geneSCOPE image installs an exact local **v1.0.2** source snapshot from
`docker/genescope/vendor/geneSCOPE-v1.0.2`. It does not install the historical
v1.0.0 tree at `docker/genescope/geneSCOPE/`, and it does not follow a mutable
remote branch during the build.

After the v1.0.2 package candidate has passed its own checks, synchronize it into
this repository once:

```bash
bash docker/genescope/sync_genescope_vendor.sh /absolute/path/to/geneSCOPE-v1.0.2
```

The Docker build then runs independent assertions for the package version, API
defaults, and the native canonical Lee S2 calculation. A build fails if any
assertion differs from the frozen contract.

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

Provide the sample and raw Xenium directory, then run:

```bash
SAMPLE_ID=P5 OUTS="/absolute/path/to/xenium_outs" bash benchmark.sh
```

`SAMPLE_ID` must be `P1`, `P2`, `P5`, or `LN`. The matching repository ROI file
is selected automatically; `ROI_CSV` can override it. Other parameters remain
configurable through environment variables:

```bash
OUTS="/path/to/xenium_outs" \
SAMPLE_ID="P5" \
THREADS=16 \
SEED=1 \
N_RUNS=1 \
bash benchmark.sh
```

Each sample is isolated under `final_out/<SAMPLE_ID>/`. The frozen edge
benchmark ranks by raw `weight`, uses positive weights, does not refill unmapped
Top-N positions, and evaluates Top N = 10, 30, 50, 100, and 1000. geneSCOPE uses
FDR < 0.05 for P1/P2 and no edge-FDR filter for P5/LN; comparator methods use no
edge-FDR filter. The module benchmark retains the published random-mixing null:
it samples module-to-background cross-pairs, with STRING-absent pairs scored as zero.

## Four-sample correction rerun

The unified entry point does not contain local raw-data paths. Supply each input
directory explicitly:

```bash
GENESCOPE_P1_OUTS="/absolute/path/to/P1/outs" \
GENESCOPE_P2_OUTS="/absolute/path/to/P2/outs" \
GENESCOPE_P5_OUTS="/absolute/path/to/P5/outs" \
GENESCOPE_LN_OUTS="/absolute/path/to/LN/outs" \
RESULT_ROOT="/absolute/path/to/correction-results" \
THREADS=16 \
bash correction-analysis/run_all_samples.sh
```

The exact parameters and module anchors are in
`correction-analysis/samples.tsv`. Delta adjustment is
eligible-universe-scaled BH for the selected Top-N, and the entry point requires
the selected set to equal the complete eligible universe. Display pairs use
`q_Delta < 0.05`, `L > 0`, `r < 0.05`, `pct1 > 20`, and `pct2 > 20`, followed by
full-precision `Delta = L - r` ranking.

To check an assembled reference-result bundle, run:

```bash
Rscript correction-analysis/verify_reference_results.R /absolute/path/to/result-bundle
```

The gate checks package/formula/permutation provenance, exact gene-to-module
membership (including a label-invariant partition digest), Top6, Wilcoxon
values, and the P5 10/30/55-um multiscale anchors.
The complete LN Delta bundle is kept outside git; provide it with
`LN_COMPLETE_DELTA_DIR=/absolute/path/to/LN_delta_complete_v102`. Its expected
artifact SHA-256 values are frozen in
`correction-analysis/reference_external_bundle_hashes.tsv`. A P5 dendrogram
audit outside the bundle can similarly be supplied with `P5_DENDRO_AUDIT_TSV`.

The module-level comparison is recomputed with the published null: for a module
of size `m`, sample `choose(m, 2)` module-to-outside-background cross-pairs
without replacement, and score missing STRING pairs as zero.

```bash
Rscript correction-analysis/run_anchor_cross_benchmark.R \
  --benchmark-root /absolute/path/to/published-benchmark-root \
  --reanalysis-root /absolute/path/to/full_reanalysis_shuffle_v102_20260726 \
  --output-root /absolute/path/to/anchor_cross_benchmark_v102
```

Supply that frozen output to the full verifier with
`ANCHOR_CROSS_RESULTS=/absolute/path/to/anchor_cross_benchmark_v102`. The gate
requires 200,000 valid null draws per geneSCOPE module and confirms that
geneSCOPE ranks first by the observed and null-adjusted module-level summaries
in all four samples.

The P5 multiscale result itself is reproduced, rather than only checked, with:

```bash
GENESCOPE_P5_OUTS="/absolute/path/to/P5/outs" \
RESULT_ROOT="/absolute/path/to/P5-multiscale-results" \
THREADS=16 \
bash correction-analysis/run_p5_multiscale.sh
```

Correction figure workflows also install and verify the same vendored package:

```bash
GENESCOPE_P5_OUTS="/absolute/path/to/P5/outs" \
GENESCOPE_P5_SCOPE_RDS="/absolute/path/to/P5_scope_shuffle_v102.rds" \
GENESCOPE_P5_TOP_PAIRS="/absolute/path/to/P5_toplvsr_all_shuffleFDR.tsv" \
bash main-text-scripts/run_frozen_workflow.sh P5

GENESCOPE_LN_OUTS="/absolute/path/to/LN/outs" \
GENESCOPE_LN_SCOPE_RDS="/absolute/path/to/LN_scope_shuffle_v102.rds" \
GENESCOPE_LN_TOP_PAIRS="/absolute/path/to/LN_top_pairs_complete_delta_v102.tsv" \
bash main-text-scripts/run_frozen_workflow.sh LN
```

These entry points require a clean paper commit, generate in a new staging
directory, record raw-input/ROI/mapping/workflow hashes, and render only from
the hash-pinned authoritative analysis objects and complete pair tables. They
verify the complete bundle and publish the figure directory only after every
freeze gate succeeds.

## Outputs

`benchmark.sh` writes everything under `final_out/<SAMPLE_ID>/`:

- Raw method outputs:
  - `final_out/<SAMPLE_ID>/genescope/repeat_###/`
  - `final_out/<SAMPLE_ID>/giotto/repeat_###/`
  - `final_out/<SAMPLE_ID>/hotspot/repeat_###/`
  - `final_out/<SAMPLE_ID>/seagal/` (container-level `stats.tsv` + `run.log`) and its `repeat_###/`
  - “All-edges” runs:
    - `final_out/<SAMPLE_ID>/giotto-alledges/repeat_001/`
    - `final_out/<SAMPLE_ID>/hotspot-alledges/repeat_001/`
    - `final_out/<SAMPLE_ID>/seagal-alledges/repeat_001/`
- Staged benchmark root for `mapping.R`:
  - `final_out/<SAMPLE_ID>/for_compare/<method>/repeat_001/edges_all.tsv`
  - `final_out/<SAMPLE_ID>/for_compare/<method>/repeat_001/modules.tsv`
- STRINGdb mapping + downstream panels:
  - `final_out/<SAMPLE_ID>/string1step/`
  - `final_out/<SAMPLE_ID>/edge-level/`
  - `final_out/<SAMPLE_ID>/module-level/`
- Runtime panels input + plots:
  - `final_out/<SAMPLE_ID>/runtime_panels_input/`
  - `final_out/<SAMPLE_ID>/runtime_panels_input/plots_runtime/`

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
