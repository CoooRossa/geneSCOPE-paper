# seagal bench container

This image runs SEAGAL (Python) on Xenium outs with ROI clip + repeat benchmarking.

## Build

docker build --no-cache -t bench/seagal:stable docker/seagal

## Smoke test

docker run --rm --entrypoint bash bench/seagal:stable -lc \
'micromamba run -n tool python -c "import seagal, scanpy, squidpy; print(\"OK\")"'

## Run (single or repeat)

docker run --rm \
  --memory=16g --memory-swap=16g --cpus=16 \
  -v "/Users/haenolabcho/Documents/Code/genescope_docker_bench:/data:ro" \
  -v "/Users/haenolabcho/Documents/Code/genescope_docker_bench/results/seagal:/out" \
  bench/seagal:stable \
  --data_dir "/data/XENIUM_OUTS" \
  --coord_file "/data/XENIUM_OUTS/roi.csv" \
  --outdir /out \
  --ncores 16 \
  --seed 1 \
  --sample_sec 1 \
  --repeat 1 \
  --max_cells 0 \
  --top_genes 0

## Outputs

Top-level (in `outdir/`):

- resource_limits.txt
- version.txt
- stats.tsv
- run.log
- replicate_01/
- replicate_02/
- ...

Each `replicate_XX/` contains:

- run.log
- metrics.tsv
- summary.tsv
- run_args.json
- meta.json
- seagal_results.tsv
- edges_all.tsv
- modules.tsv
- _seagal_input/ (count.csv + meta.csv)
- run_stats.tsv

### metrics.tsv (minimum keys)

- wall_time_sec
- n_cells
- n_genes
- n_sig
- seed
- roi_applied
- max_cells

### seagal_results.tsv

Columns follow SEAGAL global L output, including:

- gene_1, gene_2
- L
- L.p_value
- L.FDR
- pair
- -log10(FDR)
- Association

### edges_all.tsv

Standardized edge export for benchmarking (pattern genes only; full-gene all-to-all disabled).

Columns (fixed order):

- gene_a, gene_b
- weight
- p_value, fdr
- n_obs
- edge_type (always `seagal_global_L_pattern`)

### modules.tsv

Hard module assignment (module_id starts at 1; unassigned is -1):

- gene
- module_id

## CLI flags

Required:
- --data_dir
- --outdir

Optional:
- --coord_file (ROI CSV)
- --ncores (default 1)
- --seed (default 1)
- --sample_sec (default 1)
- --repeat (default 1)
- --max_cells (default 0; no downsample)
- --top_genes (default 0; no filter)
- --grid_um (unused; kept for CLI parity)
- --extra_args_json (JSON string or path)

### extra_args_json (advanced)

Supported keys:
- min_counts (default 150)
- min_cells (default 10)
- n_permutation (default 99)
- permute_ratio (default 0.2)
- fdr_cutoff (default 0.05)
- l_cutoff (default 0.1)
- indep (default true)
- use_pattern_genes (default true)
- svg_I (optional; alias: pattern_i)
- svg_topk (optional; alias: pattern_topk)
- modules_enabled (default 1)
- modules_nmax (default 6)
- modules_n (optional; fix module count)

## ROI CSV format

- Lines starting with `#` are ignored.
- Header must contain `X,Y` or `x,y` (case-insensitive).
- Non-numeric rows are dropped.
- Consecutive duplicate points are removed.
- If the first and last points are identical, the last is dropped before closing.
- Polygon is auto-closed (first point appended if needed).
- At least 3 vertices are required after cleaning.

## Acceptance checklist

1) Build

docker build --no-cache -t bench/seagal:stable docker/seagal

2) Smoke test

docker run --rm --entrypoint bash bench/seagal:stable -lc \
'micromamba run -n tool python -c "import seagal, scanpy, squidpy; print(\"OK\")"'

3) Real-run (repeat=1)

- run.log
- stats.tsv (at least header)
- replicate_01/metrics.tsv
- replicate_01/seagal_results.tsv

4) Repeat=10

- replicate_01 ... replicate_10
- each with metrics.tsv
