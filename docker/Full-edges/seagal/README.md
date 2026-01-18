# seagal (Docker)

This directory provides a Docker image recipe and runner for SEAGAL Xenium.

## Files

- `Dockerfile`: Docker build recipe (recommended)
- `env.yaml`: Conda environment specification
- `build.sh`: Convenience script to build the Docker image
- `run.sh`: Container entrypoint (invoked by `docker run`)
- `monitor_cgroup.sh`: Optional cgroup resource sampler
- `run_seagal_xenium.py`: SEAGAL Xenium runner

## Build

Run in this directory:

```bash
./build.sh
```

Override the image tag (default: `seagal:full-edges`):

```bash
IMAGE_REF="seagal:full-edges" ./build.sh
```

Extra Docker build args (optional):

```bash
DOCKER_BUILD_EXTRA_ARGS="--no-cache" ./build.sh
```

## Run

`run.sh` expects **in-container paths**, so mount host directories into the container.

Example (outs -> `/data/outs`, output -> `/out`, ROI CSV -> `/coord/roi.csv`):

```bash
mkdir -p /host/out
docker run --rm \
  -v /host/xenium_outs:/data/outs:ro \
  -v /host/out:/out \
  -v /host/roi.csv:/coord/roi.csv:ro \
  seagal:full-edges \
  --data_dir /data/outs \
  --outdir /out \
  --coord_file /coord/roi.csv \
  --ncores 64 \
  --seed 1 \
  --repeat 1
```

More args:

```bash
docker run --rm seagal:full-edges --help
```

## Outputs (full-gene edges_all + pattern clustering)

- `edges_all.tsv`: Computes Global_L for **all genes**. The implementation processes genes in blocks and calls SEAGAL
  `spatial_association()` per block to avoid building an `n_obs × n_pairs` full-gene matrix at once. By default,
  `p_value/fdr` is set to `1.0` (`seagal_edges_n_permutation=0`).
- `modules.tsv`: Clusters pattern genes selected by `spatial_pattern_genes` (silhouette-selected k; equivalent to the
  SEAGAL `genemodules` approach), then maps the result back to the full gene table (non-pattern genes use `-1`).
- Use `--grid_um` and/or `--max_cells` to control `n_obs`. Use
  `--extra_args_json '{"seagal_edges_block_genes":300,"seagal_edges_memory_budget_gib":64}'` to control the gene-blocking
  strategy for `edges_all.tsv`.
