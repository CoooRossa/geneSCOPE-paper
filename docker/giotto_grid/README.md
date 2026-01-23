# giotto_10x_grid (Docker)

This directory provides a Docker image recipe for the Giotto Xenium runner (grid-only).

## Files

- `Dockerfile`: Docker build recipe (recommended)
- `env.yaml`: Conda environment specification
- `build.sh`: Convenience script to build the Docker image
- `run.sh`: Container entrypoint (invoked by `docker run`)
- `monitor_cgroup.sh`: Optional cgroup resource sampler
- `run_giotto_xenium.R`: Giotto Xenium runner

## Build

Run in this directory:

```bash
./build.sh
```

Override the image tag (default: `giotto_10x_grid:full-edges`):

```bash
IMAGE_REF="giotto_10x_grid:full-edges" ./build.sh
```

Extra Docker build args (optional):

```bash
DOCKER_BUILD_EXTRA_ARGS="--no-cache" ./build.sh
```

## Run

`run.sh` expects **in-container paths**, so mount host directories into the container.

Example (Xenium outs -> `/data`, output -> `/out`):

```bash
docker run --rm \
  -v /host/path/to/xenium_outs:/data:ro \
  -v /host/path/to/output:/out \
  giotto_10x_grid:full-edges \
  --data_dir /data \
  --outdir /out \
  --grid_stepsize 100 \
  --threads 8
```

If the ROI CSV is not under the `xenium_outs` directory, mount an additional read-only directory to `/roi`:

```bash
docker run --rm \
  -v /host/path/to/xenium_outs:/data:ro \
  -v /host/path/to/roi_dir:/roi:ro \
  -v /host/path/to/output:/out \
  giotto_10x_grid:full-edges \
  --data_dir /data \
  --coord_file /roi/your_roi.csv \
  --outdir /out \
  --grid_stepsize 100
```

More args:

```bash
docker run --rm giotto_10x_grid:full-edges --help
```

## Outputs

The container writes only two files to `--outdir`:

- `edges_all.tsv` (columns: `gene_a`, `gene_b`, `weight`)
- `modules.tsv` (columns: `gene`, `module_id`)
