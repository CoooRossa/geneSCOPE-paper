# hotspot (Docker)

This directory provides a Docker image recipe and runner for Hotspot Xenium.

## Files

- `Dockerfile`: Docker build recipe (recommended)
- `env.yaml`: Conda environment specification
- `build.sh`: Convenience script to build the Docker image
- `run.sh`: Container entrypoint (invoked by `docker run`)
- `monitor_cgroup.sh`: Optional cgroup resource sampler
- `run_hotspot_xenium.py`: Hotspot Xenium runner

## Build

Run in this directory:

```bash
./build.sh
```

Override the image tag (default: `hotspot:full-edges`):

```bash
IMAGE_REF="hotspot:full-edges" ./build.sh
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
  hotspot:full-edges \
  --data_dir /data/outs \
  --outdir /out \
  --coord_file /coord/roi.csv \
  --ncores 64
```

More args:

```bash
docker run --rm hotspot:full-edges --help
```
