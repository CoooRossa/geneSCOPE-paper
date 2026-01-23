# hotspotsc 1.1.3 Docker runner (AnnData route)

This container pins hotspotsc==1.1.3 (PyPI latest as of release Apr 19, 2025).
Hotspot v1+ uses AnnData as the primary interface.

## Build
docker build -t hotspotsc:1.1.3 -f Dockerfile .

## Run (Xenium outs)
docker run --rm -u "$(id -u):$(id -g)" \
  -v /ABS/PATH/TO/XENIUM_OUTS:/data/outs:ro \
  -v /ABS/PATH/TO/OUTDIR:/out \
  hotspotsc:1.1.3 \
  --xenium-outs /data/outs \
  --outdir /out \
  --model bernoulli \
  --n-neighbors 50

## Run (Xenium outs + ROI)
ROI CSV is a polygon (QuPath-style) with `x`,`y` columns (or the first 2 columns).

docker run --rm -u "$(id -u):$(id -g)" \
  -v /ABS/PATH/TO/XENIUM_OUTS:/data/outs:ro \
  -v /ABS/PATH/TO/OUTDIR:/out \
  hotspotsc:1.1.3 \
  --xenium-outs /data/outs \
  --coord_file /data/outs/roi.csv \
  --outdir /out \
  --model bernoulli \
  --n-neighbors 50

Outputs:
- /out/module.tsv
- /out/edges_all.tsv

## Run (h5ad)
docker run --rm -u "$(id -u):$(id -g)" \
  -v /ABS/PATH/TO/adata.h5ad:/data/adata.h5ad:ro \
  -v /ABS/PATH/TO/OUTDIR:/out \
  hotspotsc:1.1.3 \
  --h5ad /data/adata.h5ad \
  --outdir /out \
  --model bernoulli \
  --n-neighbors 50

Notes:
- This runner also accepts some legacy flags for compatibility (e.g. `--data_dir`, `--threads`), but it always writes only `module.tsv` and `edges_all.tsv` and enforces `--repeat 1`.
