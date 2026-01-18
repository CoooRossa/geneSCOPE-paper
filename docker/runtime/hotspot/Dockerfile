FROM mambaorg/micromamba:1.5.8

ARG HOTSPOT_GIT="https://github.com/YosefLab/Hotspot.git"
ARG HOTSPOT_REF="v0.9.1"
ARG HOTSPOT_VERSION="0.9.1"

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
      git \
      time \
      procps \
      ca-certificates \
    && rm -rf /var/lib/apt/lists/*
USER ${MAMBA_USER}

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba create -y -n tool -f /tmp/env.yaml && micromamba install -n tool -c conda-forge seaborn && micromamba clean -a -y

ENV PATH=/opt/conda/envs/tool/bin:$PATH
ENV MAMBA_DOCKERFILE_ACTIVATE=1

RUN set -euo pipefail; \
  if ! micromamba run -n tool python -m pip install --no-deps "hotspot==${HOTSPOT_VERSION}"; then \
    micromamba run -n tool python -m pip install --no-deps "git+${HOTSPOT_GIT}@${HOTSPOT_REF}"; \
  fi

RUN micromamba run -n tool python - <<'PY'
import hotspot
import anndata
import numpy
import pandas
import scipy
import pyarrow
import shapely
print("OK")
PY

WORKDIR /opt/app
COPY --chown=$MAMBA_USER:$MAMBA_USER run_hotspot_xenium.py /opt/app/run_hotspot_xenium.py
COPY --chown=$MAMBA_USER:$MAMBA_USER run.sh /opt/app/run.sh
COPY --chown=$MAMBA_USER:$MAMBA_USER monitor_cgroup.sh /opt/app/monitor_cgroup.sh
RUN chmod +x /opt/app/run.sh /opt/app/monitor_cgroup.sh

ENTRYPOINT ["/opt/app/run.sh"]
