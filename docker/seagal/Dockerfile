FROM mambaorg/micromamba:1.5.8

ARG SEAGAL_GIT="https://github.com/linhuawang/SEAGAL.git"
ARG SEAGAL_REF="main"
ARG SEAGAL_VERSION="2.5.2"

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
      git \
      time \
      procps \
    && rm -rf /var/lib/apt/lists/*
USER ${MAMBA_USER}

# ---- create conda env with python + deps ----
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba create -y -n tool -f /tmp/env.yaml && micromamba clean -a -y

ENV PATH=/opt/conda/envs/tool/bin:$PATH
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# ---- install seagal (PyPI preferred; GitHub fallback) ----
RUN set -euo pipefail; \
  if ! micromamba run -n tool python -m pip install --no-deps "seagal==${SEAGAL_VERSION}"; then \
    micromamba run -n tool python -m pip install --no-deps "git+${SEAGAL_GIT}@${SEAGAL_REF}"; \
  fi

# ---- smoke test ----
RUN micromamba run -n tool python - <<'PY'
import seagal
import scanpy
import squidpy
import numpy
import pandas
import scipy
print("OK")
PY

WORKDIR /opt/app
COPY --chown=$MAMBA_USER:$MAMBA_USER run_seagal_xenium.py /opt/app/run_seagal_xenium.py
COPY --chown=$MAMBA_USER:$MAMBA_USER run.sh /opt/app/run.sh
COPY --chown=$MAMBA_USER:$MAMBA_USER monitor_cgroup.sh /opt/app/monitor_cgroup.sh
RUN chmod +x /opt/app/run.sh /opt/app/monitor_cgroup.sh

ENTRYPOINT ["/opt/app/run.sh"]
