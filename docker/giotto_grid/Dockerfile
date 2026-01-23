FROM mambaorg/micromamba:1.5.8

ARG GIOTTO_REPO="drieslab/Giotto"
ARG GIOTTO_REF="suite"

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
      git \
      time \
      procps \
      build-essential \
      gfortran \
      libfontconfig1-dev \
      libfreetype6-dev \
      libgdal-dev \
      libglpk-dev \
      libproj-dev \
      libgeos-dev \
      libudunits2-dev \
      libssl-dev \
      libcurl4-openssl-dev \
      libxml2-dev \
    && rm -rf /var/lib/apt/lists/*
USER ${MAMBA_USER}

# ---- create conda env ----
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba create -y -n tool -f /tmp/env.yaml && micromamba clean -a -y
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# ---- install Giotto Suite ----
RUN micromamba run -n tool R -q -e " \
  options(repos=c(CRAN='https://cloud.r-project.org')); \
  try(pak::pak_update(), silent = TRUE); \
  spec <- if (nzchar('${GIOTTO_REF}')) sprintf('%s@%s', '${GIOTTO_REPO}', '${GIOTTO_REF}') else '${GIOTTO_REPO}'; \
  message(sprintf('Installing %s', spec)); \
  pak::pkg_install(spec); \
  "

# ---- smoke test ----
RUN micromamba run -n tool R -q -e 'suppressPackageStartupMessages(library(Giotto)); cat("OK");'

WORKDIR /opt/app
COPY --chown=$MAMBA_USER:$MAMBA_USER run_giotto_xenium.R /opt/app/run_giotto_xenium.R
COPY --chown=$MAMBA_USER:$MAMBA_USER monitor_cgroup.sh   /opt/app/monitor_cgroup.sh
COPY --chown=$MAMBA_USER:$MAMBA_USER run.sh              /opt/app/run.sh
RUN chmod +x /opt/app/run.sh /opt/app/monitor_cgroup.sh /opt/app/run_giotto_xenium.R

ENTRYPOINT ["/opt/app/run.sh"]
