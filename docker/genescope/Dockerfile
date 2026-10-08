FROM mambaorg/micromamba:1.5.8

ARG GENESCOPE_VERSION=1.2.0
ARG GENESCOPE_COMMIT=e5e7af8a3b5129ace885dc765e72ea75dd4a6665

# ---- create conda env with R + deps ----
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba create -y -n tool -f /tmp/env.yaml && \
    micromamba env export -n tool --explicit > /tmp/geneSCOPE-conda-explicit.txt && \
    micromamba list -n tool --json > /tmp/geneSCOPE-conda-packages.json && \
    micromamba clean -a -y
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# ---- install the release used by the paper ----
RUN micromamba run -n tool R -e \
    "remotes::install_github('CoooRossa/geneSCOPE', ref = '${GENESCOPE_COMMIT}', dependencies = FALSE, upgrade = 'never'); stopifnot(as.character(utils::packageVersion('geneSCOPE')) == '${GENESCOPE_VERSION}')" \
    -e "installed_commit <- utils::packageDescription('geneSCOPE', fields = 'RemoteSha'); stopifnot(identical(installed_commit, '${GENESCOPE_COMMIT}')); dir.create('/tmp/geneSCOPE_provenance', showWarnings = FALSE); writeLines(installed_commit, '/tmp/geneSCOPE_provenance/source_commit')" && \
    mv /tmp/geneSCOPE-conda-explicit.txt /tmp/geneSCOPE_provenance/conda-explicit.txt && \
    mv /tmp/geneSCOPE-conda-packages.json /tmp/geneSCOPE_provenance/conda-packages.json

# ---- app scripts ----
WORKDIR /opt/app
COPY --chown=$MAMBA_USER:$MAMBA_USER monitor_cgroup.sh /opt/app/monitor_cgroup.sh
COPY --chown=$MAMBA_USER:$MAMBA_USER run.sh /opt/app/run.sh
COPY --chown=$MAMBA_USER:$MAMBA_USER run_genescope_xenium.R /opt/app/run_genescope_xenium.R
RUN chmod +x /opt/app/run.sh /opt/app/monitor_cgroup.sh

ENTRYPOINT ["/opt/app/run.sh"]
