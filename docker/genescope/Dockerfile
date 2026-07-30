FROM mambaorg/micromamba:1.5.8

ARG GENESCOPE_VERSION=1.2.0

# ---- create conda env with R + deps ----
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba create -y -n tool -f /tmp/env.yaml && micromamba clean -a -y
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# ---- install and verify the corrected release ----
RUN micromamba run -n tool R -e \
    "remotes::install_github('CoooRossa/geneSCOPE', ref = 'v${GENESCOPE_VERSION}', dependencies = FALSE, upgrade = 'never'); stopifnot(as.character(utils::packageVersion('geneSCOPE')) == '${GENESCOPE_VERSION}')"

# ---- app scripts ----
WORKDIR /opt/app
COPY --chown=$MAMBA_USER:$MAMBA_USER monitor_cgroup.sh /opt/app/monitor_cgroup.sh
COPY --chown=$MAMBA_USER:$MAMBA_USER run.sh /opt/app/run.sh
COPY --chown=$MAMBA_USER:$MAMBA_USER run_genescope_xenium.R /opt/app/run_genescope_xenium.R
RUN chmod +x /opt/app/run.sh /opt/app/monitor_cgroup.sh

ENTRYPOINT ["/opt/app/run.sh"]
