FROM mambaorg/micromamba:1.5.8

# ---- build args: local source path (relative to build context) ----
ARG GENESCOPE_LOCAL=geneSCOPE

# ---- create conda env with R + deps ----
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba create -y -n tool -f /tmp/env.yaml && micromamba clean -a -y
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# ---- install geneSCOPE from local source ----
COPY --chown=$MAMBA_USER:$MAMBA_USER ${GENESCOPE_LOCAL} /tmp/geneSCOPE_src
RUN find /tmp/geneSCOPE_src -type f \( -name '*.o' -o -name '*.so' -o -name '*.dll' -o -name '*.dylib' -o -name '*.a' \) -delete && \
    micromamba run -n tool R CMD INSTALL --preclean /tmp/geneSCOPE_src && \
    rm -rf /tmp/geneSCOPE_src

# ---- app scripts ----
WORKDIR /opt/app
COPY --chown=$MAMBA_USER:$MAMBA_USER monitor_cgroup.sh /opt/app/monitor_cgroup.sh
COPY --chown=$MAMBA_USER:$MAMBA_USER run.sh /opt/app/run.sh
COPY --chown=$MAMBA_USER:$MAMBA_USER run_genescope_xenium.R /opt/app/run_genescope_xenium.R
RUN chmod +x /opt/app/run.sh /opt/app/monitor_cgroup.sh

ENTRYPOINT ["/opt/app/run.sh"]
