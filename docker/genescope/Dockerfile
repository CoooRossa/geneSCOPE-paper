FROM mambaorg/micromamba:1.5.8

# The release source is vendored deliberately so that a rebuild does not depend
# on the state of a remote branch.  Populate it with sync_genescope_vendor.sh.
ARG GENESCOPE_LOCAL=vendor/geneSCOPE-v1.2.0
ARG GENESCOPE_VERSION=1.2.0

# ---- create conda env with R + deps ----
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba create -y -n tool -f /tmp/env.yaml && \
    micromamba env export -n tool --explicit > /tmp/geneSCOPE-conda-explicit.txt && \
    micromamba list -n tool --json > /tmp/geneSCOPE-conda-packages.json && \
    micromamba clean -a -y
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# ---- install and verify the exact correction candidate ----
COPY --chown=$MAMBA_USER:$MAMBA_USER ${GENESCOPE_LOCAL} /tmp/geneSCOPE_src
COPY --chown=$MAMBA_USER:$MAMBA_USER freeze_assertions.R /tmp/freeze_assertions.R
RUN test "$(sed -n 's/^Version:[[:space:]]*//p' /tmp/geneSCOPE_src/DESCRIPTION)" = "${GENESCOPE_VERSION}" && \
    test -s /tmp/geneSCOPE_src/.freeze-source-commit && \
    test -s /tmp/geneSCOPE_src/.freeze-tree-sha256 && \
    test "$(tr -d '[:space:]' < /tmp/geneSCOPE_src/.freeze-source-commit)" != "unknown" && \
    tr -d '[:space:]' < /tmp/geneSCOPE_src/.freeze-source-commit | grep -Eq '^[0-9a-f]{40}$' && \
    EXPECTED_TREE="$(tr -d '[:space:]' < /tmp/geneSCOPE_src/.freeze-tree-sha256)" && \
    OBSERVED_TREE="$(cd /tmp/geneSCOPE_src && find . -type f ! -name '.freeze-source-commit' ! -name '.freeze-tree-sha256' -print | LC_ALL=C sort | while IFS= read -r source_file; do sha256sum "${source_file}"; done | sha256sum | awk '{print $1}')" && \
    test "${OBSERVED_TREE}" = "${EXPECTED_TREE}" && \
    grep -Fq 'return arma::dot(row_sums, row_sums);' /tmp/geneSCOPE_src/src/2.LeeL.cpp && \
    find /tmp/geneSCOPE_src -type f \( -name '*.o' -o -name '*.so' -o -name '*.dll' -o -name '*.dylib' -o -name '*.a' \) -delete && \
    micromamba run -n tool R CMD INSTALL --preclean /tmp/geneSCOPE_src && \
    micromamba run -n tool Rscript /tmp/freeze_assertions.R "${GENESCOPE_VERSION}" && \
    mkdir -p /tmp/geneSCOPE_provenance && \
    cp /tmp/geneSCOPE_src/.freeze-source-commit /tmp/geneSCOPE_provenance/source_commit && \
    cp /tmp/geneSCOPE_src/.freeze-tree-sha256 /tmp/geneSCOPE_provenance/vendor_tree_sha256 && \
    mv /tmp/geneSCOPE-conda-explicit.txt /tmp/geneSCOPE_provenance/conda-explicit.txt && \
    mv /tmp/geneSCOPE-conda-packages.json /tmp/geneSCOPE_provenance/conda-packages.json && \
    rm -rf /tmp/geneSCOPE_src /tmp/freeze_assertions.R

# ---- app scripts ----
WORKDIR /opt/app
COPY --chown=$MAMBA_USER:$MAMBA_USER monitor_cgroup.sh /opt/app/monitor_cgroup.sh
COPY --chown=$MAMBA_USER:$MAMBA_USER run.sh /opt/app/run.sh
COPY --chown=$MAMBA_USER:$MAMBA_USER run_genescope_xenium.R /opt/app/run_genescope_xenium.R
RUN chmod +x /opt/app/run.sh /opt/app/monitor_cgroup.sh

ENTRYPOINT ["/opt/app/run.sh"]
