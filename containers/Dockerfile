# Container image for nf-g4rnaseq (R analysis stack + tools).
# Build (linux/amd64 = nf-core / HPC / cloud default):
#   docker buildx build --platform linux/amd64 -t ghcr.io/spegray/nf-g4rnaseq-r:0.1.0 \
#       -f containers/Dockerfile containers/
# Then set params.r_container (nextflow.config) to the pushed tag/digest.
#
# This bundles the dependencies that were painful to build from source on some
# platforms (GenomeInfoDbData, igraph, gfortran-free binaries) so users on any
# machine get identical, reproducible results.
FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="nf-g4rnaseq-r" \
      org.opencontainers.image.description="R/Bioconductor + tooling for the nf-g4rnaseq pipeline"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH
# sanity check at build time
RUN R -e 'stopifnot(all(sapply(c("DESeq2","pqsfinder","clusterProfiler","Biostrings"), requireNamespace, quietly=TRUE)))' && \
    which featureCounts samtools STAR fastp >/dev/null
