FROM mambaorg/micromamba:1.5.10-bookworm-slim

LABEL org.opencontainers.image.title="ADToolbox"
LABEL org.opencontainers.image.description="ADToolbox with metagenomics pipeline dependencies"

USER root
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/opt/conda/bin:${PATH}

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        curl \
        git \
        gzip \
        procps \
        rsync \
        tar \
        unzip \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md /tmp/adtoolbox/
COPY --chown=$MAMBA_USER:$MAMBA_USER adtoolbox /tmp/adtoolbox/adtoolbox

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 \
        pip \
        fastp \
        mmseqs2 \
        sra-tools \
        vsearch \
    && micromamba clean --all --yes

RUN micromamba run -n base pip install --no-cache-dir /tmp/adtoolbox \
    && micromamba run -n base pip check

USER root
RUN ln -sf /opt/conda/bin/ADToolbox /usr/local/bin/ADToolbox \
    && ln -sf /opt/conda/bin/fasterq-dump /usr/local/bin/fasterq-dump \
    && ln -sf /opt/conda/bin/fastp /usr/local/bin/fastp \
    && ln -sf /opt/conda/bin/mmseqs /usr/local/bin/mmseqs \
    && ln -sf /opt/conda/bin/prefetch /usr/local/bin/prefetch \
    && ln -sf /opt/conda/bin/vsearch /usr/local/bin/vsearch \
    && printf '#!/usr/bin/env bash\nexec ADToolbox "$@"\n' > /usr/local/bin/adtoolbox \
    && chmod +x /usr/local/bin/adtoolbox

USER $MAMBA_USER
WORKDIR /workspace

CMD ["adtoolbox", "--help"]
