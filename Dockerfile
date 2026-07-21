FROM mambaorg/micromamba:1.5.10-bookworm-slim

LABEL org.opencontainers.image.title="ADToolbox"
LABEL org.opencontainers.image.description="ADToolbox with metagenomics pipeline dependencies"
LABEL org.opencontainers.image.version="1.1.13"

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

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 \
        pip \
        cutadapt \
        fastp \
        mmseqs2 \
        ncbi-datasets-cli \
        r-base \
        r-digest \
        bioconductor-dada2 \
        sra-tools \
        vsearch \
    && micromamba clean --all --yes

COPY pyproject.toml README.md /tmp/adtoolbox/
COPY adtoolbox /tmp/adtoolbox/adtoolbox

RUN micromamba run -n base pip install --no-cache-dir /tmp/adtoolbox \
    && micromamba run -n base pip check

USER root
RUN ln -sf /opt/conda/bin/adtoolbox /usr/local/bin/adtoolbox \
    && ln -sf /opt/conda/bin/dataformat /usr/local/bin/dataformat \
    && ln -sf /opt/conda/bin/datasets /usr/local/bin/datasets \
    && ln -sf /opt/conda/bin/cutadapt /usr/local/bin/cutadapt \
    && ln -sf /opt/conda/bin/fasterq-dump /usr/local/bin/fasterq-dump \
    && ln -sf /opt/conda/bin/fastp /usr/local/bin/fastp \
    && ln -sf /opt/conda/bin/mmseqs /usr/local/bin/mmseqs \
    && ln -sf /opt/conda/bin/prefetch /usr/local/bin/prefetch \
    && ln -sf /opt/conda/bin/Rscript /usr/local/bin/Rscript \
    && ln -sf /opt/conda/bin/vsearch /usr/local/bin/vsearch

USER $MAMBA_USER
WORKDIR /workspace

CMD ["adtoolbox", "--help"]
