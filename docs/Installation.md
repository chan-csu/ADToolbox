# Installation

ADToolbox runs as a locally installed Python package or from a container image. Anything
that shells out to an external bioinformatics tool works either way: put the tool on your
`PATH`, or let ADToolbox call it inside Docker or Apptainer/Singularity. Pipeline steps can
additionally be submitted to Slurm, which is how the toolbox scales to HPC clusters.

!!! tip "Use a virtual environment"
    ADToolbox pins several scientific dependencies. Install it into a dedicated
    environment rather than your system Python:

    ```bash
    python -m venv .venv && source .venv/bin/activate   # or conda/mamba
    ```

**Requirements:** Python 3.11 or newer.

## Install the package

=== "conda / mamba (recommended)"

    The environment file installs ADToolbox **and every external tool** both the
    amplicon and shotgun routes need, so nothing else has to be on your `PATH`:

    ```bash
    curl -O https://raw.githubusercontent.com/chan-csu/ADToolbox/main/environment.yml
    conda env create -f environment.yml   # or: mamba env create -f environment.yml
    conda activate adtoolbox
    ```

    This is the least-fuss way to get a fully working pipeline — it pulls fastp,
    Cutadapt, DADA2, VSEARCH, MMseqs2, SRA Toolkit, and the NCBI Datasets CLI
    alongside the Python package, so you can skip the [External tools](#external-tools)
    section entirely.

    !!! note "A `conda install -c bioconda adtoolbox` one-liner is on the way"
        The [Bioconda recipe](https://github.com/bioconda/bioconda-recipes/pull/67450)
        is under review. Once it merges, `conda install -c bioconda adtoolbox`
        will install the package directly; until then, use the environment file
        above.

=== "pip"

    ```bash
    pip install adtoolbox
    ```

    Installs the Python package only. The metagenomics pipeline additionally
    needs the [external tools](#external-tools) on your `PATH`, or a container.

=== "From source"

    ```bash
    git clone https://github.com/chan-csu/ADToolbox.git
    cd ADToolbox
    pip install .
    ```

=== "Editable (development)"

    ```bash
    git clone https://github.com/chan-csu/ADToolbox.git
    cd ADToolbox
    pip install -e .
    ```

Verify the install:

```bash
adtoolbox --version
adtoolbox --help
```

## Optional extras

The base install covers databases, the metagenomics pipeline, and ADM simulation. Features
with heavier dependencies are opt-in:

| Extra | Install | Enables |
| --- | --- | --- |
| `dashboard` | `pip install "adtoolbox[dashboard]"` | Interactive Dash visualization, including Escher maps (`--report dash`). |
| `blackbox` | `pip install "adtoolbox[blackbox]"` | [`BlackBoxOptimizer`](api-optimize.md#adtoolbox.optimize.BlackBoxOptimizer) via OpenBox. |
| `genetic` | `pip install "adtoolbox[genetic]"` | [`GeneticOptimizer`](api-optimize.md#adtoolbox.optimize.GeneticOptimizer) via PyGAD. |
| `surrogate` | `pip install "adtoolbox[surrogate]"` | [`SurrogateOptimizer`](api-optimize.md#adtoolbox.optimize.SurrogateOptimizer) via PyTorch. |
| `optimize` | `pip install "adtoolbox[optimize]"` | All three optimizer backends at once. |

Combine them as needed:

```bash
pip install "adtoolbox[optimize,dashboard]"
```

## External tools

These are only required for the metagenomics pipeline. If you plan to use containers, skip
this section — the ADToolbox image already contains all of them. The `conda / mamba`
install above also pulls every one of them.

| Tool | Used for | Assay |
| --- | --- | --- |
| [fastp](https://github.com/OpenGene/fastp) | Adapter trimming and quality filtering. | amplicon + shotgun |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | Protein alignment against the ADToolbox enzyme database. | amplicon + shotgun |
| [SRA Toolkit](https://github.com/ncbi/sra-tools) | `prefetch` and `fasterq-dump` for downloading reads from SRA. | both (SRA input) |
| [Cutadapt](https://cutadapt.readthedocs.io/) | IUPAC-aware PCR primer removal. | amplicon only |
| [DADA2](https://benjjneb.github.io/dada2/) and R | Error learning, ASV inference, paired-read merging, and chimera removal. | amplicon only |
| [VSEARCH](https://github.com/torognes/vsearch) | Mapping representative ASVs to GTDB; also the legacy feature backend. | amplicon only |
| [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/) | Batched assembly downloads. | amplicon (genomes) |

!!! tip "Shotgun-only setups are lighter"
    The shotgun route needs just **fastp** and **MMseqs2** (plus SRA Toolkit for
    `--input-type sra`). If you never run the amplicon route, you can skip
    Cutadapt, DADA2/R, VSEARCH, and the NCBI Datasets CLI.

The quickest local install is conda/mamba:

```bash
mamba install -c conda-forge -c bioconda \
  fastp cutadapt r-base r-digest bioconductor-dada2 \
  vsearch mmseqs2 ncbi-datasets-cli sra-tools
```

This installs standalone DADA2 rather than QIIME2. Set `container = "None"` in the execution profile so Slurm jobs use these commands from the active environment.

## Containers

The published image bundles ADToolbox with every external tool.

=== "Docker"

    ```bash
    docker run --rm parsaghadermazi/adtoolbox:latest adtoolbox --help

    # Mount the current directory so the container can read and write your data
    docker run --rm -v "$PWD:/workspace" parsaghadermazi/adtoolbox:latest \
      adtoolbox metagenomics process --help
    ```

=== "Build locally"

    ```bash
    git clone https://github.com/chan-csu/ADToolbox.git
    cd ADToolbox
    docker build -t adtoolbox:local .
    docker run --rm -v "$PWD:/workspace" adtoolbox:local adtoolbox --help
    ```

=== "Apptainer / Singularity"

    ```bash
    apptainer pull adtoolbox.sif docker://parsaghadermazi/adtoolbox:latest
    apptainer exec adtoolbox.sif adtoolbox --help
    ```

You do not have to run ADToolbox itself inside the container. A local install can call the
container for individual pipeline steps with `--container docker` or
`--container singularity`, which is usually the more convenient arrangement:

```bash
adtoolbox metagenomics align-genome \
  --name my_genome \
  --input-file ./genomes/my_genome.fna \
  --output-dir ./alignment \
  --protein-db ./database/Protein_DB.fasta \
  --container docker
```

## HPC and Slurm

For cluster runs, describe each pipeline step in a TOML execution profile that selects the
backend, container, and resource request. ADToolbox then generates and submits the Slurm
scripts, waits for them with `sbatch --wait`, and can retry failed jobs.

```toml
backend = "local"
container = "None"

[steps.align_short_reads]
backend = "slurm"
cpus = 24
memory = "150G"
time = "12:00:00"
```

A complete reference profile ships at `reference_data/metagenomics_pipeline.toml`. See
[Execution profiles](Metagenomics_Pipeline.md#execution-profiles) for every option.

## Run without installing

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/chan-csu/ADToolbox/HEAD)

Binder launches the [example notebooks](Notebooks.md) in a browser with no local setup.
Escher map visualization is not available there.

## Next steps

Head to the [Quickstart](Quickstart.md) to download the reference databases and run your
first simulation.
