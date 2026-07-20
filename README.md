# ADToolbox

[![Documentation](https://img.shields.io/badge/docs-chan--csu.github.io-teal)](https://chan-csu.github.io/ADToolbox/)
[![PyPI version](https://badge.fury.io/py/adtoolbox.svg)](https://badge.fury.io/py/adtoolbox)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/chan-csu/ADToolbox/HEAD)

**From raw amplicon reads to a calibrated anaerobic digestion model.**

ADToolbox is developed in the Chan Lab at Colorado State University. It connects
metagenomics evidence, curated reaction and feed databases, and dynamic ADM1 / e-ADM
simulations into one reproducible Python and command-line workflow.

Anaerobic digestion models such as ADM1 lump the microbial community into a handful of
guilds whose initial biomass you are expected to guess. ADToolbox replaces that guess with
measurement: it takes 16S amplicon data, maps it through GTDB and a curated
enzyme-to-reaction database, and produces the microbial COD allocation the model needs.

## 📖 [Read the documentation](https://chan-csu.github.io/ADToolbox/)

## Install

```bash
pip install adtoolbox               # base install, Python 3.11+
pip install "adtoolbox[optimize]"   # + parameter-tuning backends
pip install "adtoolbox[dashboard]"  # + interactive Dash/Escher visualization
```

Or use the container, which bundles fastp, VSEARCH, MMseqs2, and the SRA Toolkit:

```bash
docker run --rm parsaghadermazi/adtoolbox:latest adtoolbox --help
```

See the [installation guide](https://chan-csu.github.io/ADToolbox/Installation/) for
source installs, extras, external tools, and HPC notes.

## Quick start

```bash
# 1. Download the reference databases
adtoolbox database download-all-databases --output-dir ./database

# 2. Turn amplicon samples into microbial COD allocations
adtoolbox metagenomics process \
  --input ./samples.tsv --input-type sra \
  --output-dir ./process --sra-dir ./sra \
  --database-dir ./database --execute

# 3. Simulate
adtoolbox adm e-adm --models-json reference_data/models.json --report csv
```

The full walkthrough is in the
[Quickstart](https://chan-csu.github.io/ADToolbox/Quickstart/).

## Documentation map

| Page | Contents |
| --- | --- |
| [Quickstart](https://chan-csu.github.io/ADToolbox/Quickstart/) | Install, databases, first simulation, first pipeline run. |
| [Metagenomics Pipeline](https://chan-csu.github.io/ADToolbox/Metagenomics_Pipeline/) | Sample tables, execution profiles, Slurm, output schemas. |
| [ADM Models](https://chan-csu.github.io/ADToolbox/ADM_Models/) | Input contracts, stoichiometry, rate laws, inhibition terms. |
| [Parameter Tuning](https://chan-csu.github.io/ADToolbox/Optimization/) | Fitting kinetic parameters to experimental data. |
| [CLI](https://chan-csu.github.io/ADToolbox/CLI/) | Every command and option. |
| [Python API](https://chan-csu.github.io/ADToolbox/API/) | Generated reference for all modules. |

## Try it without installing

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/chan-csu/ADToolbox/HEAD)

The [example notebooks](Examples/) run on Binder. Note that Escher map visualization is
not available there.

## Building the docs locally

```bash
pip install -r docs/requirements.txt
mkdocs serve
```

## Contact

Developed in the [Chan Lab](https://github.com/chan-csu) at Colorado State University.

| | |
| --- | --- |
| Parsa Ghadermazi | parsa.ghadermazi@colostate.edu |
| Ethan Rimelman | rimelman@colostate.edu |
| Siu Hung Joshua Chan (PI) | joshua.chan@colostate.edu |
