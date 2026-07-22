# ADToolbox Command Line Interface

ADToolbox turns metagenomics evidence into model-ready anaerobic digestion
simulations. This is the terminal summary; the full documentation lives at
https://chan-csu.github.io/ADToolbox/

## No setup step

There is no base directory to initialize and no global configuration. Every
command takes the paths it needs as options. If a required path is omitted, the
command prompts for it.

## Modules

| Module | Purpose |
| --- | --- |
| `database` | Initialize, edit, download, and build ADToolbox databases. |
| `metagenomics` | Download reads or genomes, align them, and run the processing pipeline. |
| `adm` | Run and visualize the ADM1 and e-ADM models. |
| `docs` | Print this documentation in the terminal. |

Every command and subcommand supports `-h` / `--help`:

```
adtoolbox --help
adtoolbox database --help
adtoolbox metagenomics process --help
adtoolbox adm e-adm --help
```

## Typical workflow

1. Download the reference databases:

```
adtoolbox database download-all-databases --output-dir ./database
```

2. Convert amplicon samples into microbial COD allocations. The input is a
   CSV/TSV table with either a `sample` + `accession` column pair for SRA, or
   `sample` + `read_1` (+ optional `read_2`) for local FASTQ files:

```
adtoolbox metagenomics process \
  --input ./samples.tsv \
  --input-type sra \
  --output-dir ./process \
  --sra-dir ./sra \
  --database-dir ./database \
  --execute
```

Without `--execute` the pipeline is a dry run: it parses existing files and
writes the shell scripts for the missing external steps without running fastp,
VSEARCH, MMseqs2, or the SRA tools. Each sample gets its own output folder
containing `cod_profile.csv`, `provenance.json`, `pipeline.log`, and a
`scratch/` folder of intermediates.

On a cluster, pass `--execution-profile` a TOML file to run individual steps
through Slurm, and run one `--stage` at a time (`download`, `preprocess`,
`allocate`).

3. Run a model:

```
adtoolbox adm e-adm --models-json reference_data/models.json --report csv
adtoolbox adm adm1 --parameters-dir ./ADM_Parameters --report dash
```

Use `--report dash` for the interactive dashboard (requires
`pip install "adtoolbox[dashboard]"`) or `--report csv` to write a table.
Models can be loaded from one consolidated JSON file via `--models-json`, from a
directory of conventionally named files via `--parameters-dir`, or from six
explicit `--model-parameters` / `--base-parameters` / `--initial-conditions` /
`--inlet-conditions` / `--reactions` / `--species` paths.

## Containers

Commands that call external tools accept `--container docker` or
`--container singularity` to run those tools inside the ADToolbox image instead
of requiring them on your PATH. `--container None` (the default) runs them
locally.

## Learn more

- Full documentation: https://chan-csu.github.io/ADToolbox/
- Source and issues: https://github.com/chan-csu/ADToolbox
