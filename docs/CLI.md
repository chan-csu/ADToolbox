# ADToolbox Command Line Interface

The ADToolbox command line interface is installed as `adtoolbox`.

```bash
adtoolbox --help
adtoolbox --version
```

The CLI does not create or store a global project directory. Commands that need files or directories accept those paths directly. If a required path is omitted, the command prompts for it.

## Modules

| Module | Purpose |
| --- | --- |
| `database` | Initialize, edit, download, and build ADToolbox databases. |
| `metagenomics` | Download genomes or SRA data and align genomes to protein databases. |
| `adm` | Run ADM1 and e-ADM models. |
| `docs` | Print package documentation in the terminal. |

Every command supports `-h` and `--help`.

```bash
adtoolbox database --help
adtoolbox adm adm1 --help
```

## Database

Database commands work with explicit file paths. The path can point to a local database you already maintain, a reference file copied from `reference_data`, or a new file you want ADToolbox to create.

| Command | Required path option | Purpose |
| --- | --- | --- |
| `initialize-feed-db` | `--feed-db` | Create an empty feed TSV. |
| `add-feed` | `--feed-db` | Add one feed row to a feed TSV. |
| `show-feed-db` | `--feed-db` | Print the feed TSV, optionally filtered by feed name. |
| `download-feed-db` | `--feed-db` | Download the reference feed TSV. |
| `initialize-metagenomics-studies-db` | `--studies-db` | Create an empty metagenomics studies TSV. |
| `add-metagenomics-study` | `--studies-db` | Add one metagenomics study row. |
| `initialize-protein-db` | `--protein-db` | Create an empty protein FASTA database. |
| `add-protein` | `--protein-db` | Add one protein FASTA entry. |
| `download-reaction-db` | `--reaction-db` | Download reaction metadata as CSV. |
| `download-seed-reaction-db` | `--seed-reaction-db`, `--seed-compound-db` | Download SEED reaction and compound JSON files. |
| `build-protein-db` | `--reaction-db`, `--protein-db` | Build a protein FASTA database from reaction metadata. |
| `download-protein-db` | `--protein-db` | Download the reference protein FASTA database. |
| `download-amplicon-to-genome-dbs` | `--output-dir` | Download amplicon-to-genome mapping databases. |
| `download-all-databases` | `--output-dir` | Download the standard ADToolbox database bundle. |

Examples:

```bash
adtoolbox database initialize-feed-db --feed-db ./database/feed_db.tsv

adtoolbox database add-feed \
  --feed-db ./database/feed_db.tsv \
  --name "food waste" \
  --carbohydrates 42 \
  --proteins 20 \
  --lipids 18 \
  --tss 80 \
  --si 5 \
  --xi 15 \
  --reference "example reference"

adtoolbox database show-feed-db --feed-db ./database/feed_db.tsv
adtoolbox database show-feed-db --feed-db ./database/feed_db.tsv --filter "food waste"
```

To download the full reference database bundle:

```bash
adtoolbox database download-all-databases --output-dir ./database
```

To build a protein database from reaction metadata:

```bash
adtoolbox database build-protein-db \
  --reaction-db ./database/Reaction_Metadata.csv \
  --protein-db ./database/Protein_DB.fasta
```

## Metagenomics

Metagenomics commands also take explicit inputs and output directories.

| Command | Purpose |
| --- | --- |
| `download-sra` | Download reads from SRA by sample accession. |
| `download-genome` | Download a genome from NCBI by genome accession. |
| `align-genome` | Align one genome to a protein FASTA database. |
| `align-multiple-genomes` | Align multiple genomes listed in a JSON manifest. |
| `find-representative-genomes` | Find representative genomes from a repseqs FASTA file. |
| `process` | Process a table of SRA accessions or local amplicon reads into model-ready e-ADM microbial COD allocations. |

Use `--container None` for local execution, or `--container docker` / `--container singularity` when running through a container backend.

Examples:

```bash
adtoolbox metagenomics download-sra \
  --sample-accession SRR28403133 \
  --output-dir ./metagenomics/sra \
  --container None

adtoolbox metagenomics download-genome \
  --genome-accession GCA_021152825.1 \
  --output-dir ./metagenomics/genomes \
  --container None

adtoolbox metagenomics align-genome \
  --name GCA_021152825_1 \
  --input-file ./metagenomics/genomes/GCA_021152825.1.fna \
  --output-dir ./metagenomics/alignment \
  --protein-db ./database/Protein_DB.fasta \
  --container None
```

For multiple genomes, the input JSON maps genome names to input files:

```json
{
  "genome_1": "./genomes/genome_1.fna",
  "genome_2": "./genomes/genome_2.fna"
}
```

Run the batch alignment with:

```bash
adtoolbox metagenomics align-multiple-genomes \
  --input-file ./metagenomics/genomes.json \
  --output-dir ./metagenomics/alignment \
  --protein-db ./database/Protein_DB.fasta \
  --container None
```

### Processing Pipeline

`process` is the table-driven pipeline for converting amplicon evidence into e-ADM microbial biomass/COD allocations. Each run creates one sample folder per table row under `--output-dir`. Clean result files stay at the sample-folder root, while generated commands, raw alignments, VSEARCH intermediates, and other working files are written under `scratch/`.

| File | Meaning |
| --- | --- |
| `cod_profile.csv` | Final normalized `X_*` allocation for the ADM model, as `sample`, `group`, `value`. |
| `ec_counts.csv` | EC counts when the mode produces direct EC evidence, as `sample`, `ec`, `count`. |
| `feature_abundances.csv` | Amplicon feature abundances, as `sample`, `feature_id`, `abundance`. |
| `representative_genomes.csv` | GTDB mapping from feature IDs to representative genomes. |
| `genome_abundances.csv` | Per-sample genome abundances after feature-to-genome aggregation. |
| `genome_cods.csv` | Genome-level `X_*` profiles when genomes are involved, as `sample`, `genome_id`, `group`, `value`. |
| `provenance.json` | Inputs, thresholds, databases, and artifacts used. |
| `pipeline.log` | Step-by-step log for the sample. |
| `scratch/` | Intermediate files, generated scripts, GTDB matches, and genome alignment files. |

By default, the command is a dry run for external tools: it parses existing files and writes commands for missing trimming, feature-building, and alignment steps. Add `--execute` to run fastp, VSEARCH, and MMseqs commands.

Execution behavior can be controlled with a TOML profile. The repository includes an example at `reference_data/metagenomics_pipeline.toml`.

```toml
backend = "local"
container = "apptainer"
image = "docker://parsaghadermazi/adtoolbox:latest"

[slurm]
# retries = 1
# retry_delay_seconds = 60

[steps.download_sra]
backend = "slurm"
container = "apptainer"
cpus = 4
memory = "16G"
time = "04:00:00"

[steps.trim_reads]
backend = "slurm"
container = "apptainer"
cpus = 4
memory = "8G"
time = "01:00:00"
# retries = 2

[steps.trim_reads.settings]
threads = 4
minimum_length = 100
# Optional explicit adapters. When omitted, fastp auto-detects common adapters.
# adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
# adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

[steps.build_amplicon_features]
backend = "slurm"
container = "apptainer"
cpus = 8
memory = "24G"
time = "04:00:00"

[steps.build_amplicon_features.settings]
threads = 8
identity = 0.97
maxee = 1.0

[steps.align_to_gtdb]
backend = "slurm"
container = "apptainer"
cpus = 8
memory = "32G"
time = "04:00:00"

[steps.align_to_gtdb.settings]
vsearch_threads = 8
vsearch_similarity = 0.97

[steps.align_genome]
backend = "slurm"
container = "apptainer"
cpus = 12
memory = "48G"
time = "08:00:00"

[steps.align_short_reads]
backend = "slurm"
container = "apptainer"
cpus = 24
memory = "150G"
time = "12:00:00"
```

For SRA accessions, the input table must include `sample` and `accession` columns:

```tsv
sample	accession
sample_01	SRR28403133
sample_02	SRR28403134
```

Run with:

```bash
adtoolbox metagenomics process \
  --input ./metagenomics/sra_samples.tsv \
  --input-type sra \
  --output-dir ./metagenomics/process \
  --sra-dir ./metagenomics/sra \
  --adapter-1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  --adapter-2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --amplicon-to-genome-db ./database/amplicon_to_genome \
  --genomes-dir ./metagenomics/genomes \
  --protein-db ./database/Protein_DB.fasta \
  --reaction-db ./database/Reaction_Metadata.csv \
  --execution-profile reference_data/metagenomics_pipeline.toml \
  --execute
```

For local FASTQ/FASTQ.GZ files, the input table must include `sample`, `read_1`, and optionally `read_2`:

```tsv
sample	read_1	read_2
sample_01	./fastq/sample_01_R1.fastq.gz	./fastq/sample_01_R2.fastq.gz
sample_02	./fastq/sample_02_R1.fastq.gz	./fastq/sample_02_R2.fastq.gz
```

```bash
adtoolbox metagenomics process \
  --input ./metagenomics/read_samples.tsv \
  --input-type reads \
  --output-dir ./metagenomics/process \
  --adapter-1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  --adapter-2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --amplicon-to-genome-db ./database/amplicon_to_genome \
  --genomes-dir ./metagenomics/genomes \
  --protein-db ./database/Protein_DB.fasta \
  --reaction-db ./database/Reaction_Metadata.csv \
  --execution-profile reference_data/metagenomics_pipeline.toml \
  --execute
```

On Slurm, the safer pattern is staged execution:

```bash
adtoolbox metagenomics process --input ./metagenomics/sra_samples.tsv --input-type sra --stage download --output-dir ./metagenomics/process --sra-dir ./metagenomics/sra --execution-profile reference_data/metagenomics_pipeline.toml --execute
adtoolbox metagenomics process --input ./metagenomics/sra_samples.tsv --input-type sra --stage preprocess --output-dir ./metagenomics/process --sra-dir ./metagenomics/sra --execution-profile reference_data/metagenomics_pipeline.toml --execute
adtoolbox metagenomics process --input ./metagenomics/sra_samples.tsv --input-type sra --stage allocate --output-dir ./metagenomics/process --genomes-dir ./metagenomics/genomes --reaction-db ./database/Reaction_Metadata.csv --execution-profile reference_data/metagenomics_pipeline.toml --execute
```

Each batch run writes `batch_summary.json` under `--output-dir`.

## ADM

The ADM CLI currently exposes two model families:

| Command | Model key | Purpose |
| --- | --- | --- |
| `adm1` | `adm1` | Run the original ADM1 model. |
| `e-adm` | `e_adm` | Run the extended e-ADM model. |

The recommended input format is one consolidated model JSON file keyed by model name. The repository includes a reference example at `reference_data/models.json`.

```bash
adtoolbox adm adm1 --models-json reference_data/models.json --report csv
adtoolbox adm e-adm --models-json reference_data/models.json --report csv
```

When `--report csv` is used, the CLI asks where to save the output CSV. When `--report dash` is used, or when `--report` is omitted, the CLI opens the interactive Dash visualization.

The e-ADM command can also accept a control-state JSON file:

```bash
adtoolbox adm e-adm \
  --models-json reference_data/models.json \
  --control-states ./control_states.json \
  --report csv
```

If `--control-states` is omitted, e-ADM uses:

```json
{
  "S_H_ion": 3.162277660168379e-7
}
```

The consolidated model JSON has this structure:

```json
{
  "adm1": {
    "model_parameters": {},
    "base_parameters": {},
    "initial_conditions": {},
    "inlet_conditions": {},
    "reactions": {},
    "species": {}
  },
  "e_adm": {
    "model_parameters": {},
    "base_parameters": {},
    "initial_conditions": {},
    "inlet_conditions": {},
    "reactions": {},
    "species": {}
  }
}
```

Each ADM command can also load six separate JSON files:

| Option | Contents |
| --- | --- |
| `--model-parameters` | Kinetic and model-specific parameters. |
| `--base-parameters` | Shared physical and biochemical constants. |
| `--initial-conditions` | Initial state values. |
| `--inlet-conditions` | Influent state values. |
| `--reactions` | Reaction names and ordering. |
| `--species` | Species names and ordering. |

You can pass those files directly:

```bash
adtoolbox adm adm1 \
  --model-parameters ./ADM_Parameters/adm1_model_parameters.json \
  --base-parameters ./ADM_Parameters/adm1_base_parameters.json \
  --initial-conditions ./ADM_Parameters/adm1_initial_conditions.json \
  --inlet-conditions ./ADM_Parameters/adm1_inlet_conditions.json \
  --reactions ./ADM_Parameters/adm1_reactions.json \
  --species ./ADM_Parameters/adm1_species.json \
  --report csv
```

Or pass a directory containing consistently named files:

```bash
adtoolbox adm adm1 --parameters-dir ./ADM_Parameters --report csv
adtoolbox adm e-adm --parameters-dir ./ADM_Parameters --report csv
```

## Documentation

Print the package README in the terminal:

```bash
adtoolbox docs --show
```

## Reference Data

This repository includes clean reference JSON files that mirror the current database shape:

| File | Contents |
| --- | --- |
| `reference_data/models.json` | ADM1 and e-ADM requirements keyed by model name. |
| `reference_data/feeds.json` | Feed entries keyed by normalized feed name. |
| `reference_data/experiments.json` | Experiment entries keyed by normalized experiment name. |

These files are meant as portable examples and test fixtures. Production runs can point the CLI at any equivalent local files.
