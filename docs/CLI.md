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
| `Database` | Initialize, edit, download, and build ADToolbox databases. |
| `Metagenomics` | Download genomes or SRA data and align genomes to protein databases. |
| `ADM` | Run ADM1 and e-ADM models. |
| `Documentations` | Print package documentation in the terminal. |

Every command supports `-h` and `--help`.

```bash
adtoolbox Database --help
adtoolbox ADM adm1 --help
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
adtoolbox Database initialize-feed-db --feed-db ./database/feed_db.tsv

adtoolbox Database add-feed \
  --feed-db ./database/feed_db.tsv \
  --name "food waste" \
  --carbohydrates 42 \
  --proteins 20 \
  --lipids 18 \
  --tss 80 \
  --si 5 \
  --xi 15 \
  --reference "example reference"

adtoolbox Database show-feed-db --feed-db ./database/feed_db.tsv
adtoolbox Database show-feed-db --feed-db ./database/feed_db.tsv --filter "food waste"
```

To download the full reference database bundle:

```bash
adtoolbox Database download-all-databases --output-dir ./database
```

To build a protein database from reaction metadata:

```bash
adtoolbox Database build-protein-db \
  --reaction-db ./database/Reaction_Metadata.csv \
  --protein-db ./database/Protein_DB.fasta
```

## Metagenomics

Metagenomics commands also take explicit inputs and output directories.

| Command | Purpose |
| --- | --- |
| `download_from_sra` | Download reads from SRA by sample accession. |
| `download_genome` | Download a genome from NCBI by genome accession. |
| `align-genome` | Align one genome to a protein FASTA database. |
| `align-multiple-genomes` | Align multiple genomes listed in a JSON manifest. |
| `find-representative-genomes` | Find representative genomes from a repseqs FASTA file. |

Use `--container None` for local execution, or `--container docker` / `--container singularity` when running through a container backend.

Examples:

```bash
adtoolbox Metagenomics download_from_sra \
  --sample-accession SRR28403133 \
  --output-dir ./metagenomics/sra \
  --container None

adtoolbox Metagenomics download_genome \
  --genome-accession GCA_021152825.1 \
  --output-dir ./metagenomics/genomes \
  --container None

adtoolbox Metagenomics align-genome \
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
adtoolbox Metagenomics align-multiple-genomes \
  --input-file ./metagenomics/genomes.json \
  --output-dir ./metagenomics/alignment \
  --protein-db ./database/Protein_DB.fasta \
  --container None
```

## ADM

The ADM CLI currently exposes two model families:

| Command | Model key | Purpose |
| --- | --- | --- |
| `adm1` | `adm1` | Run the original ADM1 model. |
| `e-adm` | `e_adm` | Run the extended e-ADM model. |

The recommended input format is one consolidated model JSON file keyed by model name. The repository includes a reference example at `reference_data/models.json`.

```bash
adtoolbox ADM adm1 --models-json reference_data/models.json --report csv
adtoolbox ADM e-adm --models-json reference_data/models.json --report csv
```

When `--report csv` is used, the CLI asks where to save the output CSV. When `--report dash` is used, or when `--report` is omitted, the CLI opens the interactive Dash visualization.

The e-ADM command can also accept a control-state JSON file:

```bash
adtoolbox ADM e-adm \
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

For compatibility with older database layouts, each ADM command can still load six separate JSON files:

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
adtoolbox ADM adm1 \
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
adtoolbox ADM adm1 --parameters-dir ./ADM_Parameters --report csv
adtoolbox ADM e-adm --parameters-dir ./ADM_Parameters --report csv
```

For `e-adm`, the CLI first looks for `e_adm_*.json` files and then falls back to the legacy `e_adm_2_*.json` file names.

## Documentation

Print the package README in the terminal:

```bash
adtoolbox Documentations --show
```

## Reference Data

This repository includes clean reference JSON files that mirror the current database shape:

| File | Contents |
| --- | --- |
| `reference_data/models.json` | ADM1 and e-ADM requirements keyed by model name. |
| `reference_data/feeds.json` | Feed entries keyed by normalized feed name. |
| `reference_data/experiments.json` | Experiment entries keyed by normalized experiment name. |

These files are meant as portable examples and test fixtures. Production runs can point the CLI at any equivalent local files.
