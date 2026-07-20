# Metagenomics Pipeline

ADToolbox converts metagenomics evidence into model-ready e-ADM microbial COD allocations from a sample table. The current pipeline is exposed through the `adtoolbox Metagenomics process` CLI command and supports:

- SRA accession tables
- local FASTQ/FASTQ.GZ read tables
- raw amplicon reads processed with fastp and VSEARCH

Each sample row writes a dedicated output folder containing clean result CSVs, `pipeline.log`, and `provenance.json`. Generated command scripts, raw alignment files, VSEARCH intermediates, and other working files are kept under that sample's `scratch/` folder.

## Amplicon Reads

Raw amplicon reads are handled in three stages:

1. Trim adapters and short reads with fastp. Explicit adapters can be supplied, otherwise fastp auto-detects common adapters.
2. Build a feature table and representative sequence FASTA with VSEARCH.
3. Map representative sequences to GTDB, connect genomes to ADToolbox protein alignments, and aggregate the resulting e-ADM COD allocation.

The feature-generation path is intentionally opinionated: VSEARCH UNOISE-style denoising is used to produce representative sequences and a feature table with minimal dependencies.

For local reads, use a table with `sample`, `read_1`, and optionally `read_2`:

```tsv
sample	read_1	read_2
sample_01	./fastq/sample_01_R1.fastq.gz	./fastq/sample_01_R2.fastq.gz
sample_02	./fastq/sample_02_R1.fastq.gz	./fastq/sample_02_R2.fastq.gz
```

For SRA samples, use a table with `sample` and `accession`:

```tsv
sample	accession
sample_01	SRR28403133
sample_02	SRR28403134
```

```bash
adtoolbox Metagenomics process \
  --input ./metagenomics/read_samples.tsv \
  --input-type reads \
  --output-dir ./metagenomics/process \
  --adapter-1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  --adapter-2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --amplicon-to-genome-db ./database/amplicon_to_genome \
  --genomes-dir ./metagenomics/genomes \
  --reaction-db ./database/Reaction_Metadata.csv \
  --execution-profile reference_data/metagenomics_pipeline.toml \
  --execute
```

For SRA tables, use the same command with `--input-type sra` and an SRA download directory:

```bash
adtoolbox Metagenomics process \
  --input ./metagenomics/sra_samples.tsv \
  --input-type sra \
  --output-dir ./metagenomics/process \
  --sra-dir ./metagenomics/sra \
  --amplicon-to-genome-db ./database/amplicon_to_genome \
  --genomes-dir ./metagenomics/genomes \
  --reaction-db ./database/Reaction_Metadata.csv \
  --execution-profile reference_data/metagenomics_pipeline.toml \
  --execute
```

## Execution Profiles

Use a TOML execution profile to choose local or Slurm execution per step. The reference profile is `reference_data/metagenomics_pipeline.toml`.

Important step names are:

- `download_sra`
- `trim_reads`
- `build_amplicon_features`
- `align_to_gtdb`
- `align_genome`
- `align_short_reads`

Without `--execute`, ADToolbox writes scripts and Slurm files but does not run or submit them.

## Python Step API

Each external step can also be called directly from Python. The direct step methods use the same TOML execution profile and return the same artifact structure used by the CLI.

```python
from adtoolbox import configs, core

mg = core.Metagenomics(configs.Metagenomics(database_dir="./database"))

trim = mg.run_trim_reads_step(
    sample_name="sample_01",
    output_dir="./metagenomics/process",
    read_1="./metagenomics/sample_01/R1.fastq.gz",
    read_2="./metagenomics/sample_01/R2.fastq.gz",
    forward_primer="GTGYCAGCMGCCGCGGTAA",
    reverse_primer="GGACTACNVGGGTWTCTAAT",
    execution_profile="reference_data/metagenomics_pipeline.toml",
    execute=False,
)

features = mg.run_build_amplicon_features_step(
    sample_name="sample_01",
    output_dir="./metagenomics/process",
    read_1=trim["artifacts"]["trimmed_reads"]["read_1"],
    read_2=trim["artifacts"]["trimmed_reads"]["read_2"],
    execution_profile="reference_data/metagenomics_pipeline.toml",
    execute=False,
)
```

Available direct step methods:

- `run_trim_reads_step`
- `run_build_amplicon_features_step`
- `run_sra_download_step`
- `run_gtdb_alignment_step`
- `run_genome_alignment_step`
- `run_short_read_alignment_step`
