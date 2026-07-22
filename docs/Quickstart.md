# Quickstart

This page takes you from a clean environment to a running ADM simulation, and then to a
full metagenomics-to-model run. Each step is independent, so you can stop after the
simulation if you only care about modeling.

!!! tip "Prerequisites"
    Python 3.11 or newer. External bioinformatics tools (fastp, VSEARCH, MMseqs2,
    SRA Toolkit) are only needed for the metagenomics pipeline, and can be supplied by the
    ADToolbox container instead of installing them yourself.

## 1. Install

```bash
pip install adtoolbox
adtoolbox --version
```

Confirm the command groups are available:

```bash
adtoolbox --help
```

```title="Output"
Commands:
  adm           Run and visualize ADToolbox ADM models.
  database      Build or download databases required by ADToolbox.
  docs          Print package documentation in the terminal.
  metagenomics  Import and process metagenomics data from the command line.
```

## 2. Download the reference databases

ADToolbox does not keep a global project directory. You choose a directory, and every
command points at it explicitly.

```bash
adtoolbox database download-all-databases --output-dir ./database
```

This fetches the SEED reaction and compound databases, ADM parameter sets, the curated
protein FASTA database, reaction metadata, the feed database, study metadata, and the
amplicon-to-genome (GTDB) files.

??? note "What ends up in `./database`"

    | Path | Contents |
    | --- | --- |
    | `Protein_DB.fasta` | Curated enzyme sequences keyed by EC number. |
    | `Reaction_Metadata.csv` | EC-number to e-ADM reaction mapping. |
    | `reactions.json`, `compounds.json` | SEED reaction and compound databases. |
    | `feed_db.tsv` | Feed composition database. |
    | `Amplicon2GenomeDBs/` | GTDB SSU sequences and taxonomy for amplicon mapping. |
    | `ADM_Parameters/models.json` | ADM1 and e-ADM parameter sets. |
    | `Studies/` | Metagenomics study and experimental data references. |

    The GTDB download is the large one. If you only want to run models, you can skip this
    command entirely and use the parameter files in `reference_data/` instead.

## 3. Run your first ADM simulation

The repository ships a ready-to-run parameter bundle at `reference_data/models.json`,
which contains both `adm1` and `e_adm` keys.

=== "Interactive dashboard"

    ```bash
    adtoolbox adm e-adm --models-json reference_data/models.json --report dash
    ```

    This opens a Dash app with concentration time courses. It requires the dashboard
    extra: `pip install "adtoolbox[dashboard]"`.

=== "CSV output"

    ```bash
    adtoolbox adm e-adm --models-json reference_data/models.json --report csv
    ```

    The CLI prompts for the output path and writes one row per time point.

=== "Original ADM1"

    ```bash
    adtoolbox adm adm1 --models-json reference_data/models.json --report csv
    ```

The same thing from Python, which is what you want if you plan to script anything:

```python
import numpy as np
from adtoolbox import adm, utils

payload = utils.load_model_json("reference_data/models.json", "e_adm")

model = adm.Model(
    model_parameters=payload["model_parameters"],
    base_parameters=payload["base_parameters"],
    initial_conditions=payload["initial_conditions"],
    inlet_conditions=payload["inlet_conditions"],
    reactions=payload["reactions"],
    species=payload["species"],
    feed=adm.DEFAULT_FEED,
    ode_system=adm.e_adm_ode_sys,
    build_stoichiometric_matrix=adm.build_e_adm_stoichiometric_matrix,
    control_state={"S_H_ion": 10 ** -6.5},
    name="e-ADM",
    switch="DAE",
)

solution = model.solve_model(np.linspace(0, 30, 300))
model.plot(solution).show()
```

See the [ADM model reference](ADM_Models.md) for the meaning of every parameter file and
the full stoichiometry.

## 4. Turn amplicon data into microbial COD

This is the part that distinguishes ADToolbox from a plain ADM implementation. Describe
your samples in a table, then run one command.

=== "SRA accessions"

    ```tsv title="samples.tsv"
    sample	accession
    sample_01	SRR28403133
    sample_02	SRR28403134
    ```

    ```bash
    adtoolbox metagenomics process \
      --input ./samples.tsv \
      --input-type sra \
      --assay amplicon \
      --output-dir ./process \
      --sra-dir ./sra \
      --database-dir ./database \
      --execute
    ```

=== "Local FASTQ files"

    ```tsv title="samples.tsv"
    sample	read_1	read_2
    sample_01	./fastq/sample_01_R1.fastq.gz	./fastq/sample_01_R2.fastq.gz
    sample_02	./fastq/sample_02_R1.fastq.gz	./fastq/sample_02_R2.fastq.gz
    ```

    ```bash
    adtoolbox metagenomics process \
      --input ./samples.tsv \
      --input-type reads \
      --assay amplicon \
      --output-dir ./process \
      --database-dir ./database \
      --execute
    ```

For shotgun functional profiling, use the same sample-table formats and select the shotgun assay. This route runs fastp and one MMseqs2 translated-search job per sample, then writes `ec_counts.csv` and `cod_profile.csv` without DADA2, GTDB, or genome downloads:

```bash
adtoolbox metagenomics process \
  --input ./samples.tsv \
  --input-type reads \
  --assay shotgun \
  --output-dir ./process \
  --protein-db ./database/Protein_DB.fasta \
  --reaction-db ./database/Reaction_Metadata.csv \
  --execution-profile ./reference_data/metagenomics_pipeline.toml \
  --execute
```

!!! warning "`--execute` is opt-in"
    Without `--execute` the pipeline runs as a dry run: it parses whatever files already
    exist and writes the shell scripts for the missing external steps, but does not run
    fastp, VSEARCH, MMseqs2, or the SRA tools. Start with a dry run to inspect the
    generated commands, then re-run with `--execute`.

Each sample gets its own folder under `--output-dir`. The file you feed back into the
model is `cod_profile.csv`:

```csv title="process/sample_01/cod_profile.csv"
sample,group,value
sample_01,X_su,0.184
sample_01,X_aa,0.121
sample_01,X_ac,0.093
...
```

The [pipeline guide](Metagenomics_Pipeline.md) covers execution profiles, Slurm
submission, staged runs, and every output file.

## 5. Feed the allocation back into the model

The `group` column matches e-ADM biomass state names, so the profile can be applied
directly as initial conditions:

```python
import polars as pl

profile = pl.read_csv("./process/sample_01/cod_profile.csv")
total_biomass = 2.0  # gCOD/L of active biomass in the reactor

model.update_parameters(
    initial_conditions={
        row["group"]: row["value"] * total_biomass
        for row in profile.iter_rows(named=True)
    }
)

solution = model.solve_model(np.linspace(0, 30, 300))
```

## 6. Calibrate against experimental data

If you have measurements, fit the kinetic parameters instead of accepting the defaults:

```python
from adtoolbox import optimize

optimizer = optimize.ScipyOptimizer(
    base_model=model,
    train_data=experiments,          # list of core.Experiment
    search_space={
        "k_m_su": (5.0, 60.0),
        "K_S_su": (0.1, 1.0),
        "k_m_ac": (2.0, 20.0),
    },
)

result = optimizer.optimize(maxiter=50)
print(optimizer.best_cost, optimizer.optimized_parameters)
```

The [parameter tuning guide](Optimization.md) explains the four optimizer backends, how
the search space maps onto model parameters, and how to validate the fit.

## Where to go next

<div class="grid cards" markdown>

-   :material-dna: **Run the real pipeline**

    ---

    Execution profiles, Slurm, staged runs, and output schemas.

    [:octicons-arrow-right-24: Metagenomics pipeline](Metagenomics_Pipeline.md)

-   :material-chart-bell-curve: **Understand the model**

    ---

    Stoichiometry, rate laws, inhibition terms, and state vectors.

    [:octicons-arrow-right-24: ADM models](ADM_Models.md)

-   :material-notebook: **Follow a notebook**

    ---

    End-to-end worked examples you can run locally or on Binder.

    [:octicons-arrow-right-24: Example notebooks](Notebooks.md)

-   :material-console: **Look up a command**

    ---

    Every CLI command with its options and file expectations.

    [:octicons-arrow-right-24: CLI reference](CLI.md)

</div>
