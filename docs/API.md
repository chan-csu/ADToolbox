# Python API

Everything the CLI does is available from Python. This page explains how the modules fit
together; the per-module pages hold the generated reference.

<div class="grid cards" markdown>

-   :material-cog: **[`configs`](api-configs.md)**

    ---

    Where files live and which thresholds to use.

-   :material-database: **[`core`](api-core.md)**

    ---

    Feeds, experiments, databases, and the metagenomics pipeline.

-   :material-chart-bell-curve: **[`adm`](api-adm.md)**

    ---

    ADM1 and e-ADM construction, integration, and plotting.

-   :material-tune-variant: **[`optimize`](api-optimize.md)**

    ---

    Parameter estimation and validation against experiments.

-   :material-wrench: **[`utils`](api-utils.md)**

    ---

    FASTA helpers, MMseqs2 wrappers, Slurm scripts.

-   :material-sigma: **[`stats`](api-stats.md)**

    ---

    Distances and scaling for feature tables.

</div>

## The config/core pairing

ADToolbox has **no global project directory and no hidden state**. Every configuration
object derives its file paths from a directory you pass in, and each `core` class is
constructed with the matching config object.

```python
from adtoolbox import configs, core

metagenomics_config = configs.Metagenomics(
    "./my_metagenomics_run",       # (1)!
    database_dir="./my_database",  # (2)!
)
metagenomics = core.Metagenomics(metagenomics_config)
```

1. Run directory. Genomes, alignment outputs, SRA downloads, and generated scripts are
   written under here.
2. Database directory. The protein FASTA, reaction metadata, and GTDB files are read from
   here.

Every `core.Metagenomics` method now resolves its paths under those two directories. The
same pattern applies to the other classes:

| Config | Paired with | Governs |
| --- | --- | --- |
| [`configs.Database`](api-configs.md#adtoolbox.configs.Database) | [`core.Database`](api-core.md#adtoolbox.core.Database), [`core.SeedDB`](api-core.md#adtoolbox.core.SeedDB) | Reaction, compound, protein, feed, and study databases. |
| [`configs.Metagenomics`](api-configs.md#adtoolbox.configs.Metagenomics) | [`core.Metagenomics`](api-core.md#adtoolbox.core.Metagenomics) | Pipeline directories and alignment thresholds. |
| [`configs.Annotation`](api-configs.md#adtoolbox.configs.Annotation) | [`core.Annotation`](api-core.md#adtoolbox.core.Annotation) | MetaCyc annotation. |
| [`configs.Utils`](api-configs.md#adtoolbox.configs.Utils) | [`utils`](api-utils.md) functions | Container images and Slurm defaults. |

## Overriding defaults

Pass any keyword argument to the config constructor to override a single default; the rest
still derive from the directory.

```python
metagenomics_config = configs.Metagenomics(
    "./my_metagenomics_run",
    database_dir="./my_database",
    protein_db="./custom/Protein_DB.fasta",  # explicit override
    bit_score=50,                            # stricter alignment filter
    e_value=1e-10,
    vsearch_similarity=0.99,
)
```

Container images are configured the same way:

```python
metagenomics_config = configs.Metagenomics(
    "./my_metagenomics_run",
    database_dir="./my_database",
    adtoolbox_docker="myorg/adtoolbox:dev",
    adtoolbox_singularity="docker://myorg/adtoolbox:dev",
)
```

!!! tip "Per-step control lives in the execution profile"
    For the batch pipeline, container backend, images, CPUs, memory, and Slurm settings are
    better set per step in a TOML execution profile than on the config object. See
    [Execution profiles](Metagenomics_Pipeline.md#execution-profiles).

A `configs.Metagenomics` built without an explicit `database_dir` falls back to a
`Database` config rooted at the run directory, which is convenient for self-contained
scratch runs but not what you want when several runs share one database.

## A minimal end-to-end script

```python
import numpy as np
from adtoolbox import adm, configs, core, utils

# Databases
database = core.Database(config=configs.Database(database_dir="./database"))

# Metagenomics: sample table -> microbial COD allocation.
# Works for either assay; switch the one keyword to change routes.
metagenomics = core.Metagenomics(
    configs.Metagenomics("./run", database_dir="./database")
)
result = metagenomics.batch_sample_to_cod(
    manifest="./samples.tsv",
    input_type="reads",
    assay="amplicon",           # 16S: denoise -> GTDB -> genomes -> groups
    # assay="shotgun",          # shotgun: reads -> protein DB (no GTDB/genomes)
    output_dir="./run/process",
    execute=True,
)

# Modeling: simulate e-ADM
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
)
solution = model.solve_model(np.linspace(0, 30, 300))
model.plot(solution).show()
```

## Where to go next

- [Quickstart](Quickstart.md) — the same path with more explanation.
- [Metagenomics pipeline](Metagenomics_Pipeline.md) — per-step Python API and execution
  profiles.
- [Parameter tuning](Optimization.md) — fitting the model to data.
- [Example notebooks](Notebooks.md) — executable versions of these workflows.
