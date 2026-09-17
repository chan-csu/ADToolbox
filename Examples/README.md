# Examples — notebooks and study data

Each notebook is one stage of a chain: **sequences → microbial seed → `Experiment` objects →
calibration**. They hand off through files, so you can run any stage without re-running the
one before it.

## Notebooks

| notebook | what it does | reads | writes |
|---|---|---|---|
| `metagenomics.ipynb` | 16S (SRA) → functional-group COD per sample; Bray–Curtis NMDS, PERMANOVA, per-group enrichment, publication figure | `Studies/16s_sra_day0.csv` | `../tutorial_output/ding_day0/<sample>/cod_profile.csv` |
| `build_experimental_db.ipynb` | VFAs (mg/L) → gCOD/L via the SEED DB; attaches the microbial seed and feed; assembles `Experiment` objects | `Studies/ding_raw/`, the COD profiles above | `Studies/ding_experiments.json` |
| `build_rico_validation_db.ipynb` | Same, for the held-out validation study (Rico et al.). Rico ships denoised rep-seqs, so it skips the amplicon step | `Studies/rico_raw/` | `Studies/rico_experiments.json` |
| `parameter_tuning.ipynb` | Calibrates e-ADM; §10 leave-one-condition-out; §11 microbiome-swap ablation | `Studies/ding_experiments.json`, `Studies/16s_methane_observed.csv` | `../reference_data/calibrated_model.json` |
| `modeling.ipynb` | General tutorial — building/solving a model, feed DB, study DBs. Not part of the calibration chain | `Toymap.json`, `toy_model.json`, `feed_db.tsv` | — |

Optimisation settings for every calibration live in **one cell at the top of
`parameter_tuning.ipynb`** (`MAXITER`, `POPSIZE`, `TOL`, `SOLVE_LIMIT`, `FIT_SECONDS`,
`RANDOM_STATE`).

## `Studies/`

**Calibration / validation sets** (regenerate with the build notebooks)

- `ding_experiments.json` — 6 experiments, 3 conditions (FW+TWAS, FW+AS, FW+TWAS+AS) × 2 replicates.
  Same feedstock and pH throughout; conditions differ only in inoculum.
- `rico_experiments.json` — 9 experiments, food waste, inoculum A, 35 °C, pH 5/7/9 × 3 replicates.
  Held out: never used for calibration.

**Raw study data**

- `ding_raw/` — VFA and methane tables transcribed from the paper.
- `rico_raw/` — the authors' published 16S rep-seqs, feature table, taxonomy, acid time-course and
  sample metadata. Includes feedstock-only samples, which is what let us measure how much of the
  day-0 community arrives with the substrate.
- `16s_sra_day0.csv` — SRA accessions, Ding day 0 (what `metagenomics.ipynb` consumes).
- `16s_sra_all.csv` — the same study, all timepoints. Not used yet; kept for a future time-course run.
- `16s_methane_observed.csv` — cumulative methane, used as a held-out check (never fitted).

**Reference databases** (shipped, consumed by `adtoolbox`, not study-specific)

- `experimental_data_references.json`, `metagenomics_studies.tsv`, `feed_db.tsv`
- `16s_codigestion_experiments.json` — an older co-digestion set, currently unreferenced.

## Annotation contract

Both study sets were regenerated on **2026-08-22** under one configuration, so they are directly
comparable:

- backend **mmseqs** (not HMMER)
- marker catalog **0.4.1**
- VSEARCH identity **0.95**

Each experiment records this in a `microbiome_source` field. Two corrections are baked in:

1. **Backend.** The HMMER path misses `cat`, which is `required_all` for both chain-elongation
   panels, so `X_chain_et`/`X_chain_lac` came out 0 for every Rico sample — leaving the model
   structurally unable to produce the caproate Rico measures. mmseqs detects it.
2. **Acetoclastic methanogenesis.** `mcrA` + `acs` + `cdh` are shared by *all* methanogens, so
   hydrogenotrophs were being scored as acetoclastic (and double-counted into both groups).
   Catalog 0.4.1 adds H₂-uptake hydrogenase markers and forbids them on the `acs` route, which
   keeps *Methanosarcina* and *Methanothrix* while excluding *Methanobacterium* and relatives.

Feeds come from each study's own characterisation. Rico's inert fraction (`si`/`xi` = 9.4) is
derived from their BMP (518 mL CH₄/gVS → 317 mL CH₄/gCOD against a 350 theoretical maximum) and
agrees with their measured lignin (8.8%). Protein is **not** measured in either study and is an
assumption (18% for food waste), flagged in the `feed.source` field.

## Outputs — `../tutorial_output/` (~13 MB, untracked)

Only pipeline *results* are kept here; intermediates are deleted after each run.

- `ding_day0/<sample>/` — `cod_profile.csv` (the microbial seed) plus `cod_potential.csv`,
  `cod_evidence_qc.csv`, `genome_cods.csv`, `genome_gene_annotations.csv`, pathway scores,
  feature/genome abundances, `representative_genomes.csv`, `provenance.json`
- `rico_mmseqs/profiles/J1…J9` — Rico food waste (plus `genome_cods.json`, the per-genome cache
  that lets profiles be rebuilt without re-annotating)
- `rico_manure/profiles/J18…J27` — Rico manure
- `figs/` — figures
- `ding_setup/local_profile.toml` — execution profile that runs every pipeline step locally
  instead of via Slurm

**Not kept** (regenerable, and large): downloaded reads, `scratch/` directories, mmseqs alignment
tables and `mmseqs_tmp/`, decompressed `.fna`, VSEARCH hit tables, and pipeline logs. Removing
these took the directory from 1.1 GB to 13 MB.

Genomes are cached outside the repo at `~/ADToolbox/genome_cache/`, so changing the backend or
catalog costs an annotation re-run (~20 min) rather than a fresh download (~4 h). The ADToolbox
databases (GTDB and friends) live at `~/ADToolbox/Database` — note some notebooks still default
`DATABASE_DIR` to `REPO/database`, which is wrong on this machine and needs overriding.
