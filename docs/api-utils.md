# `utils`

Small helpers used across the toolbox: sequence file handling, MMseqs2 command
construction, SRA metadata lookups, JSON loading, and Slurm job-script generation.

```python
from adtoolbox import utils
```

Most of the MMseqs2 helpers build a shell command and either return it or run it, which is
what lets the same code path work locally, inside a container, or as a submitted Slurm
job.

---

## Sequence files

### fasta_to_dict

::: adtoolbox.utils.fasta_to_dict

### dict_to_fasta

::: adtoolbox.utils.dict_to_fasta

### extract_zipped_file

::: adtoolbox.utils.extract_zipped_file

### make_json_from_genomes

::: adtoolbox.utils.make_json_from_genomes

### Sequence_Toolkit

::: adtoolbox.utils.Sequence_Toolkit

---

## MMseqs2

### create_mmseqs_database

::: adtoolbox.utils.create_mmseqs_database

### index_mmseqs_db

::: adtoolbox.utils.index_mmseqs_db

### mmseqs_search

::: adtoolbox.utils.mmseqs_search

### mmseqs_result_db_to_tsv

::: adtoolbox.utils.mmseqs_result_db_to_tsv

---

## JSON loading

### load_multiple_json_files

::: adtoolbox.utils.load_multiple_json_files

### load_json_entry

::: adtoolbox.utils.load_json_entry

### load_model_json

::: adtoolbox.utils.load_model_json

---

## Cluster execution

### wrap_for_slurm

::: adtoolbox.utils.wrap_for_slurm

### generate_batch_script

::: adtoolbox.utils.generate_batch_script

---

## Metadata

### get_sample_metadata_from_accession

::: adtoolbox.utils.get_sample_metadata_from_accession
