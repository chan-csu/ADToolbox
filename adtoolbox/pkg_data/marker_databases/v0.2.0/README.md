# ADToolbox direct gene-to-COD marker databases v0.3.0

This directory contains the two ready-to-use backends for the v0.3.0 marker
catalog:

- `Marker_Profiles.hmm`: compact HMMER database made from selected KOfam and
  dbCAN profiles.
- `Marker_Profile_Cutoffs.csv`: KOfam adaptive thresholds, including whether
  the profile uses a full-sequence or best-domain score. dbCAN families use the
  pipeline's e-value and coverage filters because dbCAN does not publish KOfam-
  style per-family score thresholds.
- `Marker_Protein_DB.fasta`: one `hmmemit --consensus` sequence per selected
  HMM, for the faster MMseqs backend. These are profile consensuses—not isolate
  proteins—and should not be interpreted as biological reference accessions.
- `Marker_Protein_DB_mmseqs*`: the ready-to-search MMseqs database built from
  that FASTA, so samples reuse one shared target rather than rebuilding it.
- `profile_manifest.tsv`: auditable marker-to-profile mapping.
- `SOURCES.json`: pinned upstream releases and source URLs.
- `SHA256SUMS`: integrity checks for the shipped assets.

The v0.3.0 assets broaden peptide/protein utilization, lipid hydrolysis,
Entner-Doudoroff metabolism, and amino-acid catabolism. The stable directory
name is retained so existing installations continue to resolve bundled paths.

The HMM backend is the more sensitive and preferred option for assembled
genomes. The MMseqs backend is faster and is useful for high-throughput
screening. Both resolve hits to the same 96 canonical marker IDs and therefore
feed the same COD-group classifier.

## Rebuild

Run `python rebuild.py --help`. A rebuild requires Python with ADToolbox,
HMMER (`hmmconvert` and `hmmemit`), MMseqs2, and enough temporary storage for the KOfam archive. The
script accepts already-downloaded source files, which is preferable on a
cluster or for reproducible offline builds.

Before redistributing this folder outside this repository, review the upstream
KOfam/KEGG and dbCAN licensing terms recorded in `SOURCES.json`.
