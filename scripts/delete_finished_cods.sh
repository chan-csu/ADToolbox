#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  delete_finished_cods.sh RESULTS_DIR [--apply]

Find completed per-sample COD outputs under RESULTS_DIR.

By default, this is a dry run and changes nothing. Pass --apply to delete:
  <sample>/cod_profile.csv
  <sample>/genome_cods.csv
  <sample>/cod_potential.csv
  <sample>/genome_pathway_scores.csv
  <sample>/genome_gene_annotations.csv
  <sample>/cod_evidence_qc.csv

Genome downloads, feature tables, representative sequences, GTDB matches,
genome abundances, and genome-alignment/annotation files are preserved.
EOF
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
    usage >&2
    exit 2
fi

results_dir=${1%/}
apply=false

if [[ $# -eq 2 ]]; then
    if [[ $2 != "--apply" ]]; then
        usage >&2
        exit 2
    fi
    apply=true
fi

if [[ ! -d $results_dir ]]; then
    echo "Results directory does not exist: $results_dir" >&2
    exit 1
fi

profiles_found=0
files_selected=0

while IFS= read -r -d '' cod_profile; do
    # ADToolbox considers a COD CSV usable when it contains at least one data row.
    if ! awk 'NR > 1 && $0 !~ /^[[:space:]]*$/ { found=1; exit } END { exit !found }' "$cod_profile"; then
        continue
    fi

    sample_dir=${cod_profile%/cod_profile.csv}
    profiles_found=$((profiles_found + 1))

    for output in \
        "$cod_profile" \
        "$sample_dir/genome_cods.csv" \
        "$sample_dir/cod_potential.csv" \
        "$sample_dir/genome_pathway_scores.csv" \
        "$sample_dir/genome_gene_annotations.csv" \
        "$sample_dir/cod_evidence_qc.csv"; do
        if [[ ! -f $output ]]; then
            continue
        fi
        files_selected=$((files_selected + 1))
        if [[ $apply == true ]]; then
            rm -- "$output"
            echo "Deleted: $output"
        else
            echo "Would delete: $output"
        fi
    done
done < <(find "$results_dir" -mindepth 2 -maxdepth 2 -type f -name cod_profile.csv -print0)

if [[ $apply == true ]]; then
    echo "Deleted $files_selected COD output file(s) from $profiles_found completed sample(s)."
else
    echo "Dry run: selected $files_selected file(s) from $profiles_found completed sample(s)."
    echo "Run again with --apply to delete them."
fi
