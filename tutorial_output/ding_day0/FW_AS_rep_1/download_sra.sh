#!/bin/bash
set -euo pipefail
set -e
mkdir -p /Users/parsaghadermarzi/Desktop/Academics/Projects/ADToolbox/tutorial_output/ding_sra/SRR19159476
if command -v prefetch >/dev/null 2>&1 && command -v fasterq-dump >/dev/null 2>&1; then
  prefetch SRR19159476 -O /Users/parsaghadermarzi/Desktop/Academics/Projects/ADToolbox/tutorial_output/ding_sra --max-size 100000000
  fasterq-dump /Users/parsaghadermarzi/Desktop/Academics/Projects/ADToolbox/tutorial_output/ding_sra/SRR19159476/SRR19159476.sra -O /Users/parsaghadermarzi/Desktop/Academics/Projects/ADToolbox/tutorial_output/ding_sra/SRR19159476 --split-3 --temp /Users/parsaghadermarzi/Desktop/Academics/Projects/ADToolbox/tutorial_output/ding_sra/SRR19159476
  rm -f /Users/parsaghadermarzi/Desktop/Academics/Projects/ADToolbox/tutorial_output/ding_sra/SRR19159476/SRR19159476.sra
elif command -v curl >/dev/null 2>&1; then
  ena_url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR19159476&result=read_run&fields=fastq_ftp&format=tsv&download=false"
  fastq_ftp=$(curl -fsSL "$ena_url" | awk -F '\t' 'NR==2 {print $NF}')
  if [ -z "$fastq_ftp" ]; then
    echo "Could not find ENA FASTQ URLs for SRR19159476" >&2
    exit 1
  fi
  old_ifs=$IFS
  IFS=';'
  for ftp_path in $fastq_ftp; do
    IFS=$old_ifs
    file_name=$(basename "$ftp_path")
    case "$ftp_path" in
      ftp://*|https://*) download_url="$ftp_path" ;;
      *) download_url="ftp://$ftp_path" ;;
    esac
    https_url=$(printf '%s' "$download_url" | sed 's#^ftp://#https://#')
    echo "Downloading $file_name from ENA"
    curl -fsSL --retry 3 --connect-timeout 30 "$https_url" -o /Users/parsaghadermarzi/Desktop/Academics/Projects/ADToolbox/tutorial_output/ding_sra/SRR19159476/"$file_name" ||       curl -fsSL --retry 3 --connect-timeout 30 "$download_url" -o /Users/parsaghadermarzi/Desktop/Academics/Projects/ADToolbox/tutorial_output/ding_sra/SRR19159476/"$file_name"
    IFS=';'
  done
  IFS=$old_ifs
else
  echo "SRA download requires either prefetch/fasterq-dump from SRA Toolkit or curl for ENA FASTQ download." >&2
  exit 127
fi
