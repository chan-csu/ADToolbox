#!/usr/bin/env python3
"""Rebuild the compact ADToolbox marker HMM and MMseqs FASTA assets."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPOSITORY_ROOT = HERE.parents[3]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from adtoolbox.markers import MarkerCatalog, build_marker_hmm_db


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kofam-profiles", type=Path, required=True, help="KOfam profiles.tar.gz")
    parser.add_argument("--dbcan-hmm", type=Path, required=True, help="Pinned dbCAN.hmm")
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()

    missing_tools = [tool for tool in ("hmmconvert", "hmmemit", "mmseqs") if shutil.which(tool) is None]
    if missing_tools:
        raise SystemExit(f"HMMER tools are required; missing: {', '.join(missing_tools)}")
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader((HERE / "profile_manifest.tsv").open(), delimiter="\t"))

    with tempfile.TemporaryDirectory(prefix="adtoolbox-markers-") as temp_name:
        temp = Path(temp_name)
        selected = temp / "selected"
        selected.mkdir()
        ko_ids = sorted({row["profile_id"] for row in rows if row["source"] == "KOfam"})
        wanted = {f"profiles/{ko}.hmm" for ko in ko_ids}
        with tarfile.open(args.kofam_profiles, "r:gz") as archive:
            members = [member for member in archive.getmembers() if member.name in wanted]
            found = {member.name for member in members}
            if missing := wanted - found:
                raise SystemExit(f"Missing KOfam profiles: {sorted(missing)}")
            archive.extractall(selected, members=members, filter="data")

        normalized = temp / "normalized"
        normalized.mkdir()
        normalized_dbcan = normalized / "dbCAN.hmm"
        with normalized_dbcan.open("w") as handle:
            subprocess.run(["hmmconvert", str(args.dbcan_hmm.resolve())], check=True, stdout=handle)
        for ko_id in ko_ids:
            source = selected / "profiles" / f"{ko_id}.hmm"
            destination = normalized / f"{ko_id}.hmm"
            with destination.open("w") as handle:
                subprocess.run(["hmmconvert", str(source)], check=True, stdout=handle)

        build_manifest = temp / "build_manifest.tsv"
        with build_manifest.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["marker_id", "source_hmm", "profile_id", "score_threshold", "score_type"],
                delimiter="\t",
            )
            writer.writeheader()
            for row in rows:
                source_hmm = normalized_dbcan if row["source"] == "dbCAN" else normalized / f"{row['profile_id']}.hmm"
                writer.writerow(
                    {
                        "marker_id": row["marker_id"],
                        "source_hmm": source_hmm,
                        "profile_id": row["profile_id"],
                        "score_threshold": row["score_threshold"],
                        "score_type": row["score_type"],
                    }
                )

        hmm_output = output_dir / "Marker_Profiles.hmm"
        cutoff_output = output_dir / "Marker_Profile_Cutoffs.csv"
        build_marker_hmm_db(
            build_manifest,
            hmm_output,
            MarkerCatalog.from_json(),
            cutoff_output=cutoff_output,
        )

        fasta_output = output_dir / "Marker_Protein_DB.fasta"
        emitted = temp / "profile_consensus.fasta"
        with emitted.open("w") as fasta:
            subprocess.run(["hmmemit", "-c", str(hmm_output)], check=True, stdout=fasta)
        with emitted.open() as source, fasta_output.open("w") as fasta:
            for line in source:
                if line.startswith(">"):
                    emitted_id = line[1:].strip().split(maxsplit=1)[0]
                    profile_id = emitted_id.removesuffix("-consensus")
                    marker_id = profile_id.split("__", 1)[0]
                    fasta.write(f">{profile_id}|{marker_id} HMM profile consensus\n")
                else:
                    fasta.write(line)

        mmseqs_prefix = output_dir / "Marker_Protein_DB_mmseqs"
        mmseqs_outputs = [
            mmseqs_prefix,
            mmseqs_prefix.with_name(mmseqs_prefix.name + ".dbtype"),
            mmseqs_prefix.with_name(mmseqs_prefix.name + ".index"),
            mmseqs_prefix.with_name(mmseqs_prefix.name + ".lookup"),
            mmseqs_prefix.with_name(mmseqs_prefix.name + ".source"),
            mmseqs_prefix.with_name(mmseqs_prefix.name + "_h"),
            mmseqs_prefix.with_name(mmseqs_prefix.name + "_h.dbtype"),
            mmseqs_prefix.with_name(mmseqs_prefix.name + "_h.index"),
        ]
        for path in mmseqs_outputs:
            path.unlink(missing_ok=True)
        subprocess.run(
            ["mmseqs", "createdb", fasta_output.name, mmseqs_prefix.name],
            check=True,
            cwd=output_dir,
        )

    checksummed = [
        "Marker_Profiles.hmm",
        "Marker_Profile_Cutoffs.csv",
        "Marker_Protein_DB.fasta",
        *[path.name for path in mmseqs_outputs],
        "profile_manifest.tsv",
        "SOURCES.json",
    ]
    with (output_dir / "SHA256SUMS").open("w") as handle:
        for name in checksummed:
            path = output_dir / name
            handle.write(f"{_sha256(path)}  {name}\n")


if __name__ == "__main__":
    main()
