import json
import pathlib
import subprocess
import threading

import polars as pl
import pytest

from adtoolbox import configs, core


def _write_alignment(path, rows):
    header = "\t".join(
        [
            "query",
            "target",
            "fident",
            "alnlen",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "tstart",
            "tend",
            "evalue",
            "bits",
        ]
    )
    with open(path, "w") as f:
        f.write(header + "\n")
        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")


def _read_profile(path):
    frame = pl.read_csv(path)
    return dict(zip(frame["group"].to_list(), frame["value"].to_list()))


def test_sample_to_cod_from_shotgun_alignment(tmp_path):
    reaction_db = tmp_path / "reactions.csv"
    pl.DataFrame(
        [
            {"EC_Numbers": "1.1.1.1", "e_adm_Reactions": "Uptake of sugars"},
            {"EC_Numbers": "2.2.2.2", "e_adm_Reactions": "Methanogenessis from acetate and h2"},
        ]
    ).write_csv(reaction_db)

    alignment = tmp_path / "reads.tsv"
    _write_alignment(
        alignment,
        [
            ["read_1", "P00001|1.1.1.1", 1, 100, 0, 0, 1, 100, 1, 100, 1e-20, 100],
            ["read_2", "P00002|2.2.2.2", 1, 100, 0, 0, 1, 100, 1, 100, 1e-20, 100],
            ["read_3", "P00003|9.9.9.9", 1, 100, 0, 0, 1, 100, 1, 100, 1e-20, 100],
        ],
    )

    metagenomics = core.Metagenomics(configs.Metagenomics(csv_reaction_db=reaction_db))
    result = metagenomics.sample_to_cod(
        sample_name="sample_a",
        output_dir=tmp_path / "out",
        mode="shotgun-alignment",
        alignment_file=alignment,
        verbose=False,
    )

    cod_profile = _read_profile(result["artifacts"]["cod_profile"])

    assert cod_profile["X_su"] == 0.5
    assert cod_profile["X_Me_ac"] == 0.5
    assert sum(cod_profile.values()) == 1
    assert result["artifacts"]["cod_profile"].endswith("cod_profile.csv")
    assert result["artifacts"]["ec_counts"].endswith("ec_counts.csv")


def test_sample_to_cod_weights_genome_alignments(tmp_path):
    reaction_db = tmp_path / "reactions.csv"
    pl.DataFrame(
        [
            {"EC_Numbers": "1.1.1.1", "e_adm_Reactions": "Uptake of sugars"},
            {"EC_Numbers": "2.2.2.2", "e_adm_Reactions": "Uptake of amino acids"},
        ]
    ).write_csv(reaction_db)

    alignments = tmp_path / "alignments"
    alignments.mkdir()
    _write_alignment(
        alignments / "Alignment_Results_mmseq_genome_a.tsv",
        [["gene_1", "P00001|1.1.1.1", 1, 100, 0, 0, 1, 100, 1, 100, 1e-20, 100]],
    )
    _write_alignment(
        alignments / "Alignment_Results_mmseq_genome_b.tsv",
        [["gene_2", "P00002|2.2.2.2", 1, 100, 0, 0, 1, 100, 1, 100, 1e-20, 100]],
    )
    abundances = tmp_path / "abundances.json"
    abundances.write_text(json.dumps({"genome_a": 0.75, "genome_b": 0.25}))

    metagenomics = core.Metagenomics(configs.Metagenomics(csv_reaction_db=reaction_db))
    result = metagenomics.sample_to_cod(
        sample_name="sample_b",
        output_dir=tmp_path / "out",
        mode="genome-alignments",
        genome_abundances=abundances,
        genome_alignments=alignments,
        verbose=False,
    )

    cod_profile = _read_profile(result["artifacts"]["cod_profile"])

    assert cod_profile["X_su"] == 0.75
    assert cod_profile["X_aa"] == 0.25
    assert result["artifacts"]["genome_cods"].endswith("genome_cods.csv")


def test_extract_relative_abundances_uses_feature_ids(tmp_path):
    feature_table = tmp_path / "feature-table.tsv"
    feature_table.write_text(
        "# Constructed from biom file\n"
        "#OTU ID\tsample_a\tsample_b\n"
        "asv_1\t10\t5\n"
        "asv_2\t30\t5\n"
    )

    metagenomics = core.Metagenomics(configs.Metagenomics())
    abundances = metagenomics.extract_relative_abundances(feature_table, sample_names=["sample_a"])

    assert abundances == {"sample_a": {"asv_2": 0.75, "asv_1": 0.25}}


def test_sample_repseq_selection_accepts_vsearch_size_headers(tmp_path):
    rep_seqs = tmp_path / "rep-seqs.fasta"
    rep_seqs.write_text(
        ">asv_1;size=20\n"
        "ACGT\n"
        ">asv_2;size=10\n"
        "TGCA\n"
    )

    output = tmp_path / "sample_repseqs.fasta"
    metagenomics = core.Metagenomics(configs.Metagenomics())
    metagenomics._write_sample_repseqs(rep_seqs, {"asv_1": 1.0}, output)

    assert output.read_text() == ">asv_1\nACGT\n"


def test_empty_gtdb_alignment_has_no_representative_genomes(tmp_path):
    matches = tmp_path / "matches.blast"
    matches.write_text("")

    metagenomics = core.Metagenomics(configs.Metagenomics())

    assert metagenomics.get_genomes_from_gtdb_alignment(matches) == {}


def test_sample_repseq_selection_normalizes_feature_table_size_suffixes(tmp_path):
    rep_seqs = tmp_path / "rep-seqs.fasta"
    rep_seqs.write_text(
        ">asv_1;size=20\n"
        "ACGT\n"
        ">asv_2;size=10\n"
        "TGCA\n"
    )

    output = tmp_path / "sample_repseqs.fasta"
    metagenomics = core.Metagenomics(configs.Metagenomics())
    metagenomics._write_sample_repseqs(rep_seqs, {"asv_1;size=20": 1.0}, output)

    assert output.read_text() == ">asv_1\nACGT\n"


def test_alignment_files_from_path_indexes_gtdb_prefix_before_contig(tmp_path):
    alignment = tmp_path / "Alignment_Results_mmseq_GCF000001.1~NZABC01000001.1.tsv"
    alignment.write_text("")

    metagenomics = core.Metagenomics(configs.Metagenomics())
    alignments = metagenomics._alignment_files_from_path(tmp_path)

    assert alignments["GCF000001.1"] == str(alignment)
    assert alignments["GCF000001.1~NZABC01000001.1"] == str(alignment)


def test_genome_files_from_dir_indexes_assembly_accession_prefix(tmp_path):
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    genome = genome_dir / "GCF_000146505.1_ASM14650v1_genomic.fna.gz"
    genome.write_text("")

    metagenomics = core.Metagenomics(configs.Metagenomics())
    genome_files = metagenomics._genome_files_from_dir(genome_dir)

    assert genome_files["GCF_000146505.1"] == str(genome.resolve())


def test_extract_genome_info_df_is_tall_and_recursive(tmp_path):
    genome_dir = tmp_path / "genomes" / "GCF_000146505.1_ASM14650v1"
    genome_dir.mkdir(parents=True)
    genome = genome_dir / "GCF_000146505.1_ASM14650v1_genomic.fna.gz"
    genome.write_text("")
    excluded = genome_dir / "GCF_000146505.1_ASM14650v1_cds_from_genomic.fna.gz"
    excluded.write_text("")

    metagenomics = core.Metagenomics(configs.Metagenomics())
    info = metagenomics.extract_genome_info_df(tmp_path / "genomes")

    assert info.columns == ["genome_id", "assembly_accession", "assembly_name", "path"]
    assert info.height == 1
    row = info.to_dicts()[0]
    assert row["assembly_accession"] == "GCF_000146505.1"
    assert row["path"] == str(genome.resolve())


def test_download_genome_uses_ncbi_https_layout(tmp_path):
    metagenomics = core.Metagenomics(configs.Metagenomics())
    script = metagenomics.download_genome("GCF_000146505.1", tmp_path / "genomes")[0]

    assert "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/505" in script
    assert "curl -fL" in script
    assert "rsync" not in script


def test_download_genomes_uses_ncbi_datasets_batch_and_bounded_workers(tmp_path):
    metagenomics = core.Metagenomics(configs.Metagenomics())
    script = metagenomics.download_genomes(
        ["GCF_000146505.1", "GCA_937889405.1", "GCF_000146505.1"],
        tmp_path / "genomes",
        tmp_path / "work",
        max_workers=6,
    )[0]

    assert "datasets download genome accession --inputfile" in script
    assert "--dehydrated" in script
    assert "datasets rehydrate" in script
    assert "--max-workers 6" in script
    assert "datasets --version" in script
    assert "falling back to NCBI HTTPS downloads" in script
    assert "curl -fL --retry 3" in script
    assert "Skipping unavailable or deprecated genome" in script
    assert 'failed_file="$batch_dir/failed_genomes.txt"' in script
    assert (tmp_path / "work" / "accessions.txt").read_text().splitlines() == [
        "GCA_937889405.1",
        "GCF_000146505.1",
    ]
    assert "GCF_000146505.1" not in script
    assert subprocess.run(["bash", "-n"], input=script, text=True).returncode == 0


def test_multiline_container_commands_use_bash_for_pipefail(tmp_path):
    metagenomics = core.Metagenomics(configs.Metagenomics())
    command = "set -euo pipefail\necho ready"

    apptainer = metagenomics._wrap_external_command(
        command,
        container="apptainer",
        mounts=[tmp_path],
        image="image.sif",
    )
    docker = metagenomics._wrap_external_command(
        command,
        container="docker",
        mounts=[tmp_path],
        image="image:latest",
    )

    assert "image.sif bash -lc" in apptainer
    assert "image:latest bash -lc" in docker
    assert " sh -lc" not in apptainer
    assert " sh -lc" not in docker


def test_amplicon_pipeline_downloads_missing_genome_before_alignment(monkeypatch, tmp_path):
    feature_table = tmp_path / "feature-table.tsv"
    feature_table.write_text("#OTU ID\tsample_a\nasv_1\t10\nasv_2\t5\n")
    rep_seqs = tmp_path / "rep-seqs.fasta"
    rep_seqs.write_text(">asv_1\nACGT\n>asv_2\nTGCA\n")
    matches = tmp_path / "matches.blast"
    matches.write_text(
        "asv_1\tGB_GCF_000146505.1~contig\t100\t4\t0\t0\t1\t4\t1\t4\t0\t10\n"
        "asv_2\tGB_GCA_000000001.1~contig\t100\t4\t0\t0\t1\t4\t1\t4\t0\t10\n"
    )
    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir()
    calls = []

    metagenomics = core.Metagenomics(configs.Metagenomics())

    def fake_execute(script, *, step_name, **kwargs):
        calls.append(step_name)
        if step_name == "download_genomes":
            assembly_dir = genomes_dir / "GCF_000146505.1_ASM14650v1"
            assembly_dir.mkdir()
            (assembly_dir / "GCF_000146505.1_ASM14650v1_genomic.fna.gz").write_text("")
        elif step_name == "align_genomes":
            alignment = (
                tmp_path / "out" / "sample_a" / "scratch" / "genome_alignments"
                / "Alignment_Results_mmseq_GCF_000146505.1.tsv"
            )
            alignment.write_text("query\ttarget\n")
        return {"step": step_name, "status": "completed", "backend": "slurm"}

    monkeypatch.setattr(metagenomics, "_execute_step", fake_execute)
    monkeypatch.setattr(metagenomics, "cod_from_alignment", lambda *args, **kwargs: {"group_a": 1.0})
    monkeypatch.setattr(metagenomics, "aggregate_genome_cod", lambda *args, **kwargs: {"group_a": 1.0})

    result = metagenomics.sample_to_cod(
        sample_name="sample_a",
        output_dir=tmp_path / "out",
        mode="amplicon",
        genomes_dir=genomes_dir,
        feature_table=feature_table,
        rep_seqs=rep_seqs,
        gtdb_matches=matches,
        execute=True,
        verbose=False,
        execution_profile={
            "backend": "local",
            "container": "None",
            "slurm": {},
            "steps": {
                "download_genomes": {"backend": "slurm", "settings": {"max_workers": 4}},
                "align_genomes": {"backend": "slurm"},
            },
        },
    )

    assert calls == ["download_genomes", "align_genomes"]
    assert result["cod_profile"] == {"group_a": 1.0}
    assert "missing_genome_fastas" not in result["artifacts"]
    assert result["artifacts"]["skipped_genome_fastas"] == ["GCA_000000001.1"]


def test_amplicon_pipeline_completes_when_all_requested_genomes_are_unavailable(monkeypatch, tmp_path):
    feature_table = tmp_path / "feature-table.tsv"
    feature_table.write_text("#OTU ID\tsample_a\nasv_1\t10\n")
    rep_seqs = tmp_path / "rep-seqs.fasta"
    rep_seqs.write_text(">asv_1\nACGT\n")
    matches = tmp_path / "matches.blast"
    matches.write_text("asv_1\tGB_GCA_000000001.1~contig\t100\t4\t0\t0\t1\t4\t1\t4\t0\t10\n")
    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir()

    metagenomics = core.Metagenomics(configs.Metagenomics())
    monkeypatch.setattr(
        metagenomics,
        "_execute_step",
        lambda script, *, step_name, **kwargs: {
            "step": step_name,
            "status": "completed",
            "backend": "slurm",
        },
    )

    result = metagenomics.sample_to_cod(
        sample_name="sample_a",
        output_dir=tmp_path / "out",
        mode="amplicon",
        genomes_dir=genomes_dir,
        feature_table=feature_table,
        rep_seqs=rep_seqs,
        gtdb_matches=matches,
        execute=True,
        verbose=False,
        execution_profile={"backend": "local", "container": "None", "slurm": {}, "steps": {}},
    )

    assert result["status"] == "completed_no_available_genomes"
    assert result["cod_profile"] == {}
    assert result["artifacts"]["skipped_genome_fastas"] == ["GCA_000000001.1"]


def test_amplicon_pipeline_submits_one_alignment_task_per_sample(monkeypatch, tmp_path):
    feature_table = tmp_path / "feature-table.tsv"
    feature_table.write_text("#OTU ID\tsample_a\nasv_1\t10\nasv_2\t5\n")
    rep_seqs = tmp_path / "rep-seqs.fasta"
    rep_seqs.write_text(">asv_1\nACGT\n>asv_2\nTGCA\n")
    matches = tmp_path / "matches.blast"
    matches.write_text(
        "asv_1\tGB_GCF_000146505.1~contig\t100\t4\t0\t0\t1\t4\t1\t4\t0\t10\n"
        "asv_2\tGB_GCA_937889405.1~contig\t100\t4\t0\t0\t1\t4\t1\t4\t0\t10\n"
    )
    genomes_dir = tmp_path / "genomes"
    for accession in ("GCF_000146505.1", "GCA_937889405.1"):
        genome_dir = genomes_dir / accession
        genome_dir.mkdir(parents=True)
        (genome_dir / f"{accession}_genomic.fna.gz").write_text("")
    calls = []

    metagenomics = core.Metagenomics(configs.Metagenomics())

    def fake_execute(script, *, step_name, **kwargs):
        calls.append((step_name, script))
        return {"step": step_name, "status": "prepared", "backend": "slurm"}

    monkeypatch.setattr(metagenomics, "_execute_step", fake_execute)

    result = metagenomics.sample_to_cod(
        sample_name="sample_a",
        output_dir=tmp_path / "out",
        mode="amplicon",
        genomes_dir=genomes_dir,
        feature_table=feature_table,
        rep_seqs=rep_seqs,
        gtdb_matches=matches,
        execute=False,
        verbose=False,
        execution_profile={
            "backend": "slurm",
            "container": "None",
            "slurm": {},
            "steps": {"align_genomes": {"backend": "slurm"}},
        },
    )

    assert [name for name, _ in calls] == ["align_genomes"]
    alignment_script = calls[0][1]
    assert alignment_script.count("mmseqs easy-search") == 2
    assert "tmpfiles" not in alignment_script
    assert "mmseqs_tmp/GCF_000146505.1" in alignment_script
    assert "mmseqs_tmp/GCA_937889405.1" in alignment_script
    assert result["status"] == "waiting_for_genome_alignment"


def test_align_to_gtdb_reports_missing_database(tmp_path):
    config = configs.Metagenomics(amplicon2genome_db=tmp_path / "missing_gtdb")
    metagenomics = core.Metagenomics(config)

    try:
        metagenomics.align_to_gtdb(tmp_path / "rep-seqs.fasta", tmp_path / "out")
    except FileNotFoundError as exc:
        assert "No GTDB/amplicon-to-genome FASTA was found" in str(exc)
    else:
        raise AssertionError("Expected missing GTDB database error")


def test_sra_download_script_uses_split_3(tmp_path):
    metagenomics = core.Metagenomics(configs.Metagenomics())
    script, reads = metagenomics._sra_download_script(
        accession="SRR14342342",
        target_dir=tmp_path / "sra",
    )

    assert "--split-3" in script
    assert "www.ebi.ac.uk/ena/portal/api/filereport" in script
    assert "print $NF" in script
    assert "print $1" not in script
    assert "command -v prefetch" in script
    assert "command -v curl" in script
    assert "curl -fsSL --retry 3" in script
    assert reads["read_1"].endswith("SRR14342342_1.fastq")
    assert reads["read_2"].endswith("SRR14342342_2.fastq")


def test_resolved_sra_reads_accepts_single_end_output(tmp_path):
    accession_dir = tmp_path / "sra" / "SRR14342342"
    accession_dir.mkdir(parents=True)
    single = accession_dir / "SRR14342342.fastq"
    single.write_text("@read_1\nACGT\n+\n!!!!\n")

    reads = core.Metagenomics._resolved_sra_reads("SRR14342342", tmp_path / "sra", paired=True)

    assert reads == {"read_1": str(single), "read_2": None}


def test_resolved_sra_reads_accepts_gzipped_paired_output(tmp_path):
    accession_dir = tmp_path / "sra" / "SRR14342342"
    accession_dir.mkdir(parents=True)
    read_1 = accession_dir / "SRR14342342_1.fastq.gz"
    read_2 = accession_dir / "SRR14342342_2.fastq.gz"
    read_1.write_text("")
    read_2.write_text("")

    reads = core.Metagenomics._resolved_sra_reads("SRR14342342", tmp_path / "sra", paired=True)

    assert reads == {"read_1": str(read_1), "read_2": str(read_2)}


def test_shotgun_reads_profile_writes_slurm_script(tmp_path):
    reads = tmp_path / "reads.fastq"
    reads.write_text("@read_1\nACGT\n+\n!!!!\n")
    protein_db_mmseqs = tmp_path / "protein_db_mmseqs"
    protein_db_mmseqs.write_text("")

    profile = tmp_path / "profile.toml"
    profile.write_text(
        """
backend = "local"
container = "None"

[steps.align_short_reads]
backend = "slurm"
container = "None"
cpus = 16
memory = "64G"
time = "06:00:00"
"""
    )

    metagenomics = core.Metagenomics(configs.Metagenomics(protein_db_mmseqs=protein_db_mmseqs))
    result = metagenomics.sample_to_cod(
        sample_name="sample_c",
        output_dir=tmp_path / "out",
        mode="shotgun-reads",
        reads=reads,
        execution_profile=profile,
        execute=False,
        verbose=False,
    )

    step = result["artifacts"]["align_short_reads"]
    sbatch = tmp_path / "out" / "sample_c" / "scratch" / "slurm" / "align_short_reads.sbatch"

    assert step["backend"] == "slurm"
    assert step["sbatch"] == str(sbatch)
    assert "#SBATCH --cpus-per-task=16" in sbatch.read_text()
    assert "#SBATCH --mem=64G" in sbatch.read_text()


def test_task_manager_writes_events_for_local_steps(tmp_path):
    metagenomics = core.Metagenomics(configs.Metagenomics())
    logger = metagenomics._sample_pipeline_logger("sample_task", tmp_path / "sample_task", verbose=False)
    artifact = metagenomics._execute_step(
        "echo hello",
        step_name="hello",
        sample_name="sample_task",
        output_dir=tmp_path / "sample_task" / "scratch",
        logger=logger,
        execute=False,
        execution_profile={"backend": "local", "container": "None", "slurm": {}, "steps": {}},
    )

    events = tmp_path / "sample_task" / "scratch" / "task_events.jsonl"

    assert artifact["status"] == "prepared"
    assert artifact["task_events"] == str(events)
    assert '"event": "prepared"' in events.read_text()


def test_old_unknown_slurm_submission_is_treated_as_resumable_stale_state(tmp_path):
    workflow = core.MetagenomicsWorkflowState(tmp_path / "out")
    workflow.record(
        "sample_old",
        "preprocess",
        "submitted",
        artifact={"job_id": "12345", "step": "preprocess"},
    )

    class MissingJobChecker:
        @staticmethod
        def slurm_job_state(job_id):
            return None

    assert not workflow.active_submission(
        "sample_old",
        "preprocess",
        slurm_checker=MissingJobChecker(),
    )
    assert workflow.stage("sample_old", "preprocess")["status"] == "stale"


def test_task_manager_keeps_pipeline_dependencies_out_of_sbatch(tmp_path):
    profile = {
        "backend": "local",
        "container": "None",
        "slurm": {},
        "steps": {"child": {"backend": "slurm", "cpus": 2, "memory": "4G", "time": "00:30:00"}},
    }
    metagenomics = core.Metagenomics(configs.Metagenomics())
    logger = metagenomics._sample_pipeline_logger("sample_task", tmp_path / "sample_task", verbose=False)
    artifact = metagenomics._execute_step(
        "echo child",
        step_name="child",
        sample_name="sample_task",
        output_dir=tmp_path / "sample_task" / "scratch",
        logger=logger,
        execute=False,
        execution_profile=profile,
        dependencies=[{"step": "parent", "job_id": "12345"}],
    )

    sbatch = tmp_path / "sample_task" / "scratch" / "slurm" / "child.sbatch"

    assert artifact["dependencies"] == ["parent"]
    assert "#SBATCH --dependency" not in sbatch.read_text()


def test_cancelled_slurm_attempt_is_retried_as_a_new_job(monkeypatch, tmp_path):
    profile = {
        "backend": "local",
        "container": "None",
        "slurm": {"retry_delay_seconds": 1},
        "steps": {"child": {"backend": "slurm", "retries": 1}},
    }
    submissions = []

    def fake_run(command, capture_output=True, text=True, **kwargs):
        if command[0] == "sbatch":
            submissions.append(command)
            if len(submissions) == 1:
                return subprocess.CompletedProcess(command, 1, stdout="100\n", stderr="CANCELLED")
            return subprocess.CompletedProcess(command, 0, stdout="101\n", stderr="")
        raise AssertionError(f"Unexpected command: {command}")

    monkeypatch.setattr(core.subprocess, "run", fake_run)
    monkeypatch.setattr(core.time, "sleep", lambda seconds: None)

    metagenomics = core.Metagenomics(configs.Metagenomics())
    logger = metagenomics._sample_pipeline_logger("sample_task", tmp_path / "sample_task", verbose=False)
    artifact = metagenomics._execute_step(
        "echo child",
        step_name="child",
        sample_name="sample_task",
        output_dir=tmp_path / "sample_task" / "scratch",
        logger=logger,
        execute=True,
        execution_profile=profile,
    )

    events = (tmp_path / "sample_task" / "scratch" / "task_events.jsonl").read_text()
    sbatch = (tmp_path / "sample_task" / "scratch" / "slurm" / "child.sbatch").read_text()

    assert len(submissions) == 2
    assert all(command[0:3] == ["sbatch", "--wait", "--parsable"] for command in submissions)
    assert artifact["status"] == "completed"
    assert artifact["job_id"] == "101"
    assert artifact["job_ids"] == ["100", "101"]
    assert artifact["attempt_count"] == 2
    assert artifact["retries"] == 1
    assert artifact["retry_mode"] == "slurm_resubmit"
    assert "max_retries" not in sbatch
    assert "while true; do" not in sbatch
    assert '"event": "attempt_failed"' in events
    assert '"event": "completed"' in events


def test_non_retryable_slurm_status_is_not_resubmitted(monkeypatch, tmp_path):
    profile = {
        "backend": "local",
        "container": "None",
        "slurm": {"retry_delay_seconds": 1},
        "steps": {"child": {"backend": "slurm", "retries": 10}},
    }
    submissions = []

    def fake_run(command, capture_output=True, text=True, **kwargs):
        if command[0] == "sbatch":
            submissions.append(command)
            return subprocess.CompletedProcess(command, 64, stdout="200\n", stderr="no ASVs")
        raise AssertionError(f"Unexpected command: {command}")

    monkeypatch.setattr(core.subprocess, "run", fake_run)
    monkeypatch.setattr(core.time, "sleep", lambda seconds: None)
    metagenomics = core.Metagenomics(configs.Metagenomics())
    logger = metagenomics._sample_pipeline_logger("sample_task", tmp_path / "sample_task", verbose=False)

    with pytest.raises(RuntimeError, match="non-retryable status 64"):
        metagenomics._execute_step(
            "exit 64",
            step_name="child",
            sample_name="sample_task",
            output_dir=tmp_path / "sample_task" / "scratch",
            logger=logger,
            execute=True,
            execution_profile=profile,
        )

    assert len(submissions) == 1


def test_slurm_failure_raises(monkeypatch, tmp_path):
    profile = {
        "backend": "local",
        "container": "None",
        "slurm": {"retry_delay_seconds": 1},
        "steps": {"child": {"backend": "slurm"}},
    }

    def fake_run(command, capture_output=True, text=True, **kwargs):
        if command[0] == "sbatch":
            return subprocess.CompletedProcess(command, 1, stdout="100\n", stderr="job failed")
        raise AssertionError(f"Unexpected command: {command}")

    monkeypatch.setattr(core.subprocess, "run", fake_run)
    monkeypatch.setattr(core.time, "sleep", lambda seconds: None)

    metagenomics = core.Metagenomics(configs.Metagenomics())
    logger = metagenomics._sample_pipeline_logger("sample_task", tmp_path / "sample_task", verbose=False)

    try:
        metagenomics._execute_step(
            "echo child",
            step_name="child",
            sample_name="sample_task",
            output_dir=tmp_path / "sample_task" / "scratch",
            logger=logger,
            execute=True,
            execution_profile=profile,
        )
    except RuntimeError as exc:
        assert "Slurm step child failed" in str(exc)
    else:
        raise AssertionError("Expected blocking Slurm failure to raise")


def test_slurm_steps_use_blocking_sbatch(monkeypatch, tmp_path):
    profile = {
        "backend": "local",
        "container": "None",
        "slurm": {"poll_seconds": 1},
        "steps": {"child": {"backend": "slurm"}},
    }

    def fake_run(command, capture_output=True, text=True, **kwargs):
        if command[0] == "sbatch":
            return subprocess.CompletedProcess(command, 0, stdout="101\n", stderr="")
        raise AssertionError(f"Unexpected command: {command}")

    monkeypatch.setattr(core.subprocess, "run", fake_run)
    monkeypatch.setattr(core.time, "sleep", lambda seconds: None)

    metagenomics = core.Metagenomics(configs.Metagenomics())
    logger = metagenomics._sample_pipeline_logger("sample_task", tmp_path / "sample_task", verbose=False)
    artifact = metagenomics._execute_step(
        "echo child",
        step_name="child",
        sample_name="sample_task",
        output_dir=tmp_path / "sample_task" / "scratch",
        logger=logger,
        execute=True,
        execution_profile=profile,
    )

    assert artifact["status"] == "completed"
    assert artifact["job_id"] == "101"


def test_amplicon_preprocessing_writes_trim_and_vsearch_steps(tmp_path):
    read_1 = tmp_path / "sample_R1.fastq"
    read_2 = tmp_path / "sample_R2.fastq"
    read_1.write_text("@read_1\nACGTACGT\n+\n!!!!!!!!\n")
    read_2.write_text("@read_1\nACGTACGT\n+\n!!!!!!!!\n")

    profile = tmp_path / "profile.toml"
    profile.write_text(
        """
backend = "local"
container = "None"

[steps.trim_reads]
backend = "slurm"
cpus = 4
memory = "8G"
time = "01:00:00"

[steps.trim_reads.settings]
threads = 4
minimum_length = 80

[steps.build_amplicon_features]
backend = "local"
cpus = 6

[steps.build_amplicon_features.settings]
threads = 6
identity = 0.99
maxee = 0.5
"""
    )

    metagenomics = core.Metagenomics(configs.Metagenomics())
    result = metagenomics.preprocess_amplicon_sample(
        sample_name="sample_d",
        output_dir=tmp_path / "out",
        read_1=read_1,
        read_2=read_2,
        execution_profile=profile,
        execute=False,
        verbose=False,
    )

    sample_dir = tmp_path / "out" / "sample_d"
    scratch_dir = sample_dir / "scratch"
    trim_script = scratch_dir / "trim_reads.sh"
    feature_script = scratch_dir / "build_amplicon_features.sh"
    trim_sbatch = scratch_dir / "slurm" / "trim_reads.sbatch"

    assert result["artifacts"]["feature_table"].endswith("feature-table.tsv")
    assert "/scratch/amplicon_preprocess/" in result["artifacts"]["feature_table"]
    assert "fastp" in trim_script.read_text()
    assert "--length_required 80" in trim_script.read_text()
    assert "--detect_adapter_for_pe" in trim_script.read_text()
    assert "vsearch --cluster_unoise" in feature_script.read_text()
    assert "vsearch --uchime3_denovo" in feature_script.read_text()
    assert "--id 0.99" in feature_script.read_text()
    assert "#SBATCH --cpus-per-task=4" in trim_sbatch.read_text()


def test_step_level_amplicon_api_returns_artifacts(tmp_path):
    read_1 = tmp_path / "sample_R1.fastq"
    read_1.write_text("@read_1\nACGTACGT\n+\n!!!!!!!!\n")

    metagenomics = core.Metagenomics(configs.Metagenomics())
    trim_result = metagenomics.run_trim_reads_step(
        sample_name="sample_e",
        output_dir=tmp_path / "out",
        read_1=read_1,
        execute=False,
        verbose=False,
    )
    feature_result = metagenomics.run_build_amplicon_features_step(
        sample_name="sample_e",
        output_dir=tmp_path / "out",
        read_1=trim_result["artifacts"]["trimmed_reads"]["read_1"],
        execute=False,
        verbose=False,
    )

    sample_dir = tmp_path / "out" / "sample_e"
    assert trim_result["step"] == "trim_reads"
    assert feature_result["step"] == "build_amplicon_features"
    assert trim_result["artifacts"]["trimmed_reads"]["read_1"].endswith("_trimmed_R1.fastq.gz")
    assert feature_result["artifacts"]["feature_table"].endswith("feature-table.tsv")
    assert (sample_dir / "trim_reads.sh").exists()
    assert (sample_dir / "build_amplicon_features.sh").exists()


def test_amplicon_feature_step_uses_vsearch_unoise(tmp_path):
    read_1 = tmp_path / "sample_R1.fastq"
    read_1.write_text("@read_1\nACGTACGT\n+\n!!!!!!!!\n")

    metagenomics = core.Metagenomics(configs.Metagenomics())
    metagenomics.run_build_amplicon_features_step(
        sample_name="sample_f",
        output_dir=tmp_path / "out",
        read_1=read_1,
        execute=False,
        verbose=False,
    )

    feature_script = tmp_path / "out" / "sample_f" / "build_amplicon_features.sh"
    script = feature_script.read_text()
    assert "vsearch --cluster_unoise" in script
    assert "vsearch --uchime3_denovo" in script


def _write_primer_test_reads(read_1, read_2, count=10):
    forward = "TAGCCCTATGGGATGCTGCAG" + "ACGT" * 35
    reverse = "AGGACTACGGGGGTATCTAAT" + "TGCA" * 35
    with open(read_1, "w") as forward_handle, open(read_2, "w") as reverse_handle:
        for index in range(count):
            forward_handle.write(
                f"@read_{index}/1\n{forward}\n+\n{'I' * len(forward)}\n"
            )
            reverse_handle.write(
                f"@read_{index}/2\n{reverse}\n+\n{'I' * len(reverse)}\n"
            )


def test_detect_amplicon_primers_accepts_leading_spacers(tmp_path):
    read_1 = tmp_path / "sample_R1.fastq"
    read_2 = tmp_path / "sample_R2.fastq"
    _write_primer_test_reads(read_1, read_2)

    detected = core.Metagenomics(configs.Metagenomics()).detect_amplicon_primers(
        read_1,
        read_2,
    )

    assert detected["name"] == "16S_341F_806R"
    assert detected["forward_match_fraction"] == 1.0
    assert detected["reverse_match_fraction"] == 1.0


def test_detect_amplicon_primers_fails_for_unknown_reads(tmp_path):
    read_1 = tmp_path / "unknown_R1.fastq"
    read_2 = tmp_path / "unknown_R2.fastq"
    sequence = "A" * 160
    quality = "I" * len(sequence)
    read_1.write_text(f"@read/1\n{sequence}\n+\n{quality}\n")
    read_2.write_text(f"@read/2\n{sequence}\n+\n{quality}\n")

    with pytest.raises(ValueError, match="No primer pair passed"):
        core.Metagenomics(configs.Metagenomics()).detect_amplicon_primers(
            read_1,
            read_2,
        )


def test_dada2_preprocessing_detects_primers_from_original_reads_in_dry_run(tmp_path):
    read_1 = tmp_path / "sample_R1.fastq"
    read_2 = tmp_path / "sample_R2.fastq"
    _write_primer_test_reads(read_1, read_2)
    profile = {
        "backend": "local",
        "container": "None",
        "steps": {
            "trim_reads": {
                "backend": "local",
                "settings": {
                    "quality_trim": "cut_right",
                    "quality_window_size": 4,
                    "quality_mean": 20,
                },
            },
            "build_amplicon_features": {
                "backend": "local",
                "settings": {
                    "denoiser": "dada2",
                    "primer_mode": "auto",
                    "threads": 2,
                },
            },
        },
    }

    result = core.Metagenomics(configs.Metagenomics()).preprocess_amplicon_sample(
        sample_name="sample_dada2",
        output_dir=tmp_path / "out",
        read_1=read_1,
        read_2=read_2,
        execution_profile=profile,
        execute=False,
        verbose=False,
    )

    scratch = tmp_path / "out" / "sample_dada2" / "scratch"
    trim_script = (scratch / "trim_reads.sh").read_text()
    feature_script = (scratch / "build_amplicon_features.sh").read_text()
    dada2_script = (
        scratch / "amplicon_preprocess" / "sample_dada2_dada2.R"
    ).read_text()
    primer_manifest = json.loads(
        (scratch / "amplicon_preprocess" / "detected_primers.json").read_text()
    )

    assert "--cut_right --cut_right_window_size 4 --cut_right_mean_quality 20" in trim_script
    assert "cutadapt" in feature_script
    assert "Rscript" in feature_script
    assert "vsearch --cluster_unoise" not in feature_script
    assert "filterAndTrim" in dada2_script
    assert "mergePairs" in dada2_script
    assert "removeBimeraDenovo" in dada2_script
    assert "digest(sequence, algo='sha1'" in dada2_script
    assert primer_manifest["selection"]["name"] == "16S_341F_806R"
    assert primer_manifest["detection_read_1"] == str(read_1)
    assert result["artifacts"]["dada2_stats"].endswith("dada2-stats.tsv")


def test_explicit_primers_override_auto_catalog_detection(tmp_path):
    read_1 = tmp_path / "sample_R1.fastq"
    read_2 = tmp_path / "sample_R2.fastq"
    _write_primer_test_reads(read_1, read_2)

    core.Metagenomics(configs.Metagenomics()).build_amplicon_features(
        read_1=read_1,
        read_2=read_2,
        output_dir=tmp_path / "features",
        sample_name="sample_explicit",
        denoiser="dada2",
        primer_mode="auto",
        forward_primer="ACGTACGT",
        reverse_primer="TGCATGCA",
    )

    primer_manifest = json.loads(
        (tmp_path / "features" / "detected_primers.json").read_text()
    )
    assert primer_manifest["mode"] == "explicit"
    assert primer_manifest["selection"]["forward_primer"] == "ACGTACGT"


def test_batch_sample_to_cod_accepts_sra_and_reads_tables(tmp_path):
    read_1 = tmp_path / "local_R1.fastq"
    read_2 = tmp_path / "local_R2.fastq"
    read_1.write_text("@read_1\nACGTACGT\n+\n!!!!!!!!\n")
    read_2.write_text("@read_1\nACGTACGT\n+\n!!!!!!!!\n")

    sra_manifest = tmp_path / "sra_samples.tsv"
    sra_manifest.write_text("sample\taccession\nsra_sample\tSRR000001\n")
    reads_manifest = tmp_path / "read_samples.tsv"
    reads_manifest.write_text(
        "sample\tread_1\tread_2\n"
        f"fastq_sample\t{read_1}\t{read_2}\n"
    )

    metagenomics = core.Metagenomics(configs.Metagenomics())
    sra_result = metagenomics.batch_sample_to_cod(
        manifest=sra_manifest,
        input_type="sra",
        output_dir=tmp_path / "out",
        stage="download",
        execute=False,
        verbose=False,
    )
    reads_result = metagenomics.batch_sample_to_cod(
        manifest=reads_manifest,
        input_type="reads",
        output_dir=tmp_path / "out",
        stage="preprocess",
        execute=False,
        verbose=False,
    )

    assert set(sra_result["samples"]) == {"sra_sample"}
    assert set(reads_result["samples"]) == {"fastq_sample"}
    assert (tmp_path / "out" / "sra_sample" / "download_sra.sh").exists()
    assert (tmp_path / "out" / "fastq_sample" / "scratch" / "trim_reads.sh").exists()
    assert sra_result["samples"]["sra_sample"]["artifacts"]["download"]["accession"] == "SRR000001"
    assert (tmp_path / "out" / "batch_summary.json").exists()


def test_sra_download_slurm_waits_and_validates_fastq_output(monkeypatch, tmp_path):
    profile = {
        "backend": "local",
        "container": "None",
        "slurm": {},
        "steps": {"download_sra": {"backend": "slurm"}},
    }

    def fake_run(command, capture_output=True, text=True, **kwargs):
        if command[0] == "sbatch":
            accession_dir = tmp_path / "sra" / "SRR000001"
            accession_dir.mkdir(parents=True)
            (accession_dir / "SRR000001_1.fastq").write_text("@r1\nA\n+\n!\n")
            (accession_dir / "SRR000001_2.fastq").write_text("@r1\nT\n+\n!\n")
            return subprocess.CompletedProcess(command, 0, stdout="222\n", stderr="")
        raise AssertionError(f"Unexpected command: {command}")

    monkeypatch.setattr(core.subprocess, "run", fake_run)

    metagenomics = core.Metagenomics(configs.Metagenomics())
    result = metagenomics.run_sra_download_step(
        sample_name="sra_sample",
        output_dir=tmp_path / "out",
        accession="SRR000001",
        sra_dir=tmp_path / "sra",
        execute=True,
        verbose=False,
        execution_profile=profile,
    )

    artifact = result["artifacts"]["download_sra"]
    reads = result["artifacts"]["reads"]

    assert artifact["status"] == "completed"
    assert artifact["job_id"] == "222"
    assert reads["read_1"] == str(tmp_path / "sra" / "SRR000001" / "SRR000001_1.fastq")
    assert reads["read_2"] == str(tmp_path / "sra" / "SRR000001" / "SRR000001_2.fastq")


def test_pending_slurm_amplicon_pipeline_does_not_write_empty_final_cod(monkeypatch, tmp_path):
    read_1 = tmp_path / "sample_R1.fastq"
    read_2 = tmp_path / "sample_R2.fastq"
    read_1.write_text("@read_1\nACGTACGT\n+\n!!!!!!!!\n")
    read_2.write_text("@read_1\nACGTACGT\n+\n!!!!!!!!\n")
    profile = {
        "backend": "local",
        "container": "None",
        "slurm": {},
        "steps": {
            "trim_reads": {"backend": "slurm"},
            "build_amplicon_features": {"backend": "slurm"},
        },
    }

    def fake_run(command, capture_output=True, text=True, **kwargs):
        if command[0] == "sbatch":
            return subprocess.CompletedProcess(command, 0, stdout="Submitted batch job 333\n", stderr="")
        raise AssertionError(f"Unexpected command: {command}")

    monkeypatch.setattr(core.subprocess, "run", fake_run)

    metagenomics = core.Metagenomics(configs.Metagenomics())
    result = metagenomics.sample_to_cod(
        sample_name="sample_pending",
        output_dir=tmp_path / "out",
        mode="amplicon-reads",
        read_1=read_1,
        read_2=read_2,
        execute=True,
        verbose=False,
        execution_profile=profile,
    )

    assert result["status"] == "waiting_for_preprocess"
    assert not (tmp_path / "out" / "sample_pending" / "cod_profile.csv").exists()


def test_batch_workflow_waits_for_download_and_reuses_validated_output(monkeypatch, tmp_path):
    manifest = tmp_path / "samples.tsv"
    manifest.write_text("sample\taccession\nsample_sra\tSRR000001\n")
    profile = {
        "backend": "local",
        "container": "None",
        "slurm": {},
        "steps": {"download_sra": {"backend": "slurm"}},
    }
    submissions = []

    def fake_run(command, capture_output=True, text=True, **kwargs):
        if command[0] == "sbatch":
            submissions.append(command)
            accession_dir = tmp_path / "sra" / "SRR000001"
            accession_dir.mkdir(parents=True)
            (accession_dir / "SRR000001_1.fastq").write_text("@r1\nA\n+\n!\n")
            (accession_dir / "SRR000001_2.fastq").write_text("@r1\nT\n+\n!\n")
            return subprocess.CompletedProcess(command, 0, stdout="444\n", stderr="")
        raise AssertionError(f"Unexpected command: {command}")

    monkeypatch.setattr(core.subprocess, "run", fake_run)

    metagenomics = core.Metagenomics(configs.Metagenomics())
    first = metagenomics.batch_sample_to_cod(
        manifest=manifest,
        input_type="sra",
        output_dir=tmp_path / "out",
        sra_dir=tmp_path / "sra",
        stage="download",
        execute=True,
        verbose=False,
        execution_profile=profile,
    )
    second = metagenomics.batch_sample_to_cod(
        manifest=manifest,
        input_type="sra",
        output_dir=tmp_path / "out",
        sra_dir=tmp_path / "sra",
        stage="download",
        execute=True,
        verbose=False,
        execution_profile=profile,
    )

    state = json.loads((tmp_path / "out" / "workflow_state.json").read_text())

    assert len(submissions) == 1
    assert first["samples"]["sample_sra"]["status"] == "completed"
    assert second["samples"]["sample_sra"]["status"] == "completed"
    assert state["samples"]["sample_sra"]["stages"]["download_sra"]["job_id"] == "444"
    assert (tmp_path / "out" / "workflow_events.jsonl").exists()


def test_batch_runs_sample_chains_concurrently(monkeypatch, tmp_path):
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(
        "sample\taccession\n"
        "sample_1\tSRR001\n"
        "sample_2\tSRR002\n"
        "sample_3\tSRR003\n"
    )
    barrier = threading.Barrier(3)

    def fake_download(*, sample_name, accession, sra_dir, **kwargs):
        barrier.wait(timeout=2)
        accession_dir = pathlib.Path(sra_dir) / accession
        accession_dir.mkdir(parents=True)
        read_1 = accession_dir / f"{accession}_1.fastq"
        read_2 = accession_dir / f"{accession}_2.fastq"
        read_1.write_text("@r1\nA\n+\n!\n")
        read_2.write_text("@r1\nT\n+\n!\n")
        artifact = {"step": "download_sra", "status": "completed", "backend": "slurm"}
        return {
            "artifacts": {
                "download_sra": artifact,
                "accession": accession,
                "reads": {"read_1": str(read_1), "read_2": str(read_2)},
            }
        }

    metagenomics = core.Metagenomics(configs.Metagenomics())
    monkeypatch.setattr(metagenomics, "run_sra_download_step", fake_download)

    result = metagenomics.batch_sample_to_cod(
        manifest=manifest,
        input_type="sra",
        output_dir=tmp_path / "out",
        sra_dir=tmp_path / "sra",
        stage="download",
        execute=True,
        verbose=False,
        sample_workers=3,
    )

    assert {sample["status"] for sample in result["samples"].values()} == {"completed"}
    state = json.loads((tmp_path / "out" / "workflow_state.json").read_text())
    assert set(state["samples"]) == {"sample_1", "sample_2", "sample_3"}


def test_batch_continues_with_next_sample_after_retries_are_exhausted(monkeypatch, tmp_path):
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(
        "sample\taccession\n"
        "bad_sample\tSRR_BAD\n"
        "good_sample\tSRR_GOOD\n"
    )
    profile = {
        "backend": "local",
        "container": "None",
        "slurm": {"retry_delay_seconds": 1},
        "steps": {"download_sra": {"backend": "slurm", "retries": 2}},
    }
    submissions = []

    def fake_run(command, capture_output=True, text=True, **kwargs):
        if command[0] != "sbatch":
            raise AssertionError(f"Unexpected command: {command}")
        submissions.append(command)
        if "bad_sample" in command[-1]:
            return subprocess.CompletedProcess(command, 1, stdout="501\n", stderr="failed after retries")
        accession_dir = tmp_path / "sra" / "SRR_GOOD"
        accession_dir.mkdir(parents=True)
        (accession_dir / "SRR_GOOD_1.fastq").write_text("@r1\nA\n+\n!\n")
        (accession_dir / "SRR_GOOD_2.fastq").write_text("@r1\nT\n+\n!\n")
        return subprocess.CompletedProcess(command, 0, stdout="502\n", stderr="")

    monkeypatch.setattr(core.subprocess, "run", fake_run)
    monkeypatch.setattr(core.time, "sleep", lambda seconds: None)

    result = core.Metagenomics(configs.Metagenomics()).batch_sample_to_cod(
        manifest=manifest,
        input_type="sra",
        output_dir=tmp_path / "out",
        sra_dir=tmp_path / "sra",
        stage="download",
        execute=True,
        verbose=False,
        execution_profile=profile,
    )

    assert len(submissions) == 4
    assert sum("bad_sample" in command[-1] for command in submissions) == 3
    assert sum("good_sample" in command[-1] for command in submissions) == 1
    assert result["samples"]["bad_sample"]["status"] == "failed"
    assert result["samples"]["bad_sample"]["error"]["stage"] == "download_sra"
    assert result["samples"]["good_sample"]["status"] == "completed"
    bad_sbatch = (
        tmp_path / "out" / "bad_sample" / "slurm" / "download_sra.sbatch"
    ).read_text()
    assert "max_retries" not in bad_sbatch
    assert "while true; do" not in bad_sbatch


def test_batch_workflow_resumes_from_cached_preprocess_without_resubmitting(tmp_path):
    reaction_db = tmp_path / "reactions.csv"
    pl.DataFrame(
        [{"EC_Numbers": "1.1.1.1", "e_adm_Reactions": "Uptake of sugars"}]
    ).write_csv(reaction_db)
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(f"sample\tread_1\tgtdb_matches\nsample_ready\treads.fastq\t{tmp_path / 'missing_matches.blast'}\n")
    preprocess = tmp_path / "out" / "sample_ready" / "scratch" / "amplicon_preprocess"
    preprocess.mkdir(parents=True)
    (preprocess / "feature-table.tsv").write_text("#OTU ID\tsample_ready\nasv_1\t10\n")
    (preprocess / "rep-seqs.fasta").write_text(">asv_1\nACGT\n")

    metagenomics = core.Metagenomics(configs.Metagenomics(csv_reaction_db=reaction_db))
    result = metagenomics.batch_sample_to_cod(
        manifest=manifest,
        input_type="reads",
        output_dir=tmp_path / "out",
        stage="all",
        execute=False,
        verbose=False,
        genome_alignments=tmp_path / "missing_alignments",
    )

    stages = result["samples"]["sample_ready"]["stages"]

    assert stages["preprocess"]["status"] == "completed"
    assert result["samples"]["sample_ready"]["status"] == "waiting_for_gtdb_alignment"
    assert not (tmp_path / "out" / "sample_ready" / "scratch" / "trim_reads.sh").exists()


def test_dada2_profile_does_not_reuse_legacy_vsearch_feature_cache(tmp_path):
    read_1 = tmp_path / "sample_R1.fastq"
    read_2 = tmp_path / "sample_R2.fastq"
    _write_primer_test_reads(read_1, read_2)
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(
        "sample\tread_1\tread_2\n"
        f"sample_ready\t{read_1}\t{read_2}\n"
    )
    preprocess = tmp_path / "out" / "sample_ready" / "scratch" / "amplicon_preprocess"
    preprocess.mkdir(parents=True)
    (preprocess / "feature-table.tsv").write_text("#OTU ID\tsample_ready\nlegacy_1\t10\n")
    (preprocess / "rep-seqs.fasta").write_text(">legacy_1\nACGT\n")
    profile = {
        "backend": "local",
        "container": "None",
        "steps": {
            "build_amplicon_features": {
                "backend": "local",
                "settings": {"denoiser": "dada2", "primer_mode": "auto"},
            }
        },
    }

    result = core.Metagenomics(configs.Metagenomics()).batch_sample_to_cod(
        manifest=manifest,
        input_type="reads",
        output_dir=tmp_path / "out",
        stage="preprocess",
        execute=False,
        verbose=False,
        execution_profile=profile,
    )

    assert result["samples"]["sample_ready"]["status"] == "waiting_for_preprocess"
    assert (tmp_path / "out" / "sample_ready" / "scratch" / "trim_reads.sh").exists()
    assert (tmp_path / "out" / "sample_ready" / "scratch" / "build_amplicon_features.sh").exists()


def test_batch_workflow_reuses_completed_gtdb_alignment(monkeypatch, tmp_path):
    reaction_db = tmp_path / "reactions.csv"
    pl.DataFrame(
        [{"EC_Numbers": "1.1.1.1", "e_adm_Reactions": "Uptake of sugars"}]
    ).write_csv(reaction_db)
    gtdb = tmp_path / "gtdb.fa"
    gtdb.write_text(">GB_GCA_000001.1~contig\nACGT\n")
    manifest = tmp_path / "samples.tsv"
    manifest.write_text("sample\tread_1\nsample_ready\treads.fastq\n")
    preprocess = tmp_path / "out" / "sample_ready" / "scratch" / "amplicon_preprocess"
    preprocess.mkdir(parents=True)
    (preprocess / "feature-table.tsv").write_text("#OTU ID\tsample_ready\nasv_1\t10\n")
    (preprocess / "rep-seqs.fasta").write_text(">asv_1\nACGT\n")
    profile = {
        "backend": "local",
        "container": "None",
        "slurm": {},
        "steps": {"align_to_gtdb": {"backend": "slurm"}},
    }
    submissions = []

    def fake_run(command, capture_output=True, text=True, **kwargs):
        if command[0] == "sbatch":
            submissions.append(command)
            (tmp_path / "out" / "sample_ready" / "scratch" / "matches.blast").write_text("")
            return subprocess.CompletedProcess(command, 0, stdout="777\n", stderr="")
        raise AssertionError(f"Unexpected command: {command}")

    monkeypatch.setattr(core.subprocess, "run", fake_run)

    config = configs.Metagenomics(csv_reaction_db=reaction_db)
    config.gtdb_dir_fasta = str(gtdb)
    metagenomics = core.Metagenomics(config)
    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir()
    first = metagenomics.batch_sample_to_cod(
        manifest=manifest,
        input_type="reads",
        output_dir=tmp_path / "out",
        stage="all",
        execute=True,
        verbose=False,
        execution_profile=profile,
        genomes_dir=genomes_dir,
    )
    second = metagenomics.batch_sample_to_cod(
        manifest=manifest,
        input_type="reads",
        output_dir=tmp_path / "out",
        stage="all",
        execute=True,
        verbose=False,
        execution_profile=profile,
        genomes_dir=genomes_dir,
    )
    state = json.loads((tmp_path / "out" / "workflow_state.json").read_text())

    assert len(submissions) == 1
    assert first["samples"]["sample_ready"]["status"] == "waiting_for_cod"
    assert second["samples"]["sample_ready"]["status"] == "waiting_for_cod"
    assert (
        tmp_path
        / "out"
        / "sample_ready"
        / "scratch"
        / "amplicon_downstream_inputs.json"
    ).exists()
    assert state["samples"]["sample_ready"]["stages"]["align_to_gtdb"]["status"] == "completed"
