import json

import polars as pl

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


def test_task_manager_adds_slurm_dependencies(tmp_path):
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

    assert artifact["dependency_job_ids"] == ["12345"]
    assert artifact["dependencies"] == ["parent"]
    assert "#SBATCH --dependency=afterok:12345" in sbatch.read_text()


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
        forward_primer="GTGYCAGCMGCCGCGGTAA",
        reverse_primer="GGACTACNVGGGTWTCTAAT",
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
    assert "cutadapt" in trim_script.read_text()
    assert "-m 80" in trim_script.read_text()
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
