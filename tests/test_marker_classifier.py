import json
from pathlib import Path

import polars as pl
import pytest
from click.testing import CliRunner

from adtoolbox import cli, configs, core
from adtoolbox.markers import MarkerCatalog, MarkerCatalogError, MarkerClassifier, alignment_marker_hits, build_marker_hmm_db, build_marker_protein_db, hmmer_marker_hits, marker_ids_from_fasta, read_hmm_score_cutoffs


@pytest.fixture()
def classifier():
    return MarkerClassifier()


def _score(result, genome, group):
    return result.filter((pl.col("genome_id") == genome) & (pl.col("cod_group") == group))["score"][0]


def test_default_catalog_preserves_all_existing_cod_groups(classifier):
    expected = set(configs.E_ADM_MICROBIAL_GROUPS_MAPPING.values())
    assert set(classifier.catalog.cod_groups) == expected
    assert set(classifier.catalog.groups) == expected


@pytest.mark.parametrize(
    ("group", "profile"),
    [
        ("X_ch", ["GH5", "CBM3", "bgl"]),
        ("X_pr", ["aprE", "pepN", "oppA", "oppB"]),
        ("X_li", ["lipase", "fadD"]),
        ("X_su", ["ptsG", "pgi", "pfkA", "gapA", "pyk", "pfor"]),
        ("X_aa", ["oppA", "grdA", "grdB", "pfor"]),
        ("X_fa", ["fadD", "fadE", "fadB", "fadA", "etfA"]),
        ("X_VFA_deg", ["thl", "hbd", "crt", "bcd", "etfA", "etfB", "hydA", "rnfB"]),
        ("X_et", ["adh", "aldh", "acs"]),
        ("X_lac", ["lctD", "lctB", "pfor"]),
    ],
)
def test_expanded_biochemistry_panels(classifier, group, profile):
    result = classifier.classify_presence({"organism": profile})
    assert _score(result, "organism", group) > 0


def test_shared_mcr_marker_does_not_classify_methanogen(classifier):
    result = classifier.classify_presence({"incomplete_archaeon": ["mcrA"]})
    assert _score(result, "incomplete_archaeon", "X_Me_ac") == 0
    assert _score(result, "incomplete_archaeon", "X_Me_CO2") == 0


def test_partial_pathway_evidence_is_retained_but_not_called_credible(classifier):
    result = classifier.classify_presence({"partial_proteolysis": ["pepN", "oppA"]})
    row = result.filter(
        (pl.col("genome_id") == "partial_proteolysis") & (pl.col("cod_group") == "X_pr")
    ).to_dicts()[0]
    assert 0 < row["score"] < 0.5
    assert row["credible"] is False
    assert row["confidence"] == "partial"
    assert row["weighted_completeness"] > 0
    assert row["marker_coverage"] > 0


def test_community_marker_pool_weights_genomes_before_pathway_scoring(classifier):
    hits = pl.DataFrame(
        {
            "genome_id": ["common", "rare", "rare"],
            "gene_id": ["pep_common", "pep_rare_1", "pep_rare_2"],
            "marker_id": ["pepN", "oppA", "oppA"],
        }
    )

    markers, panels, profile = classifier.community_profile_from_hits(
        hits,
        {"common": 0.8, "rare": 0.2},
    )

    marker_values = dict(markers.select("marker_id", "weighted_abundance").iter_rows())
    assert marker_values["pepN"] == pytest.approx(0.8)
    # Two copies in one genome remain one marker-family contribution.
    assert marker_values["oppA"] == pytest.approx(0.2)
    selected_protein = panels.filter(
        (pl.col("cod_group") == "X_pr") & pl.col("selected")
    ).to_dicts()[0]
    expected = (0.8 * 1.0 + 0.2 * 2.0) / 11.0
    assert selected_protein["panel_id"] == "proteolysis_and_peptide_import"
    assert selected_protein["potential"] == pytest.approx(expected)
    assert profile["X_pr"] == pytest.approx(expected)


def test_community_marker_pool_is_not_zeroed_by_genome_level_strict_gate(classifier):
    hits = pl.DataFrame(
        {
            "genome_id": ["archaeon"],
            "gene_id": ["mcr_gene"],
            "marker_id": ["mcrA"],
        }
    )

    _, panels, profile = classifier.community_profile_from_hits(hits, {"archaeon": 0.25})

    assert profile["X_Me_CO2"] > 0
    hydrogenotrophic = panels.filter(
        (pl.col("cod_group") == "X_Me_CO2") & pl.col("selected")
    ).to_dicts()[0]
    assert hydrogenotrophic["credible"] is False


def test_community_marker_pool_normalizes_count_like_abundances(classifier):
    hits = pl.DataFrame(
        {
            "genome_id": ["g1", "g2"],
            "gene_id": ["gene_1", "gene_2"],
            "marker_id": ["pepN", "pepN"],
        }
    )

    markers, _, _ = classifier.community_profile_from_hits(hits, {"g1": 30, "g2": 70})

    pep = markers.filter(pl.col("marker_id") == "pepN").to_dicts()[0]
    assert pep["weighted_abundance"] == pytest.approx(1.0)


def test_acetoclastic_methanogen_profile(classifier):
    # Methanosarcina-like acetate activation and CODH/ACS evidence.
    markers = ["K00399", "mcrB", "mtrA", "acsA", "cdhC", "cdhD", "cdhE"]
    result = classifier.classify_presence({"Methanosarcina_like": markers})
    assert _score(result, "Methanosarcina_like", "X_Me_ac") > 0.7
    assert _score(result, "Methanosarcina_like", "X_Me_CO2") == 0


def test_hydrogenotrophic_methanogen_profile(classifier):
    # Methanobrevibacter-like CO2-reduction pathway coverage.
    markers = ["mcrA", "mcrB", "mtrA", "fmdA", "ftr", "mch", "mtd", "mer"]
    result = classifier.classify_presence({"Methanobrevibacter_like": markers})
    assert _score(result, "Methanobrevibacter_like", "X_Me_CO2") > 0.9
    assert _score(result, "Methanobrevibacter_like", "X_Me_ac") == 0


def test_ethanol_chain_elongator_profile(classifier):
    # Clostridium kluyveri-like reverse beta oxidation plus ethanol utilization.
    markers = ["thl", "hbd", "crt", "bcd", "etfA", "etfB", "cat", "adhE", "aldh"]
    result = classifier.classify_presence({"C_kluyveri_like": markers})
    assert _score(result, "C_kluyveri_like", "X_chain_et") > 0.8
    assert _score(result, "C_kluyveri_like", "X_chain_lac") == 0


def test_lactate_chain_elongator_profile(classifier):
    # Megasphaera elsdenii-like reverse beta oxidation plus lactate oxidation.
    markers = ["atoB", "hbd", "crt", "bcd", "etfA", "etfB", "but", "lctD", "lctB", "pfor"]
    result = classifier.classify_presence({"M_elsdenii_like": markers})
    assert _score(result, "M_elsdenii_like", "X_chain_lac") > 0.8
    assert _score(result, "M_elsdenii_like", "X_chain_et") == 0


def test_hit_thresholds_and_all_group_rows(classifier):
    hits = pl.DataFrame(
        {
            "genome_id": ["g1", "g1"],
            "marker_id": ["mcrA", "ftr"],
            "bits": [100.0, 10.0],
        }
    )
    result = classifier.classify_hits(hits, min_bits=50)
    assert result.height == len(classifier.catalog.cod_groups)
    assert result.filter(pl.col("matched_markers").str.contains("ftr")).height == 0


def test_mmseqs_target_headers_resolve_marker_aliases(classifier):
    alignment = pl.DataFrame(
        {
            "query": ["gene_1", "gene_2", "gene_3"],
            "target": ["sp|P12345|mcrA", "ref|Q12345|GH5", "sp|X00000|not_a_marker"],
            "bits": [100.0, 80.0, 200.0],
            "evalue": [1e-30, 1e-20, 1e-50],
        }
    )
    hits = alignment_marker_hits(alignment, classifier.catalog, entity_id="g1", min_bits=40, max_evalue=1e-5)
    assert set(hits["marker_id"]) == {"mcrA", "cazy_GH5"}


def test_marker_protein_fasta_coverage_reader(tmp_path, classifier):
    fasta = tmp_path / "markers.fasta"
    fasta.write_text(">ref_1|K00399\nMPEPTIDE\n>ref_2|GH5 description\nMSEQUENCE\n")
    assert marker_ids_from_fasta(fasta, classifier.catalog) == {"mcrA", "cazy_GH5"}


def test_bundled_marker_databases_cover_the_full_catalog(classifier):
    config = configs.Metagenomics()
    assert marker_ids_from_fasta(config.protein_db, classifier.catalog) == set(classifier.catalog.markers)
    assert Path(config.marker_hmm_db).is_file()
    assert Path(config.marker_hmm_cutoffs).is_file()
    cutoffs = read_hmm_score_cutoffs(config.marker_hmm_cutoffs)
    assert cutoffs["mcrA__K00399"] == {"score_threshold": 775.53, "score_type": "full"}
    assert cutoffs["bgl__K01188"]["score_type"] == "domain"


def test_build_marker_protein_db_reheaders_curated_references(tmp_path, classifier):
    source = tmp_path / "methane.faa"
    source.write_text(">P12345 reviewed protein\nMPEPTIDE\n")
    manifest = tmp_path / "manifest.csv"
    pl.DataFrame({"marker_id": ["K00399"], "source_fasta": [source.name]}).write_csv(manifest)
    output = tmp_path / "Marker_Protein_DB.fasta"
    counts = build_marker_protein_db(manifest, output, classifier.catalog)
    assert counts == {"mcrA": 1}
    assert output.read_text() == ">P12345|mcrA reviewed protein\nMPEPTIDE\n"


def test_build_and_parse_marker_hmm_database(tmp_path, classifier):
    source = tmp_path / "source.hmm"
    source.write_text(
        "HMMER3/f\nNAME  K00399\nACC   K00399.1\nLENG  100\n//\n"
        "HMMER3/f\nNAME  GH5\nLENG  120\n//\n"
    )
    manifest = tmp_path / "manifest.csv"
    pl.DataFrame(
        {
            "marker_id": ["mcrA", "GH5"],
            "source_hmm": [source.name, source.name],
            "profile_id": ["K00399", "GH5"],
            "score_threshold": [50.0, None],
            "score_type": ["domain", "full"],
        }
    ).write_csv(manifest)
    output = tmp_path / "markers.hmm"
    cutoffs = tmp_path / "cutoffs.csv"
    counts = build_marker_hmm_db(manifest, output, classifier.catalog, cutoff_output=cutoffs)
    assert counts == {"mcrA": 1, "cazy_GH5": 1}
    assert "NAME  mcrA__K00399" in output.read_text()
    assert pl.read_csv(cutoffs).to_dicts() == [
        {"profile_id": "mcrA__K00399", "score_threshold": 50.0, "score_type": "domain"}
    ]

    domtbl = tmp_path / "hits.domtbl"
    domtbl.write_text(
        "gene_1 - 300 mcrA__K00399 K00399 100 1e-50 90.0 0.0 1 1 1e-45 1e-45 85.0 0.0 1 90 5 94 2 96 0.98 methane marker\n"
    )
    hits = hmmer_marker_hits(
        domtbl,
        classifier.catalog,
        entity_id="g1",
        score_cutoffs={"mcrA__K00399": {"score_threshold": 80, "score_type": "domain"}},
    )
    assert hits.select("genome_id", "gene_id", "marker_id").to_dicts() == [
        {"genome_id": "g1", "gene_id": "gene_1", "marker_id": "mcrA"}
    ]


def test_hmmer_backend_prepares_prodigal_and_hmmsearch(tmp_path):
    genome = tmp_path / "genome.fna"
    genome.write_text(">contig\nATGAAATAG\n")
    hmm_db = tmp_path / "markers.hmm"
    hmm_db.write_text("HMMER3/f\nNAME  mcrA\nLENG  10\n//\n")
    metagenomics = core.Metagenomics(
        configs.Metagenomics(marker_backend="hmmer", marker_hmm_db=hmm_db)
    )
    script, output = metagenomics.annotate_genome_with_marker_hmms(genome, tmp_path / "out", "g1", threads=4)
    assert "prodigal" in script
    assert "hmmsearch" in script
    assert "--cpu 4" in script
    assert output.endswith("Marker_Hits_hmmer_g1~catalog-0.3.0.domtbl")


def test_core_alignment_prefers_direct_marker_pipeline(tmp_path):
    alignment = tmp_path / "markers.tsv"
    pl.DataFrame(
        {
            "query": [f"gene_{i}" for i in range(8)],
            "target": [f"ref|{marker}" for marker in ("mcrA", "mcrB", "mtrA", "fwdA", "ftr", "mch", "mtd", "mer")],
            "evalue": [1e-30] * 8,
            "bits": [100.0] * 8,
        }
    ).write_csv(alignment, separator="\t")
    metagenomics = core.Metagenomics(configs.Metagenomics(bit_score=40, e_value=1e-5))
    evidence_kind, counts = metagenomics.functional_counts_from_alignment(alignment)
    profile = metagenomics.cod_from_alignment(alignment)
    assert evidence_kind == "markers"
    assert counts["mcrA"] == 1
    assert profile["X_Me_CO2"] == pytest.approx(1.0)


def test_sample_aggregation_weights_genomes_and_reports_qc(classifier):
    scores = classifier.classify_presence(
        {
            "methanogen": ["mcrA", "mcrB", "mtrA", "fwdA", "ftr", "mch", "mtd", "mer"],
            "unclassified": ["mcrA"],
        }
    )
    abundances = pl.DataFrame(
        {"sample": ["s1", "s1"], "genome_id": ["methanogen", "unclassified"], "abundance": [25.0, 75.0]}
    )
    profile, qc = classifier.aggregate_samples(scores, abundances)
    assert profile.height == len(classifier.catalog.cod_groups)
    assert _score(profile.rename({"sample": "genome_id", "abundance": "score"}), "s1", "X_Me_CO2") == 1.0
    assert qc["classified_fraction"][0] == pytest.approx(0.25)
    assert qc["groups_detected"][0] == 1


def test_invalid_catalog_rejects_unknown_panel_marker():
    payload = {
        "schema_version": 1,
        "cod_groups": ["X_su"],
        "markers": {"known": {"aliases": []}},
        "groups": {
            "X_su": {
                "panels": [
                    {"id": "bad", "marker_weights": {"unknown": 1}, "minimum_markers": 1, "minimum_score": 0}
                ]
            }
        },
    }
    with pytest.raises(MarkerCatalogError, match="unknown markers"):
        MarkerCatalog.from_dict(payload)


def test_marker_cli_classifies_and_aggregates(tmp_path):
    hits = tmp_path / "hits.csv"
    scores = tmp_path / "scores.csv"
    abundance = tmp_path / "abundance.csv"
    profile = tmp_path / "profile.csv"
    pl.DataFrame(
        {
            "genome_id": ["g1"] * 8,
            "marker_id": ["mcrA", "mcrB", "mtrA", "fwdA", "ftr", "mch", "mtd", "mer"],
        }
    ).write_csv(hits)
    pl.DataFrame({"sample": ["s1"], "genome_id": ["g1"], "abundance": [1.0]}).write_csv(abundance)

    runner = CliRunner()
    classified = runner.invoke(cli.main, ["metagenomics", "classify-markers", "--hits", str(hits), "--output", str(scores)])
    assert classified.exit_code == 0, classified.output
    aggregated = runner.invoke(
        cli.main,
        ["metagenomics", "aggregate-marker-cod", "--scores", str(scores), "--abundances", str(abundance), "--output", str(profile)],
    )
    assert aggregated.exit_code == 0, aggregated.output
    assert profile.exists()
    assert profile.with_name("profile_qc.csv").exists()
