import json
import os
import subprocess
from pathlib import Path

import click
import numpy as np
import pandas as pd
import rich
from rich import markdown
from rich.console import Console
from rich.prompt import Prompt
from rich.table import Table

from adtoolbox import __version__, adm, configs, core, utils


CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}
console = Console()


def _database(**config_overrides):
    return core.Database(config=configs.Database(**config_overrides))


def _metagenomics_config(protein_db=None, amplicon_to_genome_db=None):
    config = configs.Metagenomics()
    if protein_db:
        config.protein_db = protein_db
    if amplicon_to_genome_db:
        config.amplicon2genome_db = amplicon_to_genome_db
        matches = list(Path(amplicon_to_genome_db).rglob(config.gtdb_dir))
        config.gtdb_dir_fasta = str(matches[0]) if matches else None
    return config


def _prompt_path(value, prompt, *, exists=None, file_okay=True, dir_okay=True, writable=False):
    path_type = click.Path(exists=bool(exists), file_okay=file_okay, dir_okay=dir_okay, writable=writable)
    path = value or click.prompt(prompt, type=path_type)
    return os.path.abspath(os.path.expanduser(path_type.convert(path, None, None)))


def _database_config_from_dir(output_dir):
    return {"database_dir": output_dir}


def _model_paths_from_dir(parameters_dir, prefix):
    return {
        "model_parameters": os.path.join(parameters_dir, f"{prefix}_model_parameters.json"),
        "base_parameters": os.path.join(parameters_dir, f"{prefix}_base_parameters.json"),
        "initial_conditions": os.path.join(parameters_dir, f"{prefix}_initial_conditions.json"),
        "inlet_conditions": os.path.join(parameters_dir, f"{prefix}_inlet_conditions.json"),
        "reactions": os.path.join(parameters_dir, f"{prefix}_reactions.json"),
        "species": os.path.join(parameters_dir, f"{prefix}_species.json"),
    }


def _resolve_model_paths(parameters_dir, prefix, **overrides):
    if parameters_dir:
        parameters_dir = _prompt_path(parameters_dir, "ADM parameter directory", exists=True, file_okay=False, dir_okay=True)
    elif not any(overrides.values()):
        parameters_dir = _prompt_path(None, "ADM parameter directory", exists=True, file_okay=False, dir_okay=True)

    paths = _model_paths_from_dir(parameters_dir, prefix) if parameters_dir else {}
    resolved = {}
    for key, value in overrides.items():
        resolved[key] = value or paths.get(key)
        if not resolved[key]:
            resolved[key] = _prompt_path(None, key.replace("_", " ").title() + " JSON path", exists=True, file_okay=True, dir_okay=False)
    return resolved


def _load_json(path):
    with open(path) as f:
        return json.load(f)


def _print_feed_table(feeds):
    if not feeds:
        return

    feed_table = Table(title="Feed Database", safe_box=True, expand=True)
    for column in feeds[0].to_dict().keys():
        feed_table.add_column(column, justify="center", style="cyan", max_width=20)
    for feed in feeds:
        feed_table.add_row(*map(str, feed.to_dict().values()))
    console.print(feed_table)


def _write_representative_genomes(results, output_dir, output_format):
    if output_format == "csv":
        pd.DataFrame.from_dict(results, orient="index").to_csv(
            os.path.join(output_dir, "representative_genomes.csv")
        )
    elif output_format == "json":
        with open(os.path.join(output_dir, "representative_genomes.json"), "w") as f:
            json.dump(results, f)
    else:
        raise click.ClickException("Please provide a valid format for the output file")


def _report_adm_solution(model, solution, report):
    if report in (None, "dash"):
        model.dash_app(solution)
    elif report == "csv":
        address = Prompt.ask("\n[yellow]Where do you want to save the csv file? ")
        model.csv_report(solution, address)
    else:
        raise click.ClickException("Please provide a valid report option")


@click.group(
    name="ADToolBox",
    context_settings=CONTEXT_SETTINGS,
    help="ADToolBox, a toolbox for anaerobic digestion modeling",
    no_args_is_help=True,
    invoke_without_command=True,
)
@click.version_option(__version__, "-v", "--version", prog_name="ADToolBox")
def main():
    pass


@main.group(name="Database", help="Build or download databases required by ADToolbox.", no_args_is_help=True)
def database():
    pass


@database.command(name="initialize-feed-db", help="Initialize the Feed DB.")
@click.option("--feed-db", help="Path where the feed database TSV should be created.")
def initialize_feed_db(feed_db):
    feed_db = _prompt_path(feed_db, "Feed database TSV path", file_okay=True, dir_okay=False, writable=True)
    _database(feed_db=feed_db).initialize_feed_db()


@database.command(name="add-feed", help="Add a feed to the feed database.")
@click.option("--feed-db", help="Path to the feed database TSV.")
@click.option("-n", "--name", help="Name of the feed to be added to the database.")
@click.option("-c", "--carbohydrates", required=True, type=float, help="Carbohydrate content in percent.")
@click.option("-p", "--proteins", required=True, type=float, help="Protein content in percent.")
@click.option("-l", "--lipids", required=True, type=float, help="Lipid content in percent.")
@click.option("-t", "--tss", required=True, type=float, help="Total suspended solid content in percent.")
@click.option("-s", "--si", required=True, type=float, help="Soluble inert content in percent.")
@click.option("-x", "--xi", required=True, type=float, help="Particulate inert content in percent.")
@click.option("-r", "--reference", required=True, help="Reference where the numbers come from.")
def add_feed(feed_db, name, carbohydrates, proteins, lipids, tss, si, xi, reference):
    feed_db = _prompt_path(feed_db, "Feed database TSV path", file_okay=True, dir_okay=False, writable=True)
    feed = core.Feed(
        name=name,
        carbohydrates=carbohydrates,
        proteins=proteins,
        lipids=lipids,
        tss=tss,
        si=si,
        xi=xi,
        reference=reference,
    )
    _database(feed_db=feed_db).add_feed_to_feed_db(feed=feed)


@database.command(name="show-feed-db", help="Show the feed database.")
@click.option("--feed-db", help="Path to the feed database TSV.")
@click.option("-f", "--filter", "feed_filter", help="Filter the feed database by feed name.")
def show_feed_db(feed_db, feed_filter):
    feed_db = _prompt_path(feed_db, "Feed database TSV path", exists=True, file_okay=True, dir_okay=False)
    db = _database(feed_db=feed_db)
    if feed_filter:
        feeds = db.get_feed_from_feed_db(field_name="name", field_value=feed_filter)
    else:
        feeds = db.get_feed_from_feed_db(field_name="name", query="")
    _print_feed_table(feeds)


@database.command(name="initialize-metagenomics-studies-db", help="Initialize the Metagenomics Studies DB.")
@click.option("--studies-db", help="Path where the metagenomics studies TSV should be created.")
def initialize_metagenomics_studies_db(studies_db):
    studies_db = _prompt_path(studies_db, "Metagenomics studies TSV path", file_okay=True, dir_okay=False, writable=True)
    _database(studies_local={"metagenomics_studies": studies_db}).initialize_metagenomics_studies_db()


@database.command(name="add-metagenomics-study", help="Add a metagenomics study to the database.")
@click.option("--studies-db", help="Path to the metagenomics studies TSV.")
@click.option("-n", "--name", required=True, help="Metagenomics study name.")
@click.option("-t", "--type", "study_type", required=True, help="Metagenomics study type.")
@click.option("-m", "--microbiome", required=True, help="Microbiome where the study belongs.")
@click.option("-s", "--sample-accession", required=True, help="SRA accession ID for the sample.")
@click.option("-c", "--comments", required=True, help="Comments on the study of interest.")
@click.option("-p", "--study-accession", required=True, help="SRA accession ID for the project.")
def add_metagenomics_study(studies_db, name, study_type, microbiome, sample_accession, comments, study_accession):
    studies_db = _prompt_path(studies_db, "Metagenomics studies TSV path", file_okay=True, dir_okay=False, writable=True)
    study = core.MetagenomicsStudy(
        name=name,
        study_type=study_type,
        microbiome=microbiome,
        sample_accession=sample_accession,
        comments=comments,
        study_accession=study_accession,
    )
    _database(studies_local={"metagenomics_studies": studies_db}).add_metagenomics_study_to_metagenomics_studies_db(metagenomics_study=study)


@database.command(name="initialize-protein-db", help="Initialize the protein database.")
@click.option("--protein-db", help="Path where the protein FASTA database should be created.")
def initialize_protein_db(protein_db):
    protein_db = _prompt_path(protein_db, "Protein database FASTA path", file_okay=True, dir_okay=False, writable=True)
    _database(protein_db=protein_db).initialize_protein_db()


@database.command(name="add-protein", help="Add a protein to the protein database.")
@click.option("--protein-db", help="Path to the protein FASTA database.")
@click.option("-i", "--uniport-id", required=True, help="Uniport ID of the protein.")
@click.option("-n", "--name", required=True, help="Name attached to the protein, usually an EC number.")
def add_protein(protein_db, uniport_id, name):
    protein_db = _prompt_path(protein_db, "Protein database FASTA path", file_okay=True, dir_okay=False, writable=True)
    _database(protein_db=protein_db).add_protein_to_protein_db(protein_id=uniport_id, header_tail=name)


@database.command(name="download-reaction-db", help="Download the reaction database in CSV format.")
@click.option("--reaction-db", help="Path where the reaction metadata CSV should be saved.")
def download_reaction_db(reaction_db):
    reaction_db = _prompt_path(reaction_db, "Reaction metadata CSV path", file_okay=True, dir_okay=False, writable=True)
    _database(csv_reaction_db=reaction_db).download_reaction_database()


@database.command(name="download-seed-reaction-db", help="Download the seed reaction database in JSON format.")
@click.option("--seed-reaction-db", help="Path where the SEED reactions JSON should be saved.")
@click.option("--seed-compound-db", help="Path where the SEED compounds JSON should be saved.")
def download_seed_reaction_db(seed_reaction_db, seed_compound_db):
    seed_reaction_db = _prompt_path(seed_reaction_db, "SEED reactions JSON path", file_okay=True, dir_okay=False, writable=True)
    seed_compound_db = _prompt_path(seed_compound_db, "SEED compounds JSON path", file_okay=True, dir_okay=False, writable=True)
    _database(reaction_db=seed_reaction_db, compound_db=seed_compound_db).download_seed_databases()


@database.command(name="build-protein-db", help="Generate the protein database for ADToolbox.")
@click.option("--reaction-db", help="Path to the reaction metadata CSV.")
@click.option("--protein-db", help="Path where the protein FASTA database should be saved.")
def build_protein_db(reaction_db, protein_db):
    reaction_db = _prompt_path(reaction_db, "Reaction metadata CSV path", exists=True, file_okay=True, dir_okay=False)
    protein_db = _prompt_path(protein_db, "Protein database FASTA path", file_okay=True, dir_okay=False, writable=True)
    db = _database(csv_reaction_db=reaction_db, protein_db=protein_db)
    ecs = core.Database.ec_from_csv(reaction_db)
    db.protein_db_from_ec(ecs)
    rich.print("[green]Protein DB built successfully")


@database.command(name="download-feed-db", help="Download the feed database.")
@click.option("--feed-db", help="Path where the feed database TSV should be saved.")
def download_feed_db(feed_db):
    feed_db = _prompt_path(feed_db, "Feed database TSV path", file_okay=True, dir_okay=False, writable=True)
    _database(feed_db=feed_db).download_feed_database()


@database.command(name="download-protein-db", help="Download the protein database in FASTA format.")
@click.option("--protein-db", help="Path where the protein FASTA database should be saved.")
def download_protein_db(protein_db):
    protein_db = _prompt_path(protein_db, "Protein database FASTA path", file_okay=True, dir_okay=False, writable=True)
    _database(protein_db=protein_db).download_protein_database()


@database.command(name="download-amplicon-to-genome-dbs", help="Download amplicon-to-genome databases.")
@click.option("--output-dir", help="Directory where the amplicon-to-genome databases should be saved.")
def download_amplicon_to_genome_dbs(output_dir):
    output_dir = _prompt_path(output_dir, "Amplicon-to-genome database directory", file_okay=False, dir_okay=True, writable=True)
    _database(amplicon_to_genome_db=output_dir).download_amplicon_to_genome_db()


@database.command(name="download-all-databases", help="Download all databases required by ADToolbox.")
@click.option("--output-dir", help="Directory where all downloaded databases should be saved.")
def download_all_databases(output_dir):
    output_dir = _prompt_path(output_dir, "Database download directory", file_okay=False, dir_okay=True, writable=True)
    _database(**_database_config_from_dir(output_dir)).download_all_databases()


@main.group(
    name="Metagenomics",
    help="Import and process metagenomics data from the command line.",
    no_args_is_help=True,
)
def metagenomics():
    pass


@metagenomics.command(name="download_from_sra", help="Download metagenomics data from SRA.")
@click.option("-s", "--sample-accession", required=True, help="SRA accession ID for the sample.")
@click.option("-o", "--output-dir", help="Output directory for downloaded data.")
@click.option("-c", "--container", default="None", show_default=True, help="Container: None, docker, or singularity.")
def download_from_sra(sample_accession, output_dir, container):
    output_dir = _prompt_path(output_dir, "Downloaded SRA output directory", file_okay=False, dir_okay=True, writable=True)
    config = _metagenomics_config()
    prefetch_script, _ = core.Metagenomics(config).seqs_from_sra(
        accession=sample_accession,
        target_dir=output_dir,
        container=container,
    )
    subprocess.run(prefetch_script, shell=True)


@metagenomics.command(name="download_genome", help="Download a genome from NCBI.")
@click.option("-g", "--genome-accession", required=True, help="NCBI accession ID for the genome.")
@click.option("-o", "--output-dir", help="Output directory for downloaded data.")
@click.option("-c", "--container", default="None", show_default=True, help="Container: None, docker, or singularity.")
def download_genome(genome_accession, output_dir, container):
    output_dir = _prompt_path(output_dir, "Downloaded genome output directory", file_okay=False, dir_okay=True, writable=True)
    config = _metagenomics_config()
    script = core.Metagenomics(config).download_genome(
        identifier=genome_accession,
        output_dir=output_dir,
        container=container,
    )
    subprocess.run(script, shell=True)


@metagenomics.command(name="align-genome", help="Align one genome to the ADToolbox protein database.")
@click.option("-n", "--name", required=True, help="Name for the genome being aligned.")
@click.option("-i", "--input-file", help="Genome FASTA or JSON file containing genome information.")
@click.option("-o", "--output-dir", help="Output directory for alignment results.")
@click.option("-c", "--container", default="None", show_default=True, help="Container: None, docker, or singularity.")
@click.option("-d", "--protein-db", "--protein-db-dir", help="Path to the protein FASTA database.")
def align_genome(name, input_file, output_dir, container, protein_db):
    input_file = _prompt_path(input_file, "Genome input file", exists=True, file_okay=True, dir_okay=False)
    output_dir = _prompt_path(output_dir, "Alignment output directory", file_okay=False, dir_okay=True, writable=True)
    protein_db = _prompt_path(protein_db, "Protein database FASTA path", exists=True, file_okay=True, dir_okay=False)
    config = _metagenomics_config(protein_db)
    script = core.Metagenomics(config).align_genome_to_protein_db(
        address=input_file,
        outdir=output_dir,
        name=name,
        container=container,
    )
    subprocess.run(script, shell=True)


@metagenomics.command(name="align-multiple-genomes", help="Align multiple genomes to the ADToolbox protein database.")
@click.option("-i", "--input-file", help="JSON file containing genome names and input files.")
@click.option("-o", "--output-dir", help="Output directory for alignment results.")
@click.option("-c", "--container", default="None", show_default=True, help="Container: None, docker, or singularity.")
@click.option("-d", "--protein-db", "--protein-db-dir", help="Path to the protein FASTA database.")
def align_multiple_genomes(input_file, output_dir, container, protein_db):
    input_file = _prompt_path(input_file, "Genome manifest JSON path", exists=True, file_okay=True, dir_okay=False)
    output_dir = _prompt_path(output_dir, "Alignment output directory", file_okay=False, dir_okay=True, writable=True)
    protein_db = _prompt_path(protein_db, "Protein database FASTA path", exists=True, file_okay=True, dir_okay=False)
    config = _metagenomics_config(protein_db)
    genomes = _load_json(input_file)
    for genome, address in genomes.items():
        script = core.Metagenomics(config).align_genome_to_protein_db(
            address=address,
            outdir=output_dir,
            name=genome,
            container=container,
        )
        subprocess.run(script, shell=True)


@metagenomics.command(name="find-representative-genomes", help="Find representative genomes from a repseqs FASTA file.")
@click.option("-i", "--input-file", help="Path to the repseqs FASTA file.")
@click.option("-o", "--output-dir", help="Output directory.")
@click.option("--amplicon-to-genome-db", help="Directory containing the amplicon-to-genome database files.")
@click.option("-c", "--container", default="None", show_default=True, help="Container: None, docker, or singularity.")
@click.option("-s", "--similarity", default=0.97, show_default=True, type=float, help="Similarity cutoff for clustering.")
@click.option("-f", "--format", "output_format", default="csv", show_default=True, type=click.Choice(["csv", "json"]), help="Output format.")
def find_representative_genomes(input_file, output_dir, amplicon_to_genome_db, container, similarity, output_format):
    input_file = _prompt_path(input_file, "Repseqs FASTA path", exists=True, file_okay=True, dir_okay=False)
    output_dir = _prompt_path(output_dir, "Representative genomes output directory", file_okay=False, dir_okay=True, writable=True)
    amplicon_to_genome_db = _prompt_path(amplicon_to_genome_db, "Amplicon-to-genome database directory", exists=True, file_okay=False, dir_okay=True)
    config = _metagenomics_config(amplicon_to_genome_db=amplicon_to_genome_db)
    config.vsearch_similarity = similarity
    script = core.Metagenomics(config).align_to_gtdb(
        query_dir=input_file,
        output_dir=output_dir,
        container=container,
    )
    subprocess.run(script, shell=True)
    results = core.Metagenomics(config).get_genomes_from_gtdb_alignment(output_dir)
    _write_representative_genomes(results, output_dir, output_format)


@main.command(name="Documentations", help="Documentations for using ADToolbox.")
@click.option("-s", "--show", is_flag=True, help="Show the README documentation.")
def documentations(show):
    if not show:
        raise click.ClickException("Please use --show to display documentation.")
    with open(configs.Documentation().readme, "r") as f:
        console.print(markdown.Markdown(f.read()))


@main.group(name="ADM", help="Run and visualize ADToolbox ADM models.", no_args_is_help=True)
def adm_group():
    pass


def _adm_options(command):
    command = click.option("--report", help="Report output: dash or csv.")(command)
    command = click.option("--metagenome-report", help="JSON metagenome report for the model.")(command)
    command = click.option("--species", help="JSON species file.")(command)
    command = click.option("--reactions", help="JSON reactions file.")(command)
    command = click.option("--inlet-conditions", help="JSON inlet conditions file.")(command)
    command = click.option("--initial-conditions", help="JSON initial conditions file.")(command)
    command = click.option("--base-parameters", help="JSON base parameters file.")(command)
    command = click.option("--model-parameters", help="JSON model parameters file.")(command)
    command = click.option("--parameters-dir", help="Directory containing the ADM JSON parameter files.")(command)
    return command


@adm_group.command(name="adm1", help="Original ADM1 model.")
@_adm_options
def adm1(
    parameters_dir,
    model_parameters,
    base_parameters,
    initial_conditions,
    inlet_conditions,
    reactions,
    species,
    metagenome_report,
    report,
):
    paths = _resolve_model_paths(
        parameters_dir,
        "adm1",
        model_parameters=model_parameters,
        base_parameters=base_parameters,
        initial_conditions=initial_conditions,
        inlet_conditions=inlet_conditions,
        reactions=reactions,
        species=species,
    )
    if metagenome_report:
        _load_json(metagenome_report)

    model = adm.Model(
        model_parameters=_load_json(paths["model_parameters"]),
        base_parameters=_load_json(paths["base_parameters"]),
        initial_conditions=_load_json(paths["initial_conditions"]),
        inlet_conditions=_load_json(paths["inlet_conditions"]),
        feed=adm.DEFAULT_FEED,
        reactions=_load_json(paths["reactions"]),
        species=_load_json(paths["species"]),
        ode_system=adm.adm1_ode_sys,
        build_stoichiometric_matrix=adm.build_adm1_stoichiometric_matrix,
        name="ADM1",
        switch="DAE",
    )
    solution = model.solve_model(t_eval=np.linspace(0, 30, 10000))
    _report_adm_solution(model, solution, report)


@adm_group.command(name="e-adm", help="eADM model.")
@_adm_options
@click.option("--control-states", help="JSON file containing control states and their values.")
def e_adm(
    parameters_dir,
    model_parameters,
    base_parameters,
    initial_conditions,
    inlet_conditions,
    reactions,
    species,
    metagenome_report,
    report,
    control_states,
):
    paths = _resolve_model_paths(
        parameters_dir,
        "e_adm_2",
        model_parameters=model_parameters,
        base_parameters=base_parameters,
        initial_conditions=initial_conditions,
        inlet_conditions=inlet_conditions,
        reactions=reactions,
        species=species,
    )
    params = utils.load_multiple_json_files(paths)

    if metagenome_report:
        _load_json(metagenome_report)

    control_state = {"S_H_ion": 10 ** (-6.5)}
    if control_states:
        control_state.update(_load_json(control_states))

    model = adm.Model(
        model_parameters=params.model_parameters,
        base_parameters=params.base_parameters,
        initial_conditions=params.initial_conditions,
        inlet_conditions=params.inlet_conditions,
        feed=adm.DEFAULT_FEED,
        reactions=params.reactions,
        species=params.species,
        ode_system=adm.e_adm_2_ode_sys,
        build_stoichiometric_matrix=adm.build_e_adm_2_stoichiometric_matrix,
        control_state=control_state,
        name="Modified_ADM1",
        switch="DAE",
    )
    solution = model.solve_model(t_eval=np.linspace(0, 30, 10000), method="BDF")
    _report_adm_solution(model, solution, report)


if __name__ == "__main__":
    main()
