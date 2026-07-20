import os
import pathlib
import warnings

from adtoolbox import PKG_DATA

"""
Path builders and static configuration values for ADToolbox.

Configuration objects do not read or write global state. Pass the directory a
workflow should use, and the object derives its file paths from that directory.
"""


ADTOOLBOX_CONTAINERS = {
    "docker_x86": "parsaghadermazi/adtoolbox:latest",
    "docker_arm64": "parsaghadermazi/adtoolbox:arm64",
    "singularity_x86": "docker://parsaghadermazi/adtoolbox:latest",
    "singularity_arm64": "docker://parsaghadermazi/adtoolbox:arm64",
}

E_ADM_REMOTE = {
    "model_parameters": "https://raw.githubusercontent.com/ParsaGhadermazi/Database/main/ADToolbox/e_adm/e_adm_model_parameters.json",
    "base_parameters": "https://raw.githubusercontent.com/ParsaGhadermazi/Database/main/ADToolbox/e_adm/e_adm_base_parameters.json",
    "initial_conditions": "https://raw.githubusercontent.com/ParsaGhadermazi/Database/main/ADToolbox/e_adm/e_adm_initial_conditions.json",
    "inlet_conditions": "https://raw.githubusercontent.com/ParsaGhadermazi/Database/main/ADToolbox/e_adm/e_adm_inlet_conditions.json",
    "reactions": "https://raw.githubusercontent.com/ParsaGhadermazi/Database/main/ADToolbox/e_adm/e_adm_reactions.json",
    "species": "https://raw.githubusercontent.com/ParsaGhadermazi/Database/main/ADToolbox/e_adm/e_adm_species.json",
}

ADM1_REMOTE = {
    "model_parameters": "https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/adm1/adm1_model_parameters.json",
    "base_parameters": "https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/adm1/adm1_base_parameters.json",
    "initial_conditions": "https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/adm1/adm1_initial_conditions.json",
    "inlet_conditions": "https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/adm1/adm1_inlet_conditions.json",
    "reactions": "https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/adm1/adm1_reactions.json",
    "species": "https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/adm1/adm1_species.json",
}

EXTERNAL_LINKS = {
    "cazy_links": [
        "http://www.cazy.org/Glycoside-Hydrolases.html",
        "http://www.cazy.org/Polysaccharide-Lyases.html",
        "http://www.cazy.org/Carbohydrate-Esterases.html",
    ],
    "amplicon2genome": {
        "Version": "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/VERSION",
        "MD5SUM": "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/MD5SUM",
        "FILE_DESCRIPTIONS": "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/FILE_DESCRIPTIONS",
        "metadata_field_desc": "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/metadata_field_desc.tsv",
        "bac120_ssu": "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_all/ssu_all.fna.gz",
    },
    "seed_rxn_url": "https://github.com/modelSEED/modelSEEDDatabase/raw/master/Biochemistry/reactions.json",
    "seed_compound_url": "https://github.com/ModelSEED/ModelSEEDDatabase/raw/master/Biochemistry/compounds.json",
}

INTERNAL_LINKS = {
    "protein_db_url": "https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Protein_DB.fasta",
    "adtoolbox_rxn_db_url": "https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Reaction_Metadata.csv",
    "feed_db_url": "https://raw.githubusercontent.com/ParsaGhadermazi/Database/main/ADToolbox/feed_db.tsv",
    "metagenomics_studies": "https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Kbase/metagenomics_studies.tsv",
    "experimental_data_db": "https://raw.githubusercontent.com/ParsaGhadermazi/Database/main/ADToolbox/experimental_data_references.json",
}

E_ADM_MICROBIAL_GROUPS_MAPPING = {
    "Hydrolysis carbohydrates": "X_ch",
    "Hydrolysis proteins": "X_pr",
    "Hydrolysis lipids": "X_li",
    "Uptake of sugars": "X_su",
    "Uptake of amino acids": "X_aa",
    "Uptake of LCFA": "X_fa",
    "Uptake of acetate_et": "X_ac_et",
    "Uptake of acetate_lac": "X_ac_lac",
    "Uptake of propionate_et": "X_chain_et",
    "Uptake of propionate_lac": "X_chain_lac",
    "Uptake of butyrate_et": "X_chain_et",
    "Uptake of butyrate_lac": "X_chain_lac",
    "Uptake of valerate": "X_VFA_deg",
    "Uptake of caproate": "X_VFA_deg",
    "Methanogenessis from acetate and h2": "X_Me_ac",
    "Methanogenessis from CO2 and h2": "X_Me_CO2",
    "Uptake of ethanol": "X_et",
    "Uptake of lactate": "X_lac",
}


def _norm(path: str | os.PathLike | None, default: str | os.PathLike = ".") -> str:
    return os.path.abspath(os.path.expanduser(os.fspath(path or default)))


def _join(root: str, *parts: str) -> str:
    return os.path.join(root, *parts)


def adm_parameter_paths(parameters_dir: str | os.PathLike, prefix: str) -> dict[str, str]:
    parameters_dir = _norm(parameters_dir)
    return {
        "model_parameters": _join(parameters_dir, f"{prefix}_model_parameters.json"),
        "base_parameters": _join(parameters_dir, f"{prefix}_base_parameters.json"),
        "initial_conditions": _join(parameters_dir, f"{prefix}_initial_conditions.json"),
        "inlet_conditions": _join(parameters_dir, f"{prefix}_inlet_conditions.json"),
        "reactions": _join(parameters_dir, f"{prefix}_reactions.json"),
        "species": _join(parameters_dir, f"{prefix}_species.json"),
    }


class Database:
    "Configuration for core.Database functionality."

    def __init__(
        self,
        database_dir: str | os.PathLike = ".",
        *,
        compound_db: str | None = None,
        reaction_db: str | None = None,
        local_compound_db: str | None = None,
        local_reaction_db: str | None = None,
        csv_reaction_db: str | None = None,
        feed_db: str | None = None,
        amplicon_to_genome_db: str | None = None,
        adm_models: str | None = None,
        cazy_links: list[str] = EXTERNAL_LINKS["cazy_links"],
        amplicon_to_genome_urls: dict = EXTERNAL_LINKS["amplicon2genome"],
        adm_parameters_urls: dict = E_ADM_REMOTE,
        adm_parameters: dict | None = None,
        seed_rxn_url: str = EXTERNAL_LINKS["seed_rxn_url"],
        seed_compound_url: str = EXTERNAL_LINKS["seed_compound_url"],
        protein_db_url: str = INTERNAL_LINKS["protein_db_url"],
        adtoolbox_rxn_db_url: str = INTERNAL_LINKS["adtoolbox_rxn_db_url"],
        feed_db_url: str = INTERNAL_LINKS["feed_db_url"],
        adtoolbox_singularity: str = ADTOOLBOX_CONTAINERS["singularity_x86"],
        adtoolbox_docker: str = ADTOOLBOX_CONTAINERS["docker_x86"],
        protein_db: str | None = None,
        adm_microbial_groups_mapping: dict = E_ADM_MICROBIAL_GROUPS_MAPPING,
        metacyc_protein_db: str | None = None,
        studies_remote: dict = INTERNAL_LINKS,
        studies_local: dict | None = None,
        check_sanity: bool = False,
    ):
        self.database_dir = _norm(database_dir)
        adm_parameters_dir = _join(self.database_dir, "ADM_Parameters")
        studies_dir = _join(self.database_dir, "Studies")

        self.compound_db = _norm(compound_db or _join(self.database_dir, "compounds.json"))
        self.reaction_db = _norm(reaction_db or _join(self.database_dir, "reactions.json"))
        self.local_compound_db = _norm(local_compound_db or _join(self.database_dir, "Local_compounds.json"))
        self.local_reaction_db = _norm(local_reaction_db or _join(self.database_dir, "Local_reactions.json"))
        self.csv_reaction_db = _norm(csv_reaction_db or _join(self.database_dir, "Reaction_Metadata.csv"))
        self.feed_db = _norm(feed_db or _join(self.database_dir, "feed_db.tsv"))
        self.amplicon_to_genome_db = _norm(amplicon_to_genome_db or _join(self.database_dir, "Amplicon2GenomeDBs"))
        self.adm_models = _norm(adm_models or _join(adm_parameters_dir, "models.json"))
        self.cazy_links = cazy_links
        self.amplicon_to_genome_urls = amplicon_to_genome_urls
        self.adm_parameters_urls = adm_parameters_urls
        self.adm_parameters = adm_parameters or adm_parameter_paths(adm_parameters_dir, "e_adm")
        self.adm_parameters = {key: _norm(value) for key, value in self.adm_parameters.items()}
        self.seed_rxn_url = seed_rxn_url
        self.seed_compound_url = seed_compound_url
        self.protein_db_url = protein_db_url
        self.adtoolbox_rxn_db_url = adtoolbox_rxn_db_url
        self.feed_db_url = feed_db_url
        self.adtoolbox_singularity = adtoolbox_singularity
        self.adtoolbox_docker = adtoolbox_docker
        self.protein_db = _norm(protein_db or _join(self.database_dir, "Protein_DB.fasta"))
        self.adm_microbial_groups_mapping = adm_microbial_groups_mapping
        self.studies_remote = studies_remote
        default_studies_local = {
            "metagenomics_studies": _join(studies_dir, "metagenomics_studies.tsv"),
            "experimental_data_db": _join(studies_dir, "experimental_data_references.json"),
        }
        if studies_local:
            default_studies_local.update(studies_local)
        self.studies_local = {key: _norm(value) for key, value in default_studies_local.items()}
        self.metacyc_protein_db = _norm(metacyc_protein_db or _join(self.database_dir, "metacyc_protein_db.fasta"))
        self.protein_db_mmseqs = pathlib.Path(self.protein_db).parent.joinpath("protein_db_mmseqs")

        # Compatibility aliases for older examples/tests that treated Database as
        # the active ADM parameter config.
        self.model_parameters = self.adm_parameters["model_parameters"]
        self.base_parameters = self.adm_parameters["base_parameters"]
        self.initial_conditions = self.adm_parameters["initial_conditions"]
        self.inlet_conditions = self.adm_parameters["inlet_conditions"]
        self.reactions = self.adm_parameters["reactions"]
        self.species = self.adm_parameters["species"]

        if check_sanity:
            self.check_adm_parameters()

    def check_adm_parameters(self):
        branches = all(
            pathlib.Path(path).parent == pathlib.Path(self.adm_parameters["model_parameters"]).parent
            for path in self.adm_parameters.values()
        )
        if not branches:
            warnings.warn("The ADM parameters are not in the same directory!")


class Metagenomics:
    "Configuration for core.Metagenomics functionality."

    gtdb_dir = "*ssu*.fna"

    def __init__(
        self,
        metagenomics_dir: str | os.PathLike = ".",
        *,
        database_dir: str | os.PathLike | None = None,
        database: Database | None = None,
        amplicon2genome_k=10,
        vsearch_similarity=0.97,
        genomes_base_dir: str | None = None,
        align_to_gtdb_outputs_dir: str | None = None,
        amplicon2genome_db: str | None = None,
        genome_alignment_script: str | None = None,
        vsearch_threads: int = 4,
        rsync_download_dir: str | None = None,
        adtoolbox_singularity: str = ADTOOLBOX_CONTAINERS["singularity_x86"],
        adtoolbox_docker: str = ADTOOLBOX_CONTAINERS["docker_x86"],
        genome_alignment_output: str | None = None,
        csv_reaction_db: str | None = None,
        sra: str | None = None,
        bit_score=40,
        e_value=10**-5,
        protein_db: str | None = None,
        protein_db_mmseqs: str | None = None,
        adm_mapping=E_ADM_MICROBIAL_GROUPS_MAPPING,
    ):
        self.metagenomics_dir = _norm(metagenomics_dir)
        database = database or Database(database_dir=database_dir or self.metagenomics_dir)

        self.k = amplicon2genome_k
        self.vsearch_similarity = vsearch_similarity
        self.genomes_base_dir = genomes_base_dir or _join(self.metagenomics_dir, "Genomes")
        self.align_to_gtdb_outputs_dir = align_to_gtdb_outputs_dir or self.genomes_base_dir
        self.amplicon2genome_db = amplicon2genome_db or database.amplicon_to_genome_db
        self.protein_db = protein_db or database.protein_db
        self.protein_db_mmseqs = protein_db_mmseqs or database.protein_db_mmseqs
        self.seed_rxn_db = database.reaction_db
        self.genome_alignment_output = genome_alignment_output or _join(self.metagenomics_dir, "Outputs")
        self.bit_score = bit_score
        self.e_value = e_value
        self.vsearch_threads = vsearch_threads
        self.csv_reaction_db = csv_reaction_db or database.csv_reaction_db
        self.sra = sra or _join(self.metagenomics_dir, "SRA")
        self.gtdb_dir_fasta = None
        matches = list(pathlib.Path(self.amplicon2genome_db).rglob(Metagenomics.gtdb_dir))
        if matches:
            self.gtdb_dir_fasta = str(matches[0])
        self.genome_alignment_script = genome_alignment_script or _join(self.metagenomics_dir, "genome_alignment_script.sh")
        self.adtoolbox_singularity = adtoolbox_singularity
        self.adtoolbox_docker = adtoolbox_docker
        self.rsync_download_dir = rsync_download_dir or _join(self.genomes_base_dir, "rsync_download.sh")
        self.adm_mapping = adm_mapping


class Annotation:
    def __init__(self, annotation_dir: str | os.PathLike = ".", *, metacyc_protein_db: str | None = None):
        self.annotation_dir = _norm(annotation_dir)
        self.metacyc_protein_db = metacyc_protein_db or _join(self.annotation_dir, "metacyc_protein_db.fasta")


class Documentation:
    def __init__(self, documentation_dir: str | os.PathLike | None = None, *, readme: str | None = None):
        self.documentation_dir = _norm(documentation_dir or PKG_DATA)
        self.readme = readme or _join(self.documentation_dir, "README.md")


class Studies:
    def __init__(
        self,
        studies_dir: str | os.PathLike = ".",
        *,
        metagenomics_studies: str | None = None,
        experimental_data_db: str | None = None,
    ):
        self.studies_dir = _norm(studies_dir)
        self.metagenomics_studies = metagenomics_studies or _join(self.studies_dir, "metagenomics_studies.tsv")
        self.experimental_data_db = experimental_data_db or _join(self.studies_dir, "experimental_data_references.json")
        self.experimental_data_references = self.experimental_data_db


class Utils:
    "Configuration for utility helpers."

    def __init__(
        self,
        utils_dir: str | os.PathLike = ".",
        *,
        slurm_template: str = os.path.join(PKG_DATA, "slurm_template.txt"),
        slurm_executer: str = "",
        slurm_wall_time: str = "24:00:00",
        slurm_job_name: str = "ADToolbox",
        slurm_outlog: str = "ADToolbox.log",
        slurm_cpus: str = "12",
        slurm_memory: str = "100G",
        slurm_save_dir: str | None = None,
        adtoolbox_singularity: str = ADTOOLBOX_CONTAINERS["singularity_x86"],
        adtoolbox_docker: str = ADTOOLBOX_CONTAINERS["docker_x86"],
    ) -> None:
        self.utils_dir = _norm(utils_dir)
        self.slurm_template = slurm_template
        self.slurm_executer = slurm_executer
        self.slurm_wall_time = slurm_wall_time
        self.slurm_job_name = slurm_job_name
        self.slurm_outlog = slurm_outlog
        self.slurm_cpus = slurm_cpus
        self.slurm_save_dir = slurm_save_dir or self.utils_dir
        self.slurm_memory = slurm_memory
        self.adtoolbox_singularity = adtoolbox_singularity
        self.adtoolbox_docker = adtoolbox_docker


_DEFAULT_DATABASE = Database()
RXN_DB = _DEFAULT_DATABASE.csv_reaction_db
Seed_RXN_DB = _DEFAULT_DATABASE.reaction_db
Seed_COMPOUNDS_DB = _DEFAULT_DATABASE.compound_db
E_ADM_LOCAL = _DEFAULT_DATABASE.adm_parameters
ADM1_LOCAL = adm_parameter_paths(os.path.join(_DEFAULT_DATABASE.database_dir, "ADM_Parameters"), "adm1")
STUDIES_LOCAL = _DEFAULT_DATABASE.studies_local
