import os
from __init__ import Main_Dir,PKG_DATA
import pathlib
import rich
import json
import warnings

"""
This module contains all the paths to the files and directories used in the program.

"""


RXN_DB = os.path.join(Main_Dir, "Database", "Reaction_Metadata.csv")

Seed_RXN_DB = os.path.join(Main_Dir, "Database", "reactions.json")
Seed_COMPOUNDS_DB = os.path.join(Main_Dir, "Database", "compounds.json")

ADTOOLBOX_CONTAINERS={
	'docker_x86':"parsaghadermazi/adtoolbox:latest",
	'docker_arm64':"parsaghadermazi/adtoolbox:arm64",
	'singularity_x86':"docker://parsaghadermazi/adtoolbox:x86",
	'singularity_arm64':"docker://parsaghadermazi/adtoolbox:arm64"}

E_ADM_2_REMOTE={
	"model_parameters":"https://raw.githubusercontent.com/ParsaGhadermazi/Database/main/ADToolbox/e_adm_2/e_adm_2_model_parameters.json",
	"base_parameters":"https://raw.githubusercontent.com/ParsaGhadermazi/Database/main/ADToolbox/e_adm_2/e_adm_2_base_parameters.json",
	"initial_conditions":"https://raw.githubusercontent.com/ParsaGhadermazi/Database/main/ADToolbox/e_adm_2/e_adm_2_initial_conditions.json",
	"inlet_conditions":"https://raw.githubusercontent.com/ParsaGhadermazi/Database/main/ADToolbox/e_adm_2/e_adm_2_inlet_conditions.json",
	"reactions":"https://raw.githubusercontent.com/ParsaGhadermazi/Database/main/ADToolbox/e_adm_2/e_adm_2_reactions.json",
	"species":"https://raw.githubusercontent.com/ParsaGhadermazi/Database/main/ADToolbox/e_adm_2/e_adm_2_species.json"
			}
E_ADM_2_LOCAL={
	"model_parameters":os.path.join(Main_Dir, "Database","ADM_Parameters","e_adm_2_model_parameters.json"),
	"base_parameters":os.path.join(Main_Dir, "Database","ADM_Parameters","e_adm_2_base_parameters.json"),
	"initial_conditions":os.path.join(Main_Dir, "Database","ADM_Parameters","e_adm_2_initial_conditions.json"),
	"inlet_conditions":os.path.join(Main_Dir, "Database","ADM_Parameters","e_adm_2_inlet_conditions.json"),
	"reactions":os.path.join(Main_Dir, "Database","ADM_Parameters","e_adm_2_reactions.json"),
	"species":os.path.join(Main_Dir, "Database","ADM_Parameters","e_adm_2_species.json"),	
}

E_ADM_LOCAL={
	"model_parameters":os.path.join(Main_Dir, "Database","ADM_Parameters","e_adm_model_parameters.json"),
	"base_parameters":os.path.join(Main_Dir, "Database","ADM_Parameters","e_adm_base_parameters.json"),
	"initial_conditions":os.path.join(Main_Dir, "Database","ADM_Parameters","e_adm_initial_conditions.json"),
	"inlet_conditions":os.path.join(Main_Dir, "Database","ADM_Parameters","e_adm_inlet_conditions.json"),
	"reactions":os.path.join(Main_Dir, "Database","ADM_Parameters","e_adm_reactions.json"),
	"species":os.path.join(Main_Dir, "Database","ADM_Parameters","e_adm_species.json"),
			}

E_ADM_REMOTE={
	"model_parameters":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/e_adm/e_adm_model_parameters.json",
	"base_parameters":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/e_adm/e_adm_base_parameters.json",
	"initial_conditions":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/e_adm/e_adm_initial_conditions.json",
	"inlet_conditions":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/e_adm/e_adm_inlet_conditions.json",
	"reactions":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/e_adm/e_adm_reactions.json",
	"species":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/e_adm/e_adm_species.json"
			}

ADM1_LOCAL={
	"model_parameters":os.path.join(Main_Dir, "Database","ADM_Parameters","adm1_model_parameters.json"),
	"base_parameters":os.path.join(Main_Dir, "Database","ADM_Parameters","adm1_base_parameters.json"),
	"initial_conditions":os.path.join(Main_Dir, "Database","ADM_Parameters","adm1_initial_conditions.json"),
	"inlet_conditions":os.path.join(Main_Dir, "Database","ADM_Parameters","adm1_inlet_conditions.json"),
	"reactions":os.path.join(Main_Dir, "Database","ADM_Parameters","adm1_reactions.json"),
	"species":os.path.join(Main_Dir, "Database","ADM_Parameters","adm1_species.json"),
			}

ADM1_REMOTE={
	"model_parameters":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/adm1/adm1_model_parameters.json",
	"base_parameters":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/adm1/adm1_base_parameters.json",
	"initial_conditions":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/adm1/adm1_initial_conditions.json",
	"inlet_conditions":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/adm1/adm1_inlet_conditions.json",
	"reactions":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/adm1/adm1_reactions.json",
	"species":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/adm1/adm1_species.json"
		}
EXTERNAL_LINKS={
	"cazy_links":["http://www.cazy.org/Glycoside-Hydrolases.html",
                  "http://www.cazy.org/Polysaccharide-Lyases.html",
                  "http://www.cazy.org/Carbohydrate-Esterases.html"
                  		],
	"amplicon2genome":{'Version': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/VERSION',
               'MD5SUM': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/MD5SUM',
               'FILE_DESCRIPTIONS': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/FILE_DESCRIPTIONS',
               'metadata_field_desc': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/metadata_field_desc.tsv',
               'bac120_ssu': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_all/ssu_all.tar.gz'
               },
	
	"seed_rxn_url":"https://github.com/modelSEED/modelSEEDDatabase/raw/master/Biochemistry/reactions.json",
 	"seed_compound_url":"https://github.com/ModelSEED/ModelSEEDDatabase/raw/master/Biochemistry/compounds.json",
	
}

INTERNAL_LINKS={
	"protein_db_url":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Protein_DB.fasta",
	"adtoolbox_rxn_db_url":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Reaction_Metadata.csv",
	"feed_db_url":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Feed_DB.json",
	"qiime_classifier_db_url":"https://data.qiime2.org/2022.11/common/silva-138-99-515-806-nb-classifier.qza",
	"metagenomics_studies":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Kbase/metagenomics_studies.tsv",
    "exmpermental_data_db":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Kbase/experimental_data_references.json"
}

STUDIES_REMOTE={
		"metagenomics_studies":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Kbase/metagenomics_studies.tsv",
    "exmpermental_data_db":"https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Kbase/experimental_data_references.json"
}

STUDIES_LOCAL={
		"metagenomics_studies":os.path.join(Main_Dir,"Database","Studies","metagenomics_studies.tsv"),
	"exmpermental_data_db":os.path.join(Main_Dir,"Database","Studies","experimental_data_references.json")
}
E_ADM_MICROBIAL_GROUPS_MAPPING={
                            "Hydrolysis carbohydrates":"X_ch",
                            "Hydrolysis proteins":"X_pr",
                            "Hydrolysis lipids":"X_li",
                            "Uptake of sugars":"X_su",
                            "Uptake of amino acids":"X_aa",
                            "Uptake of LCFA":"X_fa",
                            "Uptake of acetate_et":"X_ac_et",
                            "Uptake of acetate_lac":"X_ac_lac",
                            "Uptake of propionate_et":"X_chain_et",
                            "Uptake of propionate_lac":"X_chain_lac",
                            "Uptake of butyrate_et":"X_chain_et",
                            "Uptake of butyrate_lac":"X_chain_lac",
                            "Uptake of valerate":"X_VFA_deg",
                            "Uptake of caproate":"X_VFA_deg",
                            "Methanogenessis from acetate and h2":"X_Me_ac",
                            "Methanogenessis from CO2 and h2":"X_Me_CO2",
                            "Uptake of ethanol":"X_et",
                            "Uptake of lactate":"X_lac",
                            }


class Database:
	"An instance of this class will hold all the configuration information for core.Database functionalities."
	
	def __init__(self,
		compound_db:str=Seed_COMPOUNDS_DB,
		reaction_db:str=Seed_RXN_DB,
		local_compound_db:str=os.path.join(Main_Dir, "Database", 'Local_compounds.json'),
		local_reaction_db:str=os.path.join(Main_Dir, "Database", 'Local_reactions.json'),
		csv_reaction_db:str=os.path.join(Main_Dir, "Database", 'Reaction_Metadata.csv'),
		feed_db=os.path.join(Main_Dir, "Database", 'feed_db.tsv'),
		amplicon_to_genome_db=os.path.join(Main_Dir,'Database','Amplicon2GenomeDBs'),
		cazy_links:str=EXTERNAL_LINKS["cazy_links"],
		amplicon_to_genome_urls:dict=EXTERNAL_LINKS["amplicon2genome"],
		adm_parameters_urls:dict=E_ADM_REMOTE,
		adm_parameters:dict=E_ADM_LOCAL,
		seed_rxn_url:str =EXTERNAL_LINKS["seed_rxn_url"],
		seed_compound_url:str =EXTERNAL_LINKS["seed_compound_url"],
		protein_db_url:str =INTERNAL_LINKS["protein_db_url"],
		adtoolbox_rxn_db_url:str =INTERNAL_LINKS["adtoolbox_rxn_db_url"],
		feed_db_url:str =INTERNAL_LINKS["feed_db_url"],
		qiime_classifier_db:str=os.path.join(Main_Dir, "Database","qiime2_classifier_db" ,'qiime2_classifier_db.qza'),
		qiime_classifier_db_url:str=INTERNAL_LINKS["qiime_classifier_db_url"],
  		adtoolbox_singularity=ADTOOLBOX_CONTAINERS["singularity_x86"],
		adtoolbox_docker=ADTOOLBOX_CONTAINERS["docker_x86"],
    	protein_db=os.path.join(Main_Dir, "Database", 'Protein_DB.fasta'),
		adm_microbial_groups_mapping=E_ADM_MICROBIAL_GROUPS_MAPPING,
        metagenomics_studies_db=os.path.join(Main_Dir,"Database","Studies","metagenomics_studies.tsv"),
        experimental_data_db=os.path.join(Main_Dir,"Database","Studies","experimental_data_references.json"),
        studies_remote=STUDIES_REMOTE,
        studies_local=STUDIES_LOCAL,
        check_sanity:bool=False
		):
		self.compound_db = compound_db
		self.reaction_db = reaction_db
		self.local_compound_db = local_compound_db
		self.local_reaction_db = local_reaction_db
		self.csv_reaction_db = csv_reaction_db
		self.feed_db = feed_db
		self.amplicon_to_genome_db = amplicon_to_genome_db
		self.cazy_links = cazy_links
		self.amplicon_to_genome_urls = amplicon_to_genome_urls
		self.adm_parameters_urls = adm_parameters_urls
		self.adm_parameters = adm_parameters
		self.seed_rxn_url = seed_rxn_url
		self.seed_compound_url = seed_compound_url
		self.protein_db_url = protein_db_url
		self.adtoolbox_rxn_db_url = adtoolbox_rxn_db_url
		self.feed_db_url = feed_db_url
		self.qiime_classifier_db = qiime_classifier_db
		self.qiime_classifier_db_url = qiime_classifier_db_url
		self.adtoolbox_singularity=adtoolbox_singularity
		self.adtoolbox_docker=adtoolbox_docker
		self.protein_db=protein_db
		self.adm_microbial_groups_mapping=adm_microbial_groups_mapping
		self.metagenomics_studies_db=metagenomics_studies_db
		self.experimental_data_db=experimental_data_db
		self.studies_remote=studies_remote
		self.studies_local=studies_local
		self.protein_db_mmseqs=pathlib.Path(protein_db).parent.joinpath("protein_db_mmseqs")
		if check_sanity:
			self.check_adm_parameters()
   
	def check_adm_parameters(self):
		branches=all([pathlib.Path(x).parent==pathlib.Path(self.adm_parameters["model_parameters"]).parent for x in self.adm_parameters.values()])
		if not branches:
			warnings.warn(f"The ADM parameters are not in the same directory!")

class Metagenomics:
	"""	
	An instance of this class will hold all the configuration information for core.Metagenomics functionalities.
	"""
	### Here we have some class variables that are used in the class
	gtdb_dir="ssu_all_*.fna"
	def __init__(self, 
            amplicon2genome_k=10,
            vsearch_similarity=0.97,
            genomes_base_dir=os.path.join(Main_Dir,"Genomes"),
            align_to_gtdb_outputs_dir=os.path.join(Main_Dir,"Genomes"),
            amplicon2genome_db=Database().amplicon_to_genome_db,
            qiime_outputs_dir=os.path.join(Main_Dir,'Metagenomics_Data','QIIME_Outputs'),
			genome_alignment_script=os.path.join(Main_Dir,"Metagenomics_Data","QIIME_Outputs","genome_alignment_script.sh"),
			vsearch_threads:int=4,
			rsync_download_dir=os.path.join(Main_Dir,"Genomes","rsync_download.sh"),
			adtoolbox_singularity=ADTOOLBOX_CONTAINERS["singularity_x86"],
			adtoolbox_docker=ADTOOLBOX_CONTAINERS["docker_x86"],
            genome_alignment_output=os.path.join(Main_Dir,"Outputs"),
            csv_reaction_db=Database().csv_reaction_db,
            sra=os.path.join(Main_Dir,"Metagenomics_Analysis","SRA"),
            bit_score=40,
            e_value=10**-5,
            qiime2_docker_image="quay.io/qiime2/core:2022.2",
            qiime2_singularity_image="docker://quay.io/qiime2/core:2022.2",
            qiime2_paired_end_bash_str=os.path.join(PKG_DATA,"qiime_template_paired.txt"),
            qiime2_single_end_bash_str=os.path.join(PKG_DATA,"qiime_template_single.txt"),
			qiime_classifier_db=Database().qiime_classifier_db,
			adm_mapping=Database().adm_microbial_groups_mapping,
             ):
		self.k = amplicon2genome_k
		self.vsearch_similarity = vsearch_similarity
		self.align_to_gtdb_outputs_dir = align_to_gtdb_outputs_dir
		self.amplicon2genome_db = amplicon2genome_db
		self.qiime_outputs_dir = qiime_outputs_dir
		self.protein_db=Database().protein_db
		self.protein_db_mmseqs=Database().protein_db_mmseqs
		self.seed_rxn_db=Seed_RXN_DB
		self.genome_alignment_output = genome_alignment_output
		self.bit_score = bit_score
		self.e_value = e_value
		self.vsearch_threads=vsearch_threads
		self.csv_reaction_db=csv_reaction_db
		self.sra=sra
		self.qiime2_singularity_image=qiime2_singularity_image
		self.qiime2_docker_image=qiime2_docker_image
		self.qiime2_paired_end_bash_str=qiime2_paired_end_bash_str
		self.qiime2_single_end_bash_str=qiime2_single_end_bash_str 
		self.qiime_classifier_db=qiime_classifier_db
		if list(pathlib.Path(self.amplicon2genome_db).rglob(Metagenomics.gtdb_dir)):
			self.gtdb_dir_fasta=str(list(pathlib.Path(self.amplicon2genome_db).rglob(Metagenomics.gtdb_dir))[0])
		else:
			self.gtdb_dir_fasta=None
		self.genome_alignment_script=genome_alignment_script	
		self.adtoolbox_singularity=adtoolbox_singularity
		self.adtoolbox_docker=adtoolbox_docker
		self.rsync_download_dir=rsync_download_dir
		self.genomes_base_dir=genomes_base_dir
		self.adm_mapping=adm_mapping

class Documentation:
    def __init__(self,
                 readme=os.path.join(PKG_DATA,"README.md")):
        self.readme = readme

class Utils:
	"""
	An instance of this class will hold all the configuration information for utils module functionalities."""
	def __init__(self,
	slurm_template:str=os.path.join(PKG_DATA,"slurm_template.txt"),
	docker_template_qiime:str=None,
	singularity_template_qiime:str=None,
	slurm_executer:str='',
	slurm_wall_time:str='24:00:00',
	slurm_job_name:str='ADToolbox',
	slurm_outlog:str='ADToolbox.log',
    slurm_cpus:str="12",
	slurm_memory:str="100G",
	slurm_save_dir:str=os.getcwd(),
	adtoolbox_singularity:str=ADTOOLBOX_CONTAINERS["singularity_x86"],
	adtoolbox_docker:str=ADTOOLBOX_CONTAINERS["docker_x86"]
	) -> None:
		self.slurm_template = slurm_template
		self.docker_template_qiime = docker_template_qiime
		self.singularity_template_qiime = singularity_template_qiime
		self.slurm_executer = slurm_executer
		self.slurm_wall_time = slurm_wall_time
		self.slurm_job_name = slurm_job_name
		self.slurm_outlog=slurm_outlog
		self.slurm_cpus = slurm_cpus
		self.slurm_save_dir = slurm_save_dir
		self.slurm_memory = slurm_memory
		self.adtoolbox_singularity=adtoolbox_singularity
		self.adtoolbox_docker=adtoolbox_docker
  

def get_base_dir():
	"This function returns the current base directory of the program."
	return Main_Dir

def set_base_dir(path:str):
	"""This function changes the base directory of the program"""
	global Main_Dir
	ans=input("This will change the base directory of the program. Are you sure you want to continue? (y/n)")
	if ans == "y":
		with open(os.path.join(PKG_DATA,"ADToolbox_Configs.json"),'r') as f:
			conf = json.load(f)
			try:
				os.makedirs(path,exist_ok=True)
			except Exception:
				rich.print("[red]Base directory not changed")
				return
			else:
				
				Main_Dir=path
				conf["Base_Dir"]=path
				with open(os.path.join(PKG_DATA,"ADToolbox_Configs.json"),'w') as f:
					json.dump(conf,f)
				
				rich.print("[green]Base directory changed")	
	else:
		rich.print("[red]Base directory not changed")





if __name__=="__main__":
    t=Database()
    t.adm_parameters["model_parameters"]=os.path.join(Main_Dir, "Database","somewhere_else","e_adm_model_parameters.json")
    t.check_adm_parameters()