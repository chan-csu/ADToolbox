import os
from __init__ import Main_Dir,PKG_DATA
import pathlib

"""
This module contains all the paths to the files and directories used in the program.

"""


RXN_DB = os.path.join(Main_Dir, "Database", "Reaction_Metadata.csv")

Seed_RXN_DB = os.path.join(Main_Dir, "Database", "reactions.json")

def get_base_dir():
	return Main_Dir

def set_base_dir(path:str):
	ans=input("This will change the base directory of the program. Are you sure you want to continue? (y/n)")
	if ans == "y":
		with open(os.path.join(PKG_DATA,"ADToolbox_Configs.json"),'w') as f:
			f.write(path)
	else:
		print("Base directory not changed")
	

class Alignment:
    """
    A class for aligner configurations

    """

    def __init__(self, aligner_name="mmseqs2", e_value=10**-5, bit_score=40):
        self.aligner_name = aligner_name
        self.e_value = e_value
        self.bit_score = bit_score
    
class Reaction_Toolkit:

    def __init__(self, 
        compound_db=os.path.join(Main_Dir, "Database", 'Local_compounds.json'),
        reaction_db=os.path.join(Main_Dir, "Database", 'Local_reactions.json')):
        
        self.compound_db = compound_db
        self.reaction_db = reaction_db

class Studies:
    def __init__(self,
                 base_dir=os.path.join(Main_Dir,"Database","Studies"),
                 metagenomics_studies=os.path.join(Main_Dir,"Database","Studies","metagenomics_studies.tsv"),
                 experimental_data_references=os.path.join(Main_Dir,"Database","Studies","experimental_data_references.tsv"),
                urls={'metagenomics_studies': 'https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Kbase/metagenomics_studies.tsv',
                        'exmpermental_data_references':'https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Kbase/experimental_data_references.tsv'
                        }):
        self.metagenomics_studies = metagenomics_studies
        self.experimental_data_references = experimental_data_references
        self.base_dir = base_dir
        self.urls = urls

class Database:
	"A class for database configurations"
	
	def __init__(self,
		compound_db=os.path.join(
		     Main_Dir, "Database", 'compounds.json'),
		reaction_db=os.path.join(
		     Main_Dir, "Database", 'reactions.json'),
		local_compound_db=os.path.join(
		        Main_Dir, "Database", 'Local_compounds.json'),
		local_reaction_db=os.path.join(
		        Main_Dir, "Database", 'Local_reactions.json'),
		base_dir=os.path.join(Main_Dir,"Database"),
		csv_reaction_db=os.path.join(Main_Dir, "Database", 'Reaction_Metadata.csv'),
		feed_db=os.path.join(Main_Dir, "Database", 'feed_db.tsv'),
		protein_db=os.path.join(Main_Dir, "Database", 'Protein_DB.fasta'),
		kbase_db=Studies(),
		cazy_links:str=["http://www.cazy.org/Glycoside-Hydrolases.html",
                  "http://www.cazy.org/Polysaccharide-Lyases.html",
                  "http://www.cazy.org/Carbohydrate-Esterases.html"
                  ],
		adm_parameters_urls:dict=dict(
		adm1_model_parameters="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/ADM1/ADM1_Model_Parameters.json",
        adm1_base_parameters="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/ADM1/ADM1_Base_Parameters.json",
        adm1_initial_conditions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/ADM1/ADM1_Initial_Conditions.json",
        adm1_inlet_conditions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/ADM1/ADM1_Inlet_Conditions.json",
        adm1_reactions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/ADM1/ADM1_Reactions.json",
        adm1_species="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/ADM1/ADM1_Species.json",
		model_parameters="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Model_Parameters.json",
        base_parameters="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Base_Parameters.json",
        initial_conditions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Initial_Conditions.json",
        inlet_conditions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Inlet_Conditions.json",
        reactions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Reactions.json",
        species="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Species.json",
		),
		adm_parameters_base_dir:str=os.path.join(Main_Dir, "Database","ADM_Parameters"),
		seed_rxn_url:str ="https://github.com/modelSEED/modelSEEDDatabase/raw/master/Biochemistry/reactions.json",
		protein_db_url:str ="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/protein_db.fasta",
		adtoolbox_rxn_db_url:str ="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Reaction_Metadata.csv",
		feed_db_url:str ="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Feed_DB.json",
		qiime_classifier_db:str=os.path.join(Main_Dir, "Database","qiime2_classifier_db" ,'qiime2_classifier_db.qza'),
		qiime_classifier_db_url:str= "https://data.qiime2.org/2022.11/common/silva-138-99-515-806-nb-classifier.qza"
		):

		self.compound_db = compound_db
		self.reaction_db = reaction_db
		self.local_compound_db = local_compound_db
		self.local_reaction_db = local_reaction_db
		self.base_dir = base_dir
		self.csv_reaction_db = csv_reaction_db
		self.feed_db = feed_db
		self.protein_db = protein_db
		self.kbase_db = kbase_db
		self.seed_rxn_url = seed_rxn_url
		self.protein_db_url = protein_db_url
		self.adtoolbox_rxn_db_url = adtoolbox_rxn_db_url
		self.feed_db_url = feed_db_url
		self.qiime_classifier_db = qiime_classifier_db
		self.qiime_classifier_db_url = qiime_classifier_db_url
		self.cazy_links = cazy_links
		self.adm_parameters_urls = adm_parameters_urls
		self.adm_parameters_base_dir = adm_parameters_base_dir

class Metagenomics:
	"""	
	A class for Amplicon2Genome Configs
	"""
	### Here we have some class variables that are used in the class
	gtdb_dir="bac120_ssu_reps_r207.fna"
	def __init__(self, 
            amplicon2genome_k=10,
            amplicon2genome_similarity=0.97,
            amplicon2genome_outputs_dir=os.path.join(Main_Dir,"Genomes"),
            amplicon2genome_db=os.path.join(Main_Dir,'Database','Amplicon2GenomeDBs'),
			amplicon2genome_top_repseq_dir=os.path.join(Main_Dir,"Metagenomics_Data","QIIME_Outputs","top_repseqs.fasta"),
            qiime_outputs_dir=os.path.join(Main_Dir,'Metagenomics_Data','QIIME_Outputs'),
			genome_alignment_script=os.path.join(Main_Dir,"Metagenomics_Data","QIIME_Outputs","genome_alignment_script.sh"),
            genomes_json_info=os.path.join(Main_Dir,"Genomes","Amplicon2Genome_OutInfo.json"),
			amplicon2genome_docker="vsearch",
			vsearch_threads:int=4,
			vsearch_script_dir=os.path.join(Main_Dir,"Metagenomics_Data","QIIME_Outputs","vsearch_scripts.sh"),
            feature_table_dir=os.path.join(Main_Dir,"Metagenomics_Data","QIIME_Outputs","feature-table.tsv"),
            rep_seq_fasta=os.path.join(Main_Dir,"Metagenomics_Data","QIIME_Outputs","dna-sequences.fasta"),
            taxonomy_table_dir=os.path.join(Main_Dir,"Metagenomics_Data","QIIME_Outputs","taxonomy.tsv"),
            genome_alignment_output=os.path.join(Main_Dir,"Outputs"),
			feature_to_taxa=os.path.join(Main_Dir,"Metagenomics_Data","QIIME_Outputs","feature_to_taxa.json"),
            genome_alignment_output_json=os.path.join(Main_Dir,"Outputs","Alignment_Info.json"),
            genome_adm_map_json=os.path.join(Main_Dir,"Outputs","ADM_From_Alignment_JSON_Output.json"),
            csv_reaction_db=Database().csv_reaction_db,
			genome_relative_abundances=os.path.join(Main_Dir,"Metagenomics_Data","QIIME_Outputs","genome_relative_abundances.json"),
            sra=os.path.join(Main_Dir,"Metagenomics_Analysis","SRA"),
            bit_score=40,
            e_value=10**-5,
            qiime2_docker_image="quay.io/qiime2/core:2022.2",
            qiime2_singularity_image="docker://quay.io/qiime2/core:2022.2",
            qiime2_paired_end_bash_str=os.path.join(PKG_DATA,"qiime_template_paired.txt"),
            qiime2_single_end_bash_str=os.path.join(PKG_DATA,"qiime_template_single.txt"),
			qiime_classifier_db=Database().qiime_classifier_db,
			cod_output_json=os.path.join(Main_Dir,"Metagenomics_Data","QIIME_Outputs","cod_output.json"),

             ):
		self.k = amplicon2genome_k
		self.amplicon2genome_similarity = amplicon2genome_similarity
		self.amplicon2genome_outputs_dir = amplicon2genome_outputs_dir
		self.amplicon2genome_db = amplicon2genome_db
		self.qiime_outputs_dir = qiime_outputs_dir
		self.genomes_json_info = genomes_json_info
		self.amplicon2genome_top_repseq_dir = amplicon2genome_top_repseq_dir
		self.feature_table_dir = feature_table_dir
		self.rep_seq_fasta = rep_seq_fasta
		self.taxonomy_table_dir = taxonomy_table_dir
		self.protein_db=Database().protein_db
		self.seed_rxn_db=Seed_RXN_DB
		self.genome_alignment_output = genome_alignment_output
		self.genome_alignment_output_json=genome_alignment_output_json
		self.bit_score = bit_score
		self.e_value = e_value
		self.vsearch_threads=vsearch_threads
		self.genome_adm_map_json=genome_adm_map_json
		self.feature_to_taxa=feature_to_taxa
		self.csv_reaction_db=csv_reaction_db
		self.sra=sra
		self.amplicon2genome_docker=amplicon2genome_docker
		self.qiime2_singularity_image=qiime2_singularity_image
		self.qiime2_docker_image=qiime2_docker_image
		self.qiime2_paired_end_bash_str=qiime2_paired_end_bash_str
		self.qiime2_single_end_bash_str=qiime2_single_end_bash_str 
		self.qiime_classifier_db=qiime_classifier_db
		self.gtdb_dir_fasta=os.path.join(self.amplicon2genome_db,Metagenomics.gtdb_dir)
		self.vsearch_script_dir=vsearch_script_dir
		self.genome_alignment_script=genome_alignment_script	
		self.genome_relative_abundances=genome_relative_abundances
		self.cod_output_json=cod_output_json



	def genome_save_dir(self, accession_id: str):
		""" sets the working directory for the accession id """
		return os.path.join(self.amplicon2genome_outputs_dir, accession_id)

	def sra_work_dir(self, sra_project_id: str):
		""" sets the working directory for the SRA project id """
	
		return os.path.join(self.sra, sra_project_id)
	
	def qiime2_work_dir(self, sra_project_id: str):
		""" sets the working directory for the SRA project id """
		return os.path.join(self.sra_work_dir(sra_project_id), "qiime2")






class Documentation:
    def __init__(self,
                 readme=os.path.join(PKG_DATA,"README.md")):
        self.readme = readme



class Utils:

	def __init__(self,
	slurm_template:str=os.path.join(PKG_DATA,"slurm_template.txt"),
	docker_template_qiime:str=None,
	singularity_template_qiime:str=None,
	slurm_executer:str='',
	slurm_wall_time:str='24:00:00',
	slurm_job_name:str='ADToolbox',
    slurm_cpus:str="12",
	slurm_memory:str="100G",
	slurm_save_dir:str=os.getcwd()
	) -> None:
		self.slurm_template = slurm_template
		self.docker_template_qiime = docker_template_qiime
		self.singularity_template_qiime = singularity_template_qiime
		self.slurm_executer = slurm_executer
		self.slurm_wall_time = slurm_wall_time
		self.slurm_job_name = slurm_job_name
		self.slurm_outlog=self.slurm_job_name
		self.slurm_cpus = slurm_cpus
		self.slurm_save_dir = slurm_save_dir
		self.slurm_memory = slurm_memory

