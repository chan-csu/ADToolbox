import os
from __init__ import Main_Dir

"""
This module contains all the paths to the files and directories used in the program.

"""


RXN_DB = os.path.join(Main_Dir, "Database", "Reaction_Metadata.csv")

Seed_RXN_DB = os.path.join(Main_Dir, "Database", "reactions.json")




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
                feed_db=os.path.join(Main_Dir, "Database", 'Feed_DB.json'),
                protein_db=os.path.join(Main_Dir, "Database", 'Protein_DB.fasta'),):
        self.compound_db = compound_db
        self.reaction_db = reaction_db
        self.local_compound_db = local_compound_db
        self.local_reaction_db = local_reaction_db
        self.base_dir = base_dir
        self.csv_reaction_db = csv_reaction_db
        self.feed_db = feed_db
        self.protein_db = protein_db



class Metagenomics:
    """
    A class for Amplicon2Genome Configs

    """

    def __init__(self, 
                amplicon2genome_k=10,
                amplicon2genome_similarity=0.97,
                amplicon2genome_outputs_dir=os.path.join(Main_Dir,"Genomes"),
                amplicon2genome_db=os.path.join(Main_Dir,'Database','Amplicon2GenomeDBs'),
                qiime_outputs_dir=os.path.join(Main_Dir,'Metagenomics_Data','QIIME_Outputs'),
                genomes_json_info=os.path.join(Main_Dir,"Genomes","Amplicon2Genome_OutInfo.json"),
                feature_table_dir=os.path.join(Main_Dir,"Metagenomics_Data","QIIME_Outputs","feature-table.tsv"),
                rep_seq_fasta=os.path.join(Main_Dir,"Metagenomics_Data","QIIME_Outputs","dna-sequences.fasta"),
                taxonomy_table_dir=os.path.join(Main_Dir,"Metagenomics_Data","QIIME_Outputs","taxonomy.tsv"),
                genome_alignment_output=os.path.join(Main_Dir,"Outputs"),
                genome_alignment_output_json=os.path.join(Main_Dir,"Outputs","Alignment_Info.json"),
                vsearch=os.path.join(os.path.join(os.path.dirname(os.path.realpath(__file__)),"pkg_data")),
                genome_adm_map_json=os.path.join(Main_Dir,"Outputs","ADM_From_Alignment_JSON_Output.json"),
                csv_reaction_db=Database().csv_reaction_db,
                bit_score=40,
                e_value=10**-5

                 ):
        self.k = amplicon2genome_k
        self.amplicon2genome_similarity = amplicon2genome_similarity
        self.amplicon2genome_outputs_dir = amplicon2genome_outputs_dir
        self.amplicon2genome_db = amplicon2genome_db
        self.qiime_outputs_dir = qiime_outputs_dir
        self.genomes_json_info = genomes_json_info
        self.feature_table_dir = feature_table_dir
        self.rep_seq_fasta = rep_seq_fasta
        self.taxonomy_table_dir = taxonomy_table_dir
        self.protein_db=Database().protein_db
        self.seed_rxn_db=Seed_RXN_DB
        self.genome_alignment_output = genome_alignment_output
        self.genome_alignment_output_json=genome_alignment_output_json
        self.bit_score = bit_score
        self.e_value = e_value
        self.vsearch=vsearch
        self.genome_adm_map_json=genome_adm_map_json
        self.csv_reaction_db=csv_reaction_db



class Original_ADM1:

    def __init__(self,
                base_dir=os.path.join(Main_Dir, "Database","ADM1"),
                model_parameters=os.path.join(Main_Dir, "Database","ADM1",'ADM1_Model_Parameters.json'),
                reactions=os.path.join(Main_Dir, "Database","ADM1",'ADM1_Reactions.json'),
                species=os.path.join(Main_Dir, "Database","ADM1",'ADM1_Species.json'),
                base_parameters=os.path.join(Main_Dir, "Database","ADM1", 'ADM1_Base_Parameters.json'),
                initial_conditions=os.path.join(Main_Dir,"Database","ADM1","ADM1_Initial_Conditions.json"),
                inlet_conditions=os.path.join(Main_Dir,"Database","ADM1","ADM1_Inlet_Conditions.json"),
                metagenome_report=Metagenomics().genome_adm_map_json):
        self.model_parameters = model_parameters
        self.base_parameters = base_parameters
        self.initial_conditions = initial_conditions
        self.inlet_conditions = inlet_conditions
        self.metagenome_report = metagenome_report
        self.reactions = reactions
        self.species = species
        self.base_dir = base_dir

class Modified_ADM:

    def __init__(self,
                base_dir=os.path.join(Main_Dir, "Database","Modified_ADM"),
                model_parameters=os.path.join(Main_Dir, "Database","Modified_ADM",'Modified_ADM_Model_Parameters.json'),
                base_parameters=os.path.join(Main_Dir,"Database","Modified_ADM","Modified_ADM_Base_Parameters.json"),
                initial_conditions=os.path.join(Main_Dir,"Database","Modified_ADM","Modified_ADM_Initial_Conditions.json"),
                inlet_conditions=os.path.join(Main_Dir,"Database","Modified_ADM","Modified_ADM_Inlet_Conditions.json"),
                reactions=os.path.join(Main_Dir,"Database","Modified_ADM","Modified_ADM_Reactions.json"),
                species=os.path.join(Main_Dir,"Database","Modified_ADM","Modified_ADM_Species.json"),
                metagenome_report=Metagenomics().genome_adm_map_json):

        self.model_parameters = model_parameters
        self.base_parameters = base_parameters
        self.initial_conditions = initial_conditions
        self.inlet_conditions = inlet_conditions
        self.metagenome_report = metagenome_report
        self.reactions = reactions
        self.species = species
        self.Base_dir = base_dir

class Documentation:
    def __init__(self,
                 readme=os.path.join(os.path.dirname(os.path.realpath(__file__)),"pkg_data","README.md")):
        self.readme = readme
