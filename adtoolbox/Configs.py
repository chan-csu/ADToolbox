import os
from . import Main_Dir
from .Parameters import Initial_Conditions

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

    def __init__(self, Compound_DB=os.path.join(Main_Dir, "Database", 'Local_compounds.json'),
                 Reaction_DB=os.path.join(Main_Dir, "Database", 'Local_reactions.json')):
        self.Compound_DB = Compound_DB
        self.Reaction_DB = Reaction_DB


class Database:

    def __init__(self,
                 Compound_DB=os.path.join(
                     Main_Dir, "Database", 'compounds.json'),
                 Reaction_DB=os.path.join(
                     Main_Dir, "Database", 'reactions.json'),
                Local_Compound_DB=os.path.join(
                        Main_Dir, "Database", 'Local_compounds.json'),
                Local_Reaction_DB=os.path.join(
                        Main_Dir, "Database", 'Local_reactions.json'),

                     Base_Dir=os.path.join(Main_Dir,"Database"),
                     CSV_Reaction_DB=os.path.join(Main_Dir, "Database", 'Reaction_Metadata.csv'),
                Feed_DB=os.path.join(Main_Dir, "Database", 'Feed_DB.json'),
                Protein_DB=os.path.join(Main_Dir, "Database", 'Protein_DB.fasta'),

                     ):
        self.Compound_DB = Compound_DB
        self.Reaction_DB = Reaction_DB
        self.Local_Compound_DB = Local_Compound_DB
        self.Local_Reaction_DB = Local_Reaction_DB
        self.Base_Dir = Base_Dir
        self.CSV_Reaction_DB = CSV_Reaction_DB
        self.Feed_DB = Feed_DB
        self.Protein_DB = Protein_DB

class Metagenomics:
    """
    A class for Amplicon2Genome Configs

    """

    def __init__(self, 
                Amplicon2Genome_K=10,
                 Amplicon2Genome_Similarity=0.97,
                 Amplicon2Genome_Outputs_Dir=os.path.join(Main_Dir,"Genomes"),
                 Amplicon2Genome_DBs=os.path.join(Main_Dir,'Database','Amplicon2GenomeDBs'),
                 QIIME_Outputs_Directory=os.path.join(Main_Dir,'Metagenomics_Data','QIIME_Outputs'),
                 Genomes_JSON_Info=os.path.join(Main_Dir,"Genomes","Amplicon2Genome_OutInfo.json"),
                 Feature_Table_Dir=os.path.join(Main_Dir,"Metagenomics_Data","QIIME_Outputs","feature-table.tsv"),
                 Rep_Seq_Fasta=os.path.join(Main_Dir,"Metagenomics_Data","QIIME_Outputs","dna-sequences.fasta"),
                 Taxonomy_Table_Dir=os.path.join(Main_Dir,"Metagenomics_Data","QIIME_Outputs","taxonomy.tsv"),
                 Genome_Alignment_Output=os.path.join(Main_Dir,"Outputs","Alignment_Info.json"),
                Aligner="mmseqs2",
                bit_score=40,
                e_value=10**-5

                 ):
        self.K = Amplicon2Genome_K
        self.Amplicon2Genome_Similarity = Amplicon2Genome_Similarity
        self.Amplicon2Genome_Outputs_Dir = Amplicon2Genome_Outputs_Dir
        self.Amplicon2Genome_DB = Amplicon2Genome_DBs
        self.QIIME_Outputs_Dir = QIIME_Outputs_Directory
        self.Genomes_JSON_Info = Genomes_JSON_Info
        self.Feature_Table_Dir = Feature_Table_Dir
        self.Rep_Seq_Fasta = Rep_Seq_Fasta
        self.Taxonomy_Table_Dir = Taxonomy_Table_Dir
        self.Aligner = Aligner
        self.Protein_DB=Database().Protein_DB
        self.Genome_Alignment_Output = Genome_Alignment_Output
        self.bit_score = bit_score
        self.e_value = e_value


        
class Align_Genomes:
    """
    A class for Align_Genomes Configs

    """

    def __init__(self,
                 Genomes_JSON_Info=Metagenomics().Genomes_JSON_Info,
                 Alignment_JSON_Info=Main_Dir+"/Outputs/Alignment_Info.json",
                 Aligner=Alignment().aligner_name,
                 bit_score=Alignment().bit_score,
                 e_value=Alignment().e_value,
                 ):
        self.Genomes_JSON_Info = Genomes_JSON_Info
        self.Aligner = Aligner
        self.bit_score = bit_score
        self.e_value = e_value
        self.Alignment_JSON_Info = Alignment_JSON_Info


class Report:
    """
    A class for Report Configs

    """

    def __init__(self,
    ADM_From_Alignment_JSON_Output=os.path.join(Main_Dir,"Reports","ADM1_From_Alignment_JSON_Output.json")
    ):

        
        self.ADM_From_Alignment_JSON_Output = ADM_From_Alignment_JSON_Output
        self.Seed_RXN_DB=Reaction_Toolkit().Reaction_DB
        self.RXN_DB=Database().CSV_Reaction_DB
        self.Alignment_JSON_Info=Metagenomics().Genome_Alignment_Output
        self.bit_score=Metagenomics().bit_score
        self.e_value=Metagenomics().e_value


class Original_ADM1:

    def __init__(self,
                Base_dir=os.path.join(Main_Dir, "Database","ADM1"),
                Model_Parameters=os.path.join(Main_Dir, "Database","ADM1",'ADM1_Model_Parameters.json'),
                Reactions=os.path.join(Main_Dir, "Database","ADM1",'ADM1_Reactions.json'),
                Species=os.path.join(Main_Dir, "Database","ADM1",'ADM1_Species.json'),
                Base_Parameters=os.path.join(Main_Dir, "Database","ADM1", 'ADM1_Base_Parameters.json'),
                Initial_Conditions=os.path.join(Main_Dir,"Database","ADM1","ADM1_Initial_Conditions.json"),
                Inlet_Conditions=os.path.join(Main_Dir,"Database","ADM1","ADM1_Inlet_Conditions.json"),
                Metagenome_Report=Report().ADM_From_Alignment_JSON_Output):
        self.Model_Parameters = Model_Parameters
        self.Base_Parameters = Base_Parameters
        self.Initial_Conditions = Initial_Conditions
        self.Inlet_Conditions = Inlet_Conditions
        self.Metagenome_Report = Metagenome_Report
        self.Reactions = Reactions
        self.Species = Species
        self.Base_dir = Base_dir

class Modified_ADM:

    def __init__(self,
                Base_dir=os.path.join(Main_Dir, "Database","Modified_ADM"),
                Model_Parameters=os.path.join(Main_Dir, "Database","Modified_ADM",'Modified_ADM_Model_Parameters.json'),
                Base_Parameters=os.path.join(Main_Dir,"Database","Modified_ADM","Modified_ADM_Base_Parameters.json"),
                Initial_Conditions=os.path.join(Main_Dir,"Database","Modified_ADM","Modified_ADM_Initial_Conditions.json"),
                Inlet_Conditions=os.path.join(Main_Dir,"Database","Modified_ADM","Modified_ADM_Inlet_Conditions.json"),
                Reactions=os.path.join(Main_Dir,"Database","Modified_ADM","Modified_ADM_Reactions.json"),
                Species=os.path.join(Main_Dir,"Database","Modified_ADM","Modified_ADM_Species.json"),
                Metagenome_Report=Report().ADM_From_Alignment_JSON_Output):
        self.Model_Parameters = Model_Parameters
        self.Base_Parameters = Base_Parameters
        self.Initial_Conditions = Initial_Conditions
        self.Inlet_Conditions = Inlet_Conditions
        self.Metagenome_Report = Metagenome_Report
        self.Reactions = Reactions
        self.Species = Species
        self.Base_dir = Base_dir

class Documentation:
    def __init__(self,
                 ReadMe=os.path.join(os.path.dirname(os.path.realpath(__file__)),"pkg_data","README.md")):
        self.ReadMe = ReadMe
