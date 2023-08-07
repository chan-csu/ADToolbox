from distutils.log import warn
import subprocess
import os
from collections import UserDict
import random
from matplotlib import streamplot
import pandas as pd
import json
import numpy as np
import re
import requests
import time
from requests.adapters import HTTPAdapter
import utils
import configs
from requests.packages.urllib3.util.retry import Retry
from requests.exceptions import Timeout
from bs4 import BeautifulSoup
from datetime import datetime
from pathlib import Path
from collections import Counter
import pathlib
import tarfile
import configs
from rich.progress import track,Progress
import rich
from __init__ import Main_Dir
from typing import Union
from dataclasses import dataclass
import dataclasses
from utils import (wrap_for_slurm,
                   fasta_to_dict,
                   extract_zipped_file,
                   needs_repair,
                   create_mmseqs_database,
                   index_mmseqs_db,
                   mmseqs_search,
                   mmseqs_result_db_to_tsv) 
# import doctest
# doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
class AdditiveDict(UserDict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    def __add__(self, other:dict):
        if isinstance(other,dict): 
            return AdditiveDict({i: self[i] + other.get(i,0) for i in self.keys()})
        else:
            raise TypeError("The other object is not a dict")
            
@dataclass
class Experiment:
    """
    This class creates an interface for the experimental data to be used in different places in ADToolbox.
    First you should give each experiment a name. Time must be a list of time points in days, and there must be a time 0 point assinged to each experiment.
    variables must be a list of integers that represent the variables that are the index of the ADM species that we have concentration data for.
    data must be a list of lists. Each list in the list must be a list of concentrations for each species at each time point.
    IMPORTANT: The order of the species in the data list must match the order of the species in the variables list.
    reference is an optional argument that can be used to provide a reference for the experimental data. If using the database module 
    to query for Experiment objects you can query by name or reference or model_type. So, having a descriptive reference can be useful for querying as well.
    default model name is "modified-adm". This can be changed by passing a different model name to the model_name argument. This also helps with querying.
    
    Args:
        name (str): A unique name for the experiment.
        time (list): A list of time points in days.
        variables (list): A list of integers that represent the variables that are the index of the ADM species that we have concentration data for.
        data (list): A list of lists. Each list in the list must be a list of concentrations for each species at each time point.
        reference (str, optional): A reference for the experimental data. Defaults to ''.
        model_name (str, optional): The name of the model that the experimental data is for. Defaults to "modified-adm".
    
    Examples:
        >>> from adtoolbox import configs
        >>> import json
        >>> with open(configs.Database().species,"r") as f:
        ...     species=json.load(f)
        >>> S_su_index=species.index("S_su")
        >>> S_aa_index=species.index("S_aa")
        >>> exp=Experiment(name="Test",time=[0,1,2],variables=[S_su_index,S_aa_index],data=[[1,2,3],[4,5,6]],reference="Test reference")
        
    """
    name:str
    time: list[float]
    variables: list[int]
    data: list[list[float]]
    reference: str = ""
    model_name: str = "modified-adm"
    
    def __post_init__(self):
        self.data=np.array(self.data).T
        self.validate()
    
    def validate(self):
        assert len(self.time)==self.data.shape[0], "Number of time points must match number of rows in data."
        assert len(self.variables)==self.data.shape[1] , "Number of variables must match number of columns in data."
        assert self.time[0]==0, "Time must start at 0."
        return "successful"
    
    def to_dict(self):
        return {"name":self.name,
                "initial_conditions":self.initial_conditions,
                "time":self.time,
                "variables":self.variables,
                "data":self.data.tolist(),
                "reference":self.reference}
    

    
    
@dataclass
class Feed:

    """
    The Feed class is used to store the feed information, and later use it in the modified ADM model.
    all the entered numbers must in percentages. Carbohudrates, lipids, and proteins and si must sum up to 100, 
    and they form the total dissolved solids. Carbohydrates, lipids, proteins, and xi must sum up to 100, and they form the total suspended solids.
    
    IMPORTANT: It is assumed that lipid, proteins and carbohydrates have the same fraction in soluble and insoluble fractions.
    
    Args:
        name (str): A unique name for the feed.
        carbohydrates (float): percentage of carbohydrates in the feed.
        lipids (float): percentage of lipids in the feed.
        proteins (float): percentage of proteins in the feed.
        tss (float): percentage of total COD in the form of suspended solids.
        si (float): percentage of percentage of soluble inorganics in the TDS.
        xi (float): percentage of percentage of insoluble inorganics in the TSS.
        reference (str, optional): A reference for the feed data. Defaults to ''.    
    
    Examples:
        >>> feed=Feed(name="Test",carbohydrates=20,lipids=20,proteins=20,si=20,xi=20,tss=70)
        >>> assert feed.ch_tss==feed.lip_tss==feed.prot_tss==feed.xi_tss==0.25
        
    """
    # total_cod:float Transfer to base parameters
    name:str            # A unique name for the feed
    carbohydrates:float # percentage of carbohydrates in the feed
    lipids:float        # percentage of lipids in the feed
    proteins:float      # percentage of proteins in the feed
    tss:float           # percentage of total COD in the form of suspended solids
    si:float            # percentage of percentage of soluble inorganics in the TDS
    xi:float            # percentage of percentage of insoluble inorganics in the TSS
    reference:str=''    # A reference for the feed data

    def __post_init__(self):
        if self.carbohydrates+self.lipids+self.proteins>100:
            raise ValueError("The sum of the percentages must less than 100")
        if self.carbohydrates+self.lipids+self.proteins+self.si<1:
            warn("The sum of lipids, carbohydrates, proteins are suspiciously low! Make sure youhave input the numbers in percentages!")
        tss_sum=self.carbohydrates/100+self.lipids/100+self.proteins/100+self.xi/100
        self.ch_tss=self.carbohydrates/100/tss_sum
        self.lip_tss=self.lipids/100/tss_sum
        self.prot_tss=self.proteins/100/tss_sum
        self.xi_tss=self.xi/100/tss_sum
        tds_sum=self.carbohydrates/100+self.lipids/100+self.proteins/100+self.si/100
        self.ch_tds=self.carbohydrates/100/tds_sum
        self.lip_tds=self.lipids/100/tds_sum
        self.prot_tds=self.proteins/100/tds_sum
        self.si_tds=self.si/100/tds_sum
    
    def to_dict(self)->dict:
        return {"name":self.name,
                "carbohydrates":self.carbohydrates,
                "lipids":self.lipids,
                "proteins":self.proteins,
                "tss":self.tss,
                "si":self.si,
                "xi":self.xi,
                "reference":self.reference}

    

@dataclass
class MetagenomicsStudy:
    """
    This class is used to communicate between the metagenomics studies database and the ADM model.
    
    Args:
        name (str): The name of the metagenomics study. Its okay if it is not unique.
        study_type (str): The type of the metagenomics study. It can be "amplicon" or "WGS".
        microbiome (str): The microbiome that the metagenomics study is about.
        sample_accession (str): The SRA sample accession number of the metagenomics study. This must be unique.
        comments (str): Any comments that you want to add to the metagenomics study.
        study_accession (str): The SRA study accession number of the metagenomics study.   
    
    Examples:
        >>> study=MetagenomicsStudy(name="Test",study_type="WGS",microbiome="test_microbiome",sample_accession="test_accession",comments="test_comments",study_accession="test_study_accession")
        >>> assert study.name=="Test"

    """
    name:str
    study_type:str
    microbiome:str
    sample_accession:str
    comments:str
    study_accession:str
    
    def to_dict(self)->dict:
        return {"name":self.name,
                "study_type":self.study_type,
                "microbiome":self.microbiome,
                "sample_accession":self.sample_accession,
                "comments":self.comments,
                "study_accession":self.study_accession}

class Reaction:
    """
    This class provides a simple interface between information about biochemical reactions and multiple functionalities of ADToolbox.
    In order to instantiate a reaction object, you need to pass a dictionary of the reaction information.
    This dictionary must include 'name','stoichiometry' keys. This follows the format of the seed database.
    stoichiometry must be formatted like seed database. The seed database format is as follows:
    stoichiometry: '-1:cpd00079:0:0:\"D-glucose-6-phosphate\";1:cpd00072:0:0:\"D-fructose-6-phosphate\"'

    Args:
        data (dict): A dictionary containing the reaction information. This follows the format of the seed database.


    Examples:
        >>> A={"name":'D-glucose-6-phosphate aldose-ketose-isomerase',"stoichiometry":'-1:cpd00079:0:0:\"D-glucose-6-phosphate\";1:cpd00072:0:0:\"D-fructose-6-phosphate\"'}
        >>> a=Reaction(A)
        >>> print(a)
        D-glucose-6-phosphate aldose-ketose-isomerase

    """
    def __init__(self, data:dict)->None:
        self.data = data

    def __str__(self)->str:
        return self.data['name']

    @property
    def stoichiometry(self)->dict:
        """
        Returns the stoichiometry of the reaction by the seed id of the compounds as key and the
        stoichiometric coefficient as value.
        Examples:
            >>> A={"name":'D-glucose-6-phosphate aldose-ketose-isomerase',"stoichiometry":'-1:cpd00079:0:0:\"D-glucose-6-phosphate\";1:cpd00072:0:0:\"D-fructose-6-phosphate\"'}
            >>> a=Reaction(A)
            >>> a.stoichiometry=={'cpd00079': -1, 'cpd00072': 1}
            True
        
        Args:
            self (Reaction): An instance of the Reaction.

        Returns:
            dict: The stoichiometry of the reaction 
        """
        return {compound.split(':')[1]:float(compound.split(':')[0]) for compound in self.data['stoichiometry'].split(';') }


class Metabolite:
    """
    This class provides a simple interface between information about metabolites and multiple functionalities of ADToolbox.
    In order to instantiate a metabolite object, you need to pass a dictionary of the metabolite information.
    This dictionary must include 'name','mass','formula' keys. This follows the format of the seed database.
    formula must be formatted like seed database. The seed database format is as follows:
    formula: 'C6H12O6'
    Possibly the main advantage of instantiating a metabolite object is that it provides a COD attribute that can be used to convert
    the concentration of the metabolite from g/l to gCOD/l. This is useful for comparing the experimental data with the model outputs.

    Args:
        data (dict): A dictionary containing the metabolite information. This follows the format of the seed database.


    Examples:
        >>> A={"name":"methane","mass":16,"formula":"CH4"}
        >>> a=Metabolite(A)
        >>> print(a)
        methane

    """

    def __init__(self, data):
        self.data = data
        self.cod = self.cod_calc()

    def __str__(self) -> str:
        return self.data['name']

    def cod_calc(self)->float:
        """
        Calculates the conversion rates for g/l -> gCOD/l

        Examples:
            >>> A={"name":"methane","mass":16,"formula":"CH4"}
            >>> a=Metabolite(A)
            >>> a.cod
            4.0

        Args:
            self (Metabolite): An instance of the Metabolite class: Note

        Returns:
            float: COD conversion from g/l to gCOD/l

        """
        if self.data['formula'] and self.data['mass']:
            contents = {}
            atoms = ["H", "C", "O"]
            mw = self.data['mass']
            for atom in atoms:
                if re.search(atom+'\d*', self.data['formula']):
                    if len(re.search(atom+'\d*', self.data['formula']).group()[1:]) == 0:
                        contents[atom] = 1
                    else:
                        contents[atom] = int(
                            re.search(atom+'\d*', self.data['formula']).group()[1:])
                else:
                    contents[atom] = 0
            return 1/mw*(contents['H']+4*contents['C']-2*contents['O'])/4*32

        else:
            return 'None'


class SeedDB:

    """
    This class is designed to interact with seed database. The main advantage of using this class is that it can be used to instantiate
    a reaction and metabolite object, and it provides extra functionalities that rely on information in the seed database. For example, 
    If there is a chemical formula assigned to a metabolite in the seed database, then the informattion about the COD of that metabolite
    can be computed using the chemical formula. 
    
    Args:
        compound_db (str, optional): The path to the seed compound database. Defaults to configs.SeedDB().compound_db.
        reaction_db (str, optional): The path to the seed reaction database. Defaults to configs.SeedDB().reaction_db.
    
    Examples:
        >>> seed_db=SeedDB()
        >>> assert seed_db.compound_db==configs.SeedDB().compound_db
        >>> assert seed_db.reaction_db==configs.SeedDB().reaction_db

    """

    def __init__(self, compound_db:str =configs.SeedDB().compound_db,
                       reaction_db=configs.SeedDB().reaction_db) -> None:
        
        self.reaction_db = reaction_db
        self.compound_db = compound_db

    def instantiate_rxns(self, seed_id:str)->Reaction:
        """
        This method is used to instantiate reaction objects from the seed database.
        in order to instantiate a reaction object, you need to pass the seed identifier for that reaction.
        
        Args:
            seed_id (str): The seed identifier for the reaction.
    
        Returns:
            Reaction: An instance of the Reaction class.
        
        Examples:
            >>> seed_db=SeedDB()
            >>> rxn=seed_db.instantiate_rxns("rxn00558")
            >>> assert rxn.data["name"]=="D-glucose-6-phosphate aldose-ketose-isomerase"
        """
        db=pd.read_json(self.reaction_db)
        return Reaction(data=db[db["id"]==seed_id].to_dict(orient="records")[0])

    def instantiate_metabs(self, seed_id:str)->Metabolite:
        """
        This method is used to instantiate metabolite objects from the seed database.
        In order to instantiate a metabolite object, you need to pass the seed identifier for that metabolite.

        Args:
            seed_id (str): The seed identifier for the metabolite.
        
        Returns:
            Metabolite: An instance of the Metabolite class. 
        
        Examples:
            >>> seed_db=SeedDB()
            >>> metab=seed_db.instantiate_metabs("cpd01024")
            >>> assert metab.cod==4.0
        """
        db=pd.read_json(self.compound_db)
        return Metabolite(data=db[db["id"]==seed_id].to_dict(orient="records")[0])

    def get_seed_rxn_from_ec(self, ec_number:str)->list:
        """
        This method is used to get the seed reaction identifiers for a given EC number.

        Args:
            ec_number (str): The EC number.
        
        Returns:
            list: A list of seed reaction identifiers.
        
        Examples:
            >>> seed_db=SeedDB()
            >>> seed_rxn_list=seed_db.get_seed_rxn_from_ec("1.1.1.1")
            >>> assert len(seed_rxn_list)>0
        
        """
        db=pd.read_json(self.reaction_db)
        db=db[db["ec_numbers"].apply(lambda x: ec_number in x if x else False)]
        db.drop_duplicates("id",inplace=True,keep="first")
        return db.to_dict(orient="records")
            

class Database:

    '''
    This class is designed to supply any data requirement for ADToolbox. All functionalisties for saving, loading, and querying data are implemented here.
    ADToolbox in general contains the following databases:
    
    - The seed reaction database
    
    - The seed compound database
    
    - ADToolbox's Feed database
    
    - ADToolbox's Metagenomics studies database
    
    - ADToolbox's Experimental data database
    
    - ADToolbox's Protein database
    
    - ADToolbox's Reaction database
    
    - GTDB-tk database for bacterial and archaeal 16s rRNA sequences
    
    - ADM and modified ADM model parameters
    
    This class is instantiated with a configs.Database object. This object contains the paths to all the databases that ADToolbox uses.
    Please refer to the documentation of each method for more information on the required configurations.
    
    Args:
        config (configs.Database, optional): A configs.Database object. Defaults to configs.Database().
    
    Examples:
        >>> db=Database(config=configs.Database())
        >>> assert type(db)==Database and type(db.config)==configs.Database

    '''
    def __init__(self, config:configs.Database=configs.Database())->None:
        self.config = config


    def initialize_protein_database(self)->None:
        """This function intializes ADToolbox's protein database by creating an empty fasta file.
        Be careful, this will overwrite any existing file with the same name.
        Logically, this needs method needs config.protein_db to be defined.
        
        Required Configs:
            - config.protein_db
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==False # This is just to make sure that the following lines create the file
            >>> db=Database(config=configs.Database(protein_db=os.path.join(Main_Dir,"protein_test_db.fasta"))) # point to a test non-existing file
            >>> db._initialize_protein_database() # initialize the protein database
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==True # check if the file is created
            >>> os.remove(os.path.join(Main_Dir,"protein_test_db.fasta")) # remove the file to clean up
        """

        with open(self.config.protein_db, 'w') as f:
            pass
    
    def initialize_reaction_database(self)->None:
        r"""This function intializes ADToolbox's reaction database by creating an empty tsv file.
        Be careful, this will overwrite any existing file with the same name.
        
        Required Configs:
            - config.reaction_db
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"reaction_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(reaction_db=os.path.join(Main_Dir,"reaction_test_db.tsv")))
            >>> db._initialize_reaction_database()
            >>> assert pd.read_table(os.path.join(Main_Dir,"reaction_test_db.tsv"),delimiter="\t").shape[0]==0
            >>> assert set(pd.read_csv(os.path.join(Main_Dir,"reaction_test_db.tsv"),delimiter="\t").columns)==set(["EC_Numbers","Seed Ids","Reaction Names","ADM1_Reaction","Modified_ADM_Reactions","Pathways"])
            >>> os.remove(os.path.join(Main_Dir,"reaction_test_db.tsv"))
        
        """
        pd.DataFrame(columns=["EC_Numbers","Seed Ids","Reaction Names","ADM1_Reaction","Modified_ADM_Reactions","Pathways"]).to_csv(self.config.reaction_db,index=False,sep="\t")
        
    def initialize_feed_database(self)->None:
        r"""This function intializes ADToolbox's Feed database by creating an empty tsv file.
        Be careful, this will overwrite any existing file with the same name.
        
        Required Configs:
            - config.feed_db
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(feed_db=os.path.join(Main_Dir,"feed_test_db.tsv")))
            >>> db._initialize_feed_database()
            >>> assert pd.read_table(os.path.join(Main_Dir,"feed_test_db.tsv"),delimiter='\t').shape[0]==0
            >>> assert set(pd.read_table(os.path.join(Main_Dir,"feed_test_db.tsv"),delimiter='\t').columns)==set(["Name","Carbohydrates","Lipids","Proteins","TSS","SI","XI","Reference"])
            >>> os.remove(os.path.join(Main_Dir,"feed_test_db.tsv"))
        
        """
        pd.DataFrame(columns=["Name","Carbohydrates","Lipids","Proteins","TSS","SI","XI","Reference"]).to_csv(self.config.feed_db,index=False,sep="\t")
    
    def initialize_metagenomics_studies_database(self)->None:
        r"""This function intializes ADToolbox's Metagenomics studies database by creating an empty tsv file.
        Be careful, this will overwrite any existing file with the same name.
        
        Required Configs:
            - config.metagenomics_studies_db
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(metagenomics_studies_db=os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv")))
            >>> db._initialize_metagenomics_studies_database()
            >>> assert pd.read_table(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"),delimiter="\t").shape[0]==0
            >>> assert set(pd.read_table(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"),delimiter="\t").columns)==set(["Name","Study Type","Microbiome","Sample Accession","Comments","Study Accession"])
            >>> os.remove(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"))
        
        """
        pd.DataFrame(columns=["Name","Study Type","Microbiome","Sample Accession","Comments","Study Accession"]).to_csv(self.config.metagenomics_studies_db,index=False,sep="\t")
        
    def initialize_experimental_data_database(self)->None:
        """This function intializes ADToolbox's experimental data database by creating an empty json file.
        Be careful, this will overwrite any existing file with the same name.

        Required Configs:
            - config.experimental_data_db
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"experimental_data_test_db.json"))==False
            >>> db=Database(config=configs.Database(experimental_data_db=os.path.join(Main_Dir,"experimental_data_test_db.json")))
            >>> db._initialize_experimental_data_database()
            >>> assert pd.read_json(os.path.join(Main_Dir,"experimental_data_test_db.json")).shape[0]==0
            >>> assert set(pd.read_json(os.path.join(Main_Dir,"experimental_data_test_db.json")).columns)==set(["Name","Initial Conditions","Time","Variables","Data","Reference"])
            >>> os.remove(os.path.join(Main_Dir,"experimental_data_test_db.json"))
        
        """
        pd.DataFrame(columns=["Name","Initial Conditions","Time","Variables","Data","Reference"]).to_json(self.config.experimental_data_db)
        
    
    def filter_seed_from_ec(self, ec_list,
                            reaction_db=configs.Database().reaction_db,
                            compound_db=configs.Database().compound_db,
                            local_reaction_db=configs.Database().local_reaction_db,
                            local_compound_db=configs.Database().local_compound_db,
                            save:bool=True) -> tuple:
        """

        This function takes a list of ec numbers Generates a mini-seed JSON files. This is supposed to
        make the code a lot faster, but makes no difference in terms of outputs, and won't probably need
        frequent updates.

        Args:
            ec_list (list): The list of EC numbers that you want the mini seed database to include.
            reaction_db (str, optional): The path to the main model seed reaction database. Defaults to configs.Database().reaction_db.
            compound_db (str, optional): The path to the main model seed compound database. Defaults to configs.Database().compound_db.
            local_reaction_db (str, optional): Path to where you want the mini seed reaction database to be saved. Defaults to configs.Database().local_reaction_db.
            local_compound_db (str, optional): Path to where you want the mini seed compound database to be saved. Defaults to configs.Database().local_compound_db.
        
        Returns:
            tuple: A tuple containing the mini seed reaction database and the mini seed compound database.

        """
        ## TODO: Make this function more efficient by using pandas an search using .str.contains. This would be way faster
        with open(reaction_db, 'r') as f:
            main_reaction_db = json.load(f)
        with open(compound_db, 'r') as f:
            main_compound_db = json.load(f)

        cached_compounds = []
        local_rxn_db = []
        local_comp_db = []
        for ec in track(ec_list,description="Extracting the relevant reactions:"):
            for ind, rxn in enumerate(main_reaction_db):
                if  ec in rxn.setdefault('ec_numbers', []):
                    local_rxn_db.append(rxn)
                    for Mets in rxn["compound_ids"].split(";"):
                        if Mets not in cached_compounds:
                            cached_compounds.append(Mets)

        cached_compounds = list(set(cached_compounds))

        for compound in track(cached_compounds,description="Extracting the relevant compounds:"):
            for ind, Comp in enumerate(main_compound_db):
                if compound == Comp["id"]:
                    local_comp_db.append(Comp)


        if save:
            with open(local_reaction_db, 'w') as f:
                json.dump(local_rxn_db, f)
            with open(local_compound_db, 'w') as f:
                json.dump(local_comp_db, f)

        return local_rxn_db, local_comp_db



    # def add_protein_from_uniprot(self, uniprot_ecs):

    #     pdb.set_trace()
    #     with open(self.config.protein_db, 'a') as f:
    #         for items in track(uniprot_ecs,description="[yellow] --> Fetching the protein sequences from Uniprot: "):
    #             try:
    #                 file = Database.session.get(
    #                     f"https://rest.uniprot.org/uniprotkb/{items[0]}.fasta", timeout=10)
    #             except requests.exceptions.ConnectionError:
    #                 print("Could not fetch the sequence! Trying again ...")
    #                 time.sleep(10)
    #                 file = Database.session.get(
    #                     Base_URL+items[0]+".fasta", timeout=10)
    #                 if file.ok:
    #                     print(
    #                         f"Retry was successful: {items[0]} was fetched successfully!")
    #                 else:
    #                     print(
    #                         f"Retry was unsuccessful: {items[0]} ignored!")
    #                     # I can add this uniprot id to a log file
    #                 continue
    #             if file.ok:
    #                 f.write(f'>{items[0]}|{items[1]}\n')
    #                 f.write(''.join(file.text.split('\n')[1:-1]))
    #                 f.write('\n')
    
    def protein_db_from_ec(self, ec_list:list,mode:str="a") -> dict:
        """This function takes a list of EC numbers and fetches the protein sequences mapped to those EC numbers from Uniprot.
        The protein sequences are then written to the protein database file.
        
        Args:
            ec_list (list): A list of EC numbers.
            mode (str, optional): The mode in which the protein database file is opened. Defaults to "a".
            
        Returns: 
            dict: A dictionary containing the protein sequences mapped to the EC numbers."""
        session = requests.Session()
        retry = Retry(connect=3, backoff_factor=0.5)
        adapter = HTTPAdapter(max_retries=retry)
        session.mount('http://', adapter)
        
        protein_seqs={}
        
        with open(self.config.protein_db, mode) as f:
            for ec in track(ec_list,description="Writing the protein database:"):
                try:
                    file = session.get(
                        f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28ec%3A{ec}%29%20AND%20%28reviewed%3Atrue%29%20NOT%20%28taxonomy_id%3A2759%29%29", timeout=1000)

                except requests.exceptions.HTTPError or requests.exceptions.ConnectionError:

                    print("Request Error! Trying again ...")
                    time.sleep(30)
                    file = session.get(
                        f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28ec%3A{ec}%29%20AND%20%28reviewed%3Atrue%29%20NOT%20%28taxonomy_id%3A2759%29%29", timeout=1000)
                # This alsp does a sanity check

                except Exception:
                    print('Something went wrong!')
                text = file.text
                if text:
                    for line in text.split('\n'):
                        if line.startswith('>'):
                            uniprot=line.split('|')[1]
                            f.write(f'>{uniprot}|{ec}\n')

                        else:
                            f.write(f'{line}\n')
                            protein_seqs[uniprot]=line
        
        return protein_seqs

    def add_uniprot_to_protein_db(self,uniprot_id:str,ec:str):
        Base_URL = "https://rest.uniprot.org/uniprotkb/"
        session = requests.Session()
        retry = Retry(connect=3, backoff_factor=0.5)
        adapter = HTTPAdapter(max_retries=retry)
        session.mount('http://', adapter)
        with open(self.config.protein_db,mode="a") as f:
            try:
                file = session.get(
                    f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta", timeout=10)
            except requests.exceptions.ConnectionError:
                print("Could not fetch the sequence! Trying again ...")
                time.sleep(10)
                file = session.get(
                    Base_URL+uniprot_id+".fasta", timeout=10)
                if file.ok:
                    print(
                        f"Retry was successful: {uniprot_id} was fetched successfully!")
                else:
                    print(
                        f"Retry was unsuccessful: {uniprot_id} ignored!")
                    # I can add this uniprot id to a log file
            if file.ok:
                f.write(f'>{uniprot_id}|{ec}\n')
                f.write(''.join(file.text.split('\n')[1:-1]))
                f.write('\n')
            

    @staticmethod
    def ec_from_csv(csv_file:str, Sep=',')->list:

        """This function reads a csv file and returns a list of ec numbers. The csv file should have a column named "EC_Numbers" and the delimiter 
        that is used in in your file should be defined as an argument.
        
        Args:
            csv_file (str): The path to the csv file.
            Sep (str, optional): The delimiter used in the csv file. Defaults to ','.
        
        Returns:
            list: A cleaned list of ec numbers."""
        

        try:

            ec_list = pd.read_table(csv_file, sep=Sep,dtype="str")
            ec_list.dropna(axis=0)
            ec_list.drop_duplicates(subset="EC_Numbers", inplace=True)

            output_ec_list = list(ec_list["EC_Numbers"])
            assert len(
                output_ec_list) > 0, "The ec list is empty; Check the file"

        except FileNotFoundError:

            print("CSV file not found")

        except (pd.errors.ParserError, KeyError):

            print(
                "CSV file not in the correct format; Check the delimiter and column names!")

        except AssertionError as error:

            print(error)


        return output_ec_list


    

    def cazy_ec(self):
        '''
        Scraps the ec numbers from links of the Cazy website.

        Returns:
            list: A list of ec numbers.

        '''

        ec_list = []
        for link in Database.config.cazy_links:

            page = requests.get(link)
            soup = BeautifulSoup(page.content, "html.parser")
            results = soup.find("div", class_="cadre_principal").find_all(
                "th", class_="thec")
            for ec_number in results:
                if '-' not in ec_number.text.strip() and '.' in ec_number.text.strip():

                    ec_list.append(ec_number.text.strip())

        print("ec numbers extracted from Cazy website!")
        return ec_list

    
    def seed_from_ec(self,ec_Number:str, mode:str='single_match')->list:
        """ Given EC number, this function returns the corresponding Seed IDs from the SEED reaction database.
        
        Args:
            ec_Number (str): The EC number.
            mode (str, optional): The mode in which the function should run. Defaults to 'single_match'. You can choose between 'single_match' and 'Multiple_Match'.
            
        Returns:
            list: A list of Seed IDs."""

        ##TODO This function can use pandas to speed up the process

        with open(self.config.reaction_db,'r') as f: 

            data = json.load(f)

        if mode == 'single_match':

            for i in data:

                if i['ec_numbers']:

                    if ec_Number in i['ec_numbers']:

                        return [i["id"]]

        elif mode == 'Multiple_Match':

            matched_rxn = list(set(list([

                rxn["id"] for rxn in data if rxn['ec_numbers'] == [ec_Number]])))

            return matched_rxn

        else:
            print('Mode not recognized!')

    def instantiate_rxn_from_seed(self,seed_id_list:list)->list:
        """
        Function that idenifies reaction seed ID's and instantiates them in the Reaction class.
        Examples:
            >>> seed_id_List = ['rxn00002','rxn00003','rxn00005']
            >>> dbconfigs = configs.Database()
            >>> rxnlist = Database(dbconfigs).Instantiate_Rxn_From_Seed(seed_id_List)
            >>> assert type(rxnlist[0])==Reaction

        Args:
            seed_id_List (list): A list of relevant seed ID entries [rxn#####]

        Returns:
            rxn_list: A list including reaction instances in the database class for each seed ID in input list.

        """
        ## TODO: This function can be sped up by using pandas

        rxn_list = []
        with open(self.config.reaction_db) as f:

            data = json.load(f)

            for seed_id in seed_id_list:

                for i in data:

                    if i['id']:

                        if seed_id in i['id']:

                            rxn_list.append(Reaction(i))
                            break

        return rxn_list

    def metadata_from_ec(self,ec_list:list)->pd.DataFrame:
        """
        This function returns a pandas dataframe containing relevant pathway and reactions for
        each ec number input.
        Examples:
            >>> ec_list = ['1.1.1.1','1.1.1.2'] 
            >>> metadata= Database().Metadata_From_ec(ec_list) # doctest: +ELLIPSIS 
            Finding ...
            >>> assert type(metadata)==pd.DataFrame
            >>> assert set(metadata['ec_Numbers'].to_list())==set(ec_list)
            >>> assert set(["ec_Numbers", "seed_ids","Reaction_Names", "Pathways"])==set(metadata.columns)
            >>> assert metadata.shape==(len(ec_list),4)

        Args:
            ec_list (list): A list of relevant ec numbers.

        Returns:
            pd.DataFrame: A pandas dataframe including reaction metadata or pathways for each ec number in list.

        """
        full_table = {"ec_Numbers": [], "seed_ids": [],
                      "Reaction_Names": [], "Pathways": []}
        rich.print("Finding ec Metadata ...\n")
        for ec in track(ec_list, description= "Collecting Metadata for ec numbers"):
            seed_id = self.seed_from_ec(ec, Mode='Multiple_Match')
            full_table['ec_Numbers'].append(ec)
            full_table['seed_ids'].append(seed_id)
            Temp_rxns = Database().instantiate_rxn_from_seed(seed_id)
            Temp_rxn_Names = list([reaction.dict["name"]
                                   for reaction in Temp_rxns])
            Temp_rxn_Path = list([reaction.dict["pathways"]
                                  for reaction in Temp_rxns])
            full_table["Pathways"].append(Temp_rxn_Path)
            full_table["Reaction_Names"].append(Temp_rxn_Names)

        full_table = pd.DataFrame(full_table)

        return full_table

    def init_feedstock_database(self)-> None:
        """Initializes the feedstock database as an empty table.
        """
        feed_dict = {
                        'carbohydrates':[],
                        'lipids':[],
                        'proteins':[],
                        'tss':[],
                        'si':[],
                        'xi':[],
                        'reference':[]
                    }
        feed_df = pd.DataFrame(feed_dict)
        if os.path.exists(self.config.feed_db):
            answer=input("Feedstock database already exists. Do you want to overwrite it? (y/n)").lower()
            if answer!='y':
                rich.print("[red]Feedstock database not initialized!")
                return
            else:
                feed_df.to_csv(self.config.feed_db, index=False,sep='\t')
                rich.print("[green]Feedstock database initialized!")
        else:
            feed_df.to_csv(self.config.feed_db, index=False,sep='\t') 
            rich.print("[green]Feed stock database initialized!")            


    def add_feedstock_to_database(self, feed:Feed)->None:
        
        '''
        Adds a feedstock to the feedstock database. feed must be an instance of the Feed class.

        Args:
            feed (Feed): An instance of the Feed class.
        
        Returns:
            None

        '''
        if not isinstance(feed,Feed):
            print("Feed must be an instance of the Feed class!")
            return
        
        feed_db=pd.read_table(self.config.feed_db, sep='\t')
        feed_db=feed_db.append(dataclasses.asdict(feed), ignore_index=True)
        feed_db.to_csv(self.config.feed_db, index=False,sep='\t')
        rich.print("[green]Feedstock added to the database!")
    
    def build_mmseqs_database(self,save:bool,run:bool,container:str="None")->str:
        """Builds an indexed mmseqs database from the ADToolbox's fasta protein database.
        Args:
            save (bool): If True, the script will be saved to the current directory.
            run (bool): If True, the script will be run.
            container (str, optional): The container to run the script in. Defaults to "None".
        Returns:
            str: The script to build the mmseqs database.
        """
        script=create_mmseqs_database(self.config.protein_db,
                                      self.config.protein_db_mmseqs,
                                      container=container,
                                      save=None,
                                      run=False,
                                      config=configs.Config())
        # if container=="None":
        #     pass
        
        # elif container=="singularity":
        #     script=f"singularity exec --bind {self.config.protein_db}:{self.config.protein_db},{self.config.protein_db_mmseqs}:{self.config.protein_db_mmseqs} {self.config.adtoolbox_singularity} {script}"
        
        # elif container=="docker":
        #     script=f"docker run -v {self.config.protein_db}:{self.config.protein_db} -v {self.config.protein_db_mmseqs}:{self.config.protein_db_mmseqs} {self.config.adtoolbox_docker} {script}"
        
        # else:
        #     print("Container must be one of the following: None, singularity, docker")
        #     return
        
        if save:
            with open(pathlib.Path(save),"w") as f:
                f.write(script)
        if run:
            subprocess.run(script,shell=True)
        
        return script


    def download_adm_parameters(self)->None:
        '''
        Downloads the parameters needed for running ADM models in ADToolbox.
        
        '''
        if not os.path.exists(self.config.adm_parameters_base_dir):
            os.makedirs(self.config.adm_parameters_base_dir)
        
        for url in track(self.config.adm_parameters_urls.values(),description="Downloading ADM parameters"):
            r = requests.get(url, allow_redirects=True)
            with open(os.path.join(self.config.adm_parameters_base_dir,url.split("/")[-1]), 'wb') as f:
                f.write(r.content)
        rich.print(f"[green]ADM parameters were downloaded to {self.config.adm_parameters_base_dir}")
 
        
    def download_seed_databases(self) -> None :
        "Download the modelseed rection database"
        r = requests.get(self.config.seed_rxn_url, allow_redirects=True,stream=True)
        with open(self.config.reaction_db, 'wb') as f:
            f.write(r.content)
        rich.print(f"[green]Reaction database downloaded to {self.config.reaction_db}")
        r=requests.get(self.config.seed_compound_url,allow_redirects=True,stream=True)
        with open(self.config.compound_db, 'wb') as f:
            f.write(r.content)
        rich.print(f"[green]Compound database downloaded to {self.config.compound_db}")

    def download_protein_database(self) -> None:
        r = requests.get(self.config.protein_db_url, allow_redirects=True)
        with open(self.config.protein_db, 'wb') as f:
            f.write(r.content)
        rich.print(f"[green]Protein database downloaded to {self.config.protein_db}")
        
    def download_reaction_database(self)->None:
        r = requests.get(self.config.adtoolbox_rxn_db_url, allow_redirects=True)
        with open(self.config.csv_reaction_db, 'wb') as f:
            f.write(r.content)
        rich.print(f"[green]Reaction database downloaded to {self.config.csv_reaction_db}")

    
    def download_feed_database(self)-> None:
        r = requests.get(self.config.feed_db_url, allow_redirects=True)
        with open(self.config.feed_db, 'wb') as f:
            f.write(r.content)
        rich.print(f"[green]Feed database downloaded to {self.config.feed_db}")
    
    def download_qiime_classifier_db(self)->None:
        r = requests.get(self.config.qiime_classifier_db_url, allow_redirects=True,stream=True)
        block_size = 1024
        total_size = int(r.headers.get('content-length', 0))
        if not os.path.exists(Path(self.config.qiime_classifier_db).parent):
            os.makedirs(Path(self.config.qiime_classifier_db).parent)
        with open(self.config.qiime_classifier_db, 'wb') as f:
            with Progress() as progress:
                task = progress.add_task("Downloading the qiime's classifier database:", total=total_size)
                for data in r.iter_content(block_size):
                    progress.update(task, advance=len(data))
                    f.write(data)
        
        rich.print(f"[green]Qiime's classifier database downloaded to {self.config.qiime_classifier_db}")
            


    def get_metagenomics_studies(self) -> list:
        """
        This function will return accession numbers in all metagenomics studies on the kbase.
        Requires:
            <R>Configs.Database</R>
        Satisfies:
            <S>project_accession</S>
        """
        metagenomics_studies=pd.read_table(self.config.studies.metagenomics_studies,delimiter="\t")
        return metagenomics_studies.to_dict(orient="records")

    
    def get_experimental_data_studies(self)->list:
        """
        This function will retrieve the studies as a list of dictionaries.

        Returns:
            list: A list of dictionaries containing the studies.

        """
        experimental_studies=pd.read_table(self.config.studies.experimental_data_references,delimiter="\t")
        return experimental_studies.to_dict(orient="list")
        
         

    def download_studies_database(self)->None:
        """
        This function will download the kbase database from the remote repository.
        """
        for i in self.config.studies.urls:
            r = requests.get(self.config.studies.urls[i], allow_redirects=True)
            if not os.path.exists(self.config.studies.base_dir):
                os.makedirs(self.config.studies.base_dir,exist_ok=True)
            with open(os.path.join(self.config.studies.base_dir,self.config.studies.urls[i].split("/")[-1]), 'wb') as f:
                f.write(r.content)
            rich.print(f"[bold green]Downloaded {self.config.studies.urls[i]}[/bold green]")
    
    def download_amplicon_to_genome_db(self,progress=True):
        """
        This function will automatically download the required files for Amplicon to Genome functionality.
        """
        if not os.path.exists(self.config.amplicon_to_genome_db):
            os.mkdir(self.config.amplicon_to_genome_db)

        url = {'Version': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/VERSION',
               'MD5SUM': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/MD5SUM',
               'FILE_DESCRIPTIONS': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/FILE_DESCRIPTIONS',
               'metadata_field_desc': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/metadata_field_desc.tsv',
               'bac120_metadata': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tar.gz',
               'bac120_ssu': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_all/ssu_all.tar.gz'
               }
        if progress:
            for keys in ['Version', 'MD5SUM', 'FILE_DESCRIPTIONS']:
                with requests.get(url[keys], allow_redirects=True, stream=progress) as r:
                    total_size = int(r.headers.get('content-length', 0))
                    block_size = 1024
                    with Progress() as progress:
                        task1 = progress.add_task("Downloading " + keys, total=total_size)
                        with open(os.path.join(self.config.amplicon_to_genome_db, keys), 'wb') as f:
                            for data in r.iter_content(block_size):
                                progress.update(task1, advance=len(data))
                                f.write(data)
            with requests.get(url['metadata_field_desc'], allow_redirects=True, stream=progress) as r:
                total_size = int(r.headers.get('content-length', 0))
                block_size = 1024
                with Progress() as progress:
                    task1 = progress.add_task("Downloading metadata_field_desc.tsv", total=total_size)
                    with open(os.path.join(self.config.amplicon_to_genome_db, 'metadata_field_desc.tsv'), 'wb') as f:
                        for data in r.iter_content(block_size):
                            progress.update(task1, advance=len(data))
                            f.write(data)

            for keys in ['bac120_metadata', 'bac120_ssu']:
                with requests.get(url[keys], allow_redirects=True, stream=progress) as r:
                    total_size = int(r.headers.get('content-length', 0))
                    block_size = 1024
                    with Progress() as progress:
                        task1 = progress.add_task("Downloading " + keys, total=total_size)
                        with open(os.path.join(self.config.amplicon_to_genome_db, url[keys].split("/")[-1]), 'wb') as f:
                            for data in r.iter_content(block_size):
                                progress.update(task1, advance=len(data))
                                f.write(data)
                with tarfile.open(os.path.join(self.config.amplicon_to_genome_db, url[keys].split("/")[-1])) as f_in:
                    f_in.extractall(self.config.amplicon_to_genome_db)

                
                os.remove(os.path.join(self.config.amplicon_to_genome_db, url[keys].split("/")[-1]))
        else:
            for keys in ['Version', 'MD5SUM', 'FILE_DESCRIPTIONS']:
                with requests.get(url[keys], allow_redirects=True, stream=progress) as r:
                    with open(os.path.join(self.config.amplicon_to_genome_db, keys), 'wb') as f:
                        f.write(r.content)
            with requests.get(url['metadata_field_desc'], allow_redirects=True, stream=progress) as r:
                with open(os.path.join(self.config.amplicon_to_genome_db, 'metadata_field_desc.tsv'), 'wb') as f:
                    f.write(r.content)
            for keys in ['bac120_metadata', 'bac120_ssu']:
                with requests.get(url[keys], allow_redirects=True, stream=progress) as r:
                    with open(os.path.join(self.config.amplicon_to_genome_db, url[keys].split("/")[-1]), 'wb') as f:
                        f.write(r.content)
                with tarfile.open(os.path.join(self.config.amplicon_to_genome_db, url[keys].split("/")[-1])) as f_in:
                    f_in.extractall(self.config.amplicon_to_genome_db)
        
        rich.print("[bold green]Downloaded all the required files for Amplicon to Genome functionality.[/bold green]")
                    
                        

            
    def download_all_databases(self)->None:
        """
        This function will download all the required databases for all the functionalities of ADToolbox.
        """
        if not os.path.exists(self.config.base_dir):
            os.makedirs(self.config.base_dir)
        self.download_seed_databases()
        self.download_adm_parameters()
        self.download_protein_database()
        self.download_reaction_database()
        self.download_feed_database()
        self.download_studies_database()
        self.download_amplicon_to_genome_db()
        self.download_qiime_classifier_db()
        

     

class Metagenomics:

    """
    This is the main class for Metagenomics functionality of ADToolbox. This class contains all the methods required for metagenomics analysis 
    that ADToolbox offers.
    """
    def __init__(self,config:configs.Metagenomics)->None:
        """In order to instntiate an object from this class, you need to provide a metagenomics configs object from the configs module : configs.Metagenomics.
        Information for inputs and of each method is then obtained from the corresponding configs object. The following example shows how to instantiate an object from this class
        using the default configs object:
        
        Examples:
            >>> from adtoolbox import core, configs
            >>> config=configs.Metagenomics() ### This uses default arguments. Refer to configs module for more information.
            >>> metagenomics=core.Metagenomics(config)
        
        Args:
            config (configs.Metagenomics): A metagenomics configs object from configs module.
        
        Returns:
            None
        """
        self.config=config
            
    def find_top_taxa(
        self,
        sample_name:str,
        treshold:Union[int,float],
        mode:str='top_k',
        )->dict:
        """
        This function needs three inputs from qiime:
        1. feature table: This is the abundance of each feature in each sample (TSV).
        2. taxonomy table: This is the taxonomy of each feature (TSV). 
        3. rep seqs: This is the representative sequence of each feature (fasta).
        It then finds the top k features or features that form specific percentile of the community of the sample.
        
        Required Configs:
        
            config.feature_table_dir: The path to the feature table tsv file.
            ---------
            config.taxonomy_table_dir: The path to the taxonomy table tsv file.
            ---------
            config.rep_seq_fasta: The path to the representative sequence fasta file.
            ---------
        
        Args:
            sample_name (str): The name of the sample.
            threshold (int, float): The threshold for the top k or the percentile.
            mode (str, optional): Whether to find the top k features or features that form specific percentile of the community of the sample. Defaults to 'top_k'. Options: 'top_k', 'percentile'.
        
        Returns:
            dict: A dictionary of the top k features and their taxonomy.
        """
        ### Load all the required files
        feature_table = pd.read_table(self.config.feature_table_dir, sep='\t',skiprows=1)
        taxonomy_table = pd.read_table(self.config.taxonomy_table_dir, delimiter='\t')
        repseqs=fasta_to_dict(self.config.rep_seq_fasta)
        ### End Loading
        if mode == 'top_k':
            sorted_df=feature_table.sort_values(sample_name, ascending=False)
            top_featureids=list(sorted_df['#OTU ID'].head(treshold))
            top_taxa=[taxonomy_table[taxonomy_table['Feature ID']==featureid]['Taxon'].values[0] for featureid in top_featureids]
            top_repseqs=[repseqs[featureid] for featureid in top_featureids]
            top_abundances=list(sorted_df[sample_name].head(treshold)/sorted_df[sample_name].sum())
            
        elif mode == 'percentile':
            feature_table[sample_name]=feature_table[sample_name]/feature_table[sample_name].sum()
            sorted_df=feature_table.sort_values(sample_name, ascending=False)
            sorted_df['cumsum']=sorted_df[sample_name].cumsum()*100
            sorted_df_filtered=sorted_df[sorted_df['cumsum']<=treshold]
            top_featureids=list(sorted_df_filtered['#OTU ID'])
            top_taxa=[taxonomy_table[taxonomy_table['Feature ID']==featureid]['Taxon'].values[0] for featureid in top_featureids]
            top_repseqs=[repseqs[featureid] for featureid in top_featureids]
            top_abundances=sorted_df.loc[sorted_df_filtered.index][sample_name].values.tolist()
        else:
            raise ValueError("mode must be either 'top_k' or 'percentile'")
        
        return {'top_featureids':top_featureids,'top_taxa':top_taxa,'top_repseqs':top_repseqs,'top_abundances':top_abundances}    
        
    
    def align_to_gtdb(self,container:str="None")->str:
        """This function takes the representative sequences of the top k features and generates the script to
        align these feature sequences to gtdb using VSEARCH. If you intend to run this you either
        need to have VSEARCH installed or run it with a container option. You can use either the docker or singularity
        as container options. Otherwise you can use None and run it with the assumption that VSEARCH is installed.
        If you only want the script and not to run it, set run to False.

        Required Configs:
        
            config.top_repseq_dir: The path to fasta file including the repseqs to be aligned.
            ---------
            config.gtdb_dir_fasta: The path to the gtdb fasta database.
            ---------   
            config.align_to_gtdb_outputs_dir: The path to the directory where the outputs of this function will be savedThere are two outputs:
                1. matches.blast: This is the tabularized blast output of the alignment.
                2. Alignments: This is the raw alignment file.
            ---------
            config.vsearch_similarity: The similarity threshold for the alignment to be used by VSEARCH.
            ---------
            config.vsearch_threads: The number of threads to be used by VSEARCH.
            ---------
            config.adtoolbox_docker: The name of the docker image to be used by ADToolbox (Only if using Docker as container).
            ---------
            config.adtoolbox_singularity: The name of the singularity image to be used by ADToolbox (Only if using Singularity as container).
            ---------
        
        Args:
            container (str, optional): The container to use. Defaults to "None".
        
        Returns:
            str: The script that is supposed to be running later.
        """
        ### Load all the required files
        alignment_dir = os.path.join(self.config.align_to_gtdb_outputs_dir,'Alignments')
        match_table=os.path.join(self.config.align_to_gtdb_outputs_dir,'matches.blast')
        gtdb_dir_fasta=self.config.gtdb_dir_fasta
        ### End Loading
        query=self.config.top_repseq_dir
        dirs=[self.config.align_to_gtdb_outputs_dir,
            gtdb_dir_fasta,
            query
            ]
        for dir in dirs:
            if not pathlib.Path(dir).exists():
                os.mkdir(dir)
        if container=="None":
            bash_script=('vsearch --top_hits_only --blast6out '+
                        match_table+
                        ' --usearch_global '+ query +
                        ' --db '+ gtdb_dir_fasta +
                        ' --id ' +str(self.config.vsearch_similarity) +
                        ' --threads '+str(self.config.vsearch_threads)+
                        ' --alnout '+ alignment_dir +
                        ' --top_hits_only'+'\n')
        
        if container=="docker":
            bash_script='docker run -it '
            for dir in dirs:
                bash_script+=('-v '+dir+':'+dir+' ')
            
            bash_script += (self.config.adtoolbox_docker+' vsearch --top_hits_only --blast6out '+
                        match_table+
                        ' --usearch_global '+ query +
                        ' --db '+ gtdb_dir_fasta +
                        ' --id ' +str(self.config.vsearch_similarity) +
                        ' --threads '+str(self.config.vsearch_threads)+
                        ' --alnout '+ alignment_dir +
                        ' --top_hits_only'+'\n')
        
        if container=="singularity":
            bash_script='singularity exec '
            for dir in dirs:
                bash_script+=('-B '+str(dir)+':'+str(dir)+' ')
            
            bash_script += (self.config.adtoolbox_singularity+' vsearch --top_hits_only --blast6out '+
                        match_table+
                        ' --usearch_global '+ str(query) +
                        ' --db '+ gtdb_dir_fasta +
                        ' --id ' +str(self.config.vsearch_similarity) +
                        ' --threads '+str(self.config.vsearch_threads)+
                        ' --alnout '+ alignment_dir +
                        ' --top_hits_only'+'\n')
        return bash_script
    
    
    
    def get_genomes_from_gtdb_alignment(self,save:bool=True)->dict:
        """This function takes the alignment file generated from the align_to_gtdb function and generates the the genome information
        using the GTDB-Tk. In the outputted dictionary, the keys are feature ids and the values are the representative genomes.
        if you want to save the information as a json file set save to True.

        Required Configs:
            config.align_to_gtdb_outputs_dir: The path to the directory where the outputs of the align_to_gtdb function are saved.
            ---------
            config.feature_to_taxa: The path to the json file where the json file including feature ids and the representative genomes will be saved.
        
        Args:
            save (bool, optional): Whether to save the json file or not. Defaults to True.
        """
        matches = os.path.join(self.config.align_to_gtdb_outputs_dir,'matches.blast')
        aligned=pd.read_table(matches,header=None,delimiter='\t')
        aligned.drop_duplicates(0,inplace=True)
        aligned[1]=aligned[1].apply(lambda x: "".join(x.split('_')[1:]))
        alignment_dict=dict(zip(aligned[0],aligned[1]))
        if save:
            with open(self.config.feature_to_taxa, 'w') as f:
                json.dump(alignment_dict, f)
        
        return alignment_dict
    
    
    def download_genome(self,identifier:str,container:str="None")-> str:
        """This function downloads the genomes from NCBI using the refseq/genbank identifiers.
        If you want to save the json file that is holding the information for each genome, set save to True.
        Note that this function uses rsync to download the genomes. 

        Required Configs:
            config.genomes_base_dir: The path to the base directory where the genomes will be saved.
            ---------
            config.adtoolbox_docker: The name of the docker image to be used by ADToolbox (Only if using Docker as container).
            ---------
            config.adtoolbox_singularity: The name of the singularity image to be used by ADToolbox (Only if using Singularity as container).
            ---------
        Args:
            identifier list[str]: The list of identifiers for the genomes. It can be either refseq or genbank.
            container (str, optional): The container to use. Defaults to "None". You may select from "None", "docker", "singularity".
        
        Returns:
            str: The bash script that is used to download the genomes or to be used to download the genomes.

        """
        base_ncbi_dir = 'rsync://ftp.ncbi.nlm.nih.gov/genomes/all/'
        bash_script=""

        specific_ncbi_dir = identifier[0:3]+'/'+\
                            identifier[3:6]+'/'+\
                            identifier[6:9]+'/'+\
                            identifier[9:].split('.')[0]
            
        genome_dir=pathlib.Path(self.config.genome_save_dir(identifier))
            
        if container=="None":
            bash_script+=('\nrsync -avz --progress '+base_ncbi_dir+specific_ncbi_dir+' '+self.config.genome_save_dir(identifier))
            
        if container=="docker":
            bash_script+=('docker run -it -v '+str(genome_dir.parent)+':'+str(genome_dir.parent)+ f' {self.config.adtoolbox_docker} rsync -avz --progress '+' '+base_ncbi_dir+specific_ncbi_dir+' '+str(genome_dir))
            
        if container=="singularity":
            bash_script+=('singularity exec -B '+str(genome_dir.parent)+':'+str(genome_dir.parent)+ f' {self.config.adtoolbox_singularity} rsync -avz --progress '+' '+base_ncbi_dir+specific_ncbi_dir+' '+str(genome_dir))
        
        return bash_script
    
    def extract_genome_info(self,save:bool=True)->dict[str,str]:
        """This function extracts the genome information from the genomes base directory. If you want to save the genome information as a json file, set save to True.
        
        Required Configs:
            config.genomes_base_dir: The path to the base directory where the genomes are saved.
            ---------
        Args:
            genome_info (dict[str,str]): A dictionary containing the genome information.
            save (bool, optional): Whether to save the genome information as a json file. Defaults to True.
            container (str, optional): The container to use. Defaults to "None". You may select from "None", "docker", "singularity".
        
        Returns:
            dict[str,str]: A dictionary containing the address of the genomes that are downloaded or to be downloaded.
        """
        base_dir = pathlib.Path(self.config.genomes_base_dir)
        genome_info = {}
        for genome_dir in base_dir.iterdir():
            if genome_dir.is_dir():
                candids=list(genome_dir.rglob('*genomic.fna.gz'))
                for candid in candids:
                    if "cds" not in candid.name and "rna" not in candid.name:
                        genome_info[genome_dir.name]=str(candid.absolute())
        
        if save:
            with open(self.config.genomes_json_info, 'w') as f:
                json.dump(genome_info, f)
                
        return genome_info
     
    def align_genome_to_protein_db(
            self,
            address:str,
            name:str,
            container:str="None",
            )->tuple[str,str]:
        """
        This is a function that will align a genome to the Protein Database of the ADToolbox using mmseqs2.
        If you want to save the scripts, set save to True. Note that the alignment tables will be saved in any case.
        Note that this function uses mmseqs2 to align the genomes to the protein database. So, to run this function without
        any container you need to have mmseqs2 installed on your system. However, if you want to run this function with a container,
        you need to have the container installed on your system. You may select from "None", "docker", "singularity".

        Requires:
            config.genome_alignment_output: The path to the directory where the alignment results will be saved.
            ---------
            config.protein_db: The path to the ADToolbox protein database in fasta.
            ---------
            config.adtoolbox_docker: The name of the docker image to be used by ADToolbox (Only if using Docker as container).
            ---------
            config.adtoolbox_singularity: The name of the singularity image to be used by ADToolbox (Only if using Singularity as container).
            ---------
        Args:
            address (str): The address of the genome fasta file. The file must be in fasta format.
            run (bool, optional): Whether to run the alignment. Defaults to True.
            save (bool, optional): Whether to save the alignment scripts. Defaults to True.
            container (str, optional): The container to use. Defaults to "None". You may select from "None", "docker", "singularity".

        Returns:
            str: A dictionary containing the alignment files.
            str: The bash script that is used to align the genomes or to be used to align the genomes.
        """
 
        if container=="None":
            bash_script = ""
            alignment_file=os.path.join(self.config.genome_alignment_output,"Alignment_Results_mmseq_"+name+".tsv")
            bash_script += "mmseqs easy-search " + \
                address + " " + \
                self.config.protein_db + " " + \
                alignment_file+ ' tmp --format-mode 4 '+"\n\n"
        
        if container=="docker":
            bash_script = ""
            alignment_file=os.path.join(self.config.genome_alignment_output,"Alignment_Results_mmseq_"+name+".tsv")
            bash_script +="docker run -it "+ \
            " -v "+address+":"+address+ \
            " -v "+self.config.protein_db+":"+self.config.protein_db+ \
            " -v "+self.config.genome_alignment_output+":"+self.config.genome_alignment_output+ \
            f" {self.config.adtoolbox_docker}  mmseqs easy-search " + \
                address + " " + \
                self.config.protein_db + " " + \
                alignment_file+' tmpfiles --format-mode 4 '+"\n\n"

        if container=="singularity":
            bash_script = ""
            alignment_file=os.path.join(self.config.genome_alignment_output,"Alignment_Results_mmseq_"+name+".tsv")
            bash_script +="singularity exec "+ \
            " -B "+address+":"+address+ \
            " -B "+self.config.protein_db+":"+self.config.protein_db+ \
            " -B "+self.config.genome_alignment_output+":"+self.config.genome_alignment_output+ \
            f" {self.config.adtoolbox_singularity}  mmseqs easy-search " + \
                address + " " + \
                self.config.protein_db + " " + \
                alignment_file+' tmpfiles --format-mode 4 '+"\n\n"
        
        return  bash_script,alignment_file

    def align_short_reads_to_protein_db(self,query_seq:str,alignment_file_name:str,container:str="None",run:bool=True,save:bool=True)->tuple[str,str]:
        """This function aligns shotgun short reads to the protein database of the ADToolbox using mmseqs2.
        mmseqs wrappers in utils are used to perform this task. The result of this task is an alignment table.
        
        Required Configs:
        
            protein_db_mmseqs (str): The address of the existing/to be created protein database of the ADToolbox for mmseqs.
            --------
        Args:
            query_seq (str): The address of the query sequence.
            alignment_file_name (str): The name of the alignment file.
            container (str, optional): The container to use. Defaults to "None". You may select from "None", "docker", "singularity".
            run (bool, optional): Whether to run the alignment. Defaults to True.
            save (bool, optional): Whether to save the alignment scripts. Defaults to True.

        Returns:
            str: The bash script that is used to align the genomes or to be used to align the genomes.
            str: The address of the alignment file.
        """
        if not pathlib.Path(self.config.protein_db_mmseqs).exists():
            raise FileNotFoundError("""The protein database of the ADToolbox for mmseqs is not found. Please build it first
                                    using Database.build_mmseqs_database method.""")
        path_query=pathlib.Path(query_seq)
        script = ""
        script += create_mmseqs_database(query_seq,str(path_query.parent/path_query.name.split(".")[0]),container=container,save=None,run=False)+"\n"
        script += mmseqs_search(
            query_db=str(path_query.parent/path_query.name.split(".")[0]),
            target_db=self.config.protein_db_mmseqs,
            results_db=path_query.parent/alignment_file_name,
            run=False,
            save=None,
            container=container,
        )+"\n"
        script += mmseqs_result_db_to_tsv(
            query_db=str(path_query.parent/path_query.name.split(".")[0]),
            target_db=self.config.protein_db_mmseqs,
            results_db=path_query.parent/alignment_file_name,
            tsv_file=path_query.parent/(alignment_file_name+".tsv"),
            container=container,
            save=None,
            run=False,)+"\n"
        
        if save:
            with open(path_query.parent/"alignment_script.sh","w") as f:
                f.write(script)
        if run:
            subprocess.run(script,shell=True)
        
        return script,path_query.parent/(alignment_file_name+".tsv")
    
    def extract_ec_from_alignment(self,alignment_file:str,save:bool=True)->dict:
        """
        This function extracts the number of times an EC number is found in the alignment file when aligned to ADToolbox protein database.
        
        Required Configs:
            config.e_value: The e-value threshold for the filtering the alignment table.
            ---------
            config.bit_score: The bit score threshold for the filtering the alignment table.
            ---------
            config.ec_counts_from_alignment: The address of the json file that the results will be saved in.
            ---------
        Args:
            alignment_file (str): The address of the alignment file.
            save (bool, optional): Whether to save the results. Defaults to True.
        
        Returns:
            dict: A dictionary of EC numbers and their counts.

        """
        alignment_table = pd.read_table(alignment_file,sep='\t')
        alignment_table = alignment_table[(alignment_table['evalue']<self.config.e_value)&(alignment_table['bits']>self.config.bit_score)].drop_duplicates("query",keep="first")
        alignment_table["target"]=alignment_table["target"].apply(lambda x:x.split("|")[1])
        ec_counts=alignment_table["target"].value_counts().to_dict()
        
        if save:
            with open(self.config.ec_counts_from_alignment,"w") as f:
                json.dump(ec_counts,f)
        return ec_counts
    
    def get_cod_from_ec_counts(self,ec_counts:dict,save:bool=True)->dict:
        """This function takes a json file that comtains ec counts and converts it to ADM microbial agents counts.
        Required Configs:
            config.adm_mapping : A dictionary that maps ADM reactions to ADM microbial agents.
            ---------
            config.csv_reaction_db : The address of the reaction database of ADToolbox.
            ---------
            config.adm_cod_from_ec  : The address of the json file that the results will be saved in.
            ---------
        Args:
            save (bool, optional): Whether to save the output. Defaults to True.   
            ec_counts (dict): A dictionary containing the counts for each ec number.  
        Returns:
            dict: A dictionary containing the ADM microbial agents counts.
        """
        reaction_db = pd.read_table(self.config.csv_reaction_db, sep=',').drop_duplicates("EC_Numbers")
        reaction_db.set_index("EC_Numbers",inplace=True)
        adm_reactions_agents = {k:0 for k in self.config.adm_mapping.keys()}
        for ec in ec_counts.keys():
            l=reaction_db.loc[ec,"Modified_ADM_Reactions"].split("|")
            for adm_rxn in l: 
                adm_reactions_agents[adm_rxn]+=ec_counts[ec]
        adm_microbial_agents={}
        for k,v in self.config.adm_mapping.items():
            adm_microbial_agents[v]=adm_reactions_agents[k]
        if save:
            with open(self.config.adm_cod_from_ec,"w") as f:
                json.dump(adm_microbial_agents,f)
        return adm_microbial_agents
        
        
    @needs_repair
    def adm_from_alignment_json(self,adm_rxns,model="Modified_ADM_Reactions"):
        rt = SeedDB(reaction_db=self.config.seed_rxn_db)
        with open(self.config.genome_alignment_output_json) as f:
            json_report = json.load(f)
        reaction_db = pd.read_table(self.config.csv_reaction_db, sep=',')
        json_adm_output = {}
        for adm_reaction in track([Rxn for Rxn in adm_rxns],description="Finding Alignments to ADM reactions"):
            json_adm_output[adm_reaction] = {}
            for genome_id in json_report.keys():
                if os.path.exists(json_report[genome_id]['Alignment_File']):
                    json_adm_output[adm_reaction][json_report[genome_id]
                                                  ['NCBI_Name']] = {}
                    temp_tsv = pd.read_table(
                        json_report[genome_id]['Alignment_File'], sep='\t')
                    a_filter = (temp_tsv.iloc[:, -1] > self.config.bit_score) & (
                        temp_tsv.iloc[:, -2] < self.config.e_value)
                    temp_tsv = temp_tsv[a_filter]
                    ecs = [item.split('|')[1] for item in temp_tsv.iloc[:, 1]]
                    filter = (reaction_db['EC_Numbers'].isin(ecs)) & (reaction_db[model].str.contains(adm_reaction))
                    filtered_rxns = reaction_db[filter]
                    ecs = filtered_rxns['EC_Numbers'].tolist()

                    ec_dict = {}
                    for ec in ecs:
                        ec_dict[ec] = []
                        L = rt.instantiate_rxns(ec, "Multiple_Match")
                        [ec_dict[ec].append(
                            rxn.__str__()) for rxn in L if rxn.__str__() not in ec_dict[ec]]

                    json_adm_output[adm_reaction][json_report[genome_id]
                                                  ['NCBI_Name']] = ec_dict


        with open(self.config.genome_adm_map_json, 'w') as f:
            json.dump(json_adm_output, f)

        return json_adm_output
    
    def calculate_cod_portions_from_alignment_file(self,
                                                   alignment_file:str,
                                                   microbe_reaction_map:dict=configs.Database().adm_mapping)->dict:
        
        """
        This function calculates the portion of cod for degraders in an mmseqs alignment file.
        
        Args:
            alignment_file (str): The address of the alignment file.
            microbe_reaction_map (dict, optional): The dictionary that maps microbial agents to reactions in ADM. Defaults to configs.Database().adm_mapping.
        
        Returns:
            dict: A dictionary that maps microbial agents to their portion of cod(sum to one).
    
        """
        alignment_table = pd.read_table(alignment_file,sep='\t')
        alignment_table = alignment_table[(alignment_table['evalue']<self.config.e_value)&(alignment_table['bits']>self.config.bit_score)].drop_duplicates("query",keep="first")
        alignment_table["target"]=alignment_table["target"].apply(lambda x:x.split("|")[1])
        ec_counts=alignment_table["target"].value_counts().to_dict()
        reaction_db=pd.read_table(self.config.csv_reaction_db,delimiter=",")
        reaction_db.drop_duplicates(subset=['EC_Numbers'], keep='first',inplace=True)
        reaction_db.set_index("EC_Numbers",inplace=True)
        ec_to_adm=dict(zip(reaction_db.index,reaction_db["Modified_ADM_Reactions"]))
        cods={i:0 for i in microbe_reaction_map.keys()}
        for ec  in ec_counts.keys():
            for adm_rxn in ec_to_adm[ec].split("|"):
                cods[adm_rxn]+=ec_counts[ec]
                
        cods={microbe_reaction_map[k]:v/sum(cods.values()) for k,v in cods.items()}
        return cods
        
    
    def extract_relative_abundances(self,sample_names:Union[list[str],None]=None,save:bool=True)->dict:
        
        """
        This function will extract relative abundances for top k taxa from top k taxa table
        and feature table. It will return a dictionary with the top k taxa genomes as keys and the relative abundances as values

        """
        feature_table = pd.read_table(self.config.feature_table_dir,sep='\t',skiprows=1)
        if sample_names is None:
            sample_names = feature_table.columns[1:]
        relative_abundances={sample:[] for sample in sample_names}
        with open(os.path.join(self.config.feature_to_taxa)) as f:
            feature_genome_map = json.load(f)
        features=list(feature_genome_map.keys())
        genomes=[feature_genome_map[feature] for feature in features]
        
        for sample in sample_names:
            for feature in features:
                relative_abundances[sample].append(feature_table.loc[feature_table['#OTU ID']==feature,sample].item())

        abundances=pd.DataFrame(relative_abundances,index=genomes)
        rel_abunds=abundances/abundances.sum(axis=0)
        rel_abunds=rel_abunds.T.to_dict('index')

        if save:
            with open(self.config.genome_relative_abundances,'w') as f:
                json.dump(rel_abunds,f)
        return rel_abunds
    
    @needs_repair
    def calculate_microbial_portions(self,
                                     rel_abund_genome:dict,
                                     genome_alignment_output:dict,
                                     microbe_reaction_map:dict)-> dict:
        
        """This method calculates COD fraction of each microbial term in a model
        
        Args:

            microbe_reaction_map (dict): This dictionary must determine the association between
            microbial species in the model and the reactions that these species are involved
            in. In the case of Modified-ADM, this dictionary is the following:

            microbe_reaction_map={
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
                            "Uptake of lactate":"X_lac",} 
        
            genome_info: This dictionary is the output of the align_genomes method. It holds
            the information about the directory holding the result of aligning the genomes of
            interest and other useful information

            relative_abundances: This dictionary holds the relative abundace of each genome in the 
            community. The keys of this dictionary are genome_id similar to the output of align_genomes.
            NOTE: The relative abundances must be between 0,1
            For example: 
            relative_abundances={
                "Genome_1":0.56,
                "Genome_2":0.11,
                "Genome_3":0.22,
                "Genome_4":0.11
            }

        
        Returns:
            cod (dict): This dictionary includes the fraction of COD that belongs to each microbial species in the
            model. 
        """

        cod={}
        reaction_db=pd.read_table(self.config.csv_reaction_db,delimiter=",")
        reaction_db.drop_duplicates(subset=['EC_Numbers'], keep='first',inplace=True)
        reaction_db.set_index("EC_Numbers",inplace=True)
        model_species=list(set(microbe_reaction_map.values()))
        cod_portion=AdditiveDict([(i,0) for i in model_species])
        for genome_id in rel_abund_genome:
            try:
                temp_tsv = pd.read_table(
                        genome_alignment_output[genome_id], sep='\t',header=None)
            except:
                warn(f"The alignment file for genome {genome_id} can not be found. Please check for possible errors")
                continue
            
            filter = (temp_tsv.iloc[:, -1] > self.config.bit_score) & (
                        temp_tsv.iloc[:, -2] < self.config.e_value)
            temp_tsv = temp_tsv[filter]
            unique_ecs=list(set([i[1] for i in temp_tsv[1].str.split("|")]))
            cleaned_reaction_list=[]
            [cleaned_reaction_list.extend(map(lambda x:x.strip(" "),reaction_db.loc[ec,"Modified_ADM_Reactions"].split("|"))) for ec in unique_ecs if ec in reaction_db.index]
            pathway_counts=Counter(cleaned_reaction_list)
            pathway_counts.__delitem__("")
            SUM=sum([pathway_counts[key] for key in pathway_counts])
            for i in pathway_counts:
                pathway_counts[i]=pathway_counts[i]/SUM
            microbial_species= dict([(microbe_reaction_map[key],pathway_counts[key]) for key in pathway_counts.keys() ])
            cod_portion=cod_portion+AdditiveDict(microbial_species)*rel_abund_genome[genome_id]
            SUM=sum([cod_portion[key] for key in cod_portion])
            if SUM==0:
                for i in cod_portion:
                    
                    cod_portion[i]=0
            else:
                for i in cod_portion:
                    cod_portion[i]=cod_portion[i]/SUM
        return cod
    @needs_repair
    def seqs_from_sra(self,accession:str,save:bool=True,run:bool=True,container:str="None")-> tuple[str,str,dict]:
        """ 
        This method downloads the fastq files from the SRA database using the accession number of the project or run.
        The method uses the fasterq-dump tool to download the fastq files. The method also creates two bash scripts
        that can be used to download the fastq files. The first bash script is used to download the SRA files using the
        prefetch tool. The second bash script is used to convert the SRA files to fastq files using the fasterq-dump tool.
        depending on the prefered containerization tool, the bash scripts can be converted to singularity or docker containers commands.
        In the case of docker or singularity you do not need to install the SRA tools on your machine. However, you need to have the
        containerization tool installed on your machine.

        Requires:
            <R>project_accession</R>
            <R>Configs.Metagenomics</R>
        Satisfies:
            <S>sra_project_dir</S>
            <S>sra_seq_prefetch_bash_str</S>
            <S>sra_seq_fasterq_dump_bash_str</S>
            <S>manifest_dict</S>
            <if> save=True</if><S>sra_seq_prefetch_bash_file</S>
            <if> save=True</if><S>sra_seq_fasterq_dump_bash_file</S>
            <if> run=True</if><S>.fastq_files</S>
            
        Args:
            accession (str): The accession number of the SRA project or run
            save (bool, optional): If True, the  bash scripts will be saved in the SRA work directory. Defaults to True.
            run (bool, optional): If True, the bash scripts will be executed. Defaults to True. You must have the SRA toolkit installed in your system.
            container (str, optional): The containerization tool that will be used to run the bash scripts. Defaults to "None". Options are "None","docker","singularity"
        
        Returns:
            prefetch_script (str): The bash script that will be used to download the SRA files in python string format
            fasterq_dump_script (str): The bash script that will be used to convert the SRA files to fastq files in python string format
            sample_metadata (dict): A dictionary that contains the sample metadata
    
        """   
        if container=="None":
            prefetch_script=f"""#!/bin/bash
            prefetch {accession} -O {os.path.join(self.config.sra_work_dir(accession),"seqs")}"""

            fasterq_dump_script=f"""#!/bin/bash
            for i in $(ls ./seqs/) ; do
            fasterq-dump --split-files seqs/$i/*.sra -O seqs/$i/
            rm seqs/$i/*.sra
            done
            """
        
        elif container=="docker":
            warn("Docker is not supported yet")

        project_path=pathlib.Path(self.config.sra_work_dir(accession))
        if not project_path.exists():
            project_path.mkdir(parents=True)
        
        sample_metadata=utils.get_sample_metadata_from_accession(accession)

        if save:
            
            with open(project_path.joinpath("prefetch.sh"),"w") as f:
                f.write(prefetch_script)
            
            with open(project_path.joinpath("fasterq_dump.sh"),"w") as f:
                f.write(fasterq_dump_script)
                
        
        if run:
            subprocess.run(["bash",str(project_path.joinpath("prefetch.sh"))],cwd=str(project_path))
            subprocess.run(["bash",str(project_path.joinpath("fasterq_dump.sh"))],cwd=str(project_path))
            
            
        
        return prefetch_script,fasterq_dump_script,sample_metadata
    


    
    def run_qiime2_from_sra(self,query_seq:str,save:bool=True,run:bool=True,container:str='None') -> tuple[str,str]:
        """
        Required Configs:
        

        Args:
            query_seq (str): directory where the fastq files are located
            save (bool, optional): If True, the  bash scripts will be saved in the SRA work directory. Defaults to True.
            run (bool, optional): If True, the bash scripts will be executed. Defaults to True. You must have the SRA toolkit installed in your system.
            container (str, optional): If you want to run the qiime2 commands in a container, specify the container name here. Defaults to 'None'.
        Returns:
            qiime2_bash_str (str): The bash script that will be used to run qiime2 in python string format
            manifest (dict): The manifest file that will be used to run qiime2 in python dictionary format
    

        """
        
        
        seqs=pathlib.Path(query_seq)
        manifest_single={'sample-id':[],'absolute-filepath':[]}
        manifest_paired={'sample-id':[],'forward-absolute-filepath':[],'reverse-absolute-filepath':[]}  
        if len(list(seqs.glob("*.fastq"))) == 2:
            r1=list(seqs.glob("*_1.fastq"))[0]
            r2=list(seqs.glob("*_2.fastq"))[0]
            manifest_paired['sample-id'].append(seqs.name)
            manifest_paired['forward-absolute-filepath'].append(str(r1))
            manifest_paired['reverse-absolute-filepath'].append(str(r2))
            paired_end=True
        elif len(list(seqs.glob("*.fastq"))) == 1:
            r1=list(seqs.glob("*.fastq"))[0]
            manifest_single['sample-id'].append(seqs.name)
            manifest_single['absolute-filepath'].append(str(r1))
            paired_end=False
                    
        manifest=pd.DataFrame(manifest_single) if not paired_end else pd.DataFrame(manifest_paired)
  
        if paired_end:
            with open(self.config.qiime2_paired_end_bash_str,"r") as f:
                qiime2_bash_str=f.read()
        else:
            with open(self.config.qiime2_single_end_bash_str,"r") as f:
                qiime2_bash_str=f.read()

        manifest_dir=seqs.joinpath("manifest.tsv")

        if container=="None":
            qiime2_bash_str=qiime2_bash_str.replace("<manifest>",str(manifest_dir))
            qiime2_bash_str=qiime2_bash_str.replace("<qiime2_work_dir>",str(seqs))
            qiime2_bash_str=qiime2_bash_str.replace("<classifier>",str(self.config.qiime_classifier_db))
        
        elif container=="docker":
            qiime2_bash_str=qiime2_bash_str.splitlines()
            for idx,line in enumerate(qiime2_bash_str):
                line=line.lstrip()
                if line.startswith("qiime") or line.startswith("biom"):
                    qiime2_bash_str[idx]=f"docker run --env TMPDIR=/data/tmp -v {seqs}:/data/ -v {Path(self.config.qiime_classifier_db).parent}:/data/{Path(self.config.qiime_classifier_db).parent.name} -w /data  {self.config.qiime2_docker_image}"+" "+line
            qiime2_bash_str="\n".join(qiime2_bash_str)
            qiime2_bash_str=qiime2_bash_str.replace("<manifest>",str(manifest_dir.name))
            qiime2_bash_str=qiime2_bash_str.replace("<qiime2_work_dir>",str(seqs.name))
            qiime2_bash_str=qiime2_bash_str.replace("<classifier>",os.path.join(str(Path(self.config.qiime_classifier_db).parent.name),str(Path(self.config.qiime_classifier_db).name)))
            if not paired_end:
                manifest['absolute-filepath']=[str(pathlib.Path("/data")/seqs.name/pathlib.Path(x).parent.name/pathlib.Path(x).name) for x in manifest['absolute-filepath']]
            
            else:
                manifest['forward-absolute-filepath']=[str(pathlib.Path("/data")/seqs.name/pathlib.Path(x).parent.name/pathlib.Path(x).name) for x in manifest['forward-absolute-filepath']]
                manifest['reverse-absolute-filepath']=[str(pathlib.Path("/data")/seqs.name/pathlib.Path(x).parent.name/pathlib.Path(x).name) for x in manifest['reverse-absolute-filepath']]
        
        elif container=="singularity":
            qiime2_bash_str=qiime2_bash_str.splitlines()
            for idx,line in enumerate(qiime2_bash_str):
                line=line.lstrip()
                if line.startswith("qiime") or line.startswith("biom"):
                    qiime2_bash_str[idx]=f"singularity exec --bind  {str(seqs)}:{str(seqs)},$PWD:$PWD,{str(Path(self.config.qiime_classifier_db))}:{str(Path(self.config.qiime_classifier_db))},$SINGULARITY_TMPDIR:/tmp  {self.config.qiime2_singularity_image} " +line
            qiime2_bash_str="\n".join(qiime2_bash_str)
            qiime2_bash_str=qiime2_bash_str.replace("<manifest>",str(manifest_dir))
            qiime2_bash_str=qiime2_bash_str.replace("<qiime2_work_dir>",str(seqs))
            qiime2_bash_str=qiime2_bash_str.replace("<classifier>",str(Path(self.config.qiime_classifier_db)))

            # if not paired_end:
            #     manifest['absolute-filepath']=["/"+str(Path(seqs.name)/pathlib.Path(x).parent.name/pathlib.Path(x).name) for x in manifest['absolute-filepath']]
            # else:
            #     manifest['forward-absolute-filepath']=["/"+str(Path(seqs.name)/pathlib.Path(x).parent.name/pathlib.Path(x).name) for x in manifest['forward-absolute-filepath']]
            #     manifest['reverse-absolute-filepath']=["/"+str(Path(seqs.name)/pathlib.Path(x).parent.name/pathlib.Path(x).name) for x in manifest['reverse-absolute-filepath']]

        else:
            raise ValueError("Container must be None, singularity or docker")
        
        
        if save:
            pd.DataFrame(manifest).to_csv(seqs.joinpath("manifest.tsv"),sep='\t',index=False)
            with open(seqs.joinpath("qiime2.sh"),"w") as f:
                f.write(qiime2_bash_str)
        if run:
            subprocess.run(["bash",str(seqs.joinpath("qiime2.sh"))],cwd=str(seqs))
        
        return qiime2_bash_str,manifest
    
    def extract_taxonomy_features(self):
        pass

if __name__ == "__main__":
    db=Database(config=configs.Database(feed_db=os.path.join(Main_Dir,"feed_test_db.tsv")))
    db._initialize_feed_database()
    assert pd.read_csv(os.path.join(Main_Dir,"feed_test_db.tsv")).shape[0]==0

    
    







