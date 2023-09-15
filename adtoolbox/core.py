from distutils.log import warn
import subprocess
import os
from collections import UserDict
import pandas as pd
import time
import json
import numpy as np
import re
import requests
from requests.adapters import HTTPAdapter
import utils
import configs
from requests.packages.urllib3.util.retry import Retry
from requests.exceptions import Timeout
from bs4 import BeautifulSoup
from datetime import datetime
from pathlib import Path
from collections import Counter
from collections import namedtuple
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
    if there are specific initial concentrations for the ADM species, they can be passed as a dictionary to the initial_concentrations argument.
    reference is an optional argument that can be used to provide a reference for the experimental data. If using the database module 
    to query for Experiment objects you can query by name or reference or model_type. So, having a descriptive reference can be useful for querying as well.
    default model name is "modified-adm". This can be changed by passing a different model name to the model_name argument. This also helps with querying.
    
    Args:
        name (str): A unique name for the experiment.
        time (list): A list of time points in days.
        variables (list): A list of integers that represent the variables that are the index of the ADM species that we have concentration data for.
        data (list): A list of lists. Each list in the list must be a list of concentrations for each species at each time point.
        initial_concentrations (dict, optional): A dictionary of initial concentrations for the ADM species. Defaults to {}.
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
    initial_concentrations: dict[str,float] = dataclasses.field(default_factory=dict)
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
                "time":self.time,
                "variables":self.variables,
                "data":self.data.T.tolist(),
                "initial_concentrations":self.initial_concentrations,
                "reference":self.reference,
                "model_name":self.model_name}
    

    
    
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
        li_ch_pr=self.carbohydrates+self.lipids+self.proteins
        without_xi=100-self.xi
        self.ch_tss=self.carbohydrates/li_ch_pr*without_xi/100
        self.lip_tss=self.lipids/li_ch_pr*without_xi/100
        self.prot_tss=self.proteins/li_ch_pr*without_xi/100
        self.xi_tss=self.xi/100
        without_si=100-self.si
        self.ch_tds=self.carbohydrates/li_ch_pr*without_si/100
        self.lip_tds=self.lipids/li_ch_pr*without_si/100
        self.prot_tds=self.proteins/li_ch_pr*without_si/100
        self.si_tds=self.si/100
    
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

    def cod_calc(self,add_h:float=0,add_c:float=0,add_o:float=0)->float:
        """
        Calculates the conversion rates for g/l -> gCOD/l
        In some cases we would like to add extra atoms for COD calculations
        For example, model seed biochemistry database only uses acetate instead of acetic acid.
        The 1 hydrogen difference changes the COD conversion rate. For this reason we can add extra atoms to the formula
        to calculate the COD conversion rate without changing anything else.
        
        Args:
            add_h (float): The number of extra hydrogen atoms to add to the formula for COD calculation.
            add_c (float): The number of extra carbon atoms to add to the formula for COD calculation.
            add_o (float): The number of extra oxygen atoms to add to the formula for COD calculation.

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
            mw = self.data['mass']+add_h*1+add_c*12+add_o*16
            for atom in atoms:
                if re.search(atom+'\d*', self.data['formula']):
                    if len(re.search(atom+'\d*', self.data['formula']).group()[1:]) == 0:
                        contents[atom] = 1
                    else:
                        contents[atom] = int(
                            re.search(atom+'\d*', self.data['formula']).group()[1:])
                else:
                    contents[atom] = 0
            contents['H']+=add_h
            contents['C']+=add_c
            contents['O']+=add_o
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
        config (configs.SeedDB): An instance of the SeedDB class in the configs module. This class contains the information about the seed database.
    
    Examples:
        >>> seed_db=SeedDB(configs.SeedDB())
        >>> assert seed_db.compound_db==configs.SeedDB().compound_db
        >>> assert seed_db.reaction_db==configs.SeedDB().reaction_db

    """

    def __init__(self, config:configs.SeedDB) -> None:
        
        self.reaction_db = config.reaction_db
        self.compound_db = config.compound_db

    def instantiate_rxns(self, seed_id:str)->Reaction:
        """
        This method is used to instantiate reaction objects from the seed database.
        in order to instantiate a reaction object, you need to pass the seed identifier for that reaction.
        
        Args:
            seed_id (str): The seed identifier for the reaction.
    
        Returns:
            Reaction: An instance of the Reaction class.
            
        Required Configs:
            - config.reaction_db
        
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
        
        Required Configs:
            - config.compound_db
        
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
        
        Required Configs:
            - config.reaction_db
        
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


    def initialize_protein_db(self)->None:
        """This function intializes ADToolbox's protein database by creating an empty fasta file.
        Be careful, this will overwrite any existing file with the same name.
        Logically, this needs method needs config.protein_db to be defined.
        
        Required Configs:
            - config.protein_db
            --------
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==False # This is just to make sure that the following lines create the file
            >>> db=Database(config=configs.Database(protein_db=os.path.join(Main_Dir,"protein_test_db.fasta"))) # point to a test non-existing file
            >>> db.initialize_protein_db() # initialize the protein database
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==True # check if the file is created
            >>> os.remove(os.path.join(Main_Dir,"protein_test_db.fasta")) # remove the file to clean up
        """

        with open(self.config.protein_db, 'w') as f:
            pass
    
    def initialize_reaction_db(self)->None:
        r"""This function intializes ADToolbox's reaction database by creating an empty tsv file.
        Be careful, this will overwrite any existing file with the same name.
        
        Required Configs:
            - config.reaction_db
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"reaction_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(reaction_db=os.path.join(Main_Dir,"reaction_test_db.tsv")))
            >>> db.initialize_reaction_db()
            >>> assert pd.read_table(os.path.join(Main_Dir,"reaction_test_db.tsv"),delimiter="\t").shape[0]==0
            >>> assert set(pd.read_csv(os.path.join(Main_Dir,"reaction_test_db.tsv"),delimiter="\t").columns)==set(["ec_numbers","seed_ids","reaction_names","adm1_reaction","modified_adm_reactions","pathways"])
            >>> os.remove(os.path.join(Main_Dir,"reaction_test_db.tsv"))
        
        """
        pd.DataFrame(columns=["ec_numbers","seed_ids","reaction_names","adm1_reaction","modified_adm_reactions","pathways"]).to_csv(self.config.reaction_db,index=False,sep="\t")
        
    def initialize_feed_db(self)->None:
        r"""This function intializes ADToolbox's Feed database by creating an empty tsv file.
        Be careful, this will overwrite any existing file with the same name.
        
        Required Configs:
            - config.feed_db
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(feed_db=os.path.join(Main_Dir,"feed_test_db.tsv")))
            >>> db.initialize_feed_db()
            >>> assert pd.read_table(os.path.join(Main_Dir,"feed_test_db.tsv"),delimiter='\t').shape[0]==0
            >>> assert set(pd.read_table(os.path.join(Main_Dir,"feed_test_db.tsv"),delimiter='\t').columns)==set(["name","carbohydrates","lipids","proteins","tss","si","xi","reference"])
            >>> os.remove(os.path.join(Main_Dir,"feed_test_db.tsv"))
        
        """
        pd.DataFrame(columns=["name","carbohydrates","lipids","proteins","tss","si","xi","reference"]).to_csv(self.config.feed_db,index=False,sep="\t")
    
    def initialize_metagenomics_studies_db(self)->None:
        r"""This function intializes ADToolbox's Metagenomics studies database by creating an empty tsv file.
        Be careful, this will overwrite any existing file with the same name.
        
        Required Configs:
            - config.metagenomics_studies_db
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(metagenomics_studies_db=os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv")))
            >>> db.initialize_metagenomics_studies_db()
            >>> assert pd.read_table(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"),delimiter="\t").shape[0]==0
            >>> assert set(pd.read_table(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"),delimiter="\t").columns)==set(["name","study_type","microbiome","sample_accession","comments","study_accession"])
            >>> os.remove(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"))
        
        """
        pd.DataFrame(columns=["name","study_type","microbiome","sample_accession","comments","study_accession"]).to_csv(self.config.metagenomics_studies_db,index=False,sep="\t")
        
    def initialize_experimental_data_db(self)->None:
        """This function intializes ADToolbox's experimental data database by creating an empty json file.
        Be careful, this will overwrite any existing file with the same name.

        Required Configs:
            - config.experimental_data_db
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"experimental_data_test_db.json"))==False
            >>> db=Database(config=configs.Database(experimental_data_db=os.path.join(Main_Dir,"experimental_data_test_db.json")))
            >>> db.initialize_experimental_data_db()
            >>> assert pd.read_json(os.path.join(Main_Dir,"experimental_data_test_db.json")).shape[0]==0
            >>> with open(os.path.join(Main_Dir,"experimental_data_test_db.json"),"r") as f:
            ...     assert json.load(f)==[]
            >>> os.remove(os.path.join(Main_Dir,"experimental_data_test_db.json"))
        
        """
        pd.DataFrame(columns=["name","initial_conditions","time","variables","data","reference"]).to_json(self.config.experimental_data_db,orient="records")
        
    
    def filter_seed_from_ec(self, 
                            ec_list:list[str],
                            save:bool=False) -> tuple:
        """
        This function takes a list of EC numbers and filters the seed database to find the seed reactions that have the EC numbers in their EC number list.
        This will help to trim the large seed database to a smaller one that only contains the reactions that are relevant to the AD process.

        Args:
            ec_list (list[str]): A list of EC numbers.
            save (bool, optional): Whether to save the filtered seed database or not. Defaults to False.
        
        Returns:
            tuple: A tuple containing the filtered seed reaction database and the seed compound database, respectively.
        
        Required Configs:
        
            - config.reaction_db
            --------
            - config.compound_db
            --------
            - config.local_reaction_db
            --------
            - config.local_compound_db
            --------
            
            
        Examples:
            >>> db=Database()
            >>> seed_rxn_db,seed_compound_db=db.filter_seed_from_ec(["1.1.1.1","1.1.1.2"])
            >>> assert len(seed_rxn_db)>0 and len(seed_compound_db)>0
            >>> assert pd.read_json(configs.Database().reaction_db).shape[0]>pd.DataFrame(seed_rxn_db).shape[0]
        """
        seed_rxn_db=pd.read_json(self.config.reaction_db)
        seed_compound_db=pd.read_json(self.config.compound_db)
        seed_rxn_db=seed_rxn_db[seed_rxn_db["ec_numbers"].apply(lambda x: any(ec in x for ec in ec_list) if x else False)]
        seed_compound_db=seed_compound_db[seed_compound_db["id"].apply(lambda x: True if x in seed_rxn_db["stoichiometry"].sum() else False)]
        if save:
            seed_rxn_db.to_json(self.config.local_reaction_db)
            seed_compound_db.to_json(self.config.local_compound_db)
        return seed_rxn_db.to_dict(orient="record"),seed_compound_db.to_dict(orient="record")
        
            

    def get_protein_seqs_from_uniprot(self, uniprot_id:str) -> str:
        """
        This function takes a uniprot id and fetches the protein sequence from Uniprot.

        Args:
            uniprot_id (str): The uniprot id of the protein.
        
            
        Returns:
            str: The protein sequence.
        
        Examples:
            >>> db=Database()
            >>> seq=db.get_protein_seqs_from_uniprot("P0A9P0")
            >>> assert type(seq)==str and len(seq)>0
        """
        Base_URL = "https://rest.uniprot.org/uniprotkb/"
        session = requests.Session()
        retry = Retry(connect=3, backoff_factor=0.5)
        adapter = HTTPAdapter(max_retries=retry)
        session.mount('http://', adapter)
        try:
            file = session.get(
                f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta", timeout=10)
        except:
            print("Could not fetch the sequence! Trying again ...")
            while True:
                time.sleep(5)
                file = session.get(Base_URL+uniprot_id+".fasta", timeout=10)
                if file.ok:
                    break
           
        return ''.join(file.text.split('\n')[1:-1])
   
    def proteins_from_ec(self,ec_number:str) -> dict:
        """
        This function returns a dictionary of protein sequences for a given EC number.
        The keys are the uniprot ids and ec number compatible with ADToolbox protein database
        and the values are the protein sequences. Since ADToolbox deals with microbial process,
        only bacterial and archaeal proteins are considered.

        Args:
            ec_number (str): The EC number.
        
        Returns:
            dict: A dictionary of protein sequences.
            
        Examples:
            >>> db=Database()
            >>> protein_seqs=db.proteins_from_ec("1.1.1.1")
            >>> assert len(protein_seqs)>0
            >>> assert list(protein_seqs.keys())[0].split("|")[1]=="1.1.1.1"
        """
        session = requests.Session()
        retry = Retry(connect=3, backoff_factor=0.5)
        adapter = HTTPAdapter(max_retries=retry)
        session.mount('http://', adapter)
        protein_seqs={}
        try:
            file = session.get(
                f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28ec%3A{ec_number}%29%20AND%20%28reviewed%3Atrue%29%20NOT%20%28taxonomy_id%3A2759%29%29", timeout=30)
        except requests.exceptions.HTTPError or requests.exceptions.ConnectionError:
            print("Request Error! Trying again ...")
            time.sleep(30)
            file = session.get(
                f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28ec%3A{ec_number}%29%20AND%20%28reviewed%3Atrue%29%20NOT%20%28taxonomy_id%3A2759%29%29", timeout=30)
        # This alsp does a sanity chec
        except Exception:
            print('Something went wrong!')
        text = file.text
        if text:
            text=text.split('>')
            text.remove("")
            for seq in text:
                protein_seqs.update([(seq.split("\n")[0].split("|")[1]+"|"+ec_number, "".join(seq.split("\n")[1:]))])
                
        
        return protein_seqs


    def build_protein_db_from_reactions_db(self):
        r"""
        This function builds the protein database from the reaction database.
        It takes the reaction database and finds the protein sequences for each EC number in the reaction database.
        Then it saves the protein sequences in a fasta file.

        Required Configs:
            - config.reaction_db
            --------
            - config.protein_db
            --------
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==False
            >>> assert os.path.exists(os.path.join(Main_Dir,"reaction_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(protein_db=os.path.join(Main_Dir,"protein_test_db.fasta"),reaction_db=os.path.join(Main_Dir,"reaction_test_db.tsv")))
            >>> reaction_db=pd.DataFrame(columns=["EC_Numbers","Seed Ids","Reaction Names","ADM1_Reaction","Modified_ADM_Reactions","Pathways"])
            >>> reaction_db.loc[0,"EC_Numbers"]="1.1.1.1"
            >>> reaction_db.to_csv(os.path.join(Main_Dir,"reaction_test_db.tsv"),index=False,sep="\t")
            >>> db.build_protein_db_from_reactions_db()
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==True
            >>> assert os.path.exists(os.path.join(Main_Dir,"reaction_test_db.tsv"))==True
            >>> assert os.path.getsize(os.path.join(Main_Dir,"protein_test_db.fasta"))>0
            >>> os.remove(os.path.join(Main_Dir,"protein_test_db.fasta"))
            >>> os.remove(os.path.join(Main_Dir,"reaction_test_db.tsv"))
        """
        rxn_db=pd.read_table(self.config.reaction_db,delimiter="\t")
        ec_numbers=rxn_db["EC_Numbers"]
        ec_numbers=list(set(ec_numbers))
        protein_seqs={}
        for ec in ec_numbers:
            protein_seqs.update(self.proteins_from_ec(ec))
        with open(self.config.protein_db,"w") as f:
            for key,value in protein_seqs.items():
                f.write(">"+key+"\n")
                f.write(value+"\n")

    def cazy_ec(self)->list:
        """
        This method returns a list of EC numbers that are extracted from the Cazy website.
        This method is useful for adding more carbohydrate metabolism reactions to the reaction database.
        
        Returns:
            list: A list of EC numbers for carbohydrate metabolism found on CAZy database.
        
        Examples:
            >>> db=Database()
            >>> ec_list=db.cazy_ec()
            >>> assert len(ec_list)>0
        """

        ec_list = []
        for link in self.config.cazy_links:
            page = requests.get(link)
            soup = BeautifulSoup(page.content, "html.parser")
            results = soup.find("div", class_="cadre_principal").find_all(
                "th", class_="thec")
            for ec_number in results:
                if '-' not in ec_number.text.strip() and '.' in ec_number.text.strip():
                    ec_list.append(ec_number.text.strip())
                    
        return ec_list
          
    def add_protein_to_protein_db(self, protein_id:str, ec_number:str)->None:
        """
        This funciton adds a protein sequence to the protein database. It takes a uniprot id and an EC number it is assigned to 
        and adds the corresponding protein sequence to the protein database.
        
        Required Configs:
            - config.protein_db

        Args:
            protein_id (str): The uniprot id of the protein.
            ec_number (str): The EC number of the protein.
    
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==False
            >>> db=Database(config=configs.Database(protein_db=os.path.join(Main_Dir,"protein_test_db.fasta")))
            >>> db.add_protein_to_protein_db("P0A9P0","1.2.3.4")
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==True
            >>> assert os.path.getsize(os.path.join(Main_Dir,"protein_test_db.fasta"))>0
            >>> import utils
            >>> assert len(utils.fasta_to_dict(os.path.join(Main_Dir,"protein_test_db.fasta")))>0
            >>> os.remove(os.path.join(Main_Dir,"protein_test_db.fasta"))
        """
        if not os.path.exists(self.config.protein_db):
            self.initialize_protein_db()
        with open(self.config.protein_db,"a") as f:
            f.write(">"+protein_id+"|"+ec_number+"\n")
            f.write(self.get_protein_seqs_from_uniprot(protein_id)+"\n")
            
    def add_proteins_from_ecnumbers_to_protein_db(self, ec_numbers:list)->None:
        """
        This function adds protein sequences to the protein database from a list of EC numbers.
        It takes a list of EC numbers and finds the protein sequences for each EC number in the list.
        Then it saves the protein sequences in a fasta file.
        
        Required Configs:
            - config.protein_db
        
        Args:
            ec_numbers (list): A list of EC numbers.
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==False
            >>> db=Database(config=configs.Database(protein_db=os.path.join(Main_Dir,"protein_test_db.fasta")))
            >>> db.add_proteins_from_ecnumbers_to_protein_db(["1.1.1.1","1.1.1.2"])
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==True
            >>> import utils
            >>> assert len(utils.fasta_to_dict(os.path.join(Main_Dir,"protein_test_db.fasta")))>0
            >>> os.remove(os.path.join(Main_Dir,"protein_test_db.fasta"))
        """
        if not os.path.exists(self.config.protein_db):
            self.initialize_protein_db()
        
        protein_seqs={}
        for ec in ec_numbers:
            protein_seqs.update(self.proteins_from_ec(ec))
        
        with open(self.config.protein_db,"a") as f:
            for key,value in protein_seqs.items():
                f.write(">"+key+"\n")
                f.write(value+"\n")
        
    def add_feed_to_feed_db(self,feed:Feed)->None:
        r"""
        This function adds a feed to the feed database. It takes the feed name and the feed composition and adds them to the feed database.

        Required Configs:
            - config.feed_db

        Args:
            feed (Feed): An instance of the Feed class.
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(feed_db=os.path.join(Main_Dir,"feed_test_db.tsv")))
            >>> feed=Feed(name="test_feed",carbohydrates=10,lipids=20,proteins=30,tss=80,si=10,xi=30,reference="test")
            >>> db.add_feed_to_feed_db(feed)
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==True
            >>> assert pd.read_table(os.path.join(Main_Dir,"feed_test_db.tsv"),delimiter="\t").shape[0]>0
            >>> os.remove(os.path.join(Main_Dir,"feed_test_db.tsv"))
        
        """
        if not os.path.exists(self.config.feed_db):
            self.initialize_feed_db()
            
        if feed.name in pd.read_table(self.config.feed_db,delimiter="\t")["name"].values:
            raise ValueError("Feed already exists in the database.")
        feed_db=pd.read_table(self.config.feed_db,delimiter="\t")
        feed_db=pd.concat([feed_db,pd.DataFrame([feed.to_dict()])],ignore_index=True,axis=0)
        feed_db.to_csv(self.config.feed_db,index=False,sep="\t")
    
    def remove_feed_from_feed_db(self,field_name:str,query:str)->None:
        r"""
        This function removes studyes that contain the query in the given column, field name, from the feed database.

        Required Configs:
            - config.feed_db
        
        Args:
            field_name (str): The name of the column to query.
            query (str): The query string.
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(feed_db=os.path.join(Main_Dir,"feed_test_db.tsv")))
            >>> feed=Feed(name="test_feed",carbohydrates=10,lipids=20,proteins=30,tss=80,si=10,xi=30,reference="test")
            >>> db.add_feed_to_feed_db(feed)
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==True
            >>> assert pd.read_table(os.path.join(Main_Dir,"feed_test_db.tsv"),delimiter="\t").shape[0]>0
            >>> db.remove_feed_from_feed_db("name","test_feed")
            >>> assert pd.read_table(os.path.join(Main_Dir,"feed_test_db.tsv"),delimiter="\t").shape[0]==0
            >>> os.remove(os.path.join(Main_Dir,"feed_test_db.tsv"))
        
        """
        if not os.path.exists(self.config.feed_db):
            raise FileNotFoundError("Feed database does not exist!")
        
        
        feed_db=pd.read_table(self.config.feed_db,delimiter="\t")
        feed_db=feed_db[feed_db[field_name].str.contains(query)==False]
        feed_db.to_csv(self.config.feed_db,index=False,sep="\t")
        
    def get_feed_from_feed_db(self,field_name:str,query:str)->list[Feed]:
        r"""
        This function returns a feed from the feed database. It takes the query string and the column name to query and returns the feed that contains the query string in the given column.

        Required Configs:
            - config.feed_db
        
        Args:
            field_name (str): The name of the column to query.
            query (str): The query string.
        
        Returns:
            Feed: An instance of the Feed class.
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(feed_db=os.path.join(Main_Dir,"feed_test_db.tsv")))
            >>> feed=Feed(name="test_feed",carbohydrates=10,lipids=20,proteins=30,tss=80,si=10,xi=30,reference="test")
            >>> db.add_feed_to_feed_db(feed)
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==True
            >>> assert pd.read_table(os.path.join(Main_Dir,"feed_test_db.tsv"),delimiter="\t").shape[0]>0
            >>> feed=db.get_feed_from_feed_db("name","test_feed")
            >>> assert feed[0].name=="test_feed"
            >>> os.remove(os.path.join(Main_Dir,"feed_test_db.tsv"))
        
        """
        if not os.path.exists(self.config.feed_db):
            raise FileNotFoundError("Feed database does not exist!")
        
        feed_db=pd.read_table(self.config.feed_db,delimiter="\t")
        feed_db=feed_db[feed_db[field_name].str.contains(query)]
        return [Feed(**feed.to_dict()) for _,feed in feed_db.iterrows()]
    
    def add_metagenomics_study_to_metagenomics_studies_db(self,metagenomics_study:MetagenomicsStudy)->None:
        r"""
        This function adds a metagenomics study to the metagenomics studies database. It takes a metagenomics study and adds it to the metagenomics studies database.
        
        Required Configs:
            - config.metagenomics_studies_db
        
        Args:
            metagenomics_study (MetagenomicsStudy): An instance of the MetagenomicsStudy class.

        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(metagenomics_studies_db=os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv")))
            >>> metagenomics_study=MetagenomicsStudy(name="test_study",study_type="metagenomics",microbiome="anaerobic digester",sample_accession="test",comments="test",study_accession="test")
            >>> db.add_metagenomics_study_to_metagenomics_studies_db(metagenomics_study)
            >>> assert os.path.exists(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"))==True
            >>> assert pd.read_table(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"),delimiter="\t").shape[0]>0
            >>> os.remove(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"))
        """
        if not os.path.exists(self.config.metagenomics_studies_db):
            self.initialize_metagenomics_studies_db()
        metagenomics_studies_db=pd.read_table(self.config.metagenomics_studies_db,delimiter="\t")
        metagenomics_studies_db=pd.concat([metagenomics_studies_db,pd.DataFrame([metagenomics_study.to_dict()])],ignore_index=True,axis=0)
        metagenomics_studies_db.to_csv(self.config.metagenomics_studies_db,index=False,sep="\t")
    
    def remove_metagenomics_study_from_metagenomics_studies_db(self,field_name:str,query:str)->None:
        r"""
        This function removes studies that contain the query in the given column, field name, from the metagenomics studies database.

        Required Configs:
            - config.metagenomics_studies_db

        Args:
            field_name (str): The name of the column to query.
            query (str): The query string.

        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(metagenomics_studies_db=os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv")))
            >>> metagenomics_study=MetagenomicsStudy(name="test_study",study_type="metagenomics",microbiome="anaerobic digester",sample_accession="test",comments="test",study_accession="test")
            >>> db.add_metagenomics_study_to_metagenomics_studies_db(metagenomics_study)
            >>> assert os.path.exists(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"))==True
            >>> assert pd.read_table(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"),delimiter="\t").shape[0]>0
            >>> db.remove_metagenomics_study_from_metagenomics_studies_db("name","test_study")
            >>> assert pd.read_table(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"),delimiter="\t").shape[0]==0
            >>> os.remove(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"))
        """
        if not os.path.exists(self.config.metagenomics_studies_db):
            raise FileNotFoundError("Metagenomics studies database does not exist!")

        metagenomics_studies_db=pd.read_table(self.config.metagenomics_studies_db,delimiter="\t")
        metagenomics_studies_db=metagenomics_studies_db[metagenomics_studies_db[field_name].str.contains(query)==False]
        metagenomics_studies_db.to_csv(self.config.metagenomics_studies_db,index=False,sep="\t")
    
    def get_metagenomics_study_from_metagenomics_studies_db(self,field_name:str,query:str)->list[MetagenomicsStudy]:
        r"""
        This function returns a metagenomics study from the metagenomics studies database. It takes the query string and the column name to query and returns the metagenomics study that contains the query string in the given column.

        Required Configs:
            - config.metagenomics_studies_db
        
        Args:
            field_name (str): The name of the column to query.
            query (str): The query string.
        
        Returns:
            MetagenomicsStudy: An instance of the MetagenomicsStudy class.
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(metagenomics_studies_db=os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv")))
            >>> metagenomics_study=MetagenomicsStudy(name="test_study",study_type="metagenomics",microbiome="anaerobic digester",sample_accession="test",comments="test",study_accession="test")
            >>> db.add_metagenomics_study_to_metagenomics_studies_db(metagenomics_study)
            >>> assert os.path.exists(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"))==True
            >>> assert pd.read_table(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"),delimiter="\t").shape[0]>0
            >>> metagenomics_study=db.get_metagenomics_study_from_metagenomics_studies_db("name","test_study")
            >>> assert metagenomics_study[0].name=="test_study"
            >>> os.remove(os.path.join(Main_Dir,"metagenomics_studies_test_db.tsv"))
        """
        if not os.path.exists(self.config.metagenomics_studies_db):
            raise FileNotFoundError("Metagenomics studies database does not exist!")

        metagenomics_studies_db=pd.read_table(self.config.metagenomics_studies_db,delimiter="\t")
        metagenomics_studies_db=metagenomics_studies_db[metagenomics_studies_db[field_name].str.contains(query)]
        return [MetagenomicsStudy(**metagenomics_study.to_dict()) for _,metagenomics_study in metagenomics_studies_db.iterrows()]
    
    def add_experiment_to_experiments_db(self,experiment:Experiment)->None:
        r"""
        This function adds an experiment to the experiments database. It takes an experiment and adds it to the experiments database.
        
        Required Configs:
            - config.experimental_data_db
        
        Args:
            experiment (Experiment): An instance of the Experiment class.
        
        Examples:
            >>> import os,json
            >>> assert os.path.exists(os.path.join(Main_Dir,"experiments_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(experimental_data_db=os.path.join(Main_Dir,"experiments_test_db.json")))
            >>> experiment=Experiment(name="test_study",time=[0,1,2],variables=[2,6],data= [[1,2,3],[4,5,6]],reference="test")
            >>> db.add_experiment_to_experiments_db(experiment)
            >>> assert os.path.exists(os.path.join(Main_Dir,"experiments_test_db.json"))==True
            >>> assert os.path.getsize(os.path.join(Main_Dir,"experiments_test_db.json"))>0
            >>> os.remove(os.path.join(Main_Dir,"experiments_test_db.json"))
        """
        if not os.path.exists(self.config.experimental_data_db):
            self.initialize_experimental_data_db()
        
        if experiment.name in [experiment.name for experiment in self.get_experiment_from_experiments_db("name",experiment.name)]: 
            raise ValueError("Experiment already exists in the database!")
        
        with open(self.config.experimental_data_db,"r") as f:
            experiments_db=json.load(f)
        experiments_db.append(experiment.to_dict())
        with open(self.config.experimental_data_db,"w") as f:
            json.dump(experiments_db,f)
        
    def remove_experiment_from_experiments_db(self,field_name:str,query:str)->None:
        r"""
        This function removes experiments that contain the query in the given column, field name, from the experiments database.

        Required Configs:
            - config.experimental_data_db

        Args:
            field_name (str): The name of the column to query.
            query (str): The query string.

        Examples:
            >>> import os,json
            >>> assert os.path.exists(os.path.join(Main_Dir,"experiments_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(experimental_data_db=os.path.join(Main_Dir,"experiments_test_db.json")))
            >>> experiment=Experiment(name="test_study",time=[0,1,2],variables=[2,6],data= [[1,2,3],[4,5,6]],reference="test")
            >>> db.add_experiment_to_experiments_db(experiment)
            >>> assert os.path.exists(os.path.join(Main_Dir,"experiments_test_db.json"))==True
            >>> assert os.path.getsize(os.path.join(Main_Dir,"experiments_test_db.json"))>0
            >>> db.remove_experiment_from_experiments_db("name","test_study")
            >>> assert pd.read_json(os.path.join(Main_Dir,"experiments_test_db.json")).shape[0]==0
            >>> os.remove(os.path.join(Main_Dir,"experiments_test_db.json"))
        """
        if not os.path.exists(self.config.experimental_data_db):
            raise FileNotFoundError("Experimental data database does not exist!")

        with open(self.config.experimental_data_db,"r") as f:
            experiments_db=json.load(f)
        experiments_db=[experiment for experiment in experiments_db if query not in experiment[field_name]]
        with open(self.config.experimental_data_db,"w") as f:
            json.dump(experiments_db,f)

    def get_experiment_from_experiments_db(self,field_name:str,query:str)->list[Experiment]:
        r"""
        This function returns an experiment from the experiments database. It takes the query string and the column name to query and returns the experiment that contains the query string in the given column.

        Required Configs:
            - config.experimental_data_db
        
        Args:
            field_name (str): The name of the column to query.
            query (str): The query string.
        
        Returns:
            Experiment: An instance of the Experiment class.
        
        Examples:
            >>> import os,json
            >>> assert os.path.exists(os.path.join(Main_Dir,"experiments_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(experimental_data_db=os.path.join(Main_Dir,"experiments_test_db.json")))
            >>> experiment=Experiment(name="test_study",time=[0,1,2],variables=[2,6],data= [[1,2,3],[4,5,6]],reference="test")
            >>> db.add_experiment_to_experiments_db(experiment)
            >>> assert os.path.exists(os.path.join(Main_Dir,"experiments_test_db.json"))==True
            >>> assert os.path.getsize(os.path.join(Main_Dir,"experiments_test_db.json"))>0
            >>> experiment=db.get_experiment_from_experiments_db("name","test_study")
            >>> assert experiment[0].name=="test_study"
            >>> os.remove(os.path.join(Main_Dir,"experiments_test_db.json"))
        """
        if not os.path.exists(self.config.experimental_data_db):
            raise FileNotFoundError("Experimental data database does not exist!")

        with open(self.config.experimental_data_db,"r") as f:
            experiments_db=json.load(f)
        experiments_db=[experiment for experiment in experiments_db if query in experiment[field_name]]
        return [Experiment(**experiment) for experiment in experiments_db]
        
    def build_mmseqs_database(self,container:str="None")->str:
        """Builds an indexed mmseqs database from the ADToolbox's fasta protein database.
        
        Required Configs:
            - config.protein_db
            - config.adtoolbox_singularity
            - config.adtoolbox_docker
            
        Args:
            container (str, optional): The container to run the script with. Defaults to "None".
        Returns:
            str: The script to build the mmseqs database.
            
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==False
            >>> db=Database(config=configs.Database(protein_db=os.path.join(Main_Dir,"protein_test_db.fasta")))
            >>> db.add_protein_to_protein_db("P0A9P0","x,x,x,x")
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==True
            >>> assert os.path.getsize(os.path.join(Main_Dir,"protein_test_db.fasta"))>0
            >>> script=db.build_mmseqs_database()
            >>> assert script=="mmseqs createdb "+str(os.path.join(Main_Dir,"protein_test_db.fasta"))+" "+str(db.config.protein_db_mmseqs)
            >>> os.remove(os.path.join(Main_Dir,"protein_test_db.fasta"))
        
        """
        script=create_mmseqs_database(self.config.protein_db,
                                      self.config.protein_db_mmseqs,
                                      container=container,
                                      run=False,
                                      config=self.config)
        
        if container=="None":
            pass
        
        elif container=="singularity":
            script=f"singularity exec --bind {self.config.protein_db}:{self.config.protein_db},{self.config.protein_db_mmseqs}:{self.config.protein_db_mmseqs} {self.config.adtoolbox_singularity} {script}"
        
        elif container=="docker":
            script=f"docker run -v {self.config.protein_db}:{self.config.protein_db} -v {self.config.protein_db_mmseqs}:{self.config.protein_db_mmseqs} {self.config.adtoolbox_docker} {script}"
        
        else:
            raise ValueError("Container should be either None, singularity or docker!")
    
        return script


    def download_adm_parameters(self,verbose:bool=True)->None:
        """
        Downloads the parameters needed for running ADM models in ADToolbox.
        
        Required Configs:
            - config.adm_parameters_base_dir
            - config.adm_parameters_urls
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"adm_parameters_test"))==False
            >>> db=Database(config=configs.Database(adm_parameters_base_dir=os.path.join(Main_Dir,"adm_parameters_test")))
            >>> db.download_adm_parameters(verbose=False) 
            >>> assert os.path.exists(os.path.join(Main_Dir,"adm_parameters_test"))==True
            >>> assert len(os.listdir(os.path.join(Main_Dir,"adm_parameters_test")))==12
            >>> os.system("rm -r "+os.path.join(Main_Dir,"adm_parameters_test"))
            0
        
        Args:
        
            verbose (bool, optional): Whether to print the progress or not. Defaults to True.
        
        
        """
        if not os.path.exists(self.config.adm_parameters_base_dir):
                os.makedirs(self.config.adm_parameters_base_dir)
        if verbose:
            for url in track(self.config.adm_parameters_urls.values(),description="Downloading ADM parameters"):
                r = requests.get(url, allow_redirects=True)
                with open(os.path.join(self.config.adm_parameters_base_dir,url.split("/")[-1]), 'wb') as f:
                    f.write(r.content)
            rich.print(f"[green]ADM parameters were downloaded to {self.config.adm_parameters_base_dir}")
        
        else:
            for url in self.config.adm_parameters_urls.values():
                r = requests.get(url, allow_redirects=True)
                with open(os.path.join(self.config.adm_parameters_base_dir,url.split("/")[-1]), 'wb') as f:
                    f.write(r.content)
        
    def download_seed_databases(self,verbose:bool=True) -> None :
        """This function will download the seed databases, both compound and reaction databases.

        Required Configs:
            - config.seed_rxn_url
            - config.seed_compound_url
            - config.reaction_db
            - config.compound_db

        Args:
            verbose (bool, optional): Whether to print the progress or not. Defaults to True.
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"seed_rxn.json"))==False
            >>> assert os.path.exists(os.path.join(Main_Dir,"seed_compound.json"))==False
            >>> db=Database(config=configs.Database(reaction_db=os.path.join(Main_Dir,"seed_rxn.json"),compound_db=os.path.join(Main_Dir,"seed_compound.json")))
            >>> db.download_seed_databases(verbose=False)
            >>> assert os.path.exists(os.path.join(Main_Dir,"seed_rxn.json"))==True
            >>> assert os.path.exists(os.path.join(Main_Dir,"seed_compound.json"))==True
            >>> os.remove(os.path.join(Main_Dir,"seed_rxn.json"))
            >>> os.remove(os.path.join(Main_Dir,"seed_compound.json"))
        """
        r = requests.get(self.config.seed_rxn_url, allow_redirects=True,stream=True)
        with open(self.config.reaction_db, 'wb') as f:
            f.write(r.content)
        if verbose:
            rich.print(f"[green]Reaction database downloaded to {self.config.reaction_db}")
        r=requests.get(self.config.seed_compound_url,allow_redirects=True,stream=True)
        with open(self.config.compound_db, 'wb') as f:
            f.write(r.content)
        if verbose:
            rich.print(f"[green]Compound database downloaded to {self.config.compound_db}")

    def download_protein_database(self, verbose:bool=True) -> None:
        """
        Downloads the prebuilt protein database from the remote repository.
        
        Required Configs:
            - config.protein_db_url
            - config.protein_db
        
        Args:
            verbose (bool, optional): Whether to print the progress or not. Defaults to True.
            
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==False
            >>> db=Database(config=configs.Database(protein_db=os.path.join(Main_Dir,"protein_test_db.fasta")))
            >>> db.download_protein_database(verbose=False)
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==True
            >>> assert os.path.getsize(os.path.join(Main_Dir,"protein_test_db.fasta"))>0
            >>> os.remove(os.path.join(Main_Dir,"protein_test_db.fasta"))
        """
        r = requests.get(self.config.protein_db_url, allow_redirects=True)
        with open(self.config.protein_db, 'wb') as f:
            f.write(r.content)
        if verbose:
            rich.print(f"[green]Protein database downloaded to {self.config.protein_db}")
        
    def download_reaction_database(self,verbose:bool=True)->None:
        """
        This function will download the reaction database from the remote repository.
        
        Required Configs:
            - config.adtoolbox_rxn_db_url
            - config.csv_reaction_db
        
        Args:
            verbose (bool, optional): Whether to print the progress or not. Defaults to True.

        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"reaction_test_db.csv"))==False
            >>> db=Database(config=configs.Database(csv_reaction_db=os.path.join(Main_Dir,"reaction_test_db.csv")))
            >>> db.download_reaction_database(verbose=False)
            >>> assert os.path.exists(os.path.join(Main_Dir,"reaction_test_db.csv"))==True
            >>> assert os.path.getsize(os.path.join(Main_Dir,"reaction_test_db.csv"))>0
            >>> os.remove(os.path.join(Main_Dir,"reaction_test_db.csv"))
        """
    
        r = requests.get(self.config.adtoolbox_rxn_db_url, allow_redirects=True)
        with open(self.config.csv_reaction_db, 'wb') as f:
            f.write(r.content)
        if verbose:
            rich.print(f"[green]Reaction database downloaded to {self.config.csv_reaction_db}")

    
    def download_feed_database(self,verbose:bool=True)-> None:
        """
        This function will download the feed database from the remote repository.

        Required Configs:
            - config.feed_db_url
            - config.feed_db
        
        Args:
            verbose (bool, optional): Whether to print the progress or not. Defaults to True.

        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(feed_db=os.path.join(Main_Dir,"feed_test_db.tsv")))
            >>> db.download_feed_database(verbose=False)
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==True
            >>> assert os.path.getsize(os.path.join(Main_Dir,"feed_test_db.tsv"))>0
            >>> os.remove(os.path.join(Main_Dir,"feed_test_db.tsv"))
        """
        r = requests.get(self.config.feed_db_url, allow_redirects=True)
        with open(self.config.feed_db, 'wb') as f:
            f.write(r.content)
        if verbose:
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
            
    def download_studies_database(self,verbode:bool=True)->None:
        """
        This function will download the required files for studies functionality.

        Args:
            verbode (bool, optional): Whether to print the progress or not. Defaults to True.
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"studies_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(studies_db=os.path.join(Main_Dir,"studies_test_db.tsv")))
            >>> db.download_studies_database(verbose=False)
            >>> assert os.path.exists(os.path.join(Main_Dir,"studies_test_db.tsv"))==True
            >>> assert os.path.getsize(os.path.join(Main_Dir,"studies_test_db.tsv"))>0
            >>> os.remove(os.path.join(Main_Dir,"studies_test_db.tsv"))
        """
        for i in self.config.studies.urls:
            r = requests.get(self.config.studies.urls[i], allow_redirects=True)
            if not os.path.exists(self.config.studies.base_dir):
                os.makedirs(self.config.studies.base_dir,exist_ok=True)
            with open(os.path.join(self.config.studies.base_dir,self.config.studies.urls[i].split("/")[-1]), 'wb') as f:
                f.write(r.content)
            if verbode:
                rich.print(f"[bold green]Downloaded {self.config.studies.urls[i]}[/bold green]")
    
    def download_amplicon_to_genome_db(self,verbose:bool=True):
        """
        This function will automatically download the GTDB-tk database for genome assignment.
        
        Required Configs:
            - config.amplicon_to_genome_db
            - config.amplicon_to_genome_urls
        
        Args:
            verbose (bool, optional): Whether to print the progress or not. Defaults to True.

        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"amplicon_to_genome_test_db"))==False
            >>> db=Database(config=configs.Database(amplicon_to_genome_db=os.path.join(Main_Dir,"amplicon_to_genome_test_db")))
            >>> db.download_amplicon_to_genome_db(verbose=False)
            >>> assert os.path.exists(os.path.join(Main_Dir,"amplicon_to_genome_test_db"))==True
            >>> assert len(os.listdir(os.path.join(Main_Dir,"amplicon_to_genome_test_db")))>0
            >>> os.system("rm -r "+os.path.join(Main_Dir,"amplicon_to_genome_test_db"))
            0
        """
        if not os.path.exists(self.config.amplicon_to_genome_db):
            os.mkdir(self.config.amplicon_to_genome_db)

        url = self.config.amplicon_to_genome_urls
        if verbose:
            for keys in ['Version', 'MD5SUM', 'FILE_DESCRIPTIONS']:
                with requests.get(url[keys], allow_redirects=True, stream=True) as r:
                    total_size = int(r.headers.get('content-length', 0))
                    block_size = 1024
                    with Progress() as progress:
                        task1 = progress.add_task("Downloading " + keys, total=total_size)
                        with open(os.path.join(self.config.amplicon_to_genome_db, keys), 'wb') as f:
                            for data in r.iter_content(block_size):
                                progress.update(task1, advance=len(data))
                                f.write(data)
            with requests.get(url['metadata_field_desc'], allow_redirects=True, stream=True) as r:
                total_size = int(r.headers.get('content-length', 0))
                block_size = 1024
                with Progress() as progress:
                    task1 = progress.add_task("Downloading metadata_field_desc.tsv", total=total_size)
                    with open(os.path.join(self.config.amplicon_to_genome_db, 'metadata_field_desc.tsv'), 'wb') as f:
                        for data in r.iter_content(block_size):
                            progress.update(task1, advance=len(data))
                            f.write(data)

            for keys in ['bac120_metadata', 'bac120_ssu']:
                with requests.get(url[keys], allow_redirects=True, stream=True) as r:
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
                with requests.get(url[keys], allow_redirects=True, stream=False) as r:
                    with open(os.path.join(self.config.amplicon_to_genome_db, keys), 'wb') as f:
                        f.write(r.content)
            with requests.get(url['metadata_field_desc'], allow_redirects=True, stream=False) as r:
                with open(os.path.join(self.config.amplicon_to_genome_db, 'metadata_field_desc.tsv'), 'wb') as f:
                    f.write(r.content)
            for keys in ['bac120_metadata', 'bac120_ssu']:
                with requests.get(url[keys], allow_redirects=True, stream=False) as r:
                    with open(os.path.join(self.config.amplicon_to_genome_db, url[keys].split("/")[-1]), 'wb') as f:
                        f.write(r.content)
                with tarfile.open(os.path.join(self.config.amplicon_to_genome_db, url[keys].split("/")[-1])) as f_in:
                    f_in.extractall(self.config.amplicon_to_genome_db)
        if verbose:
            rich.print("[bold green]Downloaded all the required files for Amplicon to Genome functionality.[/bold green]")
                    
                        

            
    def download_all_databases(self,verbose:bool=True)->None:
        """
        This function will download all the required databases for all the functionalities of ADToolbox.
        NOTE: each method that this function calls is individually tested so it is skipped from testing!

        Args:
            verbose (bool, optional): Whether to print the progress or not. Defaults to True.

        Required Configs:
            - config.adm_parameters_base_dir
            - config.adm_parameters_urls
            - config.seed_rxn_url
            - config.seed_compound_url
            - config.reaction_db
            - config.compound_db
            - config.protein_db_url
            - config.protein_db
            - config.adtoolbox_rxn_db_url
            - config.csv_reaction_db
            - config.feed_db_url
            - config.feed_db
            - config.amplicon_to_genome_db
            - config.amplicon_to_genome_urls
            - config.qiime_classifier_db_url
            - config.qiime_classifier_db
            - config.studies_db
            - config.studies_urls
            
        Examples:
            >>> import os # doctest: +SKIP
            >>> db=Database(config=configs.Database()) # doctest: +SKIP
            >>> db.download_all_databases(verbose=False) # doctest: +SKIP

        """

        self.download_seed_databases(verbose=verbose)
        self.download_adm_parameters(verbose=verbose)
        self.download_protein_database(verbose=verbose)
        self.download_reaction_database(verbose=verbose)
        self.download_feed_database(verbose=verbose)
        self.download_studies_database(verbose=verbose)
        self.download_amplicon_to_genome_db(verbose=verbose)
        self.download_qiime_classifier_db(verbose=verbose)
        

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
        aligned[1]=aligned[1].apply(lambda x: ("".join(x.split('_')[1:])).split("~")[0])
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
    pass

    
    







