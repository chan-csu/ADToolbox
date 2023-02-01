from distutils.log import warn
import pdb
import subprocess
import os
import random
from matplotlib import streamplot
import pandas as pd
import pdb
import json
import numpy as np
import re
import requests
import time
import warnings
from requests.adapters import HTTPAdapter
from sympy import Li
import configs
from requests.packages.urllib3.util.retry import Retry
from requests.exceptions import Timeout
from bs4 import BeautifulSoup
from datetime import datetime
from pathlib import Path
from collections import Counter
import pathlib
import gzip
import shutil
import configs
from rich.progress import track,Progress
import rich
from __init__ import Main_Dir
from typing import Union
from utils import wrap_for_slurm,fasta_to_dict
# import doctest
# doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)

class Additive_Dict(dict):
    """ A simple dictionary with additive properties. The only difference with 
    normal dictionary is that "+" operator works with the objects of this class.
    when the values are the same between two dictionaries, their values gets added,
    if not, nothing happens!
    
    """
    def __add__(self,b):
        out=self.copy()
        for key in b.keys():
            if key in self.keys():
                out[key]+=b[key]
        return Additive_Dict(out)
    
    def __mul__(self,a):
        out=self.copy()
        for key in out.keys():
            out[key]=out[key]*a
        return Additive_Dict(out)

class Pipeline:
    """A class for running a chain of commands in series with parallelism."""
    def __init__(self,commands:list[tuple],executer:str,**kwargs):
        self.commands=commands
        self.executer=executer
        self.kbase_config=kwargs.get("kbase_config",configs.Kbase())
        self.database_config=kwargs.get("database_config",configs.Database())
        self.reaction_toolkit_config=kwargs.get("reaction_toolkit_config",configs.Reaction_Toolkit())
        self.metagenomics_config=kwargs.get("metagenomics_config",configs.Metagenomics())
        self.alignment_config=kwargs.get("alignment_config",configs.Alignment())
        self.original_adm_config=kwargs.get("original_adm_config",configs.Original_ADM1())
        self.modified_adm_config=kwargs.get("modified_adm_config",configs.Modified_ADM())
        
    
    def validate(self,verbose:bool=True):
        """Validates the pipeline commands. Runs the following checks:
        1. Input to the first command is correct."""
    
        satisfied=["Configs.Metagenomics",
                    "Configs.Database",
                    "Configs.Reaction_Toolkit",
                    "Configs.Kbase",
                    "Configs.Alignment",
                    ]
        requirements=[]
        failure={
            com.__name__ :[] for com in self.commands
         }
        success=[]
        for command in self.commands:
            satisfied.extend(Pipeline.extract_reqs_sats(command)[1])
            requirements=Pipeline.extract_reqs_sats(command)[0]
            for req in requirements:
                if req not in satisfied:
                    failure[command.__name__].append(req)
            if len(failure[command.__name__])==0:
                success.append(command)
        if verbose:
            rich.print("[bold green]The following commands are valid:[/bold green]")
            for com in success:
                failure.__delitem__(com.__name__)
                print(com.__name__)
            rich.print("[bold red]The following commands are invalid:[/bold red]")
            for com in failure.keys():
                print(f"{com} : {failure[com]}")


        return success,failure

            
        






    def run(self):
        pass

    @staticmethod
    def extract_reqs_sats(command):
        """Extracts the requirements and satisfactions from a command using the tags in the docstrings"""
        
        if "Requires:" in command.__doc__ and "Satisfies:" in command.__doc__:
            reqs=re.findall("<R>.+</R>",command.__doc__)
            for ind,i in enumerate(reqs):
                reqs[ind]=i.replace("<R>","").replace("</R>","")
            sats=re.findall("<S>.+</S>",command.__doc__)
            for ind,i in enumerate(sats):
                sats[ind]=i.replace("<S>","").replace("</S>","")
        else:
            raise Exception(f"The command docstring for {command.__name__} is not properly formatted. Please use the following tags: <R> and <S> for requirements and satisfactions respectively.")

        return reqs,sats

    

class Feed:

    """

    This class maps real feed data to ADM Feed parameters.
    Later I move some these to a markdown file.

    The definitions are based in the following link:

    https://extension.okstate.edu/fact-sheets/solids-content-of-wastewater-and-manure.html


    Carbohydrates <=> the fraction of carbohydrates in the feed
    Lipids <=> the fraction of lipids in the feed
    Proteins <=> the fraction of proteins in the feed
    TS <=> the total solids in the feed %mass
    TSS <=> the total soluble solids in the feed %mass
    'The portion of TS that remains after heating at 550
    C for 1 hour is called Total Fixed Solids (TFS);
    the portion lost during heating is Total Volatile Solids (TVS).'
    TS = TDS + TSS
    TS = TFS + TVS
    Sometimes the fixed solids content,TFS, is called the Ash Content.

    """

    def __init__(self, carbohydrates, lipids, proteins, ts, tss,si,pi):
        self.carbohydrates = carbohydrates
        self.lipids = lipids
        self.proteins = proteins
        self.ts = ts
        self.tss = tss
        self.tds = ts-tss
        self.si=si
        self.pi=pi

    def Report(self):
        return "Feed Characteristics:\nCarbohydrates: "+str(self.carbohydrates)+""+"\nLipids: "+str(self.lipids)+"\nProteins: "+str(self.proteins)+"\nFast Degrading portion: "+str(1-self.tds)+"\nSouluble Inerts: "+str(self.pi)+"\nParticulate Inerts: "+str(self.pi)


class Sequence_Toolkit:

    aa_list = ['A', 'M', 'N', 'V', 'W', 'L', 'H', 'S', 'G', 'F',
               'K', 'Q', 'E', 'S', 'P', 'I', 'C', 'Y', 'R', 'N', 'D', 'T']
    n_list = ['A', 'C', 'G', 'T']

    def __init__(self, type_: str='Protein', length: int=100):
        self.type = type_
        self.length = length

    def seq_generator(self)-> None:

        if self.type == 'Protein':
            return ''.join([random.choice(self.aa_list) for i in range(self.length)])

        elif self.type == 'DNA':
            return ''.join([random.choice(self.n_list) for i in range(self.length)])

        else:
            print('Type not recognized')

    def mutate_random(self, seq, number_of_mutations=10):
        if self.type == 'Protein':
            for i in range(number_of_mutations):
                seq = list(seq)
                seq[random.randint(0, len(seq))] = random.choice(self.aa_list)
            return ''.join(seq)

        elif self.type == 'DNA':
            for i in range(number_of_mutations):
                seq = list(seq)
                seq[random.randint(0, len(seq))] = random.choice(self.n_list)
            return ''.join(seq)


class Reaction_Toolkit:

    """
    A Class for converting ec numbers to reaction objects or extracting any information from the seed database

    """

    def __init__(self, compound_db:str =configs.Reaction_Toolkit().compound_db, reaction_db=configs.Reaction_Toolkit().reaction_db) -> None:
        self.reaction_db = reaction_db
        self.compound_db = compound_db

    def instantiate_rxns(self, ec_number, mode='Single_Match'):
        f = open(self.reaction_db)
        data = json.load(f)
        if mode == 'Single_Match':
            for i in data:
                if i['ec_numbers']:
                    if ec_number in i['ec_numbers']:
                        return Reaction(i)

        elif mode == 'Multiple_Match':
            matched_rxns = [
                Reaction(rxn) for rxn in data if rxn['ec_numbers'] == [ec_number]]
            return matched_rxns

        else:
            print('Mode not recognized!')

    def instantiate_metabs(self, seed_id):
        f = open(self.compound_db)
        data = json.load(f)
        # Maybe contain is better than == below\\To be completed
        matched_compounds = [Metabolite(met)
                            for met in data if met['id'] == seed_id]
        return matched_compounds[0]




class Reaction:
    """
    This class is used to store and process the reaction information.
    In order to instantiate a reaction object, you need to pass a dictionary of the reaction information.
    This dictionary must include 'name','stoichiometry' keys. This follows the format of the seed database.
    Examples:
        >>> A={"name":'D-glucose-6-phosphate aldose-ketose-isomerase',"stoichiometry":'-1:cpd00079:0:0:\"D-glucose-6-phosphate\";1:cpd00072:0:0:\"D-fructose-6-phosphate\"'}
        >>> a=Reaction(A)
        >>> print(a)
        D-glucose-6-phosphate aldose-ketose-isomerase

    """
    def __init__(self, dict):
        self.dict = dict

    def __str__(self):
        return self.dict['name']

    @property
    def stoichiometry(self):
        """
        Returns the stoichiometry of the reaction by the seed id of the compounds as key and the
        stoichiometric coefficient as value.
        Examples:
            >>> A={"name":'D-glucose-6-phosphate aldose-ketose-isomerase',"stoichiometry":'-1:cpd00079:0:0:\"D-glucose-6-phosphate\";1:cpd00072:0:0:\"D-fructose-6-phosphate\"'}
            >>> a=Reaction(A)
            >>> a.Stoichiometry=={'cpd00079': -1, 'cpd00072': 1}
            True
        
        Args:
            self (Reaction): An instance of the Reaction.

        Returns:
            dict: The stoichiometry of the reaction by the seed id of the compounds as key and the
        stoichiometric coefficient as value.
        """

        S = {}
        for compound in self.dict['stoichiometry'].split(';'):
            S[compound.split(':')[1]] = float(compound.split(':')[0])
        return S


class Metabolite:
    """
    Any metabolite with seed id can be an instance of this class.
    In order to instantiate a Metabolite object, you first define a dictionary in seed database format.
    This dictionary must have  "name", "mass", and "formula" keys, but it is okay if it has other keys.
    Examples:
        >>> A={"name":"methane","mass":16,"formula":"CH4"}
        >>> a=Metabolite(A)
        >>> print(a)
        methane

    """

    def __init__(self, dict):
        self.dict = dict
        self.cod = self.cod_calc()

    def __str__(self) -> str:
        return self.dict['name']
    
    def __repr__(self) -> str:
        return self.dict['name']

    def cod_calc(self)->float:
        """
        Calculates the conversion rates for g/l -> gCOD/l

        Examples:
            >>> A={"name":"methane","mass":16,"formula":"CH4"}
            >>> a=Metabolite(A)
            >>> a.COD
            4.0

        Args:
            self (Metabolite): An instance of the Metabolite class: Note

        Returns:
            float: COD conversion from g/l to gCOD/l
        

        """
        if self.dict['formula'] and self.dict['mass']:
            contents = {}
            atoms = ["H", "C", "O"]
            mw = self.dict['mass']
            for atom in atoms:
                if re.search(atom+'\d*', self.dict['formula']):
                    if len(re.search(atom+'\d*', self.dict['formula']).group()[1:]) == 0:
                        contents[atom] = 1
                    else:
                        contents[atom] = int(
                            re.search(atom+'\d*', self.dict['formula']).group()[1:])
                else:
                    contents[atom] = 0
            return 1/mw*(contents['H']+4*contents['C']-2*contents['O'])/4*32

        else:
            return 'None'

class Database:

    '''
    This class will handle all of the database tasks in ADToolBox
    '''

    def __init__(self, config:configs.Database=configs.Database()):

        self.config = config


    def _initialize_database(self):

        with open(self.config.protein_db, 'w') as f:
            pass
            # The rest is for other possible formats; fastq or gzipped fasta or gzipped fastq

    
    def filter_seed_from_ec(self, ec_list,
                            reaction_db=configs.Database().reaction_db,
                            compound_db=configs.Database().compound_db,
                            local_reaction_db=configs.Database().local_reaction_db,
                            local_compound_db=configs.Database().local_compound_db) -> tuple:
        """

        This function takes a list of ec numbers Generates a mini-seed JSON files. This is supposed to
        make the code a lot faster, but makes no difference in terms of outputs, adn won't probably need
        frequent updates.

        """
        with open(reaction_db, 'r') as f:
            main_reaction_db = json.load(f)
        with open(compound_db, 'r') as f:
            main_compound_db = json.load(f)

        RT = Reaction_Toolkit()
        cached_compounds = []
        filtered_rxns_db = {}
        local_rxn_db = []
        local_comp_db = []
        counter = 0
        for ec in ec_list:
            for ind, rxn in enumerate(main_reaction_db):

                if main_reaction_db[ind]['ec_numbers'] != None and ec in main_reaction_db[ind]['ec_numbers']:
                    local_rxn_db.append(rxn)
                    for Mets in rxn["compound_ids"].split(";"):
                        if Mets not in cached_compounds:
                            cached_compounds.append(Mets)
            counter += 1
            print(" -> Percent of ecs processed: ",end=" ")
            print("%"+str(int(Counter/len(ec_list)*100)), end="\r")

        counter = 0
        for compound in cached_compounds:
            for ind, Comp in enumerate(main_compound_db):
                if compound == Comp["id"]:
                    local_comp_db.append(Comp)
            counter += 1
            print(" -> Percent of compunds processed: ",end=" ")
            print("%"+str(int(Counter/len(cached_compounds)*100)), end="\r")

        with open(local_reaction_db, 'w') as f:
            json.dump(local_rxn_db, f)
        with open(local_compound_db, 'w') as f:
            json.dump(local_comp_db, f)

        return (local_rxn_db, local_comp_db)

    # The complete pandas objec can be used as input
    session = requests.Session()
    retry = Retry(connect=3, backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)

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
    
    def protein_db_from_ec(self, ec_list,mode="a"):
        # 5.1.1.2 and 5.1.1.20 should be distinguished: So far it seems like they are distinguished automatically by Uniprot
        with open(self.config.protein_db, mode) as f:
            for ec in track(ec_list,description="[yellow]   --> Writing the protein database:"):
                try:
                    file = Database.session.get(
                        f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28ec%3A{ec}%29%20AND%20%28reviewed%3Atrue%29%20NOT%20%28taxonomy_id%3A2759%29%29", timeout=1000)

                except requests.exceptions.HTTPError or requests.exceptions.ConnectionError:

                    print("Request Error! Trying again ...")
                    time.sleep(30)
                    file = Database.session.get(
                        f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28ec%3A{ec}%29%20AND%20%28reviewed%3Atrue%29%20NOT%20%28taxonomy_id%3A2759%29%29", timeout=1000)
                # This alsp does a sanity check

                except Exception:
                    print('Something went wrong!')
                text = file.text
                for line in text.split('\n'):
                    if line.startswith('>'):
                        uniprot=line.split('|')[1]
                        f.write(f'>{uniprot}|{ec}\n')
                    else:
                        f.write(f'{line}\n')


    @staticmethod
    def ec_from_csv(csv_file:str, Sep=',')->list:

        """This function reads a csv file and returns a list of ec numbers"""
        

        try:

            ec_list = pd.read_table(csv_file, sep=Sep,dtype="str")
            ec_list.dropna(axis=0)
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

    def ec_from_uniprot(self, uniprot_id):

        Base_URL = 'https://www.uniprot.org/uniprot/'

        try:
            file = Database.session.get(
                Base_URL+uniprot_id+".txt", timeout=10)

        except requests.exceptions.ConnectionError:

            print("Request Error! Trying again ...")
            time.sleep(30)
            file = Database.session.get(
                Base_URL+uniprot_id+".txt", timeout=10)

            if file.ok:

                print("Retry was successful!")

        for line in file:
            if re.search("ec=[0-9]+", line.decode("utf-8")):
                ec = re.sub("ec=", "", re.search(
                    "ec=[.0-9]+", line.decode("utf-8")).group(0))
                break

        else:

            print("Retry was unsuccessful!")
            return None

        return ec

    cazy_links = ["http://www.cazy.org/Glycoside-Hydrolases.html",
                  "http://www.cazy.org/Polysaccharide-Lyases.html",
                  "http://www.cazy.org/Carbohydrate-Esterases.html"
                  ]

    def cazy_ec_from_link(self):
        '''
        Extracts the ec numbers from a link to the Cazy website

        '''

        ec_list = []
        for link in Database.cazy_links:

            page = requests.get(link)
            soup = BeautifulSoup(page.content, "html.parser")
            results = soup.find("div", class_="cadre_principal").find_all(
                "th", class_="thec")
            for ec_number in results:
                if '-' not in ec_number.text.strip() and '.' in ec_number.text.strip():

                    ec_list.append(ec_number.text.strip())

        print("ec numbers extracted from Cazy website!")
        return ec_list

    
    def seed_from_ec(self,ec_Number, mode='Single_Match'):

        with open(self.config.reaction_db,'r') as f: 

            data = json.load(f)

        if mode == 'Single_Match':

            for i in data:

                if i['ec_numbers']:

                    if ec_Number in i['ec_numbers']:

                        return i["id"]

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
            seed_id = self.Seed_From_ec(ec, Mode='Multiple_Match')
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

    ###-----FeedStock_Database_Functions-----###
    def init_feedstock_database(self)-> None:
        '''
        Makes an empty feedstock database json file.
        BE CAREFUL: This will overwrite the current database file.
        '''
        dec=input("Are you sure you want to create a new database? (y/n): ")
        if dec=="y":

            with open(self.config.feed_db, 'w') as f:
                json.dump([], f)
        else:
            print("Aborting!")
            return
        ## TO-DO: Error and Failure handling

    def _add_feedstock_to_database(self, feedstock):
        
        '''
        Adds a feedstock to the feedstock database.
        A Feed Stock must be a dictionary with the following keys:
        "Name": Name of the feedstock
        "TS": Totall Solids (in COD)
        "TSS": Totall Suspended Solids (in COD)
        "Lipids": Lipids (in COD)
        "Proteins": Proteins (in COD)
        "Carbohydrates": Carbohydrates (in COD)
        "PI": Particulate Inerts
        "SI": Soluble Inerts
        "Notes": Additional notes like the source of the data

        '''
        try:
            with open(self.config.feed_db, 'r') as f:
                database=json.load(f)
        except FileNotFoundError:
            print("Feedstock database not found! Creating new one...")
            self.init_feedstock_database()
        
        finally:
        
            database.append(feedstock)
            with open(self.config.feed_db, 'w') as f:
                json.dump(database, f)
        
    def add_feedstock_to_database_from_file(self, file_path):
        '''
        Adds a feedstock to the feedstock database from a csv file.
        A Feed Stock spreadsheet must be a table with the following columns:
        "Name", "TS", "TSS", "Lipids", "Proteins", "Carbohydrates", "PI", "SI", "Notes"
        File_Path does not come with a default path, since it is advised to not populate the
        project directory with feedstock spreadsheets unwanted.
        '''
        try:
            f=pd.read_table(file_path, sep=',')
        except FileNotFoundError:
            print("File not found!")
            return
        counter=0
        for i in f.index:
            feedstock=f.loc[i,:].to_dict()
            self._add_feedstock_to_database(feedstock)
            counter+=1
        print("{} Feedstocks added to the database!".format(counter))

    @staticmethod
    def download_adm1_parameters(config):
        '''
        Downloads the parameters needed for running ADM models in ADToolbox
        
        '''
        if not os.path.exists(config.base_dir):
            os.makedirs(config.base_dir)
        model_parameters_dir=config.model_parameters
        base_parameters_dir=config.base_parameters
        initial_conditions_dir=config.initial_conditions
        inlet_conditions_dir=config.inlet_conditions
        reactions_dir=config.reactions
        species_dir=config.species
        model_parameters="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/ADM1/ADM1_Model_Parameters.json"
        base_parameters="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/ADM1/ADM1_Base_Parameters.json"
        initial_conditions="https://github.com/ParsaGhadermazi/Database/blob/main/ADToolbox/ADM1/ADM1_Initial_Conditions.json"
        inlet_conditions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/ADM1/ADM1_Inlet_Conditions.json"
        reactions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/ADM1/ADM1_Reactions.json"
        species="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/ADM1/ADM1_Species.json"
        r = requests.get(model_parameters, allow_redirects=True)
        with open(model_parameters_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(base_parameters, allow_redirects=True)
        with open(base_parameters_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(initial_conditions, allow_redirects=True)
        with open(initial_conditions_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(inlet_conditions, allow_redirects=True)
        with open(inlet_conditions_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(reactions, allow_redirects=True)
        with open(reactions_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(species, allow_redirects=True)
        with open(species_dir, 'wb') as f:
            f.write(r.content)
    @staticmethod
    def download_modified_adm_parameters(config):    
        '''
        Downloads the parameters needed for running ADM models in ADToolbox
        '''
        if not os.path.exists(config.base_dir):
            os.makedirs(config.base_dir)
        model_parameters_dir=config.model_parameters
        base_parameters_dir=config.base_parameters
        initial_conditions_dir=config.initial_conditions
        inlet_conditions_dir=config.inlet_conditions
        reactions_dir=config.reactions
        species_dir=config.species
        model_parameters="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Model_Parameters.json"
        base_parameters="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Base_Parameters.json"
        initial_conditions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Initial_Conditions.json"
        inlet_conditions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Inlet_Conditions.json"
        reactions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Reactions.json"
        species="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Species.json"
        r = requests.get(model_parameters, allow_redirects=True,stream=True)
        with open(model_parameters_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(base_parameters, allow_redirects=True)
        with open(base_parameters_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(initial_conditions, allow_redirects=True)
        with open(initial_conditions_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(inlet_conditions, allow_redirects=True)
        with open(inlet_conditions_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(reactions, allow_redirects=True)
        with open(reactions_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(species, allow_redirects=True)
        with open(species_dir, 'wb') as f:
            f.write(r.content)
        
        
    
    def download_seed_databases(self) -> None :
        "Download the modelseed rection database"
        r = requests.get(self.config.seed_rxn_url, allow_redirects=True,stream=True)
        with open(self.config.reaction_db, 'wb') as f:
            f.write(r.content)

    def download_protein_database(self) -> None:
        r = requests.get(self.config.protein_db_url, allow_redirects=True)
        with open(self.config.protein_db, 'wb') as f:
            f.write(r.content)
        
    def download_reaction_database(self)->None:
        r = requests.get(self.config.adtoolbox_rxn_db_url, allow_redirects=True)
        with open(self.config.csv_reaction_db, 'wb') as f:
            f.write(r.content)

    
    def download_feed_database(self)-> None:
        r = requests.get(self.config.feed_db_url, allow_redirects=True)
        with open(self.config.feed_db, 'wb') as f:
            f.write(r.content)
    
    def download_qiime_classifier_db(self)->None:
        r = requests.get(self.config.qiime_classifier_db_url, allow_redirects=True,stream=True)
        block_size = 1024
        total_size = int(r.headers.get('content-length', 0))
        if not os.path.exists(Path(self.config.qiime_classifier_db).parent):
            os.makedirs(Path(self.config.qiime_classifier_db).parent)
        with open(self.config.qiime_classifier_db, 'wb') as f:
            with Progress() as progress:
                task = progress.add_task("Downloading the qiime's classifier database...", total=total_size)
                for data in r.iter_content(block_size):
                    progress.update(task, advance=len(data))
                    f.write(data)
        
        rich.print(f"[bold green]Download finished[/bold green]")
            


        

    # def download_escher_files(self)-> None:
    #     for i in self.config.escher_files_urls:
    #         r = requests.get(i, allow_redirects=True)
    #         with open(os.path.join(self.config.,i.split("/")[-1]), 'wb') as f:
    #             f.write(r.content)

    def get_metagenomics_studies(self) -> list:
        """
        This function will return accession numbers in all metagenomics studies on the kbase.
        Requires:
            <R>Configs.Database</R>
        Satisfies:
            <S>project_accession</S>


        """
        try:
            metagenomics_studies=pd.read_table(self.config.kbase_db.metagenomics_studies,delimiter="ᚢ",encoding="utf-8",engine="python")

        except FileNotFoundError:

            rich.print("[bold red]KBase Database is not found. Please download/build/configure the database first.[/bold red]")


        return metagenomics_studies.to_dict(orient="list")
    
    def download_kbase_database(self)->None:
        """
        This function will download the kbase database from the remote repository.
        """
        
        
        r = requests.get(self.config.kbase_db.urls['metagenomics_studies'], allow_redirects=True)
        if not os.path.exists(self.config.kbase_db.base_dir):
            os.makedirs(self.config.kbase_db.base_dir,exist_ok=True)
        with open(os.path.join(self.config.kbase_db.metagenomics_studies), 'wb') as f:
            f.write(r.content)
        rich.print(f"[bold green]Downloaded {self.config.kbase_db.urls['metagenomics_studies']}[/bold green]")
    
        r = requests.get(self.config.kbase_db.urls['exmpermental_data_references'], allow_redirects=True)
        with open(os.path.join(self.config.kbase_db.experimental_data_references), 'wb') as f:
            f.write(r.content)

        rich.print(f"[bold green]Downloaded {self.config.kbase_db.urls['exmpermental_data_references']}[/bold green]")  

    def download_all_databases(self)->None:
        """
        This function will download all the required databases for the ADToolbox.
        """
        if not os.path.exists(self.config.base_dir):
            os.makedirs(self.config.base_dir)
        self.download_seed_databases()
        self.download_protein_database()
        self.download_reaction_database()
        self.download_feed_database()
        self.download_kbase_database()
        

     

class Metagenomics:

    """
    A Core class that handles the metganomics data
    Init this with a Metagenomics class of config module
    
    
    
    It takes three diferent types of input:


    1- QIIMEII Outputs -> Refer to the documentation for the required files and format
    2- Shotgun Short reads -> Refer to the documentation for the required files and format
    3- MAGs -> Refer to the documentation for the required files and format

    This class will use other classes to generate the required files for downstream analysis.
    """
    def __init__(self,config:configs.Metagenomics):
        self.config=config
    
    def get_requirements_a2g(self,progress=True):
        """
        This function will automatically download the required files for Amplicon to Genome functionality.
        """
        if not os.path.exists(self.config.amplicon2genome_db):
            os.mkdir(self.config.amplicon2genome_db)

        url = {'Version': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/VERSION',
               'MD5SUM': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/MD5SUM',
               'FILE_DESCRIPTIONS': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/FILE_DESCRIPTIONS',
               'metadata_field_desc': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/metadata_field_desc.tsv',
               'bac120_metadata': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tar.gz',
               'bac120_ssu': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/bac120_ssu_reps.tar.gz'
               }
        if progress:
            for keys in ['Version', 'MD5SUM', 'FILE_DESCRIPTIONS']:
                with requests.get(url[keys], allow_redirects=True, stream=progress) as r:
                    total_size = int(r.headers.get('content-length', 0))
                    block_size = 1024
                    with Progress() as progress:
                        task1 = progress.add_task("Downloading " + keys, total=total_size)
                        with open(os.path.join(self.config.amplicon2genome_db, keys), 'wb') as f:
                            for data in r.iter_content(block_size):
                                progress.update(task1, advance=len(data))
                                f.write(data)
            with requests.get(url['metadata_field_desc'], allow_redirects=True, stream=progress) as r:
                total_size = int(r.headers.get('content-length', 0))
                block_size = 1024
                with Progress() as progress:
                    task1 = progress.add_task("Downloading metadata_field_desc.tsv", total=total_size)
                    with open(os.path.join(self.config.amplicon2genome_db, 'metadata_field_desc.tsv'), 'wb') as f:
                        for data in r.iter_content(block_size):
                            progress.update(task1, advance=len(data))
                            f.write(data)

            for keys in ['bac120_metadata', 'bac120_ssu']:
                with requests.get(url[keys], allow_redirects=True, stream=progress) as r:
                    total_size = int(r.headers.get('content-length', 0))
                    block_size = 1024
                    with Progress() as progress:
                        task1 = progress.add_task("Downloading " + keys, total=total_size)
                        with open(os.path.join(self.config.amplicon2genome_db, keys+'.tar.gz'), 'wb') as f:
                            for data in r.iter_content(block_size):
                                progress.update(task1, advance=len(data))
                                f.write(data)
        else:
            for keys in ['Version', 'MD5SUM', 'FILE_DESCRIPTIONS']:
                with requests.get(url[keys], allow_redirects=True, stream=progress) as r:
                    with open(os.path.join(self.config.amplicon2genome_db, keys), 'wb') as f:
                        f.write(r.content)
            with requests.get(url['metadata_field_desc'], allow_redirects=True, stream=progress) as r:
                with open(os.path.join(self.config.amplicon2genome_db, 'metadata_field_desc.tsv'), 'wb') as f:
                    f.write(r.content)
            for keys in ['bac120_metadata', 'bac120_ssu']:
                with requests.get(url[keys], allow_redirects=True, stream=progress) as r:
                    with open(os.path.join(self.config.amplicon2genome_db, keys+'.tar.gz'), 'wb') as f:
                        f.write(r.content)

    
    def find_top_k_repseqs(self,
                            sample_names:Union[list[str],None]=None
                            ,save:bool=True)->dict:

        """This function takes three inputs from qiime:
        1. feature table: This is the abundance of each feature in each sample (TSV).
        2. taxonomy table: This is the taxonomy of each feature (TSV). 
        3. rep seqs: This is the representative sequence of each feature (fasta).
        It then finds the top k features for each sample and then finds the genomes for each feature (if there is one).
        It then returns a dictionary of the genomes for the sample(s) along with the NCBI name and the
        address where they are stored.

        Requires:

        Satisfies:

        Args:
            sample_name (str, optional): The name of the sample. Defaults to None.
            if None, it will find top k genomes for all samples.
            save (bool, optional): Whether to save the top k repseqs. Defaults to True.
        
        Returns:
            dict: A dictionary of information about the top k features/taxa
        """

        feature_table = pd.read_table(self.config.feature_table_dir, sep='\t',skiprows=1)
        taxonomy_table = pd.read_table(self.config.taxonomy_table_dir, delimiter='\t')
        repseqs=fasta_to_dict(self.config.rep_seq_fasta)

        if sample_names is None:
            pass 
        else:
            feature_table = feature_table[sample_names]

        top_featureids = {}
        top_repseqs = {}
        top_k_taxa = {}


        for sample in feature_table.columns[1:]:
            top_featureids[sample] = list(feature_table.sort_values(
                sample, ascending=False)['#OTU ID'].head(self.config.k))
            top_k_taxa[sample] = [ taxonomy_table[taxonomy_table["Feature ID"]==feature]["Taxon"].item() for feature in top_featureids[sample]]
            top_repseqs |= {feature:repseqs[feature] for feature in top_featureids[sample]}
        
        if save:
            pd.DataFrame(top_featureids).to_csv(os.path.join(self.config.amplicon2genome_outputs_dir,'Top_k_FeatureIDs.csv'))
            pd.DataFrame(top_k_taxa).to_csv(os.path.join(self.config.amplicon2genome_outputs_dir,'Top_k_Taxa.csv'))
            with open(os.path.join(self.config.amplicon2genome_top_repseq_dir), 'w') as f:
                for feature in top_repseqs:
                    f.write('>'+feature+'\n'+top_repseqs[feature]+'\n')                

        return {'top_featureids':top_featureids,'top_k_taxa':top_k_taxa,'top_repseqs':top_repseqs}

        
    
    def align_to_gtdb(self,run:bool=True,save:bool=True,container:str="None")->str:
        """This function takes the representative sequences of the top k features and generates the script to
        align these feature sequences to gtdb using VSEARCH. If you intend to run this you either
        need to have VSEARCH installed or run it with a container option. You can use either the docker or singularity
        as container options. Otherwise you can use None and run it with the assumption that VSEARCH is installed.
        If you only want the script and not to run it, set run to False.

        Requires:

        Satisfies:

        Args:
            run (bool, optional): Whether to run the script. Defaults to True.
            save (bool, optional): Whether to save the script. Defaults to True.
            container (str, optional): The container to use. Defaults to "None".
        
        Returns:
            dict: A dictionary containing the generated scripts and address of the genomes downloaded or
            tobe downloaded after running the script if you are not using this function to run it.
        """
        alignment_dir = os.path.join(self.config.amplicon2genome_outputs_dir,'Alignments')
        match_table=os.path.join(self.config.amplicon2genome_outputs_dir,'matches.blast')
        gtdb_dir_fasta=self.config.gtdb_dir_fasta
        query=self.config.amplicon2genome_top_repseq_dir
        dirs=[self.config.amplicon2genome_outputs_dir,
            gtdb_dir_fasta,
            query
            ]
        if container=="None":
            
            bash_script=("#!/bin/bash\n"+'vsearch --top_hits_only --blast6out '+
                        match_table+
                        ' --usearch_global '+ query +
                        ' --db '+ gtdb_dir_fasta +
                        ' --id ' +str(self.config.amplicon2genome_similarity) +
                        ' --threads '+str(self.config.vsearch_threads)+
                        ' --alnout '+ alignment_dir +
                        ' --top_hits_only'+'\n')
        
        if container=="docker":
            alignment_dir = os.path.join(self.config.amplicon2genome_outputs_dir,'Alignments')
            bash_script="#!/bin/bash\n" +'docker run -it '
            for dir in dirs:
                bash_script+=('-v '+dir+':'+dir+' ')
            
            bash_script += (self.config.amplicon2genome_docker+' vsearch --top_hits_only --blast6out '+
                        match_table+
                        ' --usearch_global '+ query +
                        ' --db '+ gtdb_dir_fasta +
                        ' --id ' +str(self.config.amplicon2genome_similarity) +
                        ' --threads '+str(self.config.vsearch_threads)+
                        ' --alnout '+ alignment_dir +
                        ' --top_hits_only'+'\n')
        
        if container=="singularity":
            warnings.warn("Singularity is not supported yet")

        if save:
            with open(self.config.vsearch_script_dir, 'w') as f:
                f.write(bash_script)
        
        if run:
            subprocess.run(bash_script,shell=True)
        
        return bash_script
    
    def get_genomes_from_alignment(self,save:bool=True)->dict:
        """This function takes the alignment file generated from the align_to_gtdb function and generates the the genome information
        using the GTDB-Tk. if you want to save the information as a json file set save to True.

        Requires:

        Satisfies:

        Args:
            save (bool, optional): Whether to save the feature:genome_id mapping as a json file. Defaults to True.
        
        Returns:
            dict: A dictionary containing
        """
        matches = os.path.join(self.config.amplicon2genome_outputs_dir,'matches.blast')
        aligned=pd.read_table(matches,header=None,delimiter='\t')
        aligned.drop_duplicates(1,inplace=True)
        aligned[1]=aligned[1].apply(lambda x: "".join(x.split('_')[1:]))
        alignment_dict=dict(zip(aligned[0],aligned[1]))
        if save:
            with open(self.config.feature_to_taxa, 'w') as f:
                json.dump(alignment_dict, f)
        
        return alignment_dict
    
    
    def download_genomes(self,identifiers:list[str],save:bool=True,container:str="None")->dict:
        """This function downloads the genomes from NCBI using the refseq/genbank identifiers.
        If you want to save the json file that is holding the information for each genome, set save to True.
        Note that this function uses rsync to download the genomes. Also this function does not come with a run option.
    

        Requires:

        Satisfies:

        Args:
            identifier list[str]: The list of identifiers for the genomes. It can be either refseq or genbank.
            save (bool, optional): Whether to save the genome information as a json file. Defaults to True.
            container (str, optional): The container to use. Defaults to "None". You may select from "None", "docker", "singularity".
        
        Returns:
            dict: A dictionary containing the genome information.

        """
        genome_info = {}
        base_ncbi_dir = 'rsync://ftp.ncbi.nlm.nih.gov/genomes/all/'
        for identifier in identifiers:
            specific_ncbi_dir = identifier[0:3]+'/'+\
                                identifier[3:6]+'/'+\
                                identifier[6:9]+'/'+\
                                identifier[9:].split('.')[0]
            
            genome_dir=pathlib.Path(self.config.genome_save_dir(identifier))
            
            if container=="None":
                subprocess.run('rsync -avz --progress '+base_ncbi_dir+specific_ncbi_dir+' '+self.config.genome_save_dir(identifier),shell=True)

            
            if container=="docker":
                warnings.warn("Docker is not supported yet")
                subprocess.run('docker run -it -v '+str(genome_dir.parent)+':'+str(genome_dir.parent)+ ' vsearch rsync -avz --progress '+' '+base_ncbi_dir+specific_ncbi_dir+' '+str(genome_dir),shell=True)
            
            if container=="singularity":
                warnings.warn("Singularity is not supported yet")
            
            for i in pathlib.Path(self.config.genome_save_dir(identifier)).glob('**/*fna.gz'):
                if 'rna' not in i.name.lower():
                    if 'cds' not in i.name.lower():
                        if 'genomic' in i.name.lower():
                            genomics_file=str(i)

            genome_info[identifier] = genomics_file
        
        if save:
            with open(self.config.genomes_json_info, 'w') as f:
                json.dump(genome_info, f)
        
        return genome_info
    

    def align_genomes_to_protein_db(self,run:bool=True,save:bool=True,container:str="None")->dict:
        """
        This is a function that will align genomes to the Protein Database of the ADToolbox using local alignment
        and then will return a dictionary of the alignment results with the aligned genome identifiers as key and the genome as value.
        If you want to save the scripts, set save to True. Note that the alignment tables will be saved in any case.
        Note that this function uses mmseqs2 to align the genomes to the protein database. So, to run this function without 
        any container you need to have mmseqs2 installed on your system. However, if you want to run this function with a container,
        you need to have the container installed on your system. You may select from "None", "docker", "singularity".

        Requires:

        Satisfies:

        Args:
            run (bool, optional): Whether to run the alignment. Defaults to True.
            save (bool, optional): Whether to save the alignment scripts. Defaults to True.
            container (str, optional): The container to use. Defaults to "None". You may select from "None", "docker", "singularity".
        
        Returns:
            dict: A dictionary containing the alignment results.


        """
        with open(self.config.genomes_json_info, 'r') as f:
            genomes_info = json.load(f)
        valid_genomes = 0
        for genome_ids in track(genomes_info.keys(),description="Fetching Genomes ..."):
            Totall_Genomes = len(genomes_info.keys())
            if os.path.exists(genomes_info[genome_ids]):
                valid_genomes += 1
        rich.print(f"[green]{valid_genomes}/{Totall_Genomes} Genomes are valid")


        genome_alignment_files={}
        if container=="None":
            bash_script = "#!/bin/bash\n"
            for genome_id in genomes_info.keys():
                alignment_file=os.path.join(self.config.genome_alignment_output,"Alignment_Results_mmseq_"+genome_id+".tsv")
                bash_script += "mmseqs easy-search " + \
                    genomes_info[genome_id] + " " + \
                    self.config.protein_db + " " + \
                    alignment_file+ " tmp\n\n"
                genome_alignment_files[genome_id]=alignment_file
        
        if container=="docker":
            warnings.warn("Docker is not supported yet")
            bash_script = "#!/bin/bash\n"
            for genome_id in genomes_info.keys():
                alignment_file=os.path.join(self.config.genome_alignment_output,"Alignment_Results_mmseq_"+genome_id+".tsv")
                bash_script +="docker run -it "+ \
                " -v "+genomes_info[genome_id]+":"+genomes_info[genome_id]+ \
                " -v "+self.config.protein_db+":"+self.config.protein_db+ \
                " -v "+self.config.genome_alignment_output+":"+self.config.genome_alignment_output+ \
                f" {self.config.amplicon2genome_docker}  mmseqs easy-search " + \
                    genomes_info[genome_id] + " " + \
                    self.config.protein_db + " " + \
                    alignment_file+ " tmpfiles\n\n"
                genome_alignment_files[genome_id]=alignment_file
            


        with open(self.config.genome_alignment_output_json, 'w') as f:
            json.dump(genome_alignment_files, f)
        
        if save:
            with open(self.config.genome_alignment_script, 'w') as f:
                f.write(bash_script)

        if run:

            os.system(bash_script)

        return genome_alignment_files
    
    @staticmethod
    def make_json_from_genomes(input_dir:str,output_dir:str)->dict:
        """
        This function takes a csv file directory. The CSV file must include three columns:
        1- genome_id: Identifier for the genome, preferrably; NCBI ID
        2- NCBI_Name: NCBI taxonomy name for the genome: Does not need to be in a specific format
        3- Genome_Dir: Absolute path to the fasta files: NOT .gz

        """
        genomes_json={}
        genomes_table=pd.read_table(input_dir,delimiter=",")
        genomes_table.set_index("genome_id",inplace=True)
        for[gi] in genomes_table.index:
            genomes_json[gi]={}
            genomes_json[gi]["NCBI_Name"]= genomes_table.loc[gi,"NCBI_Name"]
            genomes_json[gi]["Genome_Dir"]= genomes_table.loc[gi,"Genome_Dir"]
        with open(output_dir,"w") as fp:
            json.dump(genomes_json,fp)
        return
    


    def adm_from_alignment_json(self,adm_rxns,model="Modified_ADM_Reactions"):
        rt = Reaction_Toolkit(reaction_db=self.config.seed_rxn_db)
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
        # So far it has gotten really complicated, after making sure it works, instead of adding ec numbers I'll add [SEED,STR] so
        # that it can be used for downstream analysis
    
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

    def calculate_microbial_portions(self,microbe_reaction_map:dict,save:bool=True)-> dict:
        
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
        
        with open(self.config.genome_relative_abundances,'r') as f:
            relative_abundances=json.load(f)

        with open(self.config.genome_alignment_output_json,'r') as f:    
            genome_info=json.load(f)
        
        for sample in relative_abundances.keys():

            cod_portion=Additive_Dict([(i,0) for i in model_species])

            for genome_id in relative_abundances[sample].keys():
                temp_tsv = pd.read_table(
                            genome_info[genome_id], sep='\t',header=None)
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
                cod_portion=cod_portion+Additive_Dict(microbial_species)*relative_abundances[sample][genome_id]
                SUM=sum([cod_portion[key] for key in cod_portion])
                if SUM==0:
                    for i in cod_portion:
                        
                        cod_portion[i]=0
                else:
                    for i in cod_portion:
                        cod_portion[i]=cod_portion[i]/SUM
                
            
            cod[sample]=cod_portion
        if save:
            with open(self.config.cod_output_json,'w') as f:
                json.dump(cod,f)
        

        return cod
    
    def seqs_from_sra(self,accession:str,save:bool=True,run:bool=True,container:str="None")-> tuple[str,str]:
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
        
        Returns:
            prefetch_script (str): The bash script that will be used to download the SRA files in python string format
            fasterq_dump_script (str): The bash script that will be used to convert the SRA files to fastq files in python string format
    
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

        if save:
            
            with open(project_path.joinpath("prefetch.sh"),"w") as f:
                f.write(prefetch_script)
            
            with open(project_path.joinpath("fasterq_dump.sh"),"w") as f:
                f.write(fasterq_dump_script)
                
        
        if run:
            subprocess.run(["bash",str(project_path.joinpath("prefetch.sh"))],cwd=str(project_path))
            subprocess.run(["bash",str(project_path.joinpath("fasterq_dump.sh"))],cwd=str(project_path))
            
            
        
        return prefetch_script,fasterq_dump_script
    
    def run_qiime2_from_sra(self,accession:str,save:bool=True,run:bool=True,container:str='None') -> tuple[str,str]:
        """
        Requires:
            <R>project_accession</R>
            <R>Configs.Metagenomics</R>
            <R>sra_project_dir</R>
            <R>.fastq_files</R>
        Satisfies:
            <S>qiime2_work_dir</S>
            <S>qiime2_bash_str</S>
            <S>manifest.csv</S>
            <if>save=True</if><S>qiime2_bash_file</S>
            <if>run=True</if><S>qiime2_outputs</S>
        Args:
            accession (str): The accession number of the SRA project or run
            save (bool, optional): If True, the  bash scripts will be saved in the SRA work directory. Defaults to True.
            run (bool, optional): If True, the bash scripts will be executed. Defaults to True. You must have the SRA toolkit installed in your system.
            container (str, optional): If you want to run the qiime2 commands in a container, specify the container name here. Defaults to 'None'.
        Returns:
            qiime2_bash_str (str): The bash script that will be used to run qiime2 in python string format
            manifest (dict): The manifest file that will be used to run qiime2 in python dictionary format
        
        TODO:
            add the ability to submit slurm jobs to cluster
            add singularity support

        """
        qiime2_work_dir=pathlib.Path(self.config.qiime2_work_dir(accession))
        if not qiime2_work_dir.exists():
            qiime2_work_dir.mkdir(parents=True)
        
        sra_project_dir=pathlib.Path(self.config.sra_work_dir(accession))
        seqs=sra_project_dir.joinpath("seqs")
        manifest_single={'sample-id':[],'absolute-filepath':[]}
        manifest_paired={'sample-id':[],'forward-absolute-filepath':[],'reverse-absolute-filepath':[]}
        
        for i in seqs.iterdir():
            if i.is_dir():
                
                if len(list(i.glob("*.fastq"))) == 2:
                    r1=list(i.glob("*_1.fastq"))[0]
                    r2=list(i.glob("*_2.fastq"))[0]
                    manifest_paired['sample-id'].append(i.name)
                    manifest_paired['forward-absolute-filepath'].append(str(r1))
                    manifest_paired['reverse-absolute-filepath'].append(str(r2))
                    paired_end=True
                elif len(list(i.glob("*.fastq"))) == 1:
                    r1=list(i.glob("*.fastq"))[0]
                    manifest_single['sample-id'].append(i.name)
                    manifest_single['absolute-filepath'].append(str(r1))
                    paired_end=False
                    
        manifest=pd.DataFrame(manifest_single) if not paired_end else pd.DataFrame(manifest_paired)
  
        if paired_end:
            with open(self.config.qiime2_paired_end_bash_str,"r") as f:
                qiime2_bash_str=f.read()
        else:
            with open(self.config.qiime2_single_end_bash_str,"r") as f:
                qiime2_bash_str=f.read()

        manifest_dir=sra_project_dir.joinpath("manifest.tsv")

        if container=="None":
            qiime2_bash_str=qiime2_bash_str.replace("<manifest>",str(manifest_dir))
            qiime2_bash_str=qiime2_bash_str.replace("<qiime2_work_dir>",str(qiime2_work_dir))
            qiime2_bash_str=qiime2_bash_str.replace("<classifier>",str(self.config.qiime_classifier_db))
        
        elif container=="docker":
            qiime2_bash_str=qiime2_bash_str.splitlines()
            for idx,line in enumerate(qiime2_bash_str):
                line=line.lstrip()
                if line.startswith("qiime") or line.startswith("biom"):
                    qiime2_bash_str[idx]=f"docker run --env TMPDIR=/data/tmp -v {sra_project_dir}:/data/ -v {Path(self.config.qiime_classifier_db).parent}:/data/{Path(self.config.qiime_classifier_db).parent.name} -w /data  {self.config.qiime2_docker_image}"+" "+line
            qiime2_bash_str="\n".join(qiime2_bash_str)
            qiime2_bash_str=qiime2_bash_str.replace("<manifest>",str(manifest_dir.name))
            qiime2_bash_str=qiime2_bash_str.replace("<qiime2_work_dir>",str(qiime2_work_dir.name))
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
                    qiime2_bash_str[idx]=f"singularity exec --bind  $PWD:$PWD,{str(Path(self.config.qiime_classifier_db))}:{str(Path(self.config.qiime_classifier_db))},$SINGULARITY_TMPDIR:/tmp  {self.config.qiime2_singularity_image} " +line
            qiime2_bash_str="\n".join(qiime2_bash_str)
            qiime2_bash_str=qiime2_bash_str.replace("<manifest>",str(manifest_dir.name))
            qiime2_bash_str=qiime2_bash_str.replace("<qiime2_work_dir>",str(qiime2_work_dir.name))
            qiime2_bash_str=qiime2_bash_str.replace("<classifier>",str(Path(self.config.qiime_classifier_db)))

            # if not paired_end:
            #     manifest['absolute-filepath']=["/"+str(Path(seqs.name)/pathlib.Path(x).parent.name/pathlib.Path(x).name) for x in manifest['absolute-filepath']]
            # else:
            #     manifest['forward-absolute-filepath']=["/"+str(Path(seqs.name)/pathlib.Path(x).parent.name/pathlib.Path(x).name) for x in manifest['forward-absolute-filepath']]
            #     manifest['reverse-absolute-filepath']=["/"+str(Path(seqs.name)/pathlib.Path(x).parent.name/pathlib.Path(x).name) for x in manifest['reverse-absolute-filepath']]

        else:
            raise ValueError("Container must be None, singularity or docker")
        
        
        if save:
            pd.DataFrame(manifest).to_csv(sra_project_dir.joinpath("manifest.tsv"),sep='\t',index=False)
            with open(qiime2_work_dir.joinpath("qiime2.sh"),"w") as f:
                f.write(qiime2_bash_str)
        if run:
            subprocess.run(["bash",str(qiime2_work_dir.joinpath("qiime2.sh"))],cwd=str(qiime2_work_dir))
        
        return qiime2_bash_str,manifest
    
    def extract_taxonomy_features(self):
        pass

if __name__ == "__main__":
    config_1=configs.Metagenomics(qiime_outputs_dir="/Users/parsaghadermarzi/Desktop/ADToolbox/Metagenomics_Analysis/SRA/ERR3861428/qiime2",
    feature_table_dir="/Users/parsaghadermarzi/Desktop/ADToolbox/Metagenomics_Analysis/SRA/ERR3861428/qiime2/feature-table.tsv",
    taxonomy_table_dir="/Users/parsaghadermarzi/Desktop/ADToolbox/Metagenomics_Analysis/SRA/ERR3861428/qiime2/taxonomy.tsv",
    rep_seq_fasta="/Users/parsaghadermarzi/Desktop/ADToolbox/Metagenomics_Analysis/SRA/ERR3861428/qiime2/dna-sequences.fasta",
    amplicon2genome_top_repseq_dir="/Users/parsaghadermarzi/Desktop/ADToolbox/Metagenomics_Analysis/SRA/ERR3861428/qiime2/top10_repseqs.fasta",
    vsearch_script_dir="/Users/parsaghadermarzi/Desktop/ADToolbox/Metagenomics_Analysis/SRA/ERR3861428/qiime2/vsearch.sh",
    amplicon2genome_outputs_dir="/Users/parsaghadermarzi/Desktop/ADToolbox/Metagenomics_Analysis/SRA/ERR3861428/qiime2/",
    genomes_json_info="/Users/parsaghadermarzi/Desktop/ADToolbox/Metagenomics_Analysis/SRA/ERR3861428/qiime2/genomes.json",
    genome_alignment_output="/Users/parsaghadermarzi/Desktop/ADToolbox/Metagenomics_Analysis/SRA/ERR3861428/qiime2/",
    genome_alignment_output_json="/Users/parsaghadermarzi/Desktop/ADToolbox/Metagenomics_Analysis/SRA/ERR3861428/qiime2/alignment_info.json",
    feature_to_taxa="/Users/parsaghadermarzi/Desktop/ADToolbox/Metagenomics_Analysis/SRA/ERR3861428/qiime2/feature_to_genome.json",
    cod_output_json="/Users/parsaghadermarzi/Desktop/ADToolbox/Metagenomics_Analysis/SRA/ERR3861428/qiime2/cod_output.json",
    genome_relative_abundances="/Users/parsaghadermarzi/Desktop/ADToolbox/Metagenomics_Analysis/SRA/ERR3861428/qiime2/relative_abundances.json",
    amplicon2genome_similarity=0.9
    )
    metag=Metagenomics(config_1)
    metag.align_to_gtdb(run=True,save=True,container="docker")
    rep_gens=metag.get_genomes_from_alignment(save=True)
    metag.download_genomes([val for val in rep_gens.values()],save=True,container="docker")
    metag.align_genomes_to_protein_db(run=True,save=True,container="docker")







