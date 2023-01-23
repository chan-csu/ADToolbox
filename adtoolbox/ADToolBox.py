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
from requests.adapters import HTTPAdapter
from sympy import Li
import Configs
from requests.packages.urllib3.util.retry import Retry
from requests.exceptions import Timeout
from bs4 import BeautifulSoup
from datetime import datetime
from pathlib import Path
from collections import Counter
import pathlib
import gzip
import shutil
import Configs
import Bio_seq
from rich.progress import track,Progress
import rich
from __init__ import Main_Dir
from typing import Union
from utils import wrap_for_slurm
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
        self.kbase_config=kwargs.get("kbase_config",Configs.Kbase())
        self.database_config=kwargs.get("database_config",Configs.Database())
        self.reaction_toolkit_config=kwargs.get("reaction_toolkit_config",Configs.Reaction_Toolkit())
        self.metagenomics_config=kwargs.get("metagenomics_config",Configs.Metagenomics())
        self.alignment_config=kwargs.get("alignment_config",Configs.Alignment())
        self.original_adm_config=kwargs.get("original_adm_config",Configs.Original_ADM1())
        self.modified_adm_config=kwargs.get("modified_adm_config",Configs.Modified_ADM())
        
    
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

    def __init__(self, compound_db:str =Configs.Reaction_Toolkit().compound_db, reaction_db=Configs.Reaction_Toolkit().reaction_db) -> None:
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

    def __init__(self, config:Configs.Database=Configs.Database()):

        self.config = config


    def _initialize_database(self):

        with open(self.config.protein_db, 'w') as f:
            pass
            # The rest is for other possible formats; fastq or gzipped fasta or gzipped fastq

    
    def filter_seed_from_ec(self, ec_list,
                            reaction_db=Configs.Database().reaction_db,
                            compound_db=Configs.Database().compound_db,
                            local_reaction_db=Configs.Database().local_reaction_db,
                            local_compound_db=Configs.Database().local_compound_db) -> tuple:
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
    def __init__(self,config:Configs.Metagenomics):
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


    def amplicon2genome(self,
                        sample_name:list[str]=None,
                        download_databases: bool=False,
                        )->dict:

        """This function will take the amplicon data and will fetch the genomes.
        it is created to start from QIIMEII outputs, and download the genomes from NCBI.
        
        Args:
            sample_name (str, optional): The name of the sample. Defaults to None.
            if None, it will find top k genomes for all samples.
            download_databases (bool, optional): If True, it will download the required databases. Defaults to False.
        
        Returns:
            dict: A dictionary of the genomes for the sample(s) along with the NCBI name and the
            address where they are stored.
        """    

        if not os.path.exists(self.config.qiime_outputs_dir):
            os.mkdir(self.config.qiime_outputs_dir)
            print("Please provide the QIIME output files in: "+ self.config.qiime_outputs_dir)
            return 0
        try:
            feature_table = pd.read_table(
                self.config.feature_table_dir, sep='\t')
            starting_column=["#OTU ID","featureid"]
            for i in starting_column:
                if i in list(feature_table.columns):
                    feature_table=feature_table.rename(columns={i:"#OTU ID"})
                    break
            rel_abundances = feature_table.iloc[:, list(
                feature_table.columns).index('#OTU ID') + 1:]

        except FileNotFoundError as errormsg:
            rich.print(errormsg)
            rich.print(
                '[red]Feature Table was not found in the directory. Please check the directory and try again!')
            return

        else:
           
            rich.print(
                "[green] ---> feature_table_dir was found in the directory, and was loaded successfully!")
            
            time.sleep(3)

        if download_databases:
            Metagenomics.get_requirements_a2g()

        try:

            taxconomy_table = pd.read_table(
                self.config.taxonomy_table_dir, delimiter='\t', header=0)

        except FileNotFoundError as errormsg:
            print(errormsg)
            print(
                '[red]Taxonomy Table was not found in the directory! Please check the input directory')
            return

        else:

            rich.print(
                "[green] ---> taxonomy_table_dir is found in the directory, and was loaded successfully!")
            time.sleep(3)

        rel_abundances['#OTU ID'] = feature_table['#OTU ID']
        rel_abundances = feature_table.iloc[:, list(
            feature_table.columns).index('#OTU ID') + 1:]
        rel_abundances['#OTU ID'] = feature_table['#OTU ID']
        if sample_name is None:
            Samples = list(feature_table.columns)[list(
                feature_table.columns).index('#OTU ID') + 1:]
        else:
            Samples = sample_name
        featureids = {}
        repSeqs = {}
        taxa = {}

        top_k_taxa = {}
        for i in range(len(taxconomy_table)):
            taxa[taxconomy_table.iloc[i, 0]] = taxconomy_table.iloc[i, 1]
        genome_accessions = {}
        top_genomes = []
        with open(os.path.join(self.config.amplicon2genome_outputs_dir,'Top_k_RepSeq.fasta'), 'w') as f:
        
            for Sample in Samples:

                featureids[Sample] = list(rel_abundances.sort_values(
                Sample, ascending=False)['#OTU ID'].head(self.config.k))

                top_k_taxa[Sample] = list([taxa[Feature]
                                           for Feature in featureids[Sample]])
                [top_genomes.append(RepSeq) for RepSeq in featureids[Sample]]
            top_genomes = list(set(top_genomes))
            with open(self.config.rep_seq_fasta) as D:
                repseqs = Bio_seq.Fasta(D).Fasta_To_Dict()
            for featureid in top_genomes:
                f.write('>'+featureid+'\n'+repseqs[featureid]+'\n')

        for Sample in Samples:
            genome_accessions[Sample] = ['None']*self.config.k
        alignment_dir = os.path.join(self.config.amplicon2genome_outputs_dir,'Alignments')
        subprocess.run([os.path.join(Main_Dir,self.config.vsearch,"vsearch"), '--top_hits_only', '--blast6out', os.path.join(self.config.amplicon2genome_outputs_dir,'matches.blast'), '--usearch_global',
                        os.path.join(self.config.amplicon2genome_outputs_dir,'Top_k_RepSeq.fasta'), '--db', os.path.join(self.config.amplicon2genome_db,"bac120_ssu_reps_r207.fna"), '--id', str(self.config.amplicon2genome_similarity), '--alnout', alignment_dir])
        
        with open(alignment_dir, 'r') as f:
            alignment_dict = {}
            for lines in f:
                if lines.startswith('Query >'):
                    f.readline()
                    alignment_dict[lines.split(
                        '>')[1][:-1]] = re.search("  [A-Z][A-Z]_.*", f.readline()).group(0)[2:]
            
        for Sample in Samples:
            for featureid in featureids[Sample]:
                if featureid in alignment_dict.keys():
                    genome_accessions[Sample][featureids[Sample].index(
                        featureid)] = alignment_dict[featureid]

        Base_NCBI_Dir = 'rsync://ftp.ncbi.nlm.nih.gov/genomes/all/'
        for feature_id in alignment_dict.keys():
            Genome_Out = os.path.join(self.config.amplicon2genome_outputs_dir,alignment_dict[feature_id])
            # if not os.path.exists(Outputs_Dir+alignment_dict[feature_id]):
            #    os.makedirs(Outputs_Dir+alignment_dict[feature_id])

        # I think the next line can be done much more regorously by regexp which is easy

            Specific_NCBI_Dir = alignment_dict[feature_id][3:6]+'/'+alignment_dict[feature_id][7:10]+'/' +\
                alignment_dict[feature_id][10:13]+'/' + \
                alignment_dict[feature_id][13:16]

            subprocess.run(['rsync', '--copy-links', '--times', '--verbose',
                            '--recursive', Base_NCBI_Dir+Specific_NCBI_Dir, self.config.amplicon2genome_outputs_dir])

        pd.DataFrame(featureids).to_csv(os.path.join(self.config.amplicon2genome_outputs_dir,'Top_k_featureids.csv'))
        pd.DataFrame(genome_accessions).to_csv(os.path.join(self.config.amplicon2genome_outputs_dir,'Top_k_genome_accessions.csv'))
        pd.DataFrame(top_k_taxa).to_csv(os.path.join(self.config.amplicon2genome_outputs_dir,'top_k_taxa.csv'))

        if len(list(Path(self.config.amplicon2genome_outputs_dir).rglob('*_genomic.fna.gz'))) > 0:
            for path in Path(self.config.amplicon2genome_outputs_dir).rglob('*_genomic.fna.gz'):
                if "cds" not in path.__str__() and "rna" not in path.__str__():
                    with gzip.open(path, 'rb') as f_in:
                        with open(path.__str__()[:-3], 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
        genomes_dirs = list([path.__str__()
                            for path in Path(self.config.amplicon2genome_outputs_dir).rglob('*_genomic.fna')])
        output_json = {}
        for sample in Samples:
            for genome in genome_accessions[sample]:
                if genome != 'None':
                    if genome not in output_json.keys():
                        output_json[genome] = {}
                        output_json[genome]["NCBI_Name"] = top_k_taxa[sample][genome_accessions[sample].index(
                            genome)]
                        for Items in genomes_dirs:
                            if genome[3:] in Items:
                                output_json[genome]["Genome_Dir"] = Items
        with open(os.path.join(self.config.amplicon2genome_outputs_dir,"Amplicon2Genome_OutInfo.json"), "w") as J:
            json.dump(output_json, J)

        return output_json
    
    def align_genomes(self):
        """
        This is a function that will align genomes to the Protein Database of the ADToolbox using local alignment
        and then will return a dictionary of the alignment results with the aligned genome identifiers as key and the genome as value
        """
        with open(self.config.genomes_json_info, 'r') as f:
            genomes_info = json.load(f)
        valid_genomes = 0
        for genome_ids in track(genomes_info.keys(),description="Fetching Genomes ..."):
            Totall_Genomes = len(genomes_info.keys())
            if os.path.exists(genomes_info[genome_ids]["Genome_Dir"]):
                valid_genomes += 1
        rich.print(f"[green]{valid_genomes}/{Totall_Genomes} Genomes are valid")
        time.sleep(3)

        
        for genome_id in genomes_info.keys():
            subprocess.run(['mmseqs', 'easy-search', genomes_info[genome_id]['Genome_Dir'], self.config.protein_db,
                            os.path.join(self.config.genome_alignment_output,"Alignment_Results_mmseq_"+genomes_info[genome_id]['Genome_Dir'].split("/")[-1][:-4]+".tsv"), 'tmp'])
            genomes_info[genome_id]['Alignment_File'] = os.path.join(self.config.genome_alignment_output, \
                "Alignment_Results_mmseq_" + \
                genomes_info[genome_id]['Genome_Dir'].split(
                    "/")[-1][:-4]+".tsv")

        with open(self.config.genome_alignment_output_json, 'w') as f:
            json.dump(genomes_info, f)

        return genomes_info
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
    
    def extract_relative_abundances(self,sample_names:list[str])->dict:
        
        """
        This function will extract relative abundances for top k taxa from top k taxa table
        and feature table. It will return a dictionary with the top k taxa as keys and the relative abundances as values

        """
        relative_abundances={}
        for sample in sample_names:
            with open(os.path.join(self.config.amplicon2genome_outputs_dir,'Top_k_RepSeq.fasta'),'r') as f:
                Top_k_RepSeq = Bio_seq.Fasta(f).Fasta_To_Dict()
            feature_genome_map={}
            top_k_features = pd.read_table(os.path.join(self.config.amplicon2genome_outputs_dir,'Top_k_featureids.csv'),sep=',')
            top_k_genomes = pd.read_table(os.path.join(self.config.amplicon2genome_outputs_dir,'Top_k_genome_accessions.csv'),sep=',')
            feature_table = pd.read_table(os.path.join(self.config.amplicon2genome_outputs_dir,'feature-table.tsv'),sep='\t')
            starting_column=["#OTU ID","FeatureID"]
            for i in starting_column:
                if i in list(feature_table.columns):
                    feature_table.rename(columns={i:"#OTU ID"},inplace=True)
                    break
            feature_table.set_index('#OTU ID',inplace=True)
            col_ind=top_k_genomes.columns.get_loc(sample)
            
            for row in range(top_k_features.shape[0]):

                if top_k_genomes.iloc[row,col_ind] !="None":
                    feature_genome_map[top_k_features.iloc[row,col_ind]]=top_k_genomes.iloc[row,col_ind]
                else:
                    feature_genome_map[top_k_features.iloc[row,col_ind]]="None"

            valid_genomes = [(feature,feature_genome_map[feature]) for feature in feature_genome_map.keys() if feature_genome_map[feature]!="None"]
            warn(f"{len(feature_genome_map.keys())-len(valid_genomes)} out of {len(feature_genome_map.keys())} features were not assigned to a genome")   
            abundances=dict([(feature[1],feature_table.loc[feature[0],sample]) for feature in valid_genomes])
            relative_abundances[sample]= dict([(genome,abundances[genome]/sum(abundances.values())) for genome in abundances.keys()])
        
        return relative_abundances

    def calculate_microbial_portions(self,microbe_reaction_map:dict,genome_info:dict,relative_abundances:dict)-> dict:
        
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
        for sample in relative_abundances.keys():
        
            reaction_db=pd.read_table(self.config.csv_reaction_db,delimiter=",")
            reaction_db.drop_duplicates(subset=['EC_Numbers'], keep='first',inplace=True)

            reaction_db.set_index("EC_Numbers",inplace=True)
            model_species=list(set(microbe_reaction_map.values()))
            cod_portion=Additive_Dict([(i,0) for i in model_species])

            for genome_id in relative_abundances[sample].keys():
                temp_tsv = pd.read_table(
                            genome_info[genome_id]['Alignment_File'], sep='\t',header=None)
                filter = (temp_tsv.iloc[:, -1] > self.config.bit_score) & (
                            temp_tsv.iloc[:, -2] < self.config.e_value)
                temp_tsv = temp_tsv[filter]
                unique_ecs=list(set([i[1] for i in temp_tsv[1].str.split("|")]))
                cleaned_reaction_list=[]
                [cleaned_reaction_list.extend(map(lambda x:x.strip(" "),reaction_db.loc[ec,"Modified_ADM_Reactions"].split("|"))) for ec in unique_ecs if ec in reaction_db.index]
                pathway_counts=Counter(cleaned_reaction_list)
                SUM=sum([pathway_counts[key] for key in pathway_counts])
                for i in pathway_counts:
                    pathway_counts[i]=pathway_counts[i]/SUM
                microbial_species= dict([(microbe_reaction_map[key],pathway_counts[key]) for key in pathway_counts.keys() ])
                cod_portion=cod_portion+Additive_Dict(microbial_species)*relative_abundances[sample][genome_id]
            
            cod[sample]=cod_portion
        

        return cod
    
    def seqs_from_sra(self,accession:str,save:bool=True,run:bool=True)-> tuple[str,str]:
        """ 
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
        prefetch_script=f"""#!/bin/bash
        prefetch {accession} -O {os.path.join(self.config.sra_work_dir(accession),"seqs")}"""

        fasterq_dump_script=f"""#!/bin/bash
        for i in $(ls ./seqs/) ; do
        fasterq-dump --split-files seqs/$i/*.sra -O seqs/$i/
        rm seqs/$i/*.sra
        done
        """
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
                    qiime2_bash_str[idx]=f"singularity exec --bind  $PWD:$PWD,{str(Path(self.config.qiime_classifier_db))}:$PWD/{str(Path(self.config.qiime_classifier_db).name)}  {self.config.qiime2_singularity_image} " +line
            qiime2_bash_str="\n".join(qiime2_bash_str)
            qiime2_bash_str=qiime2_bash_str.replace("<manifest>",str(manifest_dir.name))
            qiime2_bash_str=qiime2_bash_str.replace("<qiime2_work_dir>",str(qiime2_work_dir.name))
            qiime2_bash_str=qiime2_bash_str.replace("<classifier>",os.path.join(str(Path(self.config.qiime_classifier_db).parent.name),str(Path(self.config.qiime_classifier_db).name)))

            if not paired_end:
                manifest['absolute-filepath']=[str(pathlib.Path("/data")/seqs.name/pathlib.Path(x).parent.name/pathlib.Path(x).name) for x in manifest['absolute-filepath']]
            else:
                manifest['forward-absolute-filepath']=[str(seqs.name/pathlib.Path(x).parent.name/pathlib.Path(x).name) for x in manifest['forward-absolute-filepath']]
                manifest['reverse-absolute-filepath']=[str(seqs.name/pathlib.Path(x).parent.name/pathlib.Path(x).name) for x in manifest['reverse-absolute-filepath']]

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
    metag=Metagenomics(Configs.Metagenomics())
    metag.run_qiime2_from_sra("ERR3861428",container="singularity",save=True,run=False)