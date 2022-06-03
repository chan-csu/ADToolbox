import subprocess
import os
import random
import pandas as pd
import json
import numpy as np
import re
import requests
import time
from requests.adapters import HTTPAdapter
from sympy import Li
from ADToolbox.Parameters import Base_Parameters, Model_Parameters
from requests.packages.urllib3.util.retry import Retry
from requests.exceptions import Timeout
from bs4 import BeautifulSoup
from datetime import datetime
from pathlib import Path
import gzip
import shutil
from . import Configs
from . import Bio_seq
from rich.progress import track
import rich
from . import Main_Dir


class Project:

    """This class can be used for a high level control over a project.
    It is optional to use this class but it provides a nice wholistic view of an integrated project."""

    def __init__(self, Project_Name, Project_Hash, Project_Dir):

        pass


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

    My suggested Model for the ADM Feed inspired by the above link:

                        -------> FSS ----> PI
                     /
            ----> TSS                     -----> Carbs
          /          \                  /
         /             -------> VSS ---> -----> Proteins
    TS =>                                \  
         \                                 -----> Lipids   
          \
           \              -------> FDS ---> SI
            \           /
              -----> TDS                     -----> Carbohydrates
                        \                  /
                          ------> VDS ----> -----> Proteins
                                           \ 
                                             -----> Lipids
    """

    def __init__(self, Carbohydrates, Lipids, Proteins, TS, TSS):
        self.Carbohydrates = Carbohydrates
        self.Lipids = Lipids
        self.Proteins = Proteins
        self.TS = TS
        self.TSS = TSS
        self.TDS = TS-TSS

    def Report(self):
        return "Feed Characteristics:\nCarbohydrates: "+str(self.Carbohydrates)+""+"\nLipids: "+str(self.Lipids)+"\nProteins: "+str(self.Proteins)+"\nFast Degrading portion: "+str(self.X_f)+"\nSouluble Inerts: "+str(self.SI)+"\nParticulate Inerts: "+str(self.PI)


class Sequence_Toolkit:

    AA_List = ['A', 'M', 'N', 'V', 'W', 'L', 'H', 'S', 'G', 'F',
               'K', 'Q', 'E', 'S', 'P', 'I', 'C', 'Y', 'R', 'N', 'D', 'T']
    N_list = ['A', 'C', 'G', 'T']

    def __init__(self, Type='Protein', Length=100):
        self.Type = Type
        self.Length = Length

    def Seq_Generator(self):

        if self.Type == 'Protein':
            return ''.join([random.choice(self.AA_List) for i in range(self.Length)])

        elif self.Type == 'DNA':
            return ''.join([random.choice(self.N_list) for i in range(self.Length)])

        else:
            print('Type not recognized')

    def Mutate_Random(self, seq, Number_of_Mutations=10):
        if self.Type == 'Protein':
            for i in range(Number_of_Mutations):
                seq = list(seq)
                seq[random.randint(0, len(seq))] = random.choice(self.AA_List)
            return ''.join(seq)

        elif self.Type == 'DNA':
            for i in range(Number_of_Mutations):
                seq = list(seq)
                seq[random.randint(0, len(seq))] = random.choice(self.N_list)
            return ''.join(seq)


class Reaction_Toolkit:

    """
    A Class for converting EC numbers to reaction objects or extracting any information from the seed database

    """

    def __init__(self, Compound_DB=Configs.Reaction_Toolkit().Compound_DB, Reaction_DB=Configs.Reaction_Toolkit().Reaction_DB) -> None:
        self.Reaction_DB = Reaction_DB
        self.Compound_DB = Compound_DB

    def Instantiate_Rxns(self, EC_Number, Mode='Single_Match'):
        f = open(self.Reaction_DB)
        Data = json.load(f)
        if Mode == 'Single_Match':
            for i in Data:
                if i['ec_numbers']:
                    if EC_Number in i['ec_numbers']:
                        return Reaction(i)

        elif Mode == 'Multiple_Match':
            Matched_Rxns = [
                Reaction(Rxn) for Rxn in Data if Rxn['ec_numbers'] == [EC_Number]]
            return Matched_Rxns

        else:
            print('Mode not recognized!')

    def Instantiate_Metabs(self, Seed_ID):
        f = open(self.Compound_DB)
        Data = json.load(f)
        # Maybe contain is better than == below\\To be completed
        Matched_Compunds = [Metabolites(Met)
                            for Met in Data if Met['id'] == Seed_ID]
        return Matched_Compunds[0]


class Pathway:
    "A pathway is defined as a list of binary tuples:(Rank of the reaction in the pathway, and the reaction object)"
    # Under construction
    # Check mass balance could be added as a function

    def __init__(self) -> None:
        pass

    def __str__(self, Reaction_List):

        self.Reaction_List = Reaction_List


class Reaction:

    def __init__(self, Dict):
        self.Dict = Dict

    def __str__(self):
        return self.Dict['definition']

    @property
    def Stoichiometry(self):
        S = {}
        for compound in self.Dict['stoichiometry'].split(';'):
            S[compound.split(':')[1]] = compound.split(':')[0]
        return S


class Metabolites:

    def __init__(self, Dict):
        self.Dict = Dict
        self.COD = self.COD_Calc()

    def __str__(self):
        return self.Dict['name']

    def COD_Calc(self):
        # returns the compound's gCOD/g of the compound
        if self:
            Content = {}
            Atoms = ["H", "C", "O"]
            MW = self.Dict['mass']
            for atom in Atoms:
                if re.search(atom+'\d*', self.Dict['formula']):
                    if len(re.search(atom+'\d*', self.Dict['formula']).group()[1:]) == 0:
                        Content[atom] = 1
                    else:
                        Content[atom] = int(
                            re.search(atom+'\d*', self.Dict['formula']).group()[1:])
                else:
                    Content[atom] = 0
            return 1/MW*(Content['H']+4*Content['C']-2*Content['O'])/4*32

        else:
            return 'None'


class Subprocess:

    def __init__(self, Name, Reactions, Metabolites, Parameters, Input):
        self.Name = Name
        self.Reactions = Reactions
        self.Parameters = Parameters
        self.Input = Input
        self.Metabolites = Metabolites

    def Report(self):
        Report_string = f'This instance of {self.Name} includes the following \
            reactions:\n{[Rxn for Rxn in self.Reactions]}'
        return Report_string

    def __str__(self):
        return f'This is an instance of {self.Name} subprocess'

    @property
    def Stoichiometric_Matrix(self):
        S = np.zeros((len(self.Metabolites), len(self.Reactions)))
        for i in range(len(self.Metabolites)):
            for j in range(len(self.Reactions)):
                for keys in self.Reactions[j][0].Stoichiometry:
                    if self.Metabolites[i].Dict['id'] in keys:
                        S[i][j] = self.Reactions[j][0].Stoichiometry[keys]

        return (S, [rxn[0].Dict['name'] for rxn in self.Reactions], [Met.Dict['name'] for Met in self.Metabolites])

    @property
    def Nodes(self):
        return [compounds.Dict['name'] for compounds in self.Metabolites]

    @property
    def Edges(self):
        # After running many examples I think the best clean way to draw
        # A map for the system is to have pathway focus. Reactions of a pathway
        edge = []
        for reaction in self.Reactions:
            Rxn_Name = re.sub(':', ' | ', reaction[0].Dict['name'])
            for met in reaction[0].Stoichiometry.keys():
                Met = Reaction_Toolkit().Instantiate_Metabs(met)
                # This part is not efficient
                edge.append((Rxn_Name, Met.Dict['name']))

        return edge


class Process:

    def __init__(self, Name,
                 Suprocess_List,
                 Input_Seq,
                 Reaction_Metadata=Main_Dir+'/Database/Protein_DB.csv',
                 Protein_DB=Main_Dir+'/Database/Protein_DB.fasta',
                 Alignment_output=Main_Dir+"/Outputs/",
                 Create_DB=True,
                 Align_DB=True):

        self.Name = Name
        self.Input_Seq = Input_Seq
        self.Reaction_Metadata = Reaction_Metadata
        self.Protein_DB = Protein_DB
        self.Alignment_output = Alignment_output
        self.Suprocess_List = Suprocess_List
        self.__create_DB__() if Create_DB else None
        self.__run_diamond__() if Align_DB else None

    def _Create_Diamond_DB(self):
        subprocess.run(['diamond', 'makedb',
                        '--in', self.Protein_DB,
                        '-d', re.sub('.fasta', '', self.Protein_DB)])

    def _Run_Diamond(self, **kwargs):

        if len(kwargs) == 0:
            subprocess.run(['diamond', 'blastp',
                            '-d', re.sub('.fasta', '', self.Protein_DB),
                            '-q', self.Input_Seq,
                            '-o', self.Alignment_output+"Alignment_Results_Diamond_"+self.Name + ".tsv"])

        elif 'Query' in kwargs and 'DB' not in kwargs:
            subprocess.run(['diamond', 'blastp',
                            '-d', re.sub('.fasta', '', self.Protein_DB),
                            '-q', kwargs['Query'],
                            '-o', kwargs['OutDir']])

        elif 'Query' in kwargs and 'DB' in kwargs and "Mode" not in kwargs:
            subprocess.run(['diamond', 'blastp',
                            '-d', kwargs['DB'],
                            '-q', kwargs['Query'],
                            '-o', kwargs['OutDir']])

        elif "Mode" in kwargs and 'DB' in kwargs and 'Query' in kwargs:
            subprocess.run(['diamond', kwargs['Mode'],
                            '-d', kwargs['DB'],
                            '-q', kwargs['Query'],
                            '-o', kwargs['OutDir'], "--very-sensitive"])

        else:
            print('\nThe arguments can not be handled!\n')


    def __str__(self):
        return f'{self.Name} has {len(self.Subprocesses)} Subprocesses: {[ClassNames for ClassNames in self.Subprocesses]} '

    
    def Build_Cobra_Model(self, Name):
        # Finished but needs to be tested
        """This function will build a cobra model from the subprocesses"""
        Subprocesses = self.Subprocesses
        Model = cobra.Model(self.Name)
        Cobra_Metabs = {}
        for Subprocess in Subprocesses:
            Subprocess_Metabs = self.Get_Subprocess_Metabolites(
                Subprocess.Reactions)
            for Metab in Subprocess_Metabs:
                Temp_Metab = cobra.Metabolite(
                    Metab.Dict['id'], compartment=Subprocess.Name, formula=Metab.Dict['formula'], name=Metab.Dict['id'])
                Cobra_Metabs[Metab.Dict['id']] = Temp_Metab
            for Rxns in Subprocess.Reactions:
                # The following line can be changed based on gibbs free energy values
                Temp_Rxn = cobra.Reaction(
                    Rxns.Dict['id'], name=Rxns.Dict['id'], subsystem=Subprocess.Name, lower_bound=-10, upper_bound=10)
                Dict = {}
                for Metab in Rxns.Stoichiometry.keys():

                    Dict[Cobra_Metabs[Metab]] = Rxns.Stoichiometry[Metab]
                Temp_Rxn.add_metabolites(Dict)

                Model.add_reaction(Temp_Rxn)

        cobra.io.write_sbml_model(Model, "test_fbc2.xml")
        return Model


class Database:

    '''
    This class will handle all of the database tasks in ADToolBox
    '''

    def __init__(self, Config=Configs.Database()):

        self.Config = Config




    def _Initialize_Database(self):

        with open(self.Config.Protein_DB, 'w') as f:
            pass
            # The rest is for other possible formats; fastq or gzipped fasta or gzipped fastq

    
    def Filter_Seed_From_EC(self, EC_List,
                            Reaction_DB=Configs.Database().Reaction_DB,
                            Compound_DB=Configs.Database().Compound_DB,
                            Local_Reaction_DB=Configs.Database().Local_Reaction_DB,
                            Local_Compound_DB=Configs.Database().Local_Compound_DB) -> tuple:
        """

        This function takes a list of EC numbers Generates a mini-seed JSON files. This is supposed to
        make the code a lot faster, but makes no difference in terms of outputs, adn won't probably need
        frequent updates.

        """
        with open(Reaction_DB, 'r') as f:
            Main_Reaction_DB = json.load(f)
        with open(Compound_DB, 'r') as f:
            Main_Compound_DB = json.load(f)

        RT = Reaction_Toolkit()
        cached_Compounds = []
        Filtered_Rxns_DB = {}
        Local_Rxn_DB = []
        Local_Comp_DB = []
        Counter = 0
        for EC in EC_List:
            for ind, rxn in enumerate(Main_Reaction_DB):

                if Main_Reaction_DB[ind]['ec_numbers'] != None and EC in Main_Reaction_DB[ind]['ec_numbers']:
                    Local_Rxn_DB.append(rxn)
                    for Mets in rxn["compound_ids"].split(";"):
                        if Mets not in cached_Compounds:
                            cached_Compounds.append(Mets)
            Counter += 1
            print(" -> Percent of ECs processed: ",end=" ")
            print("%"+str(int(Counter/len(EC_List)*100)), end="\r")

        Counter = 0
        for Compound in cached_Compounds:
            for ind, Comp in enumerate(Main_Compound_DB):
                if Compound == Comp["id"]:
                    Local_Comp_DB.append(Comp)
            Counter += 1
            print(" -> Percent of compunds processed: ",end=" ")
            print("%"+str(int(Counter/len(cached_Compounds)*100)), end="\r")

        with open(Local_Reaction_DB, 'w') as f:
            json.dump(Local_Rxn_DB, f)
        with open(Local_Compound_DB, 'w') as f:
            json.dump(Local_Comp_DB, f)

        return (Local_Rxn_DB, Local_Comp_DB)

    # The complete pandas objec can be used as input
    session = requests.Session()
    retry = Retry(connect=3, backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)

    def Add_Protein_from_Uniprot(self, Uniprot_ECs):

        Base_URL = 'https://www.uniprot.org/uniprot/'
        
        with open(self.Config.Protein_DB, 'a') as f:
            for items in track(Uniprot_ECs,description="[yellow]    --> Fetching the protein sequences from Uniprot: "):
                try:
                    file = Database.session.get(
                        Base_URL+items[0]+".fasta", timeout=10)
                except requests.exceptions.ConnectionError:
                    print("Could not fetch the sequence! Trying again ...")
                    time.sleep(10)
                    file = Database.session.get(
                        Base_URL+items[0]+".fasta", timeout=10)
                    if file.ok:
                        print(
                            f"Retry was successful: {items[0]} was fetched successfully!")
                    else:
                        print(
                            f"Retry was unsuccessful: {items[0]} ignored!")
                        # I can add this uniprot id to a log file
                    continue
                if file.ok:
                    f.write(f'>{items[0]}|{items[1]}\n')
                    f.write(''.join(file.text.split('\n')[1:-1]))
                    f.write('\n')
    
    def Uniprots_from_EC(self, EC_list):
        # 5.1.1.2 and 5.1.1.20 should be distinguished: So far it seems like they are distinguished automatically by Uniprot
        Base_URL = 'https://www.uniprot.org/uniprot/?query=ec%3A'
        Uniprots = []

        for EC in track(EC_list,description="[yellow]   --> Fetching Uniprot IDs from ECs"):

            try:
                file = Database.session.get(
                    Base_URL+EC+'+taxonomy%3Abacteria+reviewed%3Ayes&format=list', timeout=1000)

            except requests.exceptions.HTTPError or requests.exceptions.ConnectionError:

                print("Request Error! Trying again ...")
                time.sleep(30)
                file = Database.session.get(
                    Base_URL+EC+'+taxonomy%3Abacteria+reviewed%3Ayes&format=list', timeout=1000)

            if file.ok:
                [Uniprots.append((Uniprot, EC))
                 for Uniprot in file.text.split('\n')[:-1]]  # This alsp does a sanity check

            else:
                print('Something went wrong!')

        return Uniprots

    @staticmethod
    def EC_From_CSV(CSV_File, Sep=','):

        # [ ]- One update could be to make the EC_Numbers column to be insensitive to the case

        try:

            EC_List = pd.read_table(CSV_File, sep=Sep)
            Output_EC_list = list(EC_List["EC_Numbers"])
            assert len(
                Output_EC_list) > 0, "The EC list is empty; Check the file"

        except FileNotFoundError:

            print("TSV file not found")

        except (pd.errors.ParserError, KeyError):

            print(
                "CSV file not in the correct format; Check the delimiter and column names!")

        except AssertionError as error:

            print(error)

        else:


            return Output_EC_list

    def EC_From_Uniprot(self, Uniprot_ID):

        Base_URL = 'https://www.uniprot.org/uniprot/'

        try:
            file = Database.session.get(
                Base_URL+Uniprot_ID+".txt", timeout=10)

        except requests.exceptions.ConnectionError:

            print("Request Error! Trying again ...")
            time.sleep(30)
            file = Database.session.get(
                Base_URL+Uniprot_ID+".txt", timeout=10)

            if file.ok:

                print("Retry was successful!")

        for line in file:
            if re.search("EC=[0-9]+", line.decode("utf-8")):
                EC = re.sub("EC=", "", re.search(
                    "EC=[.0-9]+", line.decode("utf-8")).group(0))
                break

        else:

            print("Retry was unsuccessful!")
            return None

        return EC

    Cazy_links = ["http://www.cazy.org/Glycoside-Hydrolases.html",
                  "http://www.cazy.org/Polysaccharide-Lyases.html",
                  "http://www.cazy.org/Carbohydrate-Esterases.html"
                  ]

    def Cazy_EC_from_link(self):
        '''
        Extracts the EC numbers from a link to the Cazy website

        '''

        EC_list = []
        for link in Database.Cazy_links:

            page = requests.get(link)
            soup = BeautifulSoup(page.content, "html.parser")
            results = soup.find("div", class_="cadre_principal").find_all(
                "th", class_="thec")
            for EC_number in results:
                if '-' not in EC_number.text.strip() and '.' in EC_number.text.strip():

                    EC_list.append(EC_number.text.strip())

        print("EC numbers extracted from Cazy website!")
        return EC_list

    @staticmethod
    def Seed_From_EC(EC_Number, Mode='Single_Match'):

        f = open(Main_Dir+'/reactions.json')

        Data = json.load(f)

        if Mode == 'Single_Match':

            for i in Data:

                if i['ec_numbers']:

                    if EC_Number in i['ec_numbers']:

                        return i["id"]

        elif Mode == 'Multiple_Match':

            Matched_Rxns = list(set(list([

                Rxn["id"] for Rxn in Data if Rxn['ec_numbers'] == [EC_Number]])))

            return Matched_Rxns

        else:
            print('Mode not recognized!')

    @staticmethod
    def Instantiate_Rxn_From_Seed(Seed_ID_List):
        Rxn_List = []
        with open(Main_Dir+'/reactions.json') as f:

            Data = json.load(f)

            for Seed_ID in Seed_ID_List:

                for i in Data:

                    if i['id']:

                        if Seed_ID in i['id']:

                            Rxn_List.append(Reaction(i))
                            break

        return Rxn_List

    @staticmethod
    def Metadata_From_EC(EC_list):
        full_table = {"EC_Numbers": [], "Seed_Ids": [],
                      "Reaction_Names": [], "Pathways": []}
        print("##---Finding EC Metadata:\n")
        for idx, ec in enumerate(EC_list):
            Seed_ID = Database.Seed_From_EC(ec, Mode='Multiple_Match')
            full_table['EC_Numbers'].append(ec)
            full_table['Seed_Ids'].append(Seed_ID)
            Temp_rxns = Database.Instantiate_Rxn_From_Seed(Seed_ID)
            Temp_rxn_Names = list([reaction.Dict["name"]
                                   for reaction in Temp_rxns])
            Temp_rxn_Path = list([reaction.Dict["pathways"]
                                  for reaction in Temp_rxns])
            full_table["Pathways"].append(Temp_rxn_Path)
            full_table["Reaction_Names"].append(Temp_rxn_Names)

        full_table = pd.DataFrame(full_table)

        return full_table

    ###-----FeedStock_Database_Functions-----###
    def Init_Feedstock_Database(self):
        '''
        Makes an empty feedstock database json file.
        BE CAREFUL: This will overwrite the current database file.
        '''
        Dec=input("Are you sure you want to create a new database? (y/n): ")
        if Dec=="y":

            with open(self.Config.Feed_DB, 'w') as f:
                json.dump([], f)
        else:
            print("Aborting!")
            return
        ## TO-DO: Error and Failure handling

    def _Add_Feedstock_To_Database(self, Feedstock):
        
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
            with open(self.Config.Feed_DB, 'r') as f:
                DataBase=json.load(f)
        except FileNotFoundError:
            print("Feedstock database not found! Creating new one...")
            self.Init_Feedstock_Database()
        
        finally:
        
            DataBase.append(Feedstock)
            with open(self.Config.Feed_DB, 'w') as f:
                json.dump(DataBase, f)
        
    def Add_Feedstock_To_Database_From_File(self, File_Path):
        '''
        Adds a feedstock to the feedstock database from a csv file.
        A Feed Stock spreadsheet must be a table with the following columns:
        "Name", "TS", "TSS", "Lipids", "Proteins", "Carbohydrates", "PI", "SI", "Notes"
        File_Path does not come with a default path, since it is advised to not populate the
        project directory with feedstock spreadsheets unwanted.
        '''
        try:
            f=pd.read_table(File_Path, sep=',')
        except FileNotFoundError:
            print("File not found!")
            return
        Counter=0
        for i in f.index:
            Feedstock=f.loc[i,:].to_dict()
            self._Add_Feedstock_To_Database(Feedstock)
            Counter+=1
        print("{} Feedstocks added to the database!".format(Counter))

    @staticmethod
    def Download_ADM1_Parameters(config):
        '''
        Downloads the parameters needed for running ADM Models in ADToolbox
        
        '''
        if not os.path.exists(config.Base_dir):
            os.makedirs(config.Base_dir)
        Model_Parameters_dir=config.Model_Parameters
        Base_Parameters_dir=config.Base_Parameters
        Initial_Conditions_dir=config.Initial_Conditions
        Inlet_Conditions_dir=config.Inlet_Conditions
        Reactions_dir=config.Reactions
        Species_dir=config.Species
        Model_Parameters="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/ADM1/ADM1_Model_Parameters.json"
        Base_Parameters="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/ADM1/ADM1_Base_Parameters.json"
        Initial_Conditions="https://github.com/ParsaGhadermazi/Database/blob/main/ADToolbox/ADM1/ADM1_Initial_Conditions.json"
        Inlet_Conditions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/ADM1/ADM1_Inlet_Conditions.json"
        Reactions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/ADM1/ADM1_Reactions.json"
        Species="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/ADM1/ADM1_Species.json"
        r = requests.get(Model_Parameters, allow_redirects=True)
        with open(Model_Parameters_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(Base_Parameters, allow_redirects=True)
        with open(Base_Parameters_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(Initial_Conditions, allow_redirects=True)
        with open(Initial_Conditions_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(Inlet_Conditions, allow_redirects=True)
        with open(Inlet_Conditions_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(Reactions, allow_redirects=True)
        with open(Reactions_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(Species, allow_redirects=True)
        with open(Species_dir, 'wb') as f:
            f.write(r.content)
    @staticmethod
    def Download_Modified_ADM_Parameters(config):    
        '''
        Downloads the parameters needed for running ADM Models in ADToolbox
        '''
        if not os.path.exists(config.Base_dir):
            os.makedirs(config.Base_dir)
        Model_Parameters_dir=config.Model_Parameters
        Base_Parameters_dir=config.Base_Parameters
        Initial_Conditions_dir=config.Initial_Conditions
        Inlet_Conditions_dir=config.Inlet_Conditions
        Reactions_dir=config.Reactions
        Species_dir=config.Species
        Model_Parameters="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Model_Parameters.json"
        Base_Parameters="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Base_Parameters.json"
        Initial_Conditions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Initial_Conditions.json"
        Inlet_Conditions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Inlet_Conditions.json"
        Reactions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Reactions.json"
        Species="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Modified_ADM/Modified_ADM_Species.json"
        r = requests.get(Model_Parameters, allow_redirects=True)
        with open(Model_Parameters_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(Base_Parameters, allow_redirects=True)
        with open(Base_Parameters_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(Initial_Conditions, allow_redirects=True)
        with open(Initial_Conditions_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(Inlet_Conditions, allow_redirects=True)
        with open(Inlet_Conditions_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(Reactions, allow_redirects=True)
        with open(Reactions_dir, 'wb') as f:
            f.write(r.content)
        r = requests.get(Species, allow_redirects=True)
        with open(Species_dir, 'wb') as f:
            f.write(r.content)
        
        
    @staticmethod
    def Download_Seed_Databases(directory):
        Reactions="https://github.com/ModelSEED/ModelSEEDDatabase/raw/master/Biochemistry/reactions.json"
        r = requests.get(Reactions, allow_redirects=True)
        with open(directory, 'wb') as f:
            f.write(r.content)        
    @staticmethod
    def Download_Protein_Database(directory):
        Protein_DB="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Protein_DB.fasta"
        r = requests.get(Protein_DB, allow_redirects=True)
        with open(directory, 'wb') as f:
            f.write(r.content)
        
    @staticmethod
    def Download_Reaction_Database(directory):
        Reactions="https://github.com/ParsaGhadermazi/Database/raw/main/ADToolbox/Reaction_Metadata.csv"
        r = requests.get(Reactions, allow_redirects=True)
        with open(directory, 'wb') as f:
            f.write(r.content)
    @staticmethod
    def Download_Feed_Database(directory):
        Feed="https://github.com/ParsaGhadermazi/Database/blob/main/ADToolbox/Feed_DB.json"
        r = requests.get(Feed, allow_redirects=True)
        with open(directory, 'wb') as f:
            f.write(r.content)
class Report:

    """
    This Class will be used for generating the final reports of any kind
    Genome_Alignment_JSON={Genome_ID1, Genome_ID2, ...}
    Genome_ID={{"NCBI_Name":"THE NCBI TAXONOMY NAME"},{"Alingment_TSV_File":"THE TSV FILE NAME"}}}}}} 

    """

    def __init__(self, Configuration):
        """
        This class is initiated with a Config module/class.
        A config module/class is a module/class that contains all the necessary information about the pipeline.

        """

        self.Configs = Configuration

    def ADM_From_Alignment_JSON(self,ADM_Rxns,Model="Modified_ADM_Reactions"):
        RT = Reaction_Toolkit(Reaction_DB=self.Configs.Seed_RXN_DB)
        with open(self.Configs.Alignment_JSON_Info) as f:
            JSON_Report = json.load(f)
        Reaction_DB = pd.read_table(self.Configs.RXN_DB, sep=',')
        JSON_ADM_Output = {}
        for ADM_Reaction in track([Rxn for Rxn in ADM_Rxns],description="Finding Alignments to ADM Reactions"):
            JSON_ADM_Output[ADM_Reaction] = {}
            for Genome_ID in JSON_Report.keys():
                if os.path.exists(JSON_Report[Genome_ID]['Alignment_File']):
                    JSON_ADM_Output[ADM_Reaction][JSON_Report[Genome_ID]
                                                  ['NCBI_Name']] = {}
                    Temp_TSV = pd.read_table(
                        JSON_Report[Genome_ID]['Alignment_File'], sep='\t')
                    A_Filter = (Temp_TSV.iloc[:, -1] > self.Configs.bit_score) & (
                        Temp_TSV.iloc[:, -2] < self.Configs.e_value)
                    Temp_TSV = Temp_TSV[A_Filter]
                    ECs = [item.split('|')[1] for item in Temp_TSV.iloc[:, 1]]
                    Filter = (Reaction_DB['EC_Numbers'].isin(ECs)) & (
                        Reaction_DB[Model] == ADM_Reaction)
                    Filtered_Rxns = Reaction_DB[Filter]
                    ECs = Filtered_Rxns['EC_Numbers'].tolist()

                    EC_Dict = {}
                    for EC in ECs:
                        EC_Dict[EC] = []
                        L = RT.Instantiate_Rxns(EC, "Multiple_Match")
                        [EC_Dict[EC].append(
                            rxn.__str__()) for rxn in L if rxn.__str__() not in EC_Dict[EC]]

                    JSON_ADM_Output[ADM_Reaction][JSON_Report[Genome_ID]
                                                  ['NCBI_Name']] = EC_Dict


        with open(self.Configs.ADM_From_Alignment_JSON_Output, 'w') as f:
            json.dump(JSON_ADM_Output, f)

        return JSON_ADM_Output
        # So far it has gotten really complicated, after making sure it works, instead of adding EC numbers I'll add [SEED,STR] so
        # that it can be used for downstream analysis


class Metagenomics:

    """
    A Core class that handles the metganomics data
    Init this with a Metagenomics class of Config module
    
    
    
    It takes three diferent types of input:


    1- QIIMEII Outputs -> Refer to the documentation for the required files and format
    2- Shotgun Short reads -> Refer to the documentation for the required files and format
    3- MAGs -> Refer to the documentation for the required files and format

    This class will use other classes to generate the required files for downstream analysis.
    """
    def __init__(self,Config):
        self.Config=Config
    
    def Get_Requirements_A2G(self):
        """
        This function will automatically download the required files for Amplicon to Genome functionality.
        """
        if not os.path.exists(self.Config.Amplicon2Genome_DB):
            os.mkdir(self.Config.Amplicon2Genome_DB)

        url = {'Version': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/VERSION',
               'MD5SUM': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/MD5SUM',
               'FILE_DESCRIPTIONS': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/FILE_DESCRIPTIONS',
               'metadata_field_desc': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/metadata_field_desc.tsv',
               'bac120_metadata': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tar.gz',
               'bac120_ssu': 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/bac120_ssu_reps.tar.gz'
               }

        for keys in ['Version', 'MD5SUM', 'FILE_DESCRIPTIONS']:
            r = requests.get(url[keys], allow_redirects=True)
            open(os.path.join(self.Config.Amplicon2Genome_DB,keys+'.txt'), 'wb').write(r.content)
        r = requests.get(url['metadata_field_desc'], allow_redirects=True)
        open(os.path.join(self.Config.Amplicon2Genome_DB,'metadata_field_desc.tsv'), 'wb').write(r.content)

        for keys in ['bac120_metadata', 'bac120_ssu']:
            r = requests.get(url[keys], allow_redirects=True)
            open(os.path.join(self.Config.Amplicon2Genome_DB,keys+'.tar.gz'), 'wb').write(r.content)


    def Amplicon2Genome(self,
                        Requirements=False,
                        Report=True):

        if not os.path.exists(self.Config.QIIME_Outputs_Dir):
            os.mkdir(self.Config.Config.QIIME_Outputs_Dir)
            print("Please provide the QIIME output files in: "+self.Config.Config.QIIME_Outputs_Dir)
            return 0
        try:
            Feature_Table = pd.read_table(
                self.Config.Feature_Table_Dir, sep='\t', header=1)
            Rel_Abundances = Feature_Table.iloc[:, list(
                Feature_Table.columns).index('#OTU ID') + 1:]

        except FileNotFoundError as errormsg:
            rich.print(errormsg)
            rich.print(
                '[red]Feature Table was not found in the directory. Please check the directory and try again!')
            return

        else:
           
            rich.print(
                "[green] ---> Feature_Table_Dir was found in the directory, and was loaded successfully!")
            
            time.sleep(3)

        if Requirements:
            Metagenomics.Get_Requirements()

        try:

            Taxconomy_Table = pd.read_table(
                self.Config.Taxonomy_Table_Dir, delimiter='\t', header=0)

        except FileNotFoundError as errormsg:
            print(errormsg)
            print(
                '[red]Taxonomy Table was not found in the directory! Please check the input directory')
            return

        else:

            rich.print(
                "[green] ---> Taxonomy_Table_Dir is found in the directory, and was loaded successfully!")
            time.sleep(3)

        Rel_Abundances['#OTU ID'] = Feature_Table['#OTU ID']
        Rel_Abundances = Feature_Table.iloc[:, list(
            Feature_Table.columns).index('#OTU ID') + 1:]
        Rel_Abundances['#OTU ID'] = Feature_Table['#OTU ID']
        Samples = list(Feature_Table.columns)[list(
            Feature_Table.columns).index('#OTU ID') + 1:]
        FeatureIDs = {}
        RepSeqs = {}
        Taxa = {}

        Top_K_Taxa = {}
        for i in range(len(Taxconomy_Table)):
            Taxa[Taxconomy_Table.iloc[i, 0]] = Taxconomy_Table.iloc[i, 1]
        Genome_Accessions = {}
        Top_genomes = []
        f = open(os.path.join(self.Config.Amplicon2Genome_Outputs_Dir,'Top_k_RepSeq.fasta'), 'w')
        for Sample in Samples:
            FeatureIDs[Sample] = list(Rel_Abundances.sort_values(
                Sample, ascending=False)['#OTU ID'].head(self.Config.K))
            Top_K_Taxa[Sample] = list([Taxa[Feature]
                                       for Feature in FeatureIDs[Sample]])
            [Top_genomes.append(RepSeq) for RepSeq in FeatureIDs[Sample]]
        Top_genomes = list(set(Top_genomes))
        with open(self.Config.Rep_Seq_Fasta) as D:
            Repseqs = Bio_seq.Fasta(D).Fasta_To_Dict()
        for FeatureID in Top_genomes:
            f.write('>'+FeatureID+'\n'+Repseqs[FeatureID]+'\n')

        for Sample in Samples:
            Genome_Accessions[Sample] = ['None']*self.Config.K
        Alignment_Dir = os.path.join(self.Config.Amplicon2Genome_Outputs_Dir,'Alignments')
        subprocess.run([os.path.join(Main_Dir,"bin","vsearch"), '--top_hits_only', '--blast6out', os.path.join(self.Config.Amplicon2Genome_Outputs_Dir,'matches.blast'), '--usearch_global',
                        self.Config.Amplicon2Genome_Outputs_Dir+'/Top_k_RepSeq.fasta', '--db', os.path.join(self.Config.Amplicon2Genome_DB,"bac120_ssu_reps_r207.fna"), '--id', str(self.Config.Amplicon2Genome_Similarity), '--alnout', Alignment_Dir])
        f.close()
        f = open(Alignment_Dir, 'r')
        Alignment_Dict = {}
        for lines in f:
            if lines.startswith('Query >'):
                f.readline()
                Alignment_Dict[lines.split(
                    '>')[1][:-1]] = re.search("  [A-Z][A-Z]_.*", f.readline()).group(0)[2:]
        f.close()
        for Sample in Samples:
            for FeatureID in FeatureIDs[Sample]:
                if FeatureID in Alignment_Dict.keys():
                    Genome_Accessions[Sample][FeatureIDs[Sample].index(
                        FeatureID)] = Alignment_Dict[FeatureID]

        Base_NCBI_Dir = 'rsync://ftp.ncbi.nlm.nih.gov/genomes/all/'
        for Feature_ID in Alignment_Dict.keys():
            Genome_Out = os.path.join(self.Config.Amplicon2Genome_Outputs_Dir,Alignment_Dict[Feature_ID])
            # if not os.path.exists(Outputs_Dir+Alignment_Dict[Feature_ID]):
            #    os.makedirs(Outputs_Dir+Alignment_Dict[Feature_ID])

        # I think the next line can be done much more regorously by regexp which is easy

            Specific_NCBI_Dir = Alignment_Dict[Feature_ID][3:6]+'/'+Alignment_Dict[Feature_ID][7:10]+'/' +\
                Alignment_Dict[Feature_ID][10:13]+'/' + \
                Alignment_Dict[Feature_ID][13:16]

            subprocess.run(['rsync', '--copy-links', '--times', '--verbose',
                            '--recursive', Base_NCBI_Dir+Specific_NCBI_Dir, self.Config.Amplicon2Genome_Outputs_Dir])

        pd.DataFrame(FeatureIDs).to_csv(
            os.path.join(Main_Dir,'Outputs','SelectedFeatures.csv'))
        pd.DataFrame(Genome_Accessions).to_csv(
            os.path.join(Main_Dir,'Outputs','GenomeAccessions.csv'))
        pd.DataFrame(Top_K_Taxa).to_csv(os.path.join(Main_Dir,'Outputs','TopKTaxa.csv'))

        if len(list(Path(self.Config.Amplicon2Genome_Outputs_Dir).rglob('*_genomic.fna.gz'))) > 0:
            for path in Path(self.Config.Amplicon2Genome_Outputs_Dir).rglob('*_genomic.fna.gz'):
                if "cds" not in path.__str__() and "rna" not in path.__str__():
                    with gzip.open(path, 'rb') as f_in:
                        with open(path.__str__()[:-3], 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
        Genomes_Dirs = list([path.__str__()
                            for path in Path(self.Config.Amplicon2Genome_Outputs_Dir).rglob('*_genomic.fna')])
        Output_JSON = {}
        for sample in Samples:
            for genome in Genome_Accessions[sample]:
                if genome != 'None':
                    if genome not in Output_JSON.keys():
                        Output_JSON[genome] = {}
                        Output_JSON[genome]["NCBI_Name"] = Top_K_Taxa[sample][Genome_Accessions[sample].index(
                            genome)]
                        for Items in Genomes_Dirs:
                            if genome[3:] in Items:
                                Output_JSON[genome]["Genome_Dir"] = Items
        with open(os.path.join(self.Config.Amplicon2Genome_Outputs_Dir,"Amplicon2Genome_OutInfo.json"), "w") as J:
            json.dump(Output_JSON, J)

        return Output_JSON
    
    def Align_Genomes(self):
        """
        #This is a function that will align genomes to the Protein Database of the ADToolbox using local alignment
        #and then will return a dictionary of the alignment results with the aligned genome identifiers as key and the genome as value
        """
        with open(self.Config.Genomes_JSON_Info, 'r') as f:
            Genomes_Info = json.load(f)
        Valid_Genomes = 0
        for Genome_IDs in track(Genomes_Info.keys(),description='[yellow]Validating the genomes'):
            Totall_Genomes = len(Genomes_Info.keys())
            if os.path.exists(Genomes_Info[Genome_IDs]["Genome_Dir"]):
                Valid_Genomes += 1
        rich.print(f"[green]{Valid_Genomes}/{Totall_Genomes} Genomes are valid")
        time.sleep(3)

        if self.Config.Aligner == 'mmseqs2':
            Aligned_Genomes = {}
            for Genome_ID in track(Genomes_Info.keys(),description='[yellow]Aligning the genomes'):

                subprocess.run(['mmseqs', 'easy-search', Genomes_Info[Genome_ID]['Genome_Dir'], self.Config.Protein_DB,
                                os.path.join(self.Config.Genome_Alignment_Output,"Alignment_Results_mmseq_"+Genomes_Info[Genome_ID]['Genome_Dir'].split("/")[-1][:-4]+".tsv"), 'tmp'])

                Genomes_Info[Genome_ID]['Alignment_File'] = os.path.join(self.Config.Genome_Alignment_Output, \
                    "Alignment_Results_mmseq_" + \
                    Genomes_Info[Genome_ID]['Genome_Dir'].split(
                        "/")[-1][:-4]+".tsv")

                Alignment_Results_Genome = pd.read_table(
                    Genomes_Info[Genome_ID]['Alignment_File'], header=None)

                Filter = (Alignment_Results_Genome.iloc[:, -1] >= self.Config.bit_score) & (
                    Alignment_Results_Genome.iloc[:, -2] <= self.Config.e_value)

                ECs_Matched = list(set(list(Alignment_Results_Genome[Filter].iloc[j, 1].split(
                    '|')[1] for j in range(len(Alignment_Results_Genome[Filter])))))

                Genomes_Info[Genome_ID]['ECs_Matched'] = ECs_Matched

                Matched_Rxns = {}

                for EC in ECs_Matched:
                    for rxn in Reaction_Toolkit().Instantiate_Rxns(EC, Mode='Multiple_Match'):

                        Matched_Rxns[rxn.Dict["id"]] = rxn.__str__()

                Genomes_Info[Genome_ID]['Rxns_Matched'] = Matched_Rxns
                # Filter the results
        with open(os.path.join(self.Config.Genome_Alignment_Output,"Alignment_Info.json"), 'w') as f:
            json.dump(Genomes_Info, f)

        return Genomes_Info



# class Visualization:
#     """This class will handle any visualization tasks"""

#     def Process_Graph(self):
#         g = gz.Digraph('G', filename='Process_Graph.gv', engine='fdp',
#                        node_attr={'color': 'lightblue2', 'style': 'filled'})
#         g.attr(compound='true')
#         Clusters = [('A'+str(i+1), i) for i in range(len(self.Subprocesses))]
#         for Cluster in Clusters:
#             with g.subgraph(name=Cluster[0]) as A:
#                 A.attr(label=self.Subprocesses[Cluster[1]].Name, color='black')
#                 A.edges(self.Subprocesses[Cluster[1]].Edges)
#                 # A.nodes=self.Subprocesses[Cluster[1]].Nodes
#         g.view()

#         # Add universal map
#         return g


class ADM_Mapping:

    """
    This class deals with mapping Metagenomics and Exp data to ADM1

    """
    Degrader_Reaction_Map = {
        'Hydrolysis carbohydrates': "X_ch",
        'Hydrolysis of proteins': "X_pr",
        'Hydrolysis of lipids': "X_li",
        'Uptake of sugars': 'X_su',
        'Uptake of amino acids': "X_aa",
        'Uptake of LCFA': "X_fa",
        'Uptake of valerate': "X_c4",
        'Uptake of butyrate': "X_c5",
        'Uptake of propionate': "X_pro",
        'Uptake of acetate': "X_ac",
        'Uptake of Hydrogen': "X_h2",
    }

    def __init__(self, Alignment_Info_JSON, ADM1_Parameters, Rel_Abundances_Dict):

        self.ADM1_Parameters = ADM1_Parameters
        self.Alignment_Info_JSON = Alignment_Info_JSON
        self.Rel_Abundances_Dict = Rel_Abundances_Dict

    def Update_Degrader_CODs(self, Inoculum_Total_COD, Degrader_Reaction_Map=Degrader_Reaction_Map):
        """
        This function updates the CODs of the degrader reactions
        """
        Degrader_CODs = {}
        for Degrader in Degrader_Reaction_Map.keys():
            Degrader_CODs[Degrader] = Inoculum_Total_COD[Degrader_Reaction_Map[Degrader]]
        return Degrader_CODs



class Universal_Model:
    pass

    # To-Do:
    # ------------------------------------------
    # 3. More systems biology stuff             | [ ]
    # ------------------------------------------
    # 4. ADM1 Improvement                       | [ ]
    # ------------------------------------------
    # 6. Compare with prokka or other tools     | [ ]
    # ------------------------------------------
    # 7. Build a map for the process            | [I think finally it would be something like graphviz rather than escher]
    # ------------------------------------------
    # 8. Metadata for reactions                 | []
    # ------------------------------------------
    # 9. Maybe for now graphviz is fine         | []
    # ------------------------------------------
    # 10. Metagenomics class for handling       |
    # Any type of metagenomics data, and the    | []
    # Downstream analysis.                      |
    # ------------------------------------------
    # 11. Visualization classMaybe for now      | []
    #  graphviz is fine                         |
    # ------------------------------------------
    # 12. Report class                          | []
    # ------------------------------------------
    # 13. Project class                         | [Very Low Priority]
    # ------------------------------------------

### TEST BLOCK ###


if __name__ == "__main__":
    # MG_C=Configs.Metagenomics()
    # MG=Metagenomics(MG_C)
    # MG.Align_Genomes()
    # R=Report(Configs.Report())
    # R.ADM_From_Alignment_JSON(Parameters_Modified.Reactions)
    pass
    

