import argparse
import Configs
import ADToolBox
from __init__ import __version__
from __init__ import Main_Dir
from rich.console import Console
import rich
from ADM import *
from rich.table import Table
from rich.prompt import Prompt
from rich import markdown
import pdb
console = Console()


console.rule("[bold red]ADToolBox")
class AParser(argparse.ArgumentParser):
    def _print_message(self, message, file=None):
        rich.print(message, file=file)

def main():
    DB_Class=Database(Config=Configs.Database())
    parser = AParser(prog="ADToolBox",
        description="[italic yellow]ADToolBox, a Toolbox for Anaerobic Digestion Modeling",
        epilog="ADToolBox-Chan Lab.")
    parser.add_argument("-v", "--version", action="version",help="Prints the version of ADToolBox") 
    parser.version = f"[green]ADToolBox v{__version__}"
    subparsers = parser.add_subparsers(dest="ADToolbox_Module", help='ADToolbox Modules:')
    subparser_Database = subparsers.add_parser('Database', help='Database Module of ADToolBox')
    db_subp=subparser_Database.add_subparsers(dest="Database_Module", help='Database Modules:')
    db_subp.add_parser("initialize-feed-db", help="Initialize the Feed DB")
    ExFDB=db_subp.add_parser("extend-feed-db", help="Extend the Feed Database using a CSV file")
    ExFDB.add_argument("-d", "--dir", help="CSV file to be used for extending the Feed DB",required=True)
    db_subp.add_parser("show-feed-db", help="Display the Current Feed Database")
    db_subp.add_parser("download-reaction-db", help="Downloads the reaction database in CSV format")
    db_subp.add_parser("build-protein-db", help="Generates the protein database for ADToolbox")
    db_subp.add_parser("download-feed-db", help="Downloads the feed database in JSON format")
    db_subp.add_parser("download-protein-db", help="Downloads the protein database in fasta format; You can alternatively build it from reaction database.")
    db_subp.add_parser("download-amplicon-to-genome-dbs", help="downloads amplicon to genome databases")
    
    
    
    
    
    ### Metagenomics Module ###
    subparser_Metagenomics= subparsers.add_parser('Metagenomics', help='Metagenomics Module of ADToolBox')
    MetaG_subp=subparser_Metagenomics.add_subparsers(dest='MetaG_Subparser',help='[yellow] Available Metagenomics Commands:')
    MetaG_subp_2=MetaG_subp.add_parser('Amplicon-to-Genome', help='Downloads the representative genome from each amplicon')
    MetaG_subp_2.add_argument("-i", "--Input-dir", action="store", help="Input the directory containing all the required [bold]QIIME output tables")
    MetaG_subp_2.add_argument("-o", "--Output-dir", action="store", help="Output the directory to store the representative genome files")
    MetaG_subp_2.add_argument("-a", "--Amplicon-to-Genome-DB", action="store", help="The Amplicon to Genome Database to use")
    MetaG_subp_2.add_argument("--k", action="store", help="Top k genomes to be selected",default=10)
    MetaG_subp_2.add_argument("--similarity", action="store", help="Similarity cutoff in the 16s V4 region for genome selection",default=97)
    MetaG_subp_3=MetaG_subp.add_parser('Align-Genomes', help='Align Genomes to the protein database of ADToolbox, or any other fasta with protein sequences')
    MetaG_subp_3.add_argument("-i", "--input-file", action="store", help="Input the address of the JSON file includeing information about the genomes to be aligned")
    MetaG_subp_3.add_argument("-d", "--protein-db-dir", action="store", help="Directory containing the protein database to be used for alignment")
    MetaG_subp_3.add_argument("-o", "--output-dir", action="store", help="Output the directory to store the alignment results")
    MetaG_subp_3.add_argument("-b", "--bit-score", action="store", help="Minimum Bit Score cutoff for alignment",default=40)
    MetaG_subp_3.add_argument("-e", "--e-value", action="store", help="Minimum e-vlaue score cutoff for alignment",default=100)
    MetaG_subp_4=MetaG_subp.add_parser("make-json-from-genomes",help="Generates JSON file required by Align-Genomes for custom genomes.")
    
    MetaG_subp_4.add_argument("-i", "--input-file", action="store", help="Input the address of the CSV file includeing information about the genomes to be aligned",required=True)
    MetaG_subp_4.add_argument("-o", "--output-file", action="store", help="Output the directory to store the JSON file.",required=True)



    Doc_Parser=subparsers.add_parser('Documentations', help='Documentations for using AD Toolbox')
    Doc_Parser.add_argument("-s", "--show", action="store_true", help="Documentation for a specific module")    
    
    ### ADM module ###
    subparser_ADM = subparsers.add_parser('ADM', help='ADM Module of ADToolBox')
    ADM_subp=subparser_ADM.add_subparsers(dest='ADM_Subparser',help='Available ADM Commands:')
    subparser_ADM1 = ADM_subp.add_parser('Original-ADM1',help='Original ADM 1:')
    subparser_ADM1.add_argument("--Model-Parameters", action="store", help="Model parameters for ADM 1")
    subparser_ADM1.add_argument("--Base-Parameters", action="store", help="Provide json file with base parameters for original ADM1")
    subparser_ADM1.add_argument("--Initial-Conditions", action="store", help="Provide json file with initial conditions for original ADM1")
    subparser_ADM1.add_argument("--Inlet-Conditions", action="store", help="Provide json file with inlet conditions for original ADM1")
    subparser_ADM1.add_argument("--Reactions", action="store", help="Provide json file with reactions for original ADM1")
    subparser_ADM1.add_argument("--Species", action="store", help="Provide json file with species for original ADM1")
    subparser_ADM1.add_argument("--Metagenome-Report", action="store", help="Provide json file with metagenome report for original ADM1")
    subparser_ADM1.add_argument("--Report", action="store", help="Describe how to report the results of original ADM1. Current options are: 'dash' and 'csv'")
    Mod_ADM_subp=ADM_subp.add_parser('Modified-ADM',help='Modified ADM:')
    Mod_ADM_subp.add_argument("--Model-Parameters", action="store", help="Model parameters for Modified ADM")
    Mod_ADM_subp.add_argument("--Base-Parameters", action="store", help="Provide json file with base parameters for modified ADM")
    Mod_ADM_subp.add_argument("--Initial-Conditions", action="store", help="Provide json file with initial conditions for modified ADM")
    Mod_ADM_subp.add_argument("--Inlet-Conditions", action="store", help="Provide json file with inlet conditions for modified ADM")
    Mod_ADM_subp.add_argument("--Reactions", action="store", help="Provide json file with reactions for modified ADM")
    Mod_ADM_subp.add_argument("--Species", action="store", help="Provide json file with species for modified ADM")
    Mod_ADM_subp.add_argument("--Metagenome-Report", action="store", help="Provide json file with metagenome report for modified ADM")
    Mod_ADM_subp.add_argument("--Report", action="store", help="Describe how to report the results of modified ADM. Current options are: 'dash' and 'csv'")
    ADM_subp.add_parser('show-escher-map',help='makes an escher map for modified ADM')
    ####

    # subparser_ADM.add_argument( "Modified-ADM-dash", help="Runs the dash app in the browser",action="store_true")

    # Mod_ADM_args=subparser_ADM.add_subparsers(help='Modified ADM Args:')


    subparser_Report = subparsers.add_parser('Report', help='Report Module of ADToolBox')
    subparser_Utility = subparsers.add_parser('Utility', help='Utility Module of ADToolBox')
    subparser_Configs = subparsers.add_parser('Configs', help='Configurations of ADToolBox')
    Conf_subp=subparser_Configs.add_subparsers(dest='Configs_Subparser',help='Available Configs Commands:')
    subparser_Configs1 = Conf_subp.add_parser('set-base-dir',help='Determine the address of the base directory for ADToolBox to work with')
    subparser_Configs2= Conf_subp.add_parser('build-folder-structure',help='Builds the folder structure for ADToolBox to work properly')
    subparser_Configs3= Conf_subp.add_parser('download-all-databases',help='Downloads all the databases for ADToolBox to work properly, and puts them in the right directory in Databases')
    subparser_Configs4= Conf_subp.add_parser('download-escher-files',help='Downloads all files required for running the escher map functionality of ADToolBox')
    subparser_Configs1.add_argument("-d", "--directory", action="store", help="Provide the address of the base directory for ADToolBox to work with",required=True)


    args=parser.parse_args()
    ### Block for running ADM Module ###
        ### Block for running Original ADM1 ###
    #### Metagenomics Module #####
    if args.ADToolbox_Module == 'Metagenomics' and "MetaG_Subparser" in args and args.MetaG_Subparser=="Align-Genomes":
        if "input_file" in args and "output_dir" in args:
            Meta_Config=Configs.Metagenomics(Genomes_JSON_Info=args.input_file,Genome_Alignment_Output=args.output_dir)
            ADToolBox.Metagenomics(Meta_Config).Align_Genomes()


    
    if args.ADToolbox_Module == 'Metagenomics' and "MetaG_Subparser" in args and args.MetaG_Subparser=="make-json-from-genomes":
        ADToolBox.Metagenomics.Make_JSON_from_Genomes(args.input_file,args.output_file)





    if args.ADToolbox_Module == 'Documentations' and args.show:
        with open(Configs.Documentation().ReadMe,'r') as f:
            T=f.read()
            console.print(markdown.Markdown(T))
    
    if args.ADToolbox_Module == 'Database' and "Database_Module" in args and args.Database_Module=="build-protein-db" :
            DB_Class
            ECs=DB_Class.EC_From_CSV(Configs.Database().CSV_Reaction_DB)
            rich.print(u"[bold green]\u2713 All EC numbers were extracted!\n")
            Uniprot_IDS=DB_Class.Uniprots_from_EC(ECs)
            rich.print(u"[bold green]\u2713 All Uniprot IDs were extracted!\n")
            DB_Class._Initialize_Database()
            DB_Class.Add_Protein_from_Uniprot(Uniprot_IDS)
            rich.print(u"[bold green]\u2713 Protein Database was built!\n")
    
    
    if args.ADToolbox_Module == 'Database' and  "Database_Module" in args and args.Database_Module=="download-feed-db":
        DB_Class.Download_Feed_Database(DB_Class.Config.Feed_DB)
        rich.print(u"[bold green]\u2713 Feed Database was downloaded!\n")
    
    if args.ADToolbox_Module == 'Database' and  "Database_Module" in args and args.Database_Module=="download-reaction-db":
        DB_Class.Download_Reaction_Database(DB_Class.Config.Reaction_DB)
        rich.print(u"[bold green]\u2713 Reaction Database was downloaded!\n")
    
    if args.ADToolbox_Module == 'Database' and  "Database_Module" in args and args.Database_Module=="download-protein-db":
        DB_Class.Download_Protein_Database(DB_Class.Config.Protein_DB)
        rich.print(u"[bold green]\u2713 Protein Database was downloaded!\n")
    
    if args.ADToolbox_Module == 'Database' and  "Database_Module" in args and args.Database_Module=="show-feed-db":
        with open(Configs.Database().Feed_DB, 'r') as f:
            Feed_DB = json.load(f)
        ColNames=Feed_DB[0].keys()
        Feed_table = Table(title="Feed Database",expand=True,safe_box=True)

        for i in ColNames:
            Feed_table.add_column(i, justify="center", style="cyan", no_wrap=True)
        for i in Feed_DB:
            Feed_table.add_row(*map(str,list(i.values())))
        console.print(Feed_table)
    
    if "Initialize_Feed_DB" in args and bool(args.Initialize_Feed_DB):
        ADToolBox.Database.Init_Feedstock_Database(DB_Class)
    if "Extend_Feed_DB" in args and bool(args.Extend_Feed_DB):
        DB_Class.Add_Feedstock_To_Database_From_File(args.Extend_Feed_DB)
    
    if "ADToolbox_Module" in args and args.ADToolbox_Module == 'ADM' and "ADM_Subparser" in args and args.ADM_Subparser == 'Original-ADM1':
        if args.Model_Parameters:
            with open(args.Model_Parameters) as f:
                ADM_Model_Parameters=json.load(f)
        else:
            with open(Configs.Original_ADM1().Model_Parameters) as f:
                ADM_Model_Parameters=json.load(f)
        if args.Base_Parameters:
            with open(args.Base_Parameters) as f:
                ADM_Base_Parameters=json.load(f)
        else:
            with open(Configs.Original_ADM1().Base_Parameters) as f:
                ADM_Base_Parameters=json.load(f)

        if args.Initial_Conditions:
            with open(args.Initial_Conditions) as f:
                ADM_Initial_Conditions=json.load(f)
        else:
            with open(Configs.Original_ADM1().Initial_Conditions) as f:
                ADM_Initial_Conditions=json.load(f)

        
        if args.Inlet_Conditions:
            with open(args.Inlet_Conditions) as f:
                ADM_Inlet_Conditions=json.load(f)
        else:
            with open(Configs.Original_ADM1().Inlet_Conditions) as f:
                ADM_Inlet_Conditions=json.load(f)

        
        if args.Reactions:
            with open(args.Reactions) as f:
                ADM_Reactions=json.load(f)
        else:
            with open(Configs.Original_ADM1().Reactions) as f:
                ADM_Reactions=json.load(f)
        
        if args.Species:
            with open(args.Species) as f:
                ADM_Species=json.load(f)
        else:
            with open(Configs.Original_ADM1().Species) as f:
                ADM_Species=json.load(f)
        
        if args.Metagenome_Report:
            with open(args.Metagenome_Report) as f:
                ADM_Metagenome_Report=json.load(f)
        else:
            with open(Configs.Original_ADM1().Metagenome_Report) as f:
                ADM_Metagenome_Report=json.load(f)
        
        ADM1 = Model(ADM_Model_Parameters, ADM_Base_Parameters, ADM_Initial_Conditions,
                 ADM_Inlet_Conditions, ADM_Reactions, ADM_Species, ADM1_ODE_Sys, Build_ADM1_Stoiciometric_Matrix, Metagenome_Report=ADM_Metagenome_Report, Name="ADM1", Switch="DAE")
        Sol = ADM1.Solve_Model(
        (0, 20), ADM1.Initial_Conditions[:, 0], np.linspace(0, 20, 100))

        if args.Report == 'dash' or args.Report==None:
            ADM1.Dash_App(Sol)
        elif args.Report == 'csv':
            Address=Prompt.ask("\n[yellow]Where do you want to save the csv file? ")
            ADM1.CSV_Report(Sol,Address)
        else:
            print('Please provide a valid report option')
    

    elif "ADToolbox_Module" in args and args.ADToolbox_Module == 'ADM' and "ADM_Subparser" in args and args.ADM_Subparser == 'Modified-ADM':
        if args.Model_Parameters:
            with open(args.Model_Parameters) as f:
                ADM_Model_Parameters=json.load(f)
        else:
            with open(Configs.Modified_ADM().Model_Parameters) as f:
                ADM_Model_Parameters=json.load(f)
        if args.Base_Parameters:
            with open(args.Base_Parameters) as f:
                ADM_Base_Parameters=json.load(f)
        else:
            with open(Configs.Modified_ADM().Base_Parameters) as f:
                ADM_Base_Parameters=json.load(f)

        if args.Initial_Conditions:
            with open(args.Initial_Conditions) as f:
                ADM_Initial_Conditions=json.load(f)
        else:
            with open(Configs.Modified_ADM().Initial_Conditions) as f:
                ADM_Initial_Conditions=json.load(f)

        
        if args.Inlet_Conditions:
            with open(args.Inlet_Conditions) as f:
                ADM_Inlet_Conditions=json.load(f)
        else:
            with open(Configs.Modified_ADM().Inlet_Conditions) as f:
                ADM_Inlet_Conditions=json.load(f)

        
        if args.Reactions:
            with open(args.Reactions) as f:
                ADM_Reactions=json.load(f)
        else:
            with open(Configs.Modified_ADM().Reactions) as f:
                ADM_Reactions=json.load(f)
        
        if args.Species:
            with open(args.Species) as f:
                ADM_Species=json.load(f)
        else:
            with open(Configs.Modified_ADM().Species) as f:
                ADM_Species=json.load(f)
        
        if args.Metagenome_Report:
            with open(args.Metagenome_Report) as f:
                ADM_Metagenome_Report=json.load(f)
        else:
            if os.path.exists(Configs.Modified_ADM().Metagenome_Report):
                with open(Configs.Modified_ADM().Metagenome_Report) as f:
                    ADM_Metagenome_Report=json.load(f)
            else:
                ADM_Metagenome_Report=None
                
        mod_adm1 = Model(ADM_Model_Parameters,ADM_Base_Parameters,ADM_Initial_Conditions,ADM_Inlet_Conditions,ADM_Reactions,
                    ADM_Species, Modified_ADM1_ODE_Sys, Build_Modified_ADM1_Stoiciometric_Matrix,Metagenome_Report=ADM_Metagenome_Report,Name="Modified_ADM1", Switch="DAE")
        Sol_mod_adm1 = mod_adm1.Solve_Model(
        (0, 100), mod_adm1.Initial_Conditions[:, 0], np.linspace(0, 100, 10000))

        if args.Report == 'dash' or args.Report==None:
            mod_adm1.Dash_App(Sol_mod_adm1)
        elif args.Report == 'csv':
            Address=Prompt.ask("\n[yellow]Where do you want to save the csv file? ")
            mod_adm1.CSV_Report(Sol_mod_adm1,Address)
        else:
            print('Please provide a valid report option')
    
    elif "ADToolbox_Module" in args and args.ADToolbox_Module == 'ADM' and "ADM_Subparser" in args and args.ADM_Subparser == 'show-escher-map':
        from http.server import HTTPServer, CGIHTTPRequestHandler
        os.chdir(os.path.join(Main_Dir,'Visualizations','escher'))
        server_object = HTTPServer(server_address=('', 80), RequestHandlerClass=CGIHTTPRequestHandler)
        rich.print('[green]Started the server:')
        rich.print('[yellow]Copy the following in your browser: http://localhost:80/')
        server_object.serve_forever()


        

    if "Modified_ADM_dash" in args and bool(args.Modified_ADM_dash):
        Report = os.path.join(Main_Dir, "..", "Reports",
                          "ADM_From_Alignment_JSON_Output.json")
        with open(Report, 'r') as j:
            Report = json.load(j)
        mod_adm1 = Model(PM.Model_Parameters, PM.Base_Parameters, PM.Initial_Conditions, PM.Inlet_Conditions, PM.Reactions,
                     PM.Species, Modified_ADM1_ODE_Sys, Build_Modified_ADM1_Stoiciometric_Matrix, Name="Modified_ADM1", Switch="DAE")
        Sol_mod_adm1 = mod_adm1.Solve_Model(
        (0, 100), mod_adm1.Initial_Conditions[:, 0], np.linspace(0, 100, 10000))
        mod_adm1.Dash_App(Sol_mod_adm1)
    

    #############################
    #### Configuration Block ####
    #############################

    if args.ADToolbox_Module == 'Configs' and args.Configs_Subparser=='set-base-dir':
        Set_Base_Dir(args)

    if args.ADToolbox_Module == 'Configs' and args.Configs_Subparser=='build-folder-structure':
        Build_Folder_Structure()

    if args.ADToolbox_Module == 'Configs' and args.Configs_Subparser=='download-all-databases':
        Download_All_Databases()
    
    if args.ADToolbox_Module == 'Configs' and args.Configs_Subparser=='download-escher-files':
        Download_Escher_Files()






###########################
#####Configs functions#####
###########################

def Set_Base_Dir(args):
    """
    This function sets the base directory for the ADToolbox.
    """
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),"ADToolbox_Configs.json"), 'r') as f:
        Conf = json.load(f)
    Conf["Base_Dir"] = args.directory
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),"ADToolbox_Configs.json"), 'w') as f:
        json.dump(Conf, f, indent=4)
    if not os.path.exists(args.directory):
        os.mkdir(args.directory)
        rich.print("[yellow]Directory did not exist. Created directory: [green]{}".format(args.directory))

def Build_Folder_Structure():
    """
    This function will build the folder structure undrestandable by the ADToolbox.
    """
    Required_Folders=["Database","Outputs","Reports","Metagenomics_Data","Genomes","Visualizations"]
    for i in Required_Folders:
        if not os.path.exists(os.path.join(Main_Dir,i)):
            os.mkdir(os.path.join(Main_Dir,i))
            rich.print(f"\n[yellow] {i} directory did not exist. Created directory: [green]{os.path.join(Main_Dir,i)}\n")

def Download_All_Databases():
    """
    This function will download all the databases required by the ADToolbox.
    """
    MetaG=ADToolBox.Metagenomics(Configs.Metagenomics())
    rich.print(u"[yellow] Downloading Amplicon to Genome databases ...\n")
    MetaG.Get_Requirements_A2G()
    rich.print(u"[bold green]\u2713 Amplicon to Genome databases were downloaded successfuly!\n")
    rich.print(u"[yellow] Downloading Seed database ...\n")
    ADToolBox.Database().Download_Seed_Databases(Configs.Seed_RXN_DB)
    rich.print(u"[bold green]\u2713 Seed database was downloaded successfuly!\n")
    rich.print(u"[yellow] Downloading Original-ADM1 database ...\n")
    ADToolBox.Database().Download_ADM1_Parameters(Configs.Original_ADM1())
    rich.print(u"[bold green]\u2713 Original-ADM1 database was downloaded successfuly!\n")
    rich.print(u"[yellow] Downloading Modified-ADM database ...\n")
    ADToolBox.Database().Download_Modified_ADM_Parameters(Configs.Modified_ADM())
    rich.print(u"[bold green]\u2713 Modified-ADM database was downloaded successfuly!\n")
    rich.print(u"[yellow] Downloading protein database ...\n")
    ADToolBox.Database().Download_Protein_Database(Configs.Database().Protein_DB)
    rich.print(u"[bold green]\u2713 Protein database was downloaded successfuly!\n")
    rich.print(u"[yellow] Downloading Reaction database ...\n")
    ADToolBox.Database().Download_Reaction_Database(Configs.Database().Reaction_DB)
    rich.print(u"[bold green]All the databases were downloaded successfuly!\n")
    rich.print(u"[yellow] Downloading the feed database ...\n")
    ADToolBox.Database().Download_Feed_Database(Configs.Database().Feed_DB)
    rich.print(u"[bold green]\u2713 Feed database was downloaded successfuly!\n")

def Download_Escher_Files():
    """
    This function will download the files necessary to load the escher map correcly.
    """
    if not os.path.exists(os.path.join(Main_Dir,"Visualizations","escher")):
        os.mkdir(os.path.join(Main_Dir,"Visualizations","escher"))
    rich.print(u"[yellow] Downloading Escher files ...\n")
    ADToolBox.Database().Download_Escher_Files(Configs.Database().Escher_Files)
    rich.print(u"[bold green]\u2713 Escher files were downloaded successfuly!\n")


if __name__ == "__main__":
    main()