import argparse
from re import M
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
console = Console()


console.rule("[bold red]ADToolBox")
class AParser(argparse.ArgumentParser):
    def _print_message(self, message, file=None):
        rich.print(message, file=file)

def main():
    db_class=Database(config=Configs.Database())
    parser = AParser(prog="ADToolBox",
        description="[italic yellow]ADToolBox, a Toolbox for Anaerobic Digestion Modeling",
        epilog="ADToolBox-Chan Lab.")
    parser.add_argument("-v", "--version", action="version",help="Prints the version of ADToolBox") 
    parser.version = f"[green]ADToolBox v{__version__}"
    subparsers = parser.add_subparsers(dest="ADToolbox_Module", help='ADToolbox Modules:')
    subparser_database = subparsers.add_parser('Database', help='Database Module of ADToolBox')
    db_subp=subparser_database.add_subparsers(dest="Database_Module", help='Database Modules:')
    db_subp.add_parser("initialize-feed-db", help="Initialize the Feed DB")
    exfdb=db_subp.add_parser("extend-feed-db", help="Extend the Feed Database using a CSV file")
    exfdb.add_argument("-d", "--dir", help="CSV file to be used for extending the Feed DB",required=True)
    db_subp.add_parser("show-feed-db", help="Display the Current Feed Database")
    db_subp.add_parser("download-reaction-db", help="Downloads the reaction database in CSV format")
    db_subp.add_parser("download-seed-reaction-db", help="Downloads the seed reaction database in JSON format")
    db_subp.add_parser("build-protein-db", help="Generates the protein database for ADToolbox")
    db_subp.add_parser("download-feed-db", help="Downloads the feed database in JSON format")
    db_subp.add_parser("download-protein-db", help="Downloads the protein database in fasta format; You can alternatively build it from reaction database.")
    db_subp.add_parser("download-amplicon-to-genome-dbs", help="downloads amplicon to genome databases")
    
    
    
    
    
    ### Metagenomics Module ###
    meta_config_defult=Configs.Metagenomics()
    subparser_metagenomics= subparsers.add_parser('Metagenomics', help='Metagenomics Module of ADToolBox')
    metag_subp=subparser_metagenomics.add_subparsers(dest='metag_Subparser',help='[yellow] Available Metagenomics Commands:')
    
    metag_subp_2=metag_subp.add_parser('amplicon-to-genome', help='Downloads the representative genome from each amplicon')
    metag_subp_2.add_argument("-q", "--qiime-outputs-dir", action="store", help="Input the directory to the QIIME outputs",default=meta_config_defult.qiime_outputs_dir)
    metag_subp_2.add_argument("-f", "--feature-table-dir", action="store", help="Input the directory to the feature table output from [bold]QIIME output tables",default=meta_config_defult.feature_table_dir)
    metag_subp_2.add_argument("-r", "--rep-Seq-dir", action="store", help="Input the directory to the repseq fasta output from [bold]QIIME output files",default=meta_config_defult.rep_seq_fasta)
    metag_subp_2.add_argument("-t", "--taxonomy-table-dir", action="store", help="Input the directory to the taxonomy table output from [bold]QIIME output files",default=meta_config_defult.taxonomy_table_dir)
    metag_subp_2.add_argument("-o", "--output-dir", action="store", help="Output the directory to store the representative genome files",default=meta_config_defult.amplicon2genome_outputs_dir)
    metag_subp_2.add_argument("-a", "--amplicon-to-genome-db", action="store", help="The Amplicon to Genome Database to use",default=meta_config_defult.amplicon2genome_db)
    metag_subp_2.add_argument("--k", action="store", help="Top k genomes to be selected",default=meta_config_defult.k,type=int)
    metag_subp_2.add_argument("--similarity", action="store", help="Similarity cutoff in the 16s V4 region for genome selection",default=meta_config_defult.amplicon2genome_similarity,type=float)
    
    metag_subp_3=metag_subp.add_parser('align-genomes', help='Align Genomes to the protein database of ADToolbox, or any other fasta with protein sequences')
    metag_subp_3.add_argument("-i", "--input-file", action="store", help="Input the address of the JSON file includeing information about the genomes to be aligned",default=meta_config_defult.genomes_json_info)
    metag_subp_3.add_argument("-d", "--protein-db-dir", action="store", help="Directory containing the protein database to be used for alignment",default=meta_config_defult.protein_db)
    metag_subp_3.add_argument("-o", "--output-dir", action="store", help="Output the directory to store the alignment results",default=meta_config_defult.genome_alignment_output)
    metag_subp_3.add_argument("-b", "--bit-score", action="store", help="Minimum Bit Score cutoff for alignment",default=meta_config_defult.bit_score)
    metag_subp_3.add_argument("-e", "--e-value", action="store", help="Minimum e-vlaue score cutoff for alignment",default=meta_config_defult.e_value)
    
    metag_subp_4=metag_subp.add_parser("make-json-from-genomes",help="Generates JSON file required by Align-Genomes for custom genomes.")
    metag_subp_4.add_argument("-i", "--input-file", action="store", help="Input the address of the CSV file includeing information about the genomes to be aligned",required=True)
    metag_subp_4.add_argument("-o", "--output-file", action="store", help="Output the directory to store the JSON file.",required=True)
    
    metag_subp_5=metag_subp.add_parser('map-genomes-to-adm', help='maps JSON file with genome infromation to ADM reactions')
    metag_subp_5.add_argument("-i","--input-file",action="store",help="Input the address of the JSON file includeing information about the alignment of the genomes to the protein database",default=meta_config_defult.genome_alignment_output_json)
    metag_subp_5.add_argument("-m","--model",action="store",help="Model determines which mapping system you'd like to use; Current options: 'Modified_adm_reactions'",default="Modified_adm_reactions")
    metag_subp_5.add_argument("-o","--output-dir",action="store",help="address to store the JSON report to be loaded with a model",default=meta_config_defult.genome_adm_map_json)


    doc_parser=subparsers.add_parser('Documentations', help='Documentations for using AD Toolbox')
    doc_parser.add_argument("-s", "--show", action="store_true", help="Documentation for a specific module")    
    
    ### ADM module ###
    subparser_adm = subparsers.add_parser('ADM', help='ADM Module of ADToolBox')
    adm_subp=subparser_adm.add_subparsers(dest='adm_Subparser',help='Available ADM Commands:')
    subparser_adm1 = adm_subp.add_parser('original-adm1',help='Original ADM1:')
    subparser_adm1.add_argument("--model-parameters", action="store", help="Model parameters for ADM 1")
    subparser_adm1.add_argument("--base-parameters", action="store", help="Provide json file with base parameters for original ADM1")
    subparser_adm1.add_argument("--initial-conditions", action="store", help="Provide json file with initial conditions for original ADM1")
    subparser_adm1.add_argument("--inlet-conditions", action="store", help="Provide json file with inlet conditions for original ADM1")
    subparser_adm1.add_argument("--reactions", action="store", help="Provide json file with reactions for original ADM1")
    subparser_adm1.add_argument("--species", action="store", help="Provide json file with species for original ADM1")
    subparser_adm1.add_argument("--metagenome-report", action="store", help="Provide json file with metagenome report for original ADM1")
    subparser_adm1.add_argument("--report", action="store", help="Describe how to report the results of original ADM1. Current options are: 'dash' and 'csv'")
    mod_adm_subp=adm_subp.add_parser('modified-adm',help='Modified ADM:')
    mod_adm_subp.add_argument("--model-parameters", action="store", help="Model parameters for Modified ADM")
    mod_adm_subp.add_argument("--base-parameters", action="store", help="Provide json file with base parameters for modified ADM")
    mod_adm_subp.add_argument("--initial-conditions", action="store", help="Provide json file with initial conditions for modified ADM")
    mod_adm_subp.add_argument("--inlet-conditions", action="store", help="Provide json file with inlet conditions for modified ADM")
    mod_adm_subp.add_argument("--reactions", action="store", help="Provide json file with reactions for modified ADM")
    mod_adm_subp.add_argument("--species", action="store", help="Provide json file with species for modified ADM")
    mod_adm_subp.add_argument("--metagenome-report", action="store", help="Provide json file with metagenome report for modified ADM")
    mod_adm_subp.add_argument("--control-states", action="store", help="Provide a json file that contains the control states, and their values ")
    mod_adm_subp.add_argument("--report", action="store", help="Describe how to report the results of modified ADM. Current options are: 'dash' and 'csv'")
    ####

    # subparser_adm.add_argument( "Modified-ADM-dash", help="Runs the dash app in the browser",action="store_true")

    # Mod_ADM_args=subparser_adm.add_subparsers(help='Modified ADM Args:')


    subparser_Report = subparsers.add_parser('Report', help='Report Module of ADToolBox')
    subparser_Utility = subparsers.add_parser('Utility', help='Utility Module of ADToolBox')
    subparser_Configs = subparsers.add_parser('Configs', help='Configurations of ADToolBox')
    conf_subp=subparser_Configs.add_subparsers(dest='Configs_Subparser',help='Available Configs Commands:')
    subparser_configs1 = conf_subp.add_parser('set-base-dir',help='Determine the address of the base directory for ADToolBox to work with')
    subparser_configs2= conf_subp.add_parser('build-folder-structure',help='Builds the folder structure for ADToolBox to work properly')
    subparser_configs3= conf_subp.add_parser('download-all-databases',help='Downloads all the databases for ADToolBox to work properly, and puts them in the right directory in Databases')
    subparser_configs1.add_argument("-d", "--directory", action="store", help="Provide the address of the base directory for ADToolBox to work with",required=True)


    args=parser.parse_args()
    ### Block for running ADM Module ###
        ### Block for running Original ADM1 ###
    #### Metagenomics Module #####
    if args.ADToolbox_Module == 'Metagenomics' and "metag_Subparser" in args and args.metag_Subparser=="amplicon-to-genome":
        meta_config_defult.feature_table_dir=args.feature_table_dir
        meta_config_defult.qiime_outputs_dir=args.qiime_outputs_dir
        meta_config_defult.rep_seq_fasta=args.rep_Seq_dir
        meta_config_defult.taxonomy_table_dir=args.taxonomy_table_dir
        meta_config_defult.amplicon2genome_outputs_dir=args.output_dir
        meta_config_defult.amplicon2genome_db=args.amplicon_to_genome_db
        meta_config_defult.k=args.k
        meta_config_defult.amplicon2genome_similarity=args.similarity
        ADToolBox.Metagenomics(meta_config_defult).amplicon2genome()

    if args.ADToolbox_Module == 'Metagenomics' and "metag_Subparser" in args and args.metag_Subparser=="map-genomes-to-adm":
        meta_config_defult.genome_alignment_output_json=args.input_file
        model_reactions=args.model
        meta_config_defult.genome_adm_map_json=args.output_dir
        if model_reactions=="Modified_adm_reactions":
            reactions=Configs.Modified_ADM().reactions
            with open(reactions,"r") as f:
                reactions=json.load(f)
            ADToolBox.Metagenomics(meta_config_defult).adm_from_alignment_json(reactions,Model=model_reactions)
        else:
            rich.print(f"[red] No Reaction mapping exists for {model_reactions}; execution aborted!")

        

    if args.ADToolbox_Module == 'Metagenomics' and "metag_Subparser" in args and args.metag_Subparser=="align-genomes":
        if "input_file" in args and "output_dir" in args:
            Meta_Config=Configs.Metagenomics(genomes_json_info=args.input_file,genome_alignment_output=args.output_dir,genome_alignment_output_json=os.path.join(args.output_dir,"Alignment_Info.json"))
            ADToolBox.Metagenomics(Meta_Config).Align_Genomes()

    
    if args.ADToolbox_Module == 'Metagenomics' and "metag_Subparser" in args and args.metag_Subparser=="make-json-from-genomes":
        ADToolBox.Metagenomics.make_json_from_genomes(args.input_file,args.output_file)





    if args.ADToolbox_Module == 'Documentations' and args.show:
        with open(Configs.Documentation().ReadMe,'r') as f:
            T=f.read()
            console.print(markdown.Markdown(T))
    
    if args.ADToolbox_Module == 'Database' and "Database_Module" in args and args.Database_Module=="build-protein-db" :
            ecs=db_class.ec_from_csv(Configs.Database().csv_reaction_db)
            rich.print(u"[bold green]\u2713 All EC numbers were extracted!\n")
            uniprot_ids=db_class.uniprots_from_ec(ecs)
            rich.print(u"[bold green]\u2713 All Uniprot IDs were extracted!\n")
            db_class._initialize_database()
            db_class.add_protein_from_uniprot(uniprot_ids)
            rich.print(u"[bold green]\u2713 Protein Database was built!\n")
    
    
    if args.ADToolbox_Module == 'Database' and  "Database_Module" in args and args.Database_Module=="download-feed-db":
        db_class.download_feed_database(db_class.config.feed_db)
        rich.print(u"[bold green]\u2713 Feed Database was downloaded!\n")
    
    if args.ADToolbox_Module == 'Database' and  "Database_Module" in args and args.Database_Module=="download-reaction-db":
        db_class.download_reaction_database(db_class.config.csv_reaction_db)
        rich.print(u"[bold green]\u2713 Reaction Database was downloaded!\n")
    
    if args.ADToolbox_Module == 'Database' and  "Database_Module" in args and args.Database_Module=="download-protein-db":
        db_class.download_protein_database(db_class.config.protein_db)
        rich.print(u"[bold green]\u2713 Protein Database was downloaded!\n")
    
    if args.ADToolbox_Module == 'Database' and  "Database_Module" in args and args.Database_Module=="download-seed-reaction-db":
        db_class.download_seed_databases(db_class.config.reaction_db)
        rich.print(u"[bold green]\u2713 SEED Reaction Database was downloaded!\n")
    
    if args.ADToolbox_Module == 'Database' and  "Database_Module" in args and args.Database_Module=="show-feed-db":
        with open(Configs.Database().feed_db, 'r') as f:
            feed_db = json.load(f)
        colnames=feed_db[0].keys()
        feed_table = Table(title="Feed Database",expand=True,safe_box=True)

        for i in colnames:
            feed_table.add_column(i, justify="center", style="cyan", no_wrap=True)
        for i in feed_db:
            feed_table.add_row(*map(str,list(i.values())))
        console.print(feed_table)
    
    if "Initialize_feed_db" in args and bool(args.initialize_feed_db):
        ADToolBox.Database.init_feedstock_database(db_class)
    if "extend_feed_db" in args and bool(args.extend_feed_db):
        db_class.add_feedstock_to_database_from_file(args.extend_feed_db)
    
    if "ADToolbox_Module" in args and args.ADToolbox_Module == 'ADM' and "adm_Subparser" in args and args.adm_Subparser == 'original-adm1':
        if args.model_parameters:
            with open(args.model_parameters) as f:
                asm_model_parameters=json.load(f)
        else:
            with open(Configs.Original_ADM1().model_parameters) as f:
                asm_model_parameters=json.load(f)
        if args.base_parameters:
            with open(args.base_parameters) as f:
                adm_base_parameters=json.load(f)
        else:
            with open(Configs.Original_ADM1().base_parameters) as f:
                adm_base_parameters=json.load(f)

        if args.initial_conditions:
            with open(args.initial_conditions) as f:
                adm_initial_conditions=json.load(f)
        else:
            with open(Configs.Original_ADM1().initial_conditions) as f:
                adm_initial_conditions=json.load(f)

        
        if args.inlet_conditions:
            with open(args.inlet_conditions) as f:
                adm_inlet_conditions=json.load(f)
        else:
            with open(Configs.Original_ADM1().inlet_conditions) as f:
                adm_inlet_conditions=json.load(f)

        
        if args.reactions:
            with open(args.reactions) as f:
                adm_reactions=json.load(f)
        else:
            with open(Configs.Original_ADM1().reactions) as f:
                adm_reactions=json.load(f)
        
        if args.species:
            with open(args.species) as f:
                adm_species=json.load(f)
        else:
            with open(Configs.Original_ADM1().species) as f:
                adm_species=json.load(f)
        
        if args.metagenome_report:
            with open(args.metagenome_report) as f:
                adm_metagenome_report=json.load(f)
        else:
            with open(Configs.Original_ADM1().metagenome_report) as f:
                adm_metagenome_report=json.load(f)
        
        ADM1 = Model(asm_model_parameters, adm_base_parameters, adm_initial_conditions,
                 adm_inlet_conditions, adm_reactions, adm_species, adm1_ode_sys, build_adm1_stoiciometric_Matrix, metagenome_report=adm_metagenome_report, Name="ADM1", Switch="DAE")
        Sol = ADM1.solve_model(
        (0, 20), ADM1.initial_conditions[:, 0], np.linspace(0, 20, 100))

        if args.report == 'dash' or args.report==None:
            ADM1.dash_app(Sol)
        elif args.report == 'csv':
            Address=Prompt.ask("\n[yellow]Where do you want to save the csv file? ")
            ADM1.csv_report(Sol,Address)
        else:
            print('Please provide a valid report option')
    

    elif "ADToolbox_Module" in args and args.ADToolbox_Module == 'ADM' and "adm_Subparser" in args and args.adm_Subparser == 'modified-adm':
        if args.model_parameters:
            with open(args.model_parameters) as f:
                adm_model_parameters=json.load(f)
        else:
            with open(Configs.Modified_ADM().model_parameters) as f:
                adm_model_parameters=json.load(f)
        if args.base_parameters:
            with open(args.base_parameters) as f:
                adm_base_parameters=json.load(f)
        else:
            with open(Configs.Modified_ADM().base_parameters) as f:
                adm_base_parameters=json.load(f)

        if args.initial_conditions:
            with open(args.initial_conditions) as f:
                adm_initial_conditions=json.load(f)
        else:
            with open(Configs.Modified_ADM().initial_conditions) as f:
                adm_initial_conditions=json.load(f)

        
        if args.inlet_conditions:
            with open(args.inlet_conditions) as f:
                adm_inlet_conditions=json.load(f)
        else:
            with open(Configs.Modified_ADM().inlet_conditions) as f:
                adm_inlet_conditions=json.load(f)

        
        if args.reactions:
            with open(args.reactions) as f:
                adm_reactions=json.load(f)
        else:
            with open(Configs.Modified_ADM().reactions) as f:
                adm_reactions=json.load(f)
        
        if args.species:
            with open(args.species) as f:
                adm_species=json.load(f)
        else:
            with open(Configs.Modified_ADM().species) as f:
                adm_species=json.load(f)
        
        if args.metagenome_report:
            with open(args.metagenome_report) as f:
                adm_metagenome_report=json.load(f)
        else:
            if os.path.exists(Configs.Modified_ADM().metagenome_report):
                with open(Configs.Modified_ADM().metagenome_report) as f:
                    adm_metagenome_report=json.load(f)
            else:
                adm_metagenome_report=None
        if args.control_states:
            with open(args.control_states) as f:
                adm_control_states=json.load(f)
        else:
            
            adm_control_states={}
                
        mod_adm1 = Model(adm_model_parameters, adm_base_parameters, adm_initial_conditions, adm_inlet_conditions, adm_reactions,
                     adm_species, modified_adm_ode_sys, build_modified_adm_stoichiometric_matrix,control_state=adm_control_states,name="Modified_ADM1", switch="DAE",metagenome_report=adm_metagenome_report)
        Sol_mod_adm1 = mod_adm1.solve_model(mod_adm1.initial_conditions[:, 0], np.linspace(0,30, 10000))

        if args.report == 'dash' or args.report==None:
            mod_adm1.dash_app(Sol_mod_adm1)
        elif args.report == 'csv':
            Address=Prompt.ask("\n[yellow]Where do you want to save the csv file? ")
            mod_adm1.csv_report(Sol_mod_adm1,Address)
        else:
            print('Please provide a valid report option')
    

        



    #############################
    #### Configuration Block ####
    #############################

    if args.ADToolbox_Module == 'Configs' and args.Configs_Subparser=='set-base-dir':
        set_base_dir(args)

    if args.ADToolbox_Module == 'Configs' and args.Configs_Subparser=='build-folder-structure':
        build_folder_structure()

    if args.ADToolbox_Module == 'Configs' and args.Configs_Subparser=='download-all-databases':
        download_all_databases()
    







###########################
#####Configs functions#####
###########################

def set_base_dir(args):
    """
    This function sets the base directory for the ADToolbox.
    """
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),"pkg_data","ADToolbox_Configs.json"), 'r') as f:
        Conf = json.load(f)
    Conf["Base_Dir"] = args.directory
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),"pkg_data","ADToolbox_Configs.json"), 'w') as f:
        json.dump(Conf, f, indent=4)
    if not os.path.exists(args.directory):
        os.mkdir(args.directory)
        rich.print("[yellow]Directory did not exist. Created directory: [green]{}".format(args.directory))

def build_folder_structure():
    """
    This function will build the folder structure undrestandable by the ADToolbox.
    """
    Required_Folders=["Database","Outputs","Reports","Metagenomics_Data","Genomes","bin"]
    for i in Required_Folders:
        if not os.path.exists(os.path.join(Main_Dir,i)):
            os.mkdir(os.path.join(Main_Dir,i))
            rich.print(f"\n[yellow] {i} directory did not exist. Created directory: [green]{os.path.join(Main_Dir,i)}\n")

def download_all_databases():
    """
    This function will download all the databases required by the ADToolbox.
    """
    MetaG=ADToolBox.Metagenomics(Configs.Metagenomics())
    rich.print(u"[yellow] Downloading Amplicon to Genome databases ...\n")
    MetaG.get_requirements_a2g()
    rich.print(u"[bold green]\u2713 Amplicon to Genome databases were downloaded successfuly!\n")
    rich.print(u"[yellow] Downloading Seed database ...\n")
    ADToolBox.Database().download_seed_databases(Configs.Seed_RXN_DB)
    rich.print(u"[bold green]\u2713 Seed database was downloaded successfuly!\n")
    rich.print(u"[yellow] Downloading Original-ADM1 database ...\n")
    ADToolBox.Database().download_adm1_parameters(Configs.Original_ADM1())
    rich.print(u"[bold green]\u2713 Original-ADM1 database was downloaded successfuly!\n")
    rich.print(u"[yellow] Downloading Modified-ADM database ...\n")
    ADToolBox.Database().download_modified_adm_parameters(Configs.Modified_ADM())
    rich.print(u"[bold green]\u2713 Modified-ADM database was downloaded successfuly!\n")
    rich.print(u"[yellow] Downloading protein database ...\n")
    ADToolBox.Database().download_protein_database(Configs.Database().protein_db)
    rich.print(u"[bold green]\u2713 Protein database was downloaded successfuly!\n")
    rich.print(u"[yellow] Downloading Reaction database ...\n")
    ADToolBox.Database().download_reaction_database(Configs.Database().csv_reaction_db)
    rich.print(u"[bold green]All the databases were downloaded successfuly!\n")
    rich.print(u"[yellow] Downloading the feed database ...\n")
    ADToolBox.Database().download_feed_database(Configs.Database().feed_db)
    rich.print(u"[bold green]\u2713 Feed database was downloaded successfuly!\n")

if __name__ == "__main__":
    main()