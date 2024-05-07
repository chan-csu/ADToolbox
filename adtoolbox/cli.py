import argparse
from re import M
import configs
import metagenomics_report
import core
from __init__ import __version__,Main_Dir
from rich.console import Console
import rich
import adm
from rich.table import Table
from rich.prompt import Prompt
from rich import markdown
import json
import os
import numpy as np
import subprocess
import pandas as pd
import utils
class AParser(argparse.ArgumentParser):
    def _print_message(self, message, file=None):
        rich.print(message, file=file)

def main():
    console = Console()
    console.rule("[bold red]ADToolBox")
    
    db_class=core.Database(config=configs.Database())

    parser = AParser(prog="ADToolBox",
        description="ADToolBox, a toolbox for anaerobic digestion modeling",
        epilog="ADToolBox-Chan Lab at Colorado State University.")
    parser.add_argument("-v", "--version", action="version",help="Prints the version of ADToolBox") 
    parser.version = f"[green]ADToolBox v{__version__}"
    subparsers = parser.add_subparsers(dest="ADToolbox_Module", help='ADToolbox Modules:')

    ### Database Module -------------------------------------------------------------
    
    
    #### FEED DATABASE ####
    subparser_database = subparsers.add_parser('Database', help='This module provides a command line interface to build or download the databases that ADToolbox requires')
    db_subp=subparser_database.add_subparsers(dest="database_module", help='Database commands:')
    
    db_subp.add_parser("initialize-feed-db", help="Initialize the Feed DB") 
    add_feed=db_subp.add_parser("add-feed", help="Add a feed to the feed database")
    add_feed.add_argument("-n", "--name", action="store", help="Name of the feed to be added to the database")
    add_feed.add_argument("-c", "--carbohydrates", action="store", help="Carbohydrate content of the feed to be added to the database in percentage",required=True,type=float)
    add_feed.add_argument("-p", "--proteins", action="store", help="Protein content of the feed to be added to the database in percentage",required=True,type=float)
    add_feed.add_argument("-l", "--lipids", action="store", help="Lipid content of the feed to be added to the database in percentage",required=True,type=float)
    add_feed.add_argument("-t","--tss", action="store", help="Total suspended solid content of the feed to be added to the database in percentage",required=True,type=float)
    add_feed.add_argument("-s","--si", action="store", help="Soluble inert content of the feed to be added to the database in percentage",required=True,type=float)
    add_feed.add_argument("-x","--xi", action="store", help="particulate inert content of the feed to be added to the database in percentage",required=True,type=float)
    add_feed.add_argument("-r","--reference", action="store", help="Reference where the numbers come from",required=True)
    
    show_feed_db=db_subp.add_parser("show-feed-db", help="Shows the feed database")
    show_feed_db.add_argument("-f","--filter", action="store", help="Filters the feed database",required=False)
    
    
### METAGENOMICS STUDIES DATABASE ###
    db_subp.add_parser("initialize-metagenomics-studies-db", help="Initialize the Metagenomics Studies DB")
    add_metagenomics_study=db_subp.add_parser("add-metagenomics-study", help="Add a metagenomics study to the Kbase")
    add_metagenomics_study.add_argument("-n", "--name", help="Metagenomics Study Name to be added to the Kbase",required=True)
    add_metagenomics_study.add_argument("-t", "--type", help="Metagenomics Study Type to be added to the Kbase",required=True)
    add_metagenomics_study.add_argument("-m", "--microbiome", help="Microbiome where the metagenomics study belongs to",required=True)
    add_metagenomics_study.add_argument("-s","--sample_accession", help="SRA accession ID for the sample",required=True)    
    add_metagenomics_study.add_argument("-c","--comments", help="Comments on the study of interest",required=True)
    add_metagenomics_study.add_argument("-p","--study_accession", help="SRA accession ID for the project",required=True)

### PROTEIN DATABASE ###
    db_subp.add_parser("initialize-protein-db", help="Generates the protein database for ADToolbox")
    add_protein=db_subp.add_parser("add-protein", help="Add a protein to the protein database")
    add_protein.add_argument("-i", "--uniport-id", action="store", help="Uniport ID of the protein to be added to the database",required=True)
    add_protein.add_argument("-n", "--name", action="store", help="The name to be attached to the protein, usally is EC number",required=True)

    
###  DOWNLOAD DATABASES ###

    db_subp.add_parser("download-reaction-db", help="Downloads the reaction database in CSV format")
    db_subp.add_parser("download-seed-reaction-db", help="Downloads the seed reaction database in JSON format")
    db_subp.add_parser("build-protein-db", help="Generates the protein database for ADToolbox")
    db_subp.add_parser("download-protein-db", help="Downloads the protein database in fasta format; You can alternatively build it from reaction database.")
    db_subp.add_parser("download-amplicon-to-genome-dbs", help="downloads amplicon to genome databases")
    db_subp.add_parser("download-all-databases", help="downloads all databases that are required by ADToolbox at once")

###

    ### Metagenomics Module ###
    meta_config_defult=configs.Metagenomics()
    subparser_metagenomics= subparsers.add_parser('Metagenomics', help="This module provides the import metagenomics functionalities of ADToolbox in the command line")
    metag_subp=subparser_metagenomics.add_subparsers(dest='metag_subparser',help='Available Metagenomics Commands:')
    
    metag_subp_1=metag_subp.add_parser('download_from_sra', help='This module provides a command line interface to download metagenomics data from SRA')
    metag_subp_1.add_argument("-s","--sample_accession",action="store",help="SRA accession ID for the sample",required=True)
    metag_subp_1.add_argument("-o","--output-dir",action="store",help="Output directory to store the downloaded data",required=True)
    metag_subp_1.add_argument("-c","--container",action="store",help="Container to use for the download: None, docker, or singualrity",default="None")
    
    metag_subp_1=metag_subp.add_parser('download_genome', help='This module provides a command line interface to download genomes from NCBI')
    metag_subp_1.add_argument("-g","--genome_accession",action="store",help="NCBI accession ID for the genome",required=True)
    metag_subp_1.add_argument("-o","--output-dir",action="store",help="Output directory to store the downloaded data",required=True)
    metag_subp_1.add_argument("-c","--container",action="store",help="Container to use for the download: None, docker, or singualrity",default="None")


    metag_subp_1=metag_subp.add_parser('align-genome' , help='Align Genomes to the protein database of ADToolbox, or any other fasta with protein sequences')
    metag_subp_1.add_argument("-n", "--name", action="store", help="An appropriate name for the genome that is to be aligned",required=True)
    metag_subp_1.add_argument("-i", "--input-file", action="store", help="Input the address of the JSON file includeing information about the genomes to be aligned",required=True)
    metag_subp_1.add_argument("-o", "--output-dir", action="store", help="Output the directory to store the alignment results",default=meta_config_defult.genome_alignment_output,required=True)
    metag_subp_1.add_argument("-c", "--container", action="store", help="Container to use for the alignment: None, docker, or singualrity",default="None")
    metag_subp_1.add_argument("-d", "--protein-db-dir", action="store", help="Directory containing the protein database to be used for alignment",default=meta_config_defult.protein_db,required=False)
    
    metag_subp_1=metag_subp.add_parser('align-multiple-genomes' , help='Align multiple Genomes to the protein database of ADToolbox, or any other fasta with protein sequences')
    metag_subp_1.add_argument("-i", "--input-file", action="store", help="The address to a JSON file that holds the information about the genomes",required=True)
    metag_subp_1.add_argument("-o", "--output-dir", action="store", help="Output the directory to store the alignment results",required=True)
    metag_subp_1.add_argument("-c", "--container", action="store", help="Container to use for the alignment: None, docker, or singualrity",default="None")
    metag_subp_1.add_argument("-d", "--protein-db-dir", action="store", help="Directory containing the protein database to be used for alignment",default=meta_config_defult.protein_db,required=False)

    metag_subp_1=metag_subp.add_parser('find-representative-genomes' , help='Finds representative genomes from the repseqs fasta file')
    metag_subp_1.add_argument("-i", "--input-file", action="store", help="The address to the repseqs fasta file",required=True)
    metag_subp_1.add_argument("-o", "--output-dir", action="store", help="The directory of the output file",required=True)
    metag_subp_1.add_argument("-c", "--container", action="store", help="Container to use for the alignment: None, docker, or singualrity",default="None")
    metag_subp_1.add_argument("-s", "--similarity", action="store", help="Similarity cutoff for clustering",default=0.97,type=float)
    metag_subp_1.add_argument("-f", "--format", action="store", help="Format of the output file. Either json or csv",default="csv")

    
    
    # metag_subp_1=metag_subp.add_parser('Metagenomics_Report', help='This module provides a command line interface to the metagenomics report web interface')
    # metag_subp_1.add_argument("-j","--json-file",help="Address to the JSON file that contains the metagenomics report",required=True)
    # metag_subp_1.add_argument("-m","--mds",help="Uses MDS plot if true else uses TSNE",action="store_true")
    # metag_subp_1.add_argument("-t","--threed",help="whether to use 3D plot instead of 2d",action="store_true")
    # metag_subp_1.add_argument("-n","--not-normalize",help="Do not normalize the data",action="store_true")
    

    
    # metag_subp_3=metag_subp.add_parser('align-genomes', help='Align Genomes to the protein database of ADToolbox, or any other fasta with protein sequences')
    # metag_subp_3.add_argument("-i", "--input-file", action="store", help="Input the address of the JSON file includeing information about the genomes to be aligned")
    # metag_subp_3.add_argument("-d", "--protein-db-dir", action="store", help="Directory containing the protein database to be used for alignment",default=meta_config_defult.protein_db)
    # metag_subp_3.add_argument("-o", "--output-dir", action="store", help="Output the directory to store the alignment results",default=meta_config_defult.genome_alignment_output)
    # metag_subp_3.add_argument("-b", "--bit-score", action="store", help="Minimum Bit Score cutoff for alignment",default=meta_config_defult.bit_score)
    # metag_subp_3.add_argument("-e", "--e-value", action="store", help="Minimum e-vlaue score cutoff for alignment",default=meta_config_defult.e_value)
    
    # metag_subp_4=metag_subp.add_parser("make-json-from-genomes",help="Generates JSON file required by Align-Genomes for custom genomes.")
    # metag_subp_4.add_argument("-i", "--input-file", action="store", help="Input the address of the CSV file includeing information about the genomes to be aligned",required=True)
    # metag_subp_4.add_argument("-o", "--output-file", action="store", help="Output the directory to store the JSON file.",required=True)
    
    # metag_subp_5=metag_subp.add_parser('map-genomes-to-adm', help='maps JSON file with genome infromation to ADM reactions')
    # metag_subp_5.add_argument("-i","--input-file",action="store",help="Input the address of the JSON file includeing information about the alignment of the genomes to the protein database")
    # metag_subp_5.add_argument("-m","--model",action="store",help="Model determines which mapping system you'd like to use; Current options: 'Modified_adm_reactions'",default="Modified_adm_reactions")
    # metag_subp_5.add_argument("-o","--output-dir",action="store",help="address to store the JSON report to be loaded with a model")


    doc_parser=subparsers.add_parser('Documentations', help='Documentations for using AD Toolbox')
    doc_parser.add_argument("-s", "--show", action="store_true", help="Documentation for a specific module")    
    
    ### ADM module ###
    subparser_adm = subparsers.add_parser('ADM', help='This module provides a command line interface for running and visualizing the ADToolbox ADM and modified ADM models')
    adm_subp=subparser_adm.add_subparsers(dest='adm_Subparser',help='Available ADM Commands:')
    subparser_adm1 = adm_subp.add_parser('adm1',help='Original ADM1 Model')
    subparser_adm1.add_argument("--model-parameters", action="store", help="Model parameters for ADM 1")
    subparser_adm1.add_argument("--base-parameters", action="store", help="Provide json file with base parameters for original ADM1")
    subparser_adm1.add_argument("--initial-conditions", action="store", help="Provide json file with initial conditions for original ADM1")
    subparser_adm1.add_argument("--inlet-conditions", action="store", help="Provide json file with inlet conditions for original ADM1")
    subparser_adm1.add_argument("--reactions", action="store", help="Provide json file with reactions for original ADM1")
    subparser_adm1.add_argument("--species", action="store", help="Provide json file with species for original ADM1")
    subparser_adm1.add_argument("--metagenome-report", action="store", help="Provide json file with metagenome report for original ADM1")
    subparser_adm1.add_argument("--report", action="store", help="Describe how to report the results of original ADM1. Current options are: 'dash' and 'csv'")
    
    mod_adm_subp=adm_subp.add_parser('e-adm',help='eADM Model')
    mod_adm_subp.add_argument("--model-parameters", action="store", help="Model parameters for Modified ADM")
    mod_adm_subp.add_argument("--base-parameters", action="store", help="Provide json file with base parameters for modified ADM")
    mod_adm_subp.add_argument("--initial-conditions", action="store", help="Provide json file with initial conditions for modified ADM")
    mod_adm_subp.add_argument("--inlet-conditions", action="store", help="Provide json file with inlet conditions for modified ADM")
    mod_adm_subp.add_argument("--reactions", action="store", help="Provide json file with reactions for modified ADM")
    mod_adm_subp.add_argument("--species", action="store", help="Provide json file with species for modified ADM")
    mod_adm_subp.add_argument("--metagenome-report", action="store", help="Provide json file with metagenome report for modified ADM")
    mod_adm_subp.add_argument("--control-states", action="store", help="Provide a json file that contains the control states, and their values ")
    mod_adm_subp.add_argument("--report", action="store", help="Describe how to report the results of modified ADM. Current options are: 'dash' and 'csv'")

    # Mod_ADM_args=subparser_adm.add_subparsers(help='Modified ADM Args:')

    subparser_Configs = subparsers.add_parser('Configs', help='Configurations of ADToolBox')
    subparser_Configs.add_argument("-s", "--set-base-dir", action="store", help="Set the base directory for ADToolBox to work with")
    subparser_Configs.add_argument("-g", "--get-base-dir", action="store_true", help="Get the current base directory for ADToolBox")
    
    ### PARSE ARGS ###

    args=parser.parse_args()

    ### CLI config block ###
    if args.ADToolbox_Module == 'Configs':

        if args.set_base_dir:
            configs.set_base_dir(args.set_base_dir)
        elif args.get_base_dir:
            print(configs.get_base_dir())  
    ### Database ###

    if args.ADToolbox_Module == 'Database' and "database_module" in args:
        if args.database_module=="initialize-feed-db":
            db_class.initialize_feed_db()
        elif args.database_module=="add-feed":
            feed=core.Feed(name=args.name,
                           carbohydrates=args.carbohydrates,
                           proteins=args.proteins,
                           lipids=args.lipids,
                           tss=args.tss,
                           si=args.si,
                           xi=args.xi,
                           reference=args.reference)
            db_class.add_feed_to_feed_db(feed=feed)
        
        elif args.database_module=="show-feed-db":
            if args.filter:
                t=db_class.get_feed_from_feed_db(field_name="name",field_value=args.filter)
            else:
                t=db_class.get_feed_from_feed_db(field_name="name",query="")
            if t:
                feed_table = Table(title="Feed Database",expand=True,safe_box=True)
                for i in t[0].to_dict().keys():
                    feed_table.add_column(i, justify="center", style="cyan", no_wrap=True)
                for i in t:
                    feed_table.add_row(*map(str,list(i.to_dict().values())))
                console.print(feed_table)
    if args.ADToolbox_Module == 'Database' and "database_module" in args:
        if args.database_module=="initialize-metagenomics-studies-db":
            db_class.initialize_metagenomics_studies_db()
        elif args.database_module=="add-metagenomics-study":
            metagenomics_study=core.MetagenomicsStudy(name=args.name,
                                                      study_type=args.type,
                                                      microbiome=args.microbiome,
                                                      sample_accession=args.sample_accession,
                                                      comments=args.comments,
                                                      study_accession=args.study_accession)
            db_class.add_metagenomics_study_to_metagenomics_studies_db(metagenomics_study=metagenomics_study)
        


    if args.ADToolbox_Module == 'Database' and "database_module" in args:
        if args.database_module=="initialize-protein-db":
            db_class.initialize_protein_db()
        elif args.database_module=="add-protein":
            db_class.add_protein_to_protein_db(protein_id=args.uniport_id,header_tail=args.name)
        
        elif args.database_module=="download-reaction-db":
            db_class.download_reaction_database()
        
        elif args.database_module=="download-seed-reaction-db":
            db_class.download_seed_databases()
        
        elif args.database_module=="build-protein-db":
            ecs=core.Database.ec_from_csv(configs.Database().csv_reaction_db)
            db_class.protein_db_from_ec(ecs)
            rich.print("[green]Protein DB built successfully")
        
        elif args.database_module=="download-feed-db":
            db_class.download_feed_database()
        
        elif args.database_module=="download-protein-db":
            db_class.download_protein_database()
        
        elif args.database_module=="download-amplicon-to-genome-dbs":
            db_class.download_amplicon_to_genome_db()
        
        elif args.database_module=="download-all-databases":
            db_class.download_all_databases()







    #### Metagenomics Module #####
    if args.ADToolbox_Module == 'Metagenomics' and "metag_subparser" in args and args.metag_subparser=="download_from_SRA":
        prefetch_script,sample_metadata=core.Metagenomics(meta_config_defult).seqs_from_sra(accession=args.sample_accession,target_dir=args.output_dir,container=args.container)
        subprocess.run(f"{prefetch_script}",shell=True)
        # core.Metagenomics(meta_config_defult).amplicon2genome()
        ### UNDER CONSTRUCTION ###
    
    if args.ADToolbox_Module == 'Metagenomics' and "metag_subparser" in args and args.metag_subparser=="download_genome":
        subprocess.run(core.Metagenomics(meta_config_defult).download_genome(identifier=args.genome_accession,
                                                              output_dir=args.output_dir,
                                                              container=args.container),shell=True)
    

    if args.ADToolbox_Module == 'Metagenomics' and "metag_subparser" in args and args.metag_subparser=="align-genome":
        if args.protein_db_dir:
            meta_config_defult.protein_db=args.protein_db_dir
        mg=core.Metagenomics(meta_config_defult).align_genome_to_protein_db(
            address=args.input_file,
            outdir=args.output_dir,
            name=args.name,
            container=args.container
        )
        subprocess.run(mg,shell=True)
    
    if args.ADToolbox_Module == 'Metagenomics' and "metag_subparser" in args and args.metag_subparser=="align-multiple-genomes":
        if args.protein_db_dir:
            meta_config_defult.protein_db=args.protein_db_dir
        with open(args.input_file) as f:
            genomes=json.load(f)
        for genome in genomes:
            mg=core.Metagenomics(meta_config_defult).align_genome_to_protein_db(
                address=genomes[genome],
                outdir=args.output_dir,
                name=genome,
                container=args.container
            )
            subprocess.run(mg,shell=True)

    if args.ADToolbox_Module == 'Metagenomics' and "metag_subparser" in args and args.metag_subparser=="find-representative-genomes":
        meta_config_defult.vsearch_similarity=args.similarity
        print(args.output_dir)
        mg=core.Metagenomics(meta_config_defult).align_to_gtdb(
            query_dir=args.input_file,
            output_dir=args.output_dir,
            container=args.container,
        )
        subprocess.run(mg,shell=True)
        res=core.Metagenomics(meta_config_defult).get_genomes_from_gtdb_alignment(args.output_dir)
        if args.format=="csv":
            pd.DataFrame.from_dict(res,orient="index").to_csv(os.path.join(args.output_dir,"representative_genomes.csv"))
        elif args.format=="json":
            with open(os.path.join(args.output_dir,"representative_genomes.json"),"w") as f:
                json.dump(res,f)
        else:
            raise ValueError("Please provide a valid format for the output file")
        
        

        

    # if args.ADToolbox_Module == 'Metagenomics' and "metag_Subparser" in args and args.metag_Subparser=="align-genomes":
    #     if "input_file" in args and "output_dir" in args:
    #         Meta_Config=configs.Metagenomics(genomes_json_info=args.input_file,genome_alignment_output=args.output_dir,genome_alignment_output_json=os.path.join(args.output_dir,"Alignment_Info.json"))
    #         core.Metagenomics(Meta_Config).Align_Genomes()

    
    # if args.ADToolbox_Module == 'Metagenomics' and "metag_subparser" in args:
    #     print(args)
    #     ncomp=3 if args.threed else 2
    #     method="MDS" if args.mds else "TSNE"
    #     normalize=False if args.not_normalize else True
    #     metagenomics_report.main(args.json_file,method=method,normalize=normalize,n_components=ncomp)





    if args.ADToolbox_Module == 'Documentations' and args.show:
        with open(configs.Documentation().readme,'r') as f:
            T=f.read()
            console.print(markdown.Markdown(T))
    
    if args.ADToolbox_Module == 'Database' and "Database_Module" in args and args.Database_Module=="build-protein-db" :
            ecs=db_class.ec_from_csv(configs.Database().csv_reaction_db)
            rich.print(u"[bold green]\u2713 All EC numbers were extracted!\n")
            db_class._initialize_database()
            db_class.protein_db_from_ec(ecs)
            rich.print(u"\n[bold green]\u2713 Protein Database was built!\n")
    
    
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
        with open(configs.Database().feed_db, 'r') as f:
            feed_db = json.load(f)
        colnames=feed_db[0].keys()
        feed_table = Table(title="Feed Database",expand=True,safe_box=True)

        for i in colnames:
            feed_table.add_column(i, justify="center", style="cyan", no_wrap=True)
        for i in feed_db:
            feed_table.add_row(*map(str,list(i.values())))
        console.print(feed_table)
    

    
    if "ADToolbox_Module" in args and args.ADToolbox_Module == 'ADM' and "adm_Subparser" in args and args.adm_Subparser == 'adm1':
        
        if args.model_parameters:
            with open(args.model_parameters) as f:
                adm_model_parameters=json.load(f)
        else:
            with open(configs.ADM1_LOCAL["model_parameters"]) as f:
                adm_model_parameters=json.load(f)
        
        if args.base_parameters:
            with open(args.base_parameters) as f:
                adm_base_parameters=json.load(f)
        else:
            with open(configs.ADM1_LOCAL["base_parameters"]) as f:
                adm_base_parameters=json.load(f)

        if args.initial_conditions:
            with open(args.initial_conditions) as f:
                adm_initial_conditions=json.load(f)
        else:
            with open(configs.ADM1_LOCAL["initial_conditions"]) as f:
                adm_initial_conditions=json.load(f)

        if args.inlet_conditions:
            with open(args.inlet_conditions) as f:
                adm_inlet_conditions=json.load(f)
        else:
            with open(configs.ADM1_LOCAL["inlet_conditions"]) as f:
                adm_inlet_conditions=json.load(f)

        
        if args.reactions:
            with open(args.reactions) as f:
                adm_reactions=json.load(f)
        else:
            with open(configs.ADM1_LOCAL["reactions"]) as f:
                adm_reactions=json.load(f)
        
        if args.species:
            with open(args.species) as f:
                adm_species=json.load(f)
        else:
            with open(configs.ADM1_LOCAL["species"]) as f:
                adm_species=json.load(f)
        
        if args.metagenome_report:
            with open(args.metagenome_report) as f:
                adm_metagenome_report=json.load(f)
        else:
            pass
        
        ADM1 = adm.Model(
            model_parameters=adm_model_parameters,
            base_parameters=adm_base_parameters,
            initial_conditions= adm_initial_conditions,
            inlet_conditions= adm_inlet_conditions,
            feed=adm.DEFAULT_FEED,
            reactions=adm_reactions, 
            species=adm_species, 
            ode_system=adm.adm1_ode_sys, 
            build_stoichiometric_matrix=adm.build_adm1_stoiciometric_matrix,
            name="ADM1",
            switch="DAE")
        sol = ADM1.solve_model(t_eval=np.linspace(0,30, 10000))

        if args.report == 'dash' or args.report==None:
            ADM1.dash_app(sol)
        elif args.report == 'csv':
            Address=Prompt.ask("\n[yellow]Where do you want to save the csv file? ")
            ADM1.csv_report(sol,Address)
        else:
            print('Please provide a valid report option')
    

    elif "ADToolbox_Module" in args and args.ADToolbox_Module == 'ADM' and "adm_Subparser" in args and args.adm_Subparser == 'e-adm':
        params=utils.load_multiple_json_files(configs.E_ADM_2_LOCAL)
        control_state={"S_H_ion":10**(-6.5)}
        if args.control_states:
            with open(args.control_states) as f:
                control_state.update(json.load(f))
                
        mod_adm1 = adm.Model(model_parameters=params.model_parameters,
                             base_parameters=params.base_parameters,
                             initial_conditions=params.initial_conditions, 
                             inlet_conditions=params.inlet_conditions,
                             feed=adm.DEFAULT_FEED,
                             reactions=params.reactions,
                            species=params.species,
                            ode_system=adm.e_adm_2_ode_sys,
                            build_stoichiometric_matrix=adm.build_e_adm_2_stoichiometric_matrix,
                            control_state=control_state,
                            name="Modified_ADM1",
                            switch="DAE",
                            )
        sol_mod_adm1 = mod_adm1.solve_model(t_eval=np.linspace(0,30, 10000),method="BDF")

        if args.report == 'dash' or args.report==None:
            mod_adm1.dash_app(sol_mod_adm1)
        elif args.report == 'csv':
            Address=Prompt.ask("\n[yellow]Where do you want to save the csv file? ")
            mod_adm1.csv_report(sol_mod_adm1,Address)
        else:
            print('Please provide a valid report option')
    

if __name__ == "__main__":
    main()