# ADToolbox Commandline Interface

-------------
## Initialization
First, we need to initialize the CLI:
After installing ADToolbox, type and execute the following in your terminal to initialize the base directory for ADToolbox files:

```
ADToolbox

```
You should see the following if you are running this for the first time:

```
No Base Directory Found: 
Where do you want to store your ADToolbox Data?:

```
Type the **absolute path** directory of interest. Don't worry if you mess this part up. You can change this later as well. you can type '.' for now and change this later.

You can access all the commands along with their brief explanation by:

```
ADToolbox --help

```

-------------
## ADToolbox Modules

This toolbox is comprised of different modules:

1. Configs Module

2. Database Module

3. Metagenomics Module

4. ADM Module

5. Documentations Module

-------------
### 1. Configs Module

After installation, the base working directory must be specified:

```
ADToolbox Configs --help

────────────────────────────────── ADToolBox───────────────────────────────────

usage: ADToolBox Configs [-h] [-s SET_BASE_DIR] [-g] [--get-base-dir]

optional arguments:
  -h, --help            show this help message and exit
  -s SET_BASE_DIR, --set-base-dir SET_BASE_DIR
                        Set the base directory for ADToolBox to work with
  -g, --get-base-dir    Get the current base directory for ADToolBox


```
This will give you a list of all functionalities that are related to the configurations of the toolbox. Here we go one by one in the correct order:

- set-base-dir:  The first configuration command will allow you to set the base directory for ADToolbox to work. This could be an existing folder somewhere in your files or a directory that you are willing to create. If the directory does not already exit, it will be automatically created after this command. For example if I want to set the base directory to be ADToolbox directory on my desktop the command would be, in MacOS, something like this:

```

ADToolbox Configs --set-base-dir ~/Desktop/ADToolbox

```

Anything that you will do from now on, will be saved in this directory.


-------------
### 2. Database Module

Any database that is used by ADToolbox can be modified from this module. Type the following in your commandline to find all of the database module's commands:

```
ADToolbox Database --help
──────────────────────────── ADToolBox ────────────────────────────
usage: ADToolBox Database [-h]
                          {initialize-feed-db,extend-feed-db,show-f
eed-db,download-reaction-db,download-seed-reaction-db,build-protein
-db,download-feed-db,download-protein-db,download-amplicon-to-genom
e-dbs}
                          ...

positional arguments:
  {initialize-feed-db,extend-feed-db,show-feed-db,download-reaction
-db,download-seed-reaction-db,build-protein-db,download-feed-db,dow
nload-protein-db,download-amplicon-to-genome-dbs}
                        Database Modules:
    initialize-feed-db  Initialize the Feed DB
    add-feed            Add a feed to the feed database
    show-feed-db        Shows the feed table
    download-reaction-db
                        Downloads the reaction database in CSV
                        format
    download-seed-reaction-db
                        Downloads the seed reaction database in
                        JSON format
    build-protein-db    Generates the protein database for
                        ADToolbox
    download-feed-db    Downloads the feed database in JSON
                        format
    download-protein-db
                        Downloads the protein database in fasta
                        format; You can alternatively build it
                        from reaction database.
    download-amplicon-to-genome-dbs
                        downloads amplicon to genome databases
    download-all-databases
                        downloads all databases that are required by ADToolbox at once
    show-tables         Show the list of all studies in Kbase
    add-metagenomics-study
                        Add a metagenomics study to the Kbase
    add-experimental-data-study
                        Add a study with experimental data to the Kbase
options:
  -h, --help            show this help message and exit
```

We now go over these commands one by one:

- initialize-feed-db: This will create an empty JSON file in the Database sub-directory in your base directory that will hold all the future feed information that you add. You can run this command  by:

```
ADToolbox Database initialize-feed-db

```

- add-feed: This will add feed data to the database. Such data includes: NAME (the name of the feed), CARBOHYDRATES (carbohydrate content of the feed in a percentage), PROTEINS (protein content of the feed in a percecntage), LIPIDS (lipid content of the feed in a percentage), TSS (total suspended solid content of the feed in a percentage), SI (soluable inert content of feed in a percentage), XI (particulate inert content of feed in a percentage), and REFERENCE (reference where numbers came from).

```
ADToolbox Database add-feed

```

- show-feed-db: As the name implies, this will show the user the feed database along with any values they have added to it, in the command window. 

```
ADToolbox Database show-feed-db

```

*NOTE*: Skip the following download commands if you have run ```ADToolbox Configs download-all-databases```
```

- download-reaction-db: As the name implies, this will download the ADToolbox reaction database. This is required for many important modules of the toolbox

```
ADToolbox Database download-reaction-db

```

- download-feed-db: Downloads the default feed database in JSON format

```

ADToolbox Database download-feed-db

```

- download-protein-db: Downloads the protein database in fasta format; You can alternatively build it from reaction database if you have downloaded it; Check below.

```

ADToolbox Database download-protein-db

```

- build-protein-db: Generates the protein database for ADToolbox from the reaction database:

```

ADToolbox Database build-protein-db

```

- download-amplicon-to-genome-dbs: If you need to use the 16s mapping to the protein database and ADM, you will need to download the required databases using this command:


```

ADToolbox Database download-amplicon-to-genome-dbs

```

- download-seed-reaction-db: This will download the SEED reaction database in JSON format.

```
ADToolbox Database download-seed-reaction-db

```
```
-------------
### 3. Metagenomics Module

Metagenomics module of ADToolbox is designed to input metagenomics data into consideration when designing an AD process.

You can observe all the functionalities by:

```

ADToolbox Metagenomics --help   

──────────────────────────── ADToolBox ────────────────────────────
usage: ADToolBox Metagenomics [-h]
                              {amplicon-to-genome,align-genomes,mak
e-json-from-genomes,map-genomes-to-adm}
                              ...

positional arguments:
  {amplicon-to-genome,align-genomes,make-json-from-genomes,map-geno
mes-to-adm}
                         Available Metagenomics Commands:
    amplicon-to-genome  Downloads the representative genome from
                        each amplicon
    align-genomes       Align Genomes to the protein database of
                        ADToolbox, or any other fasta with
                        protein sequences
    make-json-from-genomes
                        Generates JSON file required by Align-
                        Genomes for custom genomes.
    map-genomes-to-adm  maps JSON file with genome infromation to
                        ADM reactions

options:
  -h, --help            show this help message and exit

```

As of right now, there are 3 main submodules exist in the Metagenomics module:

- amplicon-to-genome: If you have 16s Data from QIIME, you can, hopefully, find the representative genomes for each replicon in an automated way
using this functionality by selecting different parameters:


```
ADToolbox Metagenomics amplicon-to-genome --help

────────────────────────────────── ADToolBox ───────────────────────────────────
usage: ADToolBox Metagenomics amplicon-to-genome [-h] [-q QIIME_OUTPUTS_DIR]
                                                 [-f FEATURE_TABLE_DIR]
                                                 [-r REP_SEQ_DIR]
                                                 [-t TAXONOMY_TABLE_DIR]
                                                 [-o OUTPUT_DIR]
                                                 [-a AMPLICON_TO_GENOME_DB]
                                                 [--k K]
                                                 [--similarity SIMILARITY]

options:
  -h, --help            show this help message and exit
  -q QIIME_OUTPUTS_DIR, --qiime-outputs-dir QIIME_OUTPUTS_DIR
                        Input the directory to the QIIME outputs
  -f FEATURE_TABLE_DIR, --feature-table-dir FEATURE_TABLE_DIR
                        Input the directory to the feature table output from
                        QIIME output tables
  -r REP_SEQ_DIR, --rep-Seq-dir REP_SEQ_DIR
                        Input the directory to the repseq fasta output from
                        QIIME output files
  -t TAXONOMY_TABLE_DIR, --taxonomy-table-dir TAXONOMY_TABLE_DIR
                        Input the directory to the taxonomy table output from
                        QIIME output files
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output the directory to store the representative
                        genome files
  -a AMPLICON_TO_GENOME_DB, --amplicon-to-genome-db AMPLICON_TO_GENOME_DB
                        The Amplicon to Genome Database to use
  --k K                 Top k genomes to be selected
  --similarity SIMILARITY
                        Similarity cutoff in the 16s V4 region for genome
                        selection

```

If you do not provide any arguents, the default directories and values will be used. By default, ADToolbox looks at
your base directory, that you set in Configs, in :```Metagenomics_Data/QIIME_Outputs```. 

- If you want to use your default arguments, you need to provide the directory to Feature table, taxonomy table, repseq fasta file from QIIME.

- If you do not want to use the amplicon to genome database that you downloaded in Configs or Database modules, you can point to the database directory of your interest as well by ```--amplicon-to-genome-db``` or ```-a```

- From each sample you can choose the top k abundant taxa to be selected for downloading genome, for instance ```--k 10``` will select the top 10 taxa from each sample

- When selecting genomes you need a precentage of similarity cutoff to select a representaive genome from GTDB. You do this by ```--similarity 96```. This will set 96% as a similarity cutoff.

- Finally, you need to provide the directory where you want the representative genomes to be saved. Additionally, a few extra files 
providing information about the fetched genomes will be saved here:
  
  - ```GenomeAccessions.csv``` provides the NCBI accession IDs of the found genomes.
  - ```SelectedFeatures.csv``` provides the feature IDs, hashes, from QIIME outputs for the fetched genomes.
  - ```TopKTaxa.csv``` taxonomy name of the replicons for which a genome has been found.
  - ```Amplicon2Genome_OutInfo.json``` ***important*** This JSON file is a metadata about the genomes that were found. This is later used to align these genomes to the protein database.

A complete amplicon-to-genome command will look like:

```
ADToolbox Metagenomics amplicon-to-genome -f ~/Desktop/test/feature-table.tsv \
-r ~/Desktop/test/dna-sequences.fasta \
-t ~/Desktop/test/taxonomy.tsv \
-o ~/Desktop/test/ --k 10 \
-q ~/Desktop/test/ \
-a ~/Desktop/ADToolbox/Database/Amplicon2GenomeDBs/

```


- align-genomes This submodule uses MMseqs to align a list of genomes to the protein database of ADToolbox for functional analysis.

You can run this submodule by the following arguments:

```

ADToolbox Metagenomics align-genomes --help    
────────────────────────────────── ADToolBox ───────────────────────────────────
usage: ADToolBox Metagenomics align-genomes [-h] [-i INPUT_FILE]
                                            [-d PROTEIN_DB_DIR]
                                            [-o OUTPUT_DIR] [-b BIT_SCORE]
                                            [-e E_VALUE]

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file INPUT_FILE
                        Input the address of the JSON file includeing
                        information about the genomes to be aligned
  -d PROTEIN_DB_DIR, --protein-db-dir PROTEIN_DB_DIR
                        Directory containing the protein database to be used
                        for alignment
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output the directory to store the alignment results
  -b BIT_SCORE, --bit-score BIT_SCORE
                        Minimum Bit Score cutoff for alignment
  -e E_VALUE, --e-value E_VALUE
                        Minimum e-vlaue score cutoff for alignment

```

- input-file: you need to give a JSON file similar to what you create in amplicon-to-genome, see previous submodule, so that ADToolbox finds
all the information it needs about the genomes

- protein-db-dir: you can provide ADToolboxes protein database fasta file, or any protein database that you would like to align your genomes with.
- output-dir: Describes the location to save the alignment results. 
- bits-core: Minimum bit score to filter out the alignment results.
- e_value: Minimum bit score to filter out the alignment results.

The outputs of this step includes one more import file:

```Alignment_Info.json```  -> This file is used by map-genomes-to-adm

A complete command for this submodule would looklike:


```

ADToolBox Metagenomics align-genomes \
--input-file ~/Desktop/ADToolbox/Genomes/Amplicon2Genome_OutInfo.json
--protein-db-dir ~/Desktop/ADToolbox/Database/Protein_DB.fasta
--output-dir ~/Desktop/ADToolbox/Outputs/ \
--bit-score 40 \
--e-value 0.000001

```

----------------------

- make-json-from-genomes Sometimes you have genomes either from assembly or downloading it manually. In this case you can import your genomes 
to the pipeline by making a JSON file similar to  ```Amplicon2Genome_OutInfo.json```. To this you need to make a CSV file for your genomes:

```

ADToolbox Metagenomics make-json-from-genomes  --help

────────────────────────────────── ADToolBox ───────────────────────────────────
usage: ADToolBox Metagenomics make-json-from-genomes [-h] -i INPUT_FILE -o
                                                     OUTPUT_FILE

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file INPUT_FILE
                        Input the address of the CSV file includeing
                        information about the genomes to be aligned
  -o OUTPUT_FILE, --output-file OUTPUT_FILE
                        Output the directory to store the JSON file.

```

- input-file: This input file should be the address to a CSV file exactly in the following column format:

|Genome_ID|NCBI_Name|Genome_Dir   |
|---------|---------|-------------|
|Genome1  |xyz      |~/Desktop/...|

1- Genome_ID: Identifier for the genome, preferrably; NCBI ID
2- NCBI_Name: NCBI taxonomy name for the genome: Does not need to be in a specific format
3- Genome_Dir: Absolute path to the fasta files: NOT .gz


- output-file: output directory for the generated JSON to be saved. 

An example of this command would be:

```

ADToolBox Metagenomics make-json-from-genomes --input-file ~/Desktop/MyGenomes.CSV --output-files ~/Desktop/Genome_info.json

```


--------------------

- map-genomes-to-adm: This command will take a JSON file that has Alignment info for all of the genomes, and will output the mapping to ADM models. This is how you run this command:

```

ADToolbox Metagenomics map-genomes-to-adm --help
────────────────────────────────── ADToolBox ───────────────────────────────────
usage: ADToolBox Metagenomics map-genomes-to-adm [-h] [-i INPUT_FILE]
                                                 [-m MODEL] [-o OUTPUT_DIR]

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file INPUT_FILE
                        Input the address of the JSON file includeing
                        information about the alignment of the genomes to the
                        protein database
  -m MODEL, --model MODEL
                        Model determines which mapping system you'd like to
                        use; Current options: 'Modified_ADM_Reactions'
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        address to store the JSON report to be loaded with a
                        model



```
- input-file: This should be the alignment info JSON file created in align-genome step.
- model: this must be a string defining the model that you want to map the aligned reactions to. Current option is only
"Modified_ADM_Reactions"
- output-dir: Directory to save the Genome alignment report, Metagenome report, in JSON format. This can be used later in ADM module to include the role of different genomes in the process.

An example of a complete command for this submodule is:

```
ADToolbox Metagenomics map-genomes-to-adm -i ~/Desktop/alignment_info.json -m Modified_ADM_Reactions -o ~/Desktop/ADM_Mapping_Report.json

```



-------------------------------
### 4. ADM Module

ADM module provides all the tools needed to run instances of ADM Model. This include the originsl ADM, Batstone et al., and the Modified-ADM suggested by the Authors of ADToolbox. In order to find out about all the functionalities in this module, you can run:

```

ADToolbox ADM --help
────────────────────────────────── ADToolBox ───────────────────────────────────
usage: ADToolBox ADM [-h] {original-adm1,modified-adm,show-escher-map} ...

positional arguments:
  {original-adm1,modified-adm,show-escher-map}
                        Available ADM Commands:
    original-adm1       Original ADM1:
    modified-adm        Modified ADM:

options:
  -h, --help            show this help message and exit

```

- original-adm1: If you want to run the original ADM, batstone et al, in your browser you can run this command with the required parameters in JSON format:

```
ADToolbox ADM original-adm1 --help
────────────────────────────────── ADToolBox ───────────────────────────────────
usage: ADToolBox ADM original-adm1 [-h] [--model-parameters MODEL_PARAMETERS]
                                   [--base-parameters BASE_PARAMETERS]
                                   [--initial-conditions INITIAL_CONDITIONS]
                                   [--inlet-conditions INLET_CONDITIONS]
                                   [--reactions REACTIONS] [--species SPECIES]
                                   [--metagenome-report METAGENOME_REPORT]
                                   [--report REPORT]

options:
  -h, --help            show this help message and exit
  --model-parameters MODEL_PARAMETERS
                        Model parameters for ADM 1
  --base-parameters BASE_PARAMETERS
                        Provide json file with base parameters for original
                        ADM1
  --initial-conditions INITIAL_CONDITIONS
                        Provide json file with initial conditions for original
                        ADM1
  --inlet-conditions INLET_CONDITIONS
                        Provide json file with inlet conditions for original
                        ADM1
  --reactions REACTIONS
                        Provide json file with reactions for original ADM1
  --species SPECIES     Provide json file with species for original ADM1
  --metagenome-report METAGENOME_REPORT
                        Provide json file with metagenome report for original
                        ADM1
  --report REPORT       Describe how to report the results of original ADM1.
                        Current options are: 'dash' and 'csv'


```

Every argument is optional, and their role is clear from the comments in front of them, So we just provide a full example of this command:


```
ADToolbox ADM original-adm1 \
--model-parameters ~/Desktop/Model_Parameters.json \
--base-parameters ~/Desktop/Base_Parameters.json \
--initial-conditions ~/Desktop/Initial_Conditions.json \
--inlet-conditions ~/Desktop/Inlet-Conditions.json \
--reactions ~/Desktop/Reactions.json
--species  ~/Desktop/Species.json
--metagenome-report ~/Desktop/ADM_Mapping_Report.json
--repor dash

```
**NOTE** if you choose dash for your report, the CLI will prompt you to open your browser in the instructed address, if you choose csv, it will generate a CSV file that includes concentration profiles simulated over time.

---------------

- modified-adm: This command is exactly similar to the previous one, except that it requires parameters taylored for modified ADM:

```

ADToolbox ADM modified-adm --help 
────────────────────────────────── ADToolBox ───────────────────────────────────
usage: ADToolBox ADM modified-adm [-h] [--model-parameters MODEL_PARAMETERS]
                                  [--base-parameters BASE_PARAMETERS]
                                  [--initial-conditions INITIAL_CONDITIONS]
                                  [--inlet-conditions INLET_CONDITIONS]
                                  [--reactions REACTIONS] [--species SPECIES]
                                  [--metagenome-report METAGENOME_REPORT]
                                  [--report REPORT]

options:
  -h, --help            show this help message and exit
  --model-parameters MODEL_PARAMETERS
                        Model parameters for Modified ADM
  --base-parameters BASE_PARAMETERS
                        Provide json file with base parameters for modified
                        ADM
  --initial-conditions INITIAL_CONDITIONS
                        Provide json file with initial conditions for modified
                        ADM
  --inlet-conditions INLET_CONDITIONS
                        Provide json file with inlet conditions for modified
                        ADM
  --reactions REACTIONS
                        Provide json file with reactions for modified ADM
  --species SPECIES     Provide json file with species for modified ADM
  --metagenome-report METAGENOME_REPORT
                        Provide json file with metagenome report for modified
                        ADM
  --report REPORT       Describe how to report the results of modified ADM.
                        Current options are: 'dash' and 'csv'

```

The usage is exactly the same as the original-adm

-------------------------

- show-escher-map: This command will prompt you to open an escher map for the modified-adm model in your browser with the instructed address:

```
ADToolbox ADM show-escher-map

```

-------------
### 5. Documentations Module

You can view the documentaion in your CLI using rich's markdown render. You can do this by:

```
ADToolbox Documentations --show 

```
