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
                          {initialize-feed-db,add-feed,sh
ow-feed-db,initialize-metagenomics-studies-db,add-metagen
omics-study,initialize-protein-db,add-protein,download-re
action-db,download-seed-reaction-db,build-protein-db,down
load-protein-db,download-amplicon-to-genome-dbs,download-
all-databases}
                          ...

positional arguments:
  {initialize-feed-db,add-feed,show-feed-db,initialize-me
tagenomics-studies-db,add-metagenomics-study,initialize-p
rotein-db,add-protein,download-reaction-db,download-seed-
reaction-db,build-protein-db,download-protein-db,download
-amplicon-to-genome-dbs,download-all-databases}
                        Database commands:
    initialize-feed-db  Initialize the Feed DB
    add-feed            Add a feed to the feed database
    show-feed-db        Shows the feed database
    initialize-metagenomics-studies-db
                        Initialize the Metagenomics Studies DB
    add-metagenomics-study
                        Add a metagenomics study to the Kbase
    initialize-protein-db
                        Generates the protein database for ADToolbox
    add-protein         Add a protein to the protein database           
    download-reaction-db
                        Downloads the reaction database in CSV
                        format
    download-seed-reaction-db
                        Downloads the seed reaction database in
                        JSON format
    build-protein-db    Generates the protein database for
                        ADToolbox
    download-protein-db
                        Downloads the protein database in fasta
                        format; You can alternatively build it
                        from reaction database.
    download-amplicon-to-genome-dbs
                        downloads amplicon to genome databases
    download-all-databases
                        downloads all databases that are required by ADToolbox at once
    
options:
  -h, --help            show this help message and exit
```

We will now go over these commands one by one:

- initialize-feed-db: This will create an empty JSON file in the Database sub-directory in your base directory that will hold all the future feed information that you add. You can run this command  by:

```
ADToolbox Database initialize-feed-db

```

- add-feed: This will add feed data to the database. Such data includes: the name of the feed (-n, --name), carbohydrate content of the feed in a percentage (-c, --carbohydrates), protein content of the feed in a percecntage (-p, --proteins), lipid content of the feed in a percentage (-l, --lipids), total suspended solid content of the feed in a percentage (-t, --tss), soluable inert content of feed in a percentage (-s, --si), particulate inert content of feed in a percentage (-x, --xi), and a reference where numbers came from (-r, --reference).  This command is run by:

```
ADToolbox Database add-feed

```
An example of this would look like:

```
ADToolbox Database add-feed -n "test feed" -c 20 -p 20 -l 20 -t 20 -s 20 -x 20 -r "test reference"

```

- show-feed-db: As the name implies, this will show the user the feed database along with any values they have added to it, in the command window. This command is run by:

```
ADToolbox Database show-feed-db

```
- initialize-metagenomics-studies-db: This will create an empty TSV file in the Database sub-directory in your base directory that will hold all the future information about various metagenomics studies that you add. You can run this command  by:

```
ADToolbox Database initialize-metagenomics-studies-db

```
- add-metagenomics-study: This command will add a metagenomics study to the Kbase and will require the study name (-n,--name), study type (-t, --type), microbiome where the metagenomics study belongs to (-m, --microbiome), SRA accession ID for the sample (-s, --sample_accesion), SRA accession ID for the project (-p, --study_accesion), and comments on the study of interest (-c, --comments). This command is run by:

```
ADToolbox Database add-metagenomics-study

```
An example of this would look like:

```
ADToolbox Database add-metagenomics-study  -n test_study -t 16s -m "anaerobic digestion"  -s 11111111 -c "this is just a test" -p 222222

```
- initialize-protein-db: This will create an empty JSON file in the Database sub-directory in your base directory that will hold all the future protein information that you add. You can run this command  by:

```
ADToolbox Database initialize-protein-db

```
- add-protein: As the name implies, this will add information about a protein to the empty protein database. Information about such protein includes its UniProt ID (-i, --uniprot-id), and the name attached to the protein which is usually the EC number (-n, --name). You can run this command by:

```
ADToolbox Database add-protein

```
An example of this would look like:

```
ADToolbox Database add-protein -i ATEST1 -n 1.1.1.1

```

*NOTE*: Skip the following download commands if you have run ```ADToolbox Configs download-all-databases```
```

- download-reaction-db: As the name implies, this will download the ADToolbox reaction database. This is required for many important modules of the toolbox

```
ADToolbox Database download-reaction-db


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
                              {download_from_sra}
                              ...

positional arguments:
  {download_from_sra,download_genome}
    download_from_sra   This module provides a command line interface to download
                        metagenomics data from SRA
    download_genome     This module provides a command line interface to download
                        genomes from NCBI      
    align-genome        Align genomes to the protein database
                        of ADToolbox, or any other fasta with
                        protein sequences
    align-multiple-genomes
                        Align multiple Genomes to the protein
                        database of ADToolbox, or any other
                        fasta with protein sequences                                           
    find-representative-genomes
                        Finds representative genomes from the
                        repseqs fasta file
options:
  -h, --help            show this help message and exit

```
- download_from_sra: This command takes a sample accesion ID (-s, --sample_accesion) for a sample, downloads it, and places it into a directory provided by the you (-o, --output-dir). It also requires you to state a container you are using. If you are downloading locally, put "None". Otherwise, you can use the containers docker or singularity. You can run this command by:

```
ADToolbox Metagenomics download_from_sra

```
An example of this command would look like:

```
ADToolbox Metagenomics download_from_sra -s SRR28403133 -o OUTPUT/DIRECTORY/PATHNAME -c None

```
- download_genome: This command requires you to provide a NCBI accesion ID for a genome (-g, --genome_accesion), and output directory (-o,--output-dir), and a container (-c, --container). It will then take the NCBI accesion ID for a genome and download it in the directory provided by you. If you are downloading locally, put "None" as your container option. Otherwise, you can use the containers docker or singularity. You can run this command by: 

```
ADToolbox Metagenomics download_genome

```
An example of this command would look like:

```
ADToolbox Metagenomics download_genome -g GCA021152825.1 -o OUTPUT/DIRECTORY/PATHNAME -c None

```
- align-genome: This command requires that you to give a name for the genome (-n,--name),the address of the JSON file that includes information about the genome to be aligned (-i,--input-file), and output directory to store alignment results (-o,--output-dir), a container to use for the alignment (-c,--container), and the directory containing the protein database to be used for the alignment (-d, --protein-db-dir).  If you are downloading locally, put "None" as your container option. Otherwise, you can use the containers docker or singularity. Overall, this command takes a genome and aligns it to a protein sequence. You can run this command by: 

```
ADToolBox Metagenomics align-genome

```
An example of this code would look like:

```
ADToolbox Metagenomics align-genome -n "test genome" -i INPUT/PATHNAME/OF/GENOME -o OUTPUT/PATHNAME/DIRECTORY -c None -d PATHNAME/OF/PROTEIN

```
- align-multiple-genomes: This command allows you to align multiple genomes to the protein database of ADToolbox, or any other fasta file with protein sequences. It requires to user to input the address to a JSON file that holds the information about the genomes (-i,--input-file), an output directory to store the alignment results (-o,--output-dir), a container to use for the alignment (-c,--container), and the directory containing the protein database to be used for alignment (-d,--protein-db-dir). If you are downloading locally, put "None" as your container option. Otherwise, you can use the containers docker or singularity. This command can be run by: 

```
ADToolbox Metagenomics align-multiple-genomes

```
An example of this command looks like:

```
ADToolbox Metagenomics align-multiple-genomes -i PATHNAME/TO/FILE/OF/GENOMES -o OUTPUT/DIRECTORY -c None -d DIRECTORY/OF/PROTEIN/DATABSE

```
- find-represenative-genomes: This command maps represenative amplicon sequences to a representative genome in GTDB database. It requires the user to provide the address to the repseqs fasta file (-i,--input-file), the directory of the output file (-o, --output-dir), a container used for the alignment (-c,--container), and the format of the output file which can be json or csv (-f,--format). Something optional that you can provide is the similarity cutoff for clustering; though, the default is 0.97 (-s,--similarity). If you are downloading locally, put "None" as your container option. Otherwise, you can use the containers docker or singularity. You can run this code by:

```
ADToolbox Metagenomics find-representative-genomes

```
An example of this code will look like: 

```
ADToolbox Metagenomics find-representative-genomes -i PATHNAME/TO/REPSEQS/FASTA/FILE -o PATHNAME/TO/OUTPUT/DIRECTORY -c None -f csv

```

-------------------------------
### 4. ADM Module

ADM module provides all the tools needed to run instances of ADM Model. This includes the original ADM, Batstone et al., and the Modified-ADM suggested by the Authors of ADToolbox. In order to find out about all the functionalities in this module, you can run:

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
