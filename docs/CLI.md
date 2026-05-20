# ADToolbox Commandline Interface

-------------
## Initialization
After installing ADToolbox, you can access all commands and their brief explanations by:

```
adtoolbox --help

```
The CLI does not initialize or store a global base directory. Commands that need a file or directory ask for it directly or accept it as an option.

-------------
## ADToolbox Modules

This toolbox is comprised of different modules:

1. Database Module

2. Metagenomics Module

3. ADM Module

4. Documentations Module

-------------
### 1. Database Module

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

- initialize-feed-db: This will create an empty feed database at the database path you provide. You can run this command by:

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
- initialize-metagenomics-studies-db: This will create an empty TSV file at the studies database path you provide. You can run this command by:

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
- initialize-protein-db: This will create an empty protein database at the protein database path you provide. You can run this command by:

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

*NOTE*: Skip the following download commands if you have already downloaded the required databases.
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
### 2. Metagenomics Module

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
### 3. ADM Module

ADM module provides the CLI entry points for the active ADToolbox anaerobic digestion models: ADM1 and e-ADM. In order to find out about all the functionalities in this module, you can run:

```

adtoolbox ADM --help

Usage: adtoolbox ADM [OPTIONS] COMMAND [ARGS]...

  Run and visualize ADToolbox ADM models.

Commands:
  adm1   Original ADM1 model.
  e-adm  eADM model.

```

- `adm1`: run the original ADM1 model with parameters from a directory or from explicit JSON file paths.

```
adtoolbox ADM adm1 --parameters-dir /path/to/ADToolbox/adm1 --report csv
```

- `e-adm`: run the current e-ADM model.

```
adtoolbox ADM e-adm --parameters-dir /path/to/ADToolbox/e_adm --report csv
```

Both commands accept the same parameter file options:

```
--parameters-dir
--models-json
--model-parameters
--base-parameters
--initial-conditions
--inlet-conditions
--reactions
--species
--metagenome-report
--report
```

The preferred input is `--models-json`, a single JSON file containing all ADM models keyed by model name. The e-ADM command also accepts `--control-states`, which should point to a JSON object of states that should be held constant. If you choose `dash` for your report, the CLI opens the interactive Dash interface. If you choose `csv`, it writes concentration profiles over time.

-------------
### 4. Documentations Module

You can view the documentaion in your CLI using rich's markdown render. You can do this by:

```
ADToolbox Documentations --show 

```
