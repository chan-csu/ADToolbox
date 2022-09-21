# ADToolbox Commandline Interface

Here we go over using the commandline interface, CLI, of ADToolbox.
First, we need to initialize the CLI:

## Initialization

After installing ADToolbox, type and execute the following in your terminal to initialize the base directory for ADToolbox files:

```
ADToolbox

```
You should see the following if you are running this for the first time:

```
No Base Directory Found: 
Where do you want to store your ADToolBox Data?:

```
Type the **absolute path** directory of interest. Don't worry if you mess this part up. You can change this later as well. you can type '.' for now and change this later.

You can access all the commands along with their brief explanation by:

```
ADToolbox --help

```

## ADToolbox Modules

This toolbox is comprised of different modules:

1. Configs Module

2. Database Module

3. Metagenomics Module

4. Documentations Module

5. ADM Module

6. Report Module

7. Utility Module



-------------

### 1. Configs Module

After installation, you have to download all the required files to run ADToolbox properly. To do this, first go to the Configs module by:

```

ADToolbox Configs --help

```
This will give you a list of all functionalities that are related to the configurations of the toolbox. Here we go one by one in the correct order:

- set-base-dir:  The first configuration command will allow you to set the base directory for ADToolbox to work. This could be an existing folder somewhere in your files or a directory that you are willing to create. If the directory does not already exit, it will be automatically created after this command. For example if I want to set the base directory to be ADToolbox directory on my desktop the command would be, in MacOS, something like this:

```

ADToolbox Configs set-base-dir -d ~/Desktop/ADToolbox

```

Note that you ***must*** include -d after ```set-base-dir```.

Anything that you will do from now on, will be saved in this directory.

- build-folder-structure: Now you need to build the folder structure that is understandable by ADToolbox. You will do that by:

```

ADToolbox Configs build-folder-structure

```

- download-all-databases: The next step is to download all the necessary database files for ADToolbox to work properly. You will achieve this by: 

```
ADToolbox Configs download-all-databases

```

- download-escher-files: If you want to access the esher map functionalities of ADToolbox you need to download the required files by:

```

ADToolbox Configs download-escher-files

```

An overal view of the Configs module of ADToolbox can be obtained by:

```

ADToolbox Configs --help

────────────────────────────────── ADToolBox───────────────────────────────────
usage: ADToolBox Configs [-h]
                         {set-base-dir,build-folder-structure,download-all-datab
ases,download-escher-files}
                         ...

positional arguments:
  {set-base-dir,build-folder-structure,download-all-databases,download-escher-fi
les}
                        Available Configs Commands:
    set-base-dir        Determine the address of the base directory for
                        ADToolBox to work with
    build-folder-structure
                        Builds the folder structure for ADToolBox to work
                        properly
    download-all-databases
                        Downloads all the databases for ADToolBox to work
                        properly, and puts them in the right directory in
                        Databases
    download-escher-files
                        Downloads all files required for running the escher
                        map functionality of ADToolBox

options:
  -h, --help            show this help message and exit


```
### 2. Database Module

Any database that is used by ADToolbox can be modified from this module. Type the following in your commandline to find all of the database module's commands:

```
ADToolbox Database --help

────────────────────────────────── ADToolBox───────────────────────────────────
usage: ADToolBox Database [-h]
                          {initialize-feed-db,extend-feed-db,show-feed-db,show-r
eaction-db,build-protein-db,download-feed-db,download-reaction-db,download-prote
in-db,download-amplicon-to-genome-dbs}
                          ...

positional arguments:
  {initialize-feed-db,extend-feed-db,show-feed-db,show-reaction-db,build-protein
-db,download-feed-db,download-reaction-db,download-protein-db,download-amplicon-
to-genome-dbs}
                        Database Modules:
    initialize-feed-db  Initialize the Feed DB
    extend-feed-db      Extend the Feed Database using a CSV file
    show-feed-db        Display the Current Feed Database
    build-protein-db    Generates the protein database for ADToolbox
    download-feed-db    Downloads the feed database in JSON format
    download-reaction-db
                        Downloads the reaction database in CSV format
    download-protein-db
                        Downloads the protein database in fasta format; You
                        can alternatively build it from reaction database.
    download-amplicon-to-genome-dbs
                        downloads amplicon to genome databases

options:
  -h, --help            show this help message and exit


```

We now go over these commands one by one:

- initialize-feed-db: This will create an empty JSON file in the Database sub-directory in your base directory that will hold all the future feed information that you add. You can run this command  by:

```

ADToolbox Database initialize-feed-db

```

- extend-feed-db: you can add to your current feed database from a csv file with this command. The CSV file **must** follow this column configuration:

|Name|TS|TSS|Lipids|Proteins|Carbohydrates|PI|SI|Notes|
|----|--|---|------|--------|-------------|--|--|-----|


You can add your CSV file by:

```

ADToolbox Database extend-feed-db -d [PATH TO YOUR CSV FILE]

```

- show-feed-db: You can pretty print your current feed table by

```

ADToolbox Database show-feed-db

```

This command will print a pretty table of your current feed database, which is the JSON file you initialized earlier and extended with a CSV file. This is made possible because of the great Rich library.


*NOTE*: Skipp the following download commands if you ran ```ADToolbox Configs download-all-databases```
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