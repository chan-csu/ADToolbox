# Toolbox Overview

AD Toolbox is developed in Chan Lab at Colorado State University. The main goal of this toolbox is to provide all of the tools required for modeling and optimization of anaerobic digestion process. This pipeline is comprised of different modules:

1. Database Module

2. Metagenomics Module

3. Documentations Module

4. ADM Module

5. Report Module

6. Utility Module

7. Configs Module

---------------------------------------------------

## 1- Database Module

-------------


Database module handles any task that is related to the ADToolbox. In order to find out what functions this module provides, type the following in your terminal **after installation**:

```
python ADToolbox Database --help

```

this will print the following in the console:

```
optional arguments:
  -h, --help            show this help message and exit
  --Initialize-Feed-DB  Initialize the Feed DB
  --Extend-Feed-DB EXTEND_FEED_DB
                        Extend the Feed Database using a CSV file
  --Show-Feed-DB        Display the Current Feed Database
  --Show-Reaction-DB    Display the Current Compound Database
  --Build-Protein-DB    Generates the protein database for ADToolbox

```