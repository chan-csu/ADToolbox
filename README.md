# Toolbox Overview

Parsa Ghadermazi 
parsa96@colostate.edu

AD Toolbox is developed in Chan Lab at Colorado State University. The main goal of this toolbox is to provide the tools that are useful for modeling and optimization of anaerobic digestion process. 

## Installation

Although it is not mandatory, It is greatly advised to create a python environment for this project and install ADToolbox in that environment. Check for how to build and use python's virtual environments. There are severall methods that allow you to install or use ADToolbox:


1- Install by cloning this repository:

```
git clone https://github.com/chan-csu/ADToolbox.git
cd ADToolbox
pip install .

```
2- Install using pip

```
Not available yet

```

3- Use docker

```
Not available yet

```

4- Click bellow:

```
My binder place holder 

```
-------------------
# Using ADToolbox

This Toolbox is comprised of different modules:

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
