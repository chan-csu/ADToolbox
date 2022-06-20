# Toolbox Overview

Parsa Ghadermazi 
parsa96@colostate.edu

AD Toolbox is developed in Chan Lab at Colorado State University. The main goal of this toolbox is to provide the tools that are useful for modeling and optimization of anaerobic digestion process.

** [Full Documentation Here](https://chan-csu.github.io/ADToolbox/) **


## Installation

Although it is not mandatory, It is greatly advised to create a python environment for this project and install ADToolbox in that environment. Check for how to build and use python's virtual environments. There are severall methods that allow you to install or use ADToolbox:


1- Install by cloning this repository:

```
git clone https://github.com/chan-csu/ADToolbox.git
cd ADToolbox
pip install .

```

Or if you wish to install in editable mode (contribute to the developemnt), run the following:

```
git clone https://github.com/chan-csu/ADToolbox.git
cd ADToolbox
pip install -e .

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

## Installation

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

## Using the ADToolox Commandline Interface

This toolbox is comprised of different modules:

1. Database Module

2. Metagenomics Module

3. Documentations Module

4. ADM Module

5. Report Module

6. Utility Module

7. Configs Module
-------------

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

----------


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
