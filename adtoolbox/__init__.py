import json
import os
import rich
from rich.prompt import Prompt
import sys

"""Project Setup for ADToolBox."""

__version__ = "1.0.0"
__all__=["adm","configs","__main__","cli","core","optimize","pipeline","utils","Main_Dir","PKG_DATA"]

sys.path.append(os.path.join(os.path.dirname(__file__)))

PKG_DATA=os.path.join(os.path.dirname(os.path.realpath(__file__)),"pkg_data")

if not os.path.exists(os.path.join(PKG_DATA,"ADToolbox_Configs.json")):
    with open(os.path.join(PKG_DATA,"ADToolbox_Configs.json"),"w") as f:
        json.dump({"Base_Dir":""},f)

with open(os.path.join(PKG_DATA,"ADToolbox_Configs.json"),"r") as f:
    conf = json.load(f)
    Main_Dir=conf["Base_Dir"]
if Main_Dir and os.path.exists(Main_Dir):
    pass
elif Main_Dir and not os.path.exists(Main_Dir):
    Main_Dir=input(f"Base directory is not configured properly.\nPlease input the correct path for the base directory:")
    if not os.path.exists(Main_Dir):
        os.makedirs(Main_Dir)
    rich.print(f"[green]Base directory is set to {Main_Dir}")
else:
    Main_Dir=Prompt.ask("No Base Directory Found: \nWhere do you want to store your ADToolBox Data?")  

if not os.path.exists(Main_Dir):
    os.mkdir(Main_Dir)
    rich.print(f"\nDirectory did not exist. Created directory: {Main_Dir}")

with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),"pkg_data","ADToolbox_Configs.json"),"w") as f:
    conf["Base_Dir"]=Main_Dir
    json.dump(conf,f)
