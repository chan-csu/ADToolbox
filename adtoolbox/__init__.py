import json
import os
import rich
from rich.prompt import Prompt
import sys

"""Project Setup for ADToolBox."""
__version__ = "0.1.0"
sys.path.append(os.path.join(os.path.dirname(__file__)))

with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),"ADToolbox_Configs.json"),"r") as f:
    Conf = json.load(f)
    Main_Dir=Conf["Base_Dir"]
if Main_Dir:
    pass
else:
    Main_Dir=Prompt.ask("[yellow]No Base Directory Found: \nWhere do you want to store your ADToolBox Data?")  

if not os.path.exists(Main_Dir):
    os.mkdir(Main_Dir)
    rich.print("\n[yellow]Directory did not exist. Created directory: [yellow]{}".format(Main_Dir))

with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),"ADToolbox_Configs.json"),"w") as f:
    Conf["Base_Dir"]=Main_Dir
    json.dump(Conf,f)
