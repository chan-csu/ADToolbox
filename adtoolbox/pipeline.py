# class Pipeline:
#     """A class for running a chain of commands in series with parallelism."""
#     def __init__(self,commands:list[tuple],executer:str,**kwargs):
#         self.commands=commands
#         self.executer=executer
#         self.kbase_config=kwargs.get("kbase_config",configs.Kbase())
#         self.database_config=kwargs.get("database_config",configs.Database())
#         self.reaction_toolkit_config=kwargs.get("reaction_toolkit_config",configs.Reaction_Toolkit())
#         self.metagenomics_config=kwargs.get("metagenomics_config",configs.Metagenomics())
#         self.alignment_config=kwargs.get("alignment_config",configs.Alignment())
#         self.original_adm_config=kwargs.get("original_adm_config",configs.Original_ADM1())
#         self.modified_adm_config=kwargs.get("modified_adm_config",configs.Modified_ADM())
        
    
#     def validate(self,verbose:bool=True):
#         """Validates the pipeline commands. Runs the following checks:
#         1. Input to the first command is correct."""
    
#         satisfied=["Configs.Metagenomics",
#                     "Configs.Database",
#                     "Configs.Reaction_Toolkit",
#                     "Configs.Kbase",
#                     "Configs.Alignment",
#                     ]
#         requirements=[]
#         failure={
#             com.__name__ :[] for com in self.commands
#          }
#         success=[]
#         for command in self.commands:
#             satisfied.extend(Pipeline.extract_reqs_sats(command)[1])
#             requirements=Pipeline.extract_reqs_sats(command)[0]
#             for req in requirements:
#                 if req not in satisfied:
#                     failure[command.__name__].append(req)
#             if len(failure[command.__name__])==0:
#                 success.append(command)
#         if verbose:
#             rich.print("[bold green]The following commands are valid:[/bold green]")
#             for com in success:
#                 failure.__delitem__(com.__name__)
#                 print(com.__name__)
#             rich.print("[bold red]The following commands are invalid:[/bold red]")
#             for com in failure.keys():
#                 print(f"{com} : {failure[com]}")


#         return success,failure

            
        






#     def run(self):
#         pass

#     @staticmethod
#     def extract_reqs_sats(command):
#         """Extracts the requirements and satisfactions from a command using the tags in the docstrings"""
        
#         if "Requires:" in command.__doc__ and "Satisfies:" in command.__doc__:
#             reqs=re.findall("<R>.+</R>",command.__doc__)
#             for ind,i in enumerate(reqs):
#                 reqs[ind]=i.replace("<R>","").replace("</R>","")
#             sats=re.findall("<S>.+</S>",command.__doc__)
#             for ind,i in enumerate(sats):
#                 sats[ind]=i.replace("<S>","").replace("</S>","")
#         else:
#             raise Exception(f"The command docstring for {command.__name__} is not properly formatted. Please use the following tags: <R> and <S> for requirements and satisfactions respectively.")

#         return reqs,sats