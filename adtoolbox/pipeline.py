import configs
import inspect
import rich
from typing import Callable

class Pipeline:
    """A class for running a chain of commands in series with parallelism."""
    def __init__(self,commands:list[Callable],arguments:dict,validate:bool=True):
        self.commands=commands
        self.arguments=arguments
        self._required_args=sum(list(inspect.signature(command).parameters.keys()) for command in self.commands)
        if validate:
            self.validate()    
    
    def validate(self):
        pass
    
    
    
        

            
        






    def run(self):
        pass

    @staticmethod
    def extract_reqs_sats(command):
        """Extracts the requirements and satisfactions from a command using the tags in the docstrings"""
        
        if "Requires:" in command.__doc__ and "Satisfies:" in command.__doc__:
            reqs=re.findall("<R>.+</R>",command.__doc__)
            for ind,i in enumerate(reqs):
                reqs[ind]=i.replace("<R>","").replace("</R>","")
            sats=re.findall("<S>.+</S>",command.__doc__)
            for ind,i in enumerate(sats):
                sats[ind]=i.replace("<S>","").replace("</S>","")
        else:
            raise Exception(f"The command docstring for {command.__name__} is not properly formatted. Please use the following tags: <R> and <S> for requirements and satisfactions respectively.")

        return reqs,sats