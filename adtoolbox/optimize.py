from openbox import Optimizer
import core
from openbox import space as sp
import adm
import numpy as np
from dataclasses import dataclass
import json


        
class Tuner:
    """
    This class is a wrapper around openbox to optimize a model.
    """
    def __init__(self,
                 base_model: adm.Model,
                 train_data: list[core.Experiment],
                 tuneables:dict,
                 fitness_mode:str="equalized",
                 var_type:str="model_parameters",
                 **kwargs)->None:
        """
        base_model: The model to optimize.
        train_data: The data to train the model on.
        fitness_mode: The fitness mode to use.
        kwargs: Additional arguments to pass to openbox.
        """
        self.base_model = base_model
        self.tunables = tuneables
        self.train_data = train_data
        self.fitness_mode = fitness_mode
        self.var_type = var_type
        self.kwargs = kwargs
    
    def _get_space(self)->sp.Space:
        """
        This function creates the search space for openbox.
        :return: The search space.
        """
        space = sp.Space()
        space.add_variable([space.add_variable(sp.RealVariable(name, low, high)) for name, (low, high) in self.tunables.items()])
        return space
        
    def fitness(self, parameters: dict)->float:
        """
        This function is called by openbox to evaluate a configuration.
        :param config: The configuration to evaluate.
        :return: The fitness of the configuration.
        """
        if self.fitness_mode == "equalized":
            res=0
            for experiment in self.train_data:
                self._model= self.base_model.copy()
                self._model.update_parameters(**{self.var_type:parameters})
                self._model.update_parameters(initial_conditions=experiment.initial_conditions)
                solution=self._model.solve_model(np.array(experiment.time))[:,experiment.variables]
                res+=np.sum(np.square(solution-experiment.data))
            
            return res
        else:
            raise NotImplementedError("Fitness mode not implemented.")


    def optimize(self, **kwargs)->dict:
        """
        This function optimizes the model.
        kwargs: Additional arguments to pass to openbox.
        :return: The best configuration.
        """
        opt=Optimizer(
                self.fitness,
                self._get_space(),
                **self.kwargs,

                )
        history=opt.run()
        return history
        
        