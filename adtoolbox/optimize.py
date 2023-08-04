import openbox
import adm
import numpy as np
from dataclasses import dataclass
import json

@dataclass
class Experiment:
    name:str
    initial_conditions: dict
    time: list[float]
    variables: list[int]
    data: list[list[float]]
    reference: str = None
    
    def __post_init__(self):
        self.data=np.array(self.data)
        self.validate()
    
    def validate(self):
        assert len(self.time)==self.data.shape[0], "Number of time points must match number of rows in data."
        assert len(self.variables)==self.data.shape[1] , "Number of variables must match number of columns in data."
        assert self.time[0]==0, "Time must start at 0."
    
    def to_dict(self):
        return {"name":self.name,
                "initial_conditions":self.initial_conditions,
                "time":self.time,
                "variables":self.variables,
                "data":self.data.tolist(),
                "reference":self.reference}
    
    def write_experiment(self, path:str):
        with open(path, "w") as f:
            json.dump(self.to_dict(), f)
    
    @classmethod
    def load_experiment(cls, path:str):
        with open(path, "r") as f:
            data=json.load(f)
        return cls(**data)
        
        
    
    
    
class Tuner:
    """
    This class is a wrapper around openbox to optimize a model.
    """
    def __init__(self,
                 base_model: adm.Model,
                 train_data: list[Experiment],
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
        self.kwargs.update(kwargs)
        optimizer = openbox.optimizer.PSO(self.fitness, **self.kwargs)
        optimizer.run()
        return optimizer.get_incumbent()
        
        