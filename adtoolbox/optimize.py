from openbox import Optimizer
import core,configs
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
        if set(tuneables.keys())-set(base_model.__getattribute__(var_type).keys()):
            raise ValueError("Tuneable parameters not in model parameters.")
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
        space.add_variables([sp.Real(name, low, high,default_value=(high+low)/2) for name, (low, high) in self.tunables.items()])
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
                ic={self.base_model.species[k]:experiment.data[0,idx] for idx,k in enumerate(experiment.variables) }
                ic.update(experiment.initial_concentrations)
                self._model= self.base_model.copy()
                self._model.update_parameters(**{self.var_type:parameters})
                self._model.update_parameters(initial_conditions=ic)
                solution=self._model.solve_model(np.array(experiment.time)).y[experiment.variables,:]
                res+=np.sum(np.square(solution.T-experiment.data))
            
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
                **kwargs,

                )
        history=opt.run()
        return history
        

if __name__ == "__main__":
    with open(configs.Database().initial_conditions) as file:
        initial_conditions=json.load(file)
    study=core.Experiment(
        "test_study",
        time=[0,1,2,3,4,5],
        variables=[10,12],
        data=[[1,2,3,4,5,6],[2,3,4,5,6,7]],
        reference="test_reference"
    )
    with open(configs.Database().initial_conditions) as file:
        initial_conditions=json.load(file)
    
    with open(configs.Database().base_parameters) as file:
        base_parameters=json.load(file)
    
    with open(configs.Database().model_parameters) as file:
        model_parameters=json.load(file)
    
    with open(configs.Database().inlet_conditions) as file:
        inlet_conditions=json.load(file)
    
    with open(configs.Database().reactions) as file:
        reactions=json.load(file)
    
    with open(configs.Database().species) as file:
        species=json.load(file)

    tuner=Tuner(
        base_model=adm.Model(initial_conditions=initial_conditions,
                             base_parameters=base_parameters,
                             model_parameters=model_parameters,
                            inlet_conditions=inlet_conditions,
                            feed=adm.DEFAULT_FEED,
                            reactions=reactions,
                            species=species,
                            ode_system=adm.modified_adm_ode_sys,
                            build_stoichiometric_matrix=adm.build_modified_adm_stoichiometric_matrix,
                            name="test_model"),
                        
                    train_data=[study],
                    tuneables={"Y_Me_h2":(0,1),},
                    fitness_mode="equalized",
                    var_type="model_parameters"
             )
    hist=tuner.optimize(max_runs=2)
   
