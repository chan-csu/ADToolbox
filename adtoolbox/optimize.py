from openbox import Optimizer,ParallelOptimizer
import core,configs
from openbox import space as sp
import pygad
import adm
import pandas as pd
import numpy as np
from dataclasses import dataclass
import json
import plotly.graph_objects as go
import plotly.express as px


        
class BBTuner:
    """
    This class is a wrapper around openbox to optimize a model.
    """
    def __init__(self,
                 base_model: adm.Model,
                 train_data: list[core.Experiment],
                 tuneables:dict,
                 fitness_mode:str="equalized",
                 var_type:str="model_parameters",
                 parallel:bool=False,
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
        self.parallel=parallel
    
    def _get_space(self)->sp.Space:
        """
        This function creates the search space for openbox.
        :return: The search space.
        """
        space = sp.Space()
        space.add_variables([sp.Real(name, low, high,default_value=(high+low)/2) for name, (low, high) in self.tunables.items()])
        return space
        
    def cost(self, parameters: dict)->float:
        """
        This function is called by openbox to evaluate a configuration.
        :param config: The configuration to evaluate.
        :return: The cost of the configuration.
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
        if self.parallel:
            opt=ParallelOptimizer(
                self.cost,
                self._get_space(),
                **kwargs,
                )
        else:
            opt=Optimizer(
                    self.cost,
                    self._get_space(),
                    **kwargs,
                    )
        self.history=opt.run()
        self.optimized_model=self.base_model.copy()
        self.optimized_model.update_parameters(**{self.var_type:self.history.get_incumbents()})
        return self.history
        
class GATuner:
    """
    This class is a wrapper around pyGAD to optimize a model.

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
        kwargs: Additional arguments to pass to pyGAD.
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
        Returns the search space for pyGAD.
        """
        return [{"low":low,"high":high} for name, (low, high) in self.tunables.items()]
        
    def fitness(self, ga_instance,solution,solution_idx)->float:
        """
        This function is called by pyGAD to evaluate a configuration.
        :param config: The configuration to evaluate.
        :return: The cost of the configuration.
        """
        solution_=dict(zip(self.tunables.keys(),solution))
        if self.fitness_mode == "equalized":
            res=0
            for experiment in self.train_data:
                ic={self.base_model.species[k]:experiment.data[0,idx] for idx,k in enumerate(experiment.variables) }
                ic.update(experiment.initial_concentrations)
                self._model= self.base_model.copy()
                
                self._model.update_parameters(**{self.var_type:solution_})
                self._model.update_parameters(initial_conditions=ic)
                solution=self._model.solve_model(np.array(experiment.time)).y[experiment.variables,:]
                res+=np.sum(np.square(solution.T-experiment.data))
            return 1/res
        else:
            raise NotImplementedError("Fitness mode not implemented.")


    def optimize(self, **kwargs)->dict:
        """
        This function optimizes the model.
        kwargs: Additional arguments to pass to openbox.
        :return: The best configuration.
        """
        opt=pygad.GA(
                     fitness_func=self.fitness,
                     num_genes=len(self.tunables),
                     gene_space=self._get_space(),

                     **kwargs)
        opt.run()
        
        return opt

def validate_model(model:adm.Model,data:core.Experiment,plot:bool=False)->dict[str,pd.DataFrame]:
    """
    This function can be used to compare the model's predictions to the experimental data of interest.
    
    Args:
        model: The model to validate.
        data: The experiment to validate the model on in the form of a core.Experiment object.
        plot: Whether to plot the results.
    
    Returns:
        dict[str,pd.DataFrame]: A dictionary containing the model's predictions and the experimental data.
    
    """
    pallet=px.colors.qualitative.Plotly
    ic={model.species[k]:data.data[0,idx] for idx,k in enumerate(data.variables) }
    ic.update(data.initial_concentrations)
    model.update_parameters(initial_conditions=ic)
    solution=model.solve_model(np.array(data.time)).y[data.variables,:]
    out={"model":pd.DataFrame(solution.T,index=data.time,columns=data.variables),
            "data":pd.DataFrame(data.data,index=data.time,columns=data.variables)}
    if plot:
        fig=go.Figure()
        for idx,variable in enumerate(data.variables):
            fig.add_trace(go.Scatter(x=out["model"].index,y=out["model"][variable],name=model.species[variable],mode="lines",line=dict(
                color=pallet[idx])))
            fig.add_trace(go.Scatter(x=out["data"].index,y=out["data"][variable],name=model.species[variable]+" observed",mode="markers",marker=dict(
                color=pallet[idx])))
        fig.show()
        
    return 





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
   
