from adtoolbox import core,configs
from adtoolbox.core import Experiment
from adtoolbox import Main_Dir
from adtoolbox import adm
import pandas as pd
import numpy as np
from dataclasses import dataclass
import plotly
import json
from typing import Iterable
import plotly.graph_objects as go
import plotly.express as px
import torch
import pickle
import pathlib
from collections import namedtuple
import ray
import time
import multiprocessing as mp
import logging
import os
from functools import lru_cache
start_time=time.strftime("%Y-%m-%d-%H-%M-%S")
log_dir=os.path.join(Main_Dir,"logs",f"optimize_{start_time}")

if not os.path.exists(os.path.join(Main_Dir,"logs")):
    os.makedirs(os.path.join(Main_Dir,"logs"))

logging.basicConfig(
     filename=f"{log_dir}.log",
     level=logging.INFO,
     encoding="utf-8",
     filemode="a",
     format="{asctime} - {levelname} - {message}",
     style="{",
     datefmt="%Y-%m-%d %H:%M",
 )
NUM_CORES=mp.cpu_count()
Validation=namedtuple("Validation",("r_squared","rmse"))
        
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
    
    def _get_space(self):
        """
        This function creates the search space for openbox.
        :return: The search space.
        """
        from openbox import space as sp
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
        from openbox import Optimizer,ParallelOptimizer
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
        
    def _get_space(self):
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
        import pygad
        opt=pygad.GA(
                     fitness_func=self.fitness,
                     num_genes=len(self.tunables),
                     gene_space=self._get_space(),

                     **kwargs)
        opt.run()
        
        return opt


class NNSurrogateTuner:
    def __init__(self,
             base_model: adm.Model,
             train_data: list[core.Experiment],
             tuneables:dict,
             fitness_mode:str="equalized",
             var_type:str="model_parameters",
             n_process:int=4,
             grad_steps:int=10,
             n_steps:int=100,
             initial_points:int=5,
             save_every:int=10,
             history_file_path:pathlib.Path=pathlib.Path("./history.pkl"),
             exp_std:float=5,
             )->None:

        self.base_model = base_model
        self.tunables = tuneables
        self.train_data = train_data
        self.fitness_mode = fitness_mode
        self.var_type = var_type
        self.n_process=n_process
        self.grad_steps=grad_steps
        self.n_steps=n_steps
        self.initial_points=initial_points
        self._get_space()
        self._generate_initial_population()
        self._aquired={}
        self._network=torch.nn.Sequential(
            torch.nn.Linear(len(self.tunables), 30),
            torch.nn.Tanh(),
            torch.nn.Linear(30, 30),
            torch.nn.Tanh(),
            torch.nn.Linear(30, 30),
            torch.nn.Tanh(),
            torch.nn.Linear(30, 30),
            torch.nn.Tanh(),
            torch.nn.Linear(30, 30),
            torch.nn.Tanh(),
            torch.nn.Linear(30, 30),
            torch.nn.Tanh(),
            torch.nn.Linear(30, 1),   
        )
        self.save_every=save_every
        self.history_file_path=history_file_path
        self.exp_std=exp_std
        torch.set_num_threads(1)
        self._initialized=False
        
    def _get_space(self)->np.ndarray:
        self._param_space=np.array(list(zip(*self.tunables.values()))).T
        return self._param_space.copy()
    
    def _generate_initial_population(self):
        self._popuplation=np.array([np.random.uniform(low=self._param_space[:,0],high=self._param_space[:,1]) for i in range(self.initial_points)])

    def _cost(self, parameters: dict,ode_method:str)->float:
        """
        This function is called by openbox to evaluate a configuration.
        :param config: The configuration to evaluate.
        :return: The cost of the configuration.
        """
        if self.fitness_mode == "equalized":

            res=0
            for experiment in self.train_data:
                ic=experiment.initial_concentrations.copy()
                ic.update({k:experiment.data[0,idx] for idx,k in enumerate(experiment.variables)})
                _model= self.base_model.copy()

                for k,v in _model._ic.items():
                    if k in parameters:
                        ic[k]=parameters.get(k,v)
                
                for k,v in _model.model_parameters.items():
                    if k in parameters:
                        _model.model_parameters[k]=parameters.get(k,v)
                
                for k,v in _model.base_parameters.items():
                    if k in parameters:
                        _model.base_parameters[k]=parameters.get(k,v)
                
                _model.update_parameters(initial_conditions=ic)
                _model.control_state={k:experiment.initial_concentrations[k] for k in experiment.constants}
                _model.feed=experiment.feed

                solution=_model.solve_model(np.array(experiment.time),method=ode_method).y[[_model.species.index(i) for i in experiment.variables],:]
                res+=np.sum(np.square(solution.T-experiment.data))

            return res
        else:
            raise NotImplementedError("Fitness mode not implemented.")
    

    def _suggest_parameters(self):
        if len(self._aquired)==0:
            raise ValueError("No sample is aquired yet.")
        else:
            input_data = torch.tensor(self._best_tensor, requires_grad=True,dtype=torch.float32)
            self._optimizer_inputs = torch.optim.Adam([input_data], lr=0.001)
            self._optimizer_model = torch.optim.Adam(self._network.parameters(), lr=0.001)
            for step in range(self.grad_steps):
                output = self._network(input_data)
                loss = torch.mean(output)
                loss.backward()
                self._optimizer_inputs.step()
                self._optimizer_inputs.zero_grad()
                self._optimizer_inputs.zero_grad()
        return input_data.detach().numpy()
    def _train_surrogate(self):
        net_inputs=torch.tensor(self._aquired["parameters"],dtype=torch.float32)
        net_labels=torch.tensor(self._aquired["cost"],dtype=torch.float32)
        self._optimizer_model = torch.optim.Adam(self._network.parameters(), lr=0.001)
        net_labels[net_labels>10000]=1000
        for step in range(1000):
            output=self._network(net_inputs)
            loss=torch.mean(torch.square(output-net_labels))
            loss.backward()
            self._optimizer_model.step()
            self._optimizer_model.zero_grad()
        return loss.detach().numpy()

    def optimize(self, perturbation_method:str="random",ode_method="BDF",parallel:bool=False, **kwargs)->dict:
        if not parallel:
            costs=[]
            if not self._initialized:
                logging.info("Generating initial population:")
                for num,pop in enumerate(self._popuplation):
                    costs.append(self._cost(tuple(pop),ode_method=ode_method))
                    logging.info(f"Initial population: {num+1}/{len(self._popuplation)}")

                self._aquired={"parameters":self._popuplation,"cost":[c for c in costs]}
            self._best_tensor=self._aquired["parameters"][np.argmin(self._aquired["cost"])]
            self._best_cost=np.min(self._aquired["cost"])
            for i in range(self.n_steps):
                loss=self._train_surrogate()
                print(f"Training Loss: {loss}")
                new_params=self._suggest_parameters()
                new_params[new_params<self._param_space[:,0]]=self._param_space[:,0][new_params<self._param_space[:,0]]
                new_params[new_params>self._param_space[:,1]]=self._param_space[:,1][new_params>self._param_space[:,1]]
                if (len(self._aquired["cost"])>1) and (abs(self._aquired["cost"][-1]-self._aquired["cost"][-2])<1e-3):
                    print("Local optima reached: Perturbing the best solution and continuing.")

                    if perturbation_method=="random":
                        new_params=np.random.normal(self._best_tensor,np.abs(self._best_tensor/self.exp_std))


                    elif perturbation_method=="estimate_gradient_directions":
                        diff=self._aquired["parameters"]-self._best_tensor
                        cost_diff=(np.array(self._aquired["cost"])-self._best_cost).reshape(1,-1).repeat(diff.shape[1],axis=0).T
                        diff=-np.sign(np.mean(np.sign(np.multiply(diff,cost_diff)),axis=0))
                        new_params=self._best_tensor+np.random.normal(np.multiply(self._best_tensor,diff)/self.exp_std,np.abs(self._best_tensor/self.exp_std),size=diff.shape[0])


                    elif perturbation_method=="predict_pattern":
                        pass
                    
                    
                new_params[new_params<self._param_space[:,0]]=self._param_space[:,0][new_params<self._param_space[:,0]]
                new_params[new_params>self._param_space[:,1]]=self._param_space[:,1][new_params>self._param_space[:,1]]
                ## make sure not in the local optima
                new_cost=self._cost(tuple(new_params),ode_method=ode_method)
                self._aquired["parameters"]=np.vstack((self._aquired["parameters"],new_params))
                self._aquired["cost"].append(new_cost)
                self._best_tensor=self._aquired["parameters"][np.argmin(self._aquired["cost"])]
                self._best_cost=np.min(self._aquired["cost"])
                print(f"Step {i+1}/{self.n_steps} completed.Current cost:{new_cost} Best cost: {self._best_cost}")
                if i%self.save_every==0:
                    self.history={"parameters":self._aquired["parameters"],"cost":self._aquired["cost"],"tunable_parameters":list(self.tunables.keys())}
                    with open(f"{str(self.history_file_path.absolute())}","wb") as file:
                        pickle.dump(self.history,file)
                if i>51:
                    self._aquired["parameters"]=self._aquired["parameters"][-50:,:]
                    self._aquired["cost"]=self._aquired["cost"][-50:]   


            return self.history
        
        else:
            if kwargs.get("parallel_framework","ray")=="ray":
                ray.init()
                self.base_model_id = ray.put(self.base_model)
                self.train_data_ids = [ray.put(exp) for exp in self.train_data]
                self.tunables_id = ray.put(self.tunables)
                if not self._initialized:
                    logging.info("Generating initial population:")
                    costs=ray.get([_single_cost_ray.remote(self.base_model_id,
                                                           dict(zip(self.tunables.keys(),pop)),
                                                           exp_id,
                                                           ode_method=ode_method) for exp_id in self.train_data_ids for pop in self._popuplation])
                    costs=list(np.array(costs).reshape(-1,len(self.train_data)).sum(axis=1))
                    self._aquired={"parameters":self._popuplation,"cost":[c for c in costs]}
                self._best_tensor=self._aquired["parameters"][np.argmin(self._aquired["cost"])]
                self._best_cost=np.min(self._aquired["cost"])
                for i in range(self.n_steps):
                    loss=self._train_surrogate()
                    print(f"Training Loss: {loss}")
                    new_params=self._suggest_parameters()
                    new_params[new_params<self._param_space[:,0]]=self._param_space[:,0][new_params<self._param_space[:,0]]
                    new_params[new_params>self._param_space[:,1]]=self._param_space[:,1][new_params>self._param_space[:,1]]

                    if (len(self._aquired["cost"])>1) and (abs(self._aquired["cost"][-1]-self._aquired["cost"][-2])<1e-3):
                        print("Local optima reached: Perturbing the best solution and continuing.")
                        if perturbation_method=="random":
                            new_params=np.random.normal(self._best_tensor,np.abs(self._best_tensor/self.exp_std))
                        elif perturbation_method=="estimate_gradient_directions":
                            diff=self._aquired["parameters"]-self._best_tensor
                            cost_diff=(np.array(self._aquired["cost"])-self._best_cost).reshape(1,-1).repeat(diff.shape[1],axis=0).T
                            diff=-np.sign(np.mean(np.sign(np.multiply(diff,cost_diff)),axis=0))
                            new_params=self._best_tensor+np.random.normal(np.multiply(self._best_tensor,diff)/self.exp_std,np.abs(self._best_tensor/self.exp_std),size=diff.shape[0])


                        elif perturbation_method=="predict_pattern":
                            pass
                        
                        
                        
                    new_params[new_params<self._param_space[:,0]]=self._param_space[:,0][new_params<self._param_space[:,0]]
                    new_params[new_params>self._param_space[:,1]]=self._param_space[:,1][new_params>self._param_space[:,1]]
                    ## make sure not in the local optima
                    new_cost=np.sum(ray.get([_single_cost_ray.remote(self.base_model,dict(zip(self.tunables.keys(),new_params)),exp,ode_method=ode_method) for exp in self.train_data]))
                    self._aquired["parameters"]=np.vstack((self._aquired["parameters"],new_params))
                    self._aquired["cost"].append(new_cost)
                    self._best_tensor=self._aquired["parameters"][np.argmin(self._aquired["cost"])]
                    self._best_cost=np.min(self._aquired["cost"])
                    logging.info(f"Step {i+1}/{self.n_steps} completed.Current cost:{new_cost} Best cost: {self._best_cost}")
                    if i%self.save_every==0:
                        self.history={"parameters":self._aquired["parameters"],"cost":self._aquired["cost"],"tunable_parameters":list(self.tunables.keys())}
                        with open(f"{str(self.history_file_path.absolute())}","wb") as file:
                            pickle.dump(self.history,file)
                    if i>51:
                        self._aquired["parameters"]=self._aquired["parameters"][-50:,:]
                        self._aquired["cost"]=self._aquired["cost"][-50:]
                return self.history
            
            elif kwargs.get("parallel_framework","ray")=="native":
                if not self._initialized:
                    logging.info("Generating initial population:")
                    with mp.Pool(NUM_CORES) as pool:
                        costs=[pool.apply_async(_single_cost,args=(self.base_model,dict(zip(self.tunables.keys(),pop)),exp,self.var_type,ode_method)) for exp in self.train_data for pop in self._popuplation]
                        costs=[c.get() for c in costs]
                    logging.info("All initial population costs calculated.")
                    costs=list(np.array(costs).reshape(-1,len(self.train_data)).sum(axis=1))
                    self._aquired={"parameters":self._popuplation,"cost":[c for c in costs]}
                    
                self._best_tensor=self._aquired["parameters"][np.argmin(self._aquired["cost"])]
                self._best_cost=np.min(self._aquired["cost"])
                with mp.Pool(NUM_CORES) as pool:
                    for i in range(self.n_steps):
                        loss=self._train_surrogate()
                        logging.info(f"Training Loss: {loss}")
                        new_params=self._suggest_parameters()
                        new_params[new_params<self._param_space[:,0]]=self._param_space[:,0][new_params<self._param_space[:,0]]
                        new_params[new_params>self._param_space[:,1]]=self._param_space[:,1][new_params>self._param_space[:,1]]
                        if (len(self._aquired["cost"])>1) and (abs(self._aquired["cost"][-1]-self._aquired["cost"][-2])<1e-3):
                            logging.info("Local optima reached: Perturbing the best solution and continuing.")
                            if perturbation_method=="random":
                                new_params=np.random.normal(self._best_tensor,np.abs(self._best_tensor/self.exp_std))
                            elif perturbation_method=="estimate_gradient_directions":
                                diff=self._aquired["parameters"]-self._best_tensor
                                cost_diff=(np.array(self._aquired["cost"])-self._best_cost).reshape(1,-1).repeat(diff.shape[1],axis=0).T
                                diff=-np.sign(np.mean(np.sign(np.multiply(diff,cost_diff)),axis=0))
                                new_params=self._best_tensor+np.random.normal(np.multiply(self._best_tensor,diff)/self.exp_std,np.abs(self._best_tensor/self.exp_std),size=diff.shape[0])
                    
                    
                    
                        new_params[new_params<self._param_space[:,0]]=self._param_space[:,0][new_params<self._param_space[:,0]]
                        new_params[new_params>self._param_space[:,1]]=self._param_space[:,1][new_params>self._param_space[:,1]]
                        ## make sure not in the local optima
                        new_cost=[pool.apply_async(_single_cost,args=(self.base_model,dict(zip(self.tunables.keys(),new_params)),exp,self.var_type,ode_method)) for exp in self.train_data]
                        new_cost=np.sum([c.get() for c in new_cost])
                        self._aquired["parameters"]=np.vstack((self._aquired["parameters"],new_params))
                        self._aquired["cost"].append(new_cost)
                        self._best_tensor=self._aquired["parameters"][np.argmin(self._aquired["cost"])]
                        self._best_cost=np.min(self._aquired["cost"])
                        logging.info(f"Step {i+1}/{self.n_steps} completed.Current cost:{new_cost} Best cost: {self._best_cost}")
                        if i%self.save_every==0:
                            self.history={"parameters":self._aquired["parameters"],"cost":self._aquired["cost"],"tunable_parameters":list(self.tunables.keys())}
                            with open(f"{str(self.history_file_path.absolute())}","wb") as file:
                                pickle.dump(self.history,file)
                            logging.info(f"optimization results saved to {self.history_file_path.absolute()}")
                        if i>51:
                            self._aquired["parameters"]=self._aquired["parameters"][-50:,:]
                            self._aquired["cost"]=self._aquired["cost"][-50:]
                return self.history
            else:
                raise ValueError("Parallel framework not supported.")
                

    
    def load_history_file(self,address:str)->None:
        with open(address,"rb") as f:
            info=pickle.load(f)
            
        parameters=info["parameters"][np.argmin(info["cost"])]
        cost=np.min(info["cost"])
        parameters=dict(zip(info["tunable_parameters"],parameters))
        p_=np.array([parameters[i] for i in self.tunables.keys()]).reshape(1,-1)
        self._aquired={"cost":[cost],"parameters":p_,"tunable_parameters":list(self.tunables.keys())}
        self.base_model.update_parameters(**{self.var_type:parameters})
        # self._aquired=info
        self._initialized=True
    
    def _calculate_distributed_costs(
        self,
        parameters:np.ndarray,
        experiments:Iterable[core.Experiment],
        ode_method:str="BDF",
        parallel_framework:str="ray")->np.ndarray:
        
        if parallel_framework=="ray":
            costs=ray.get([_single_cost_ray.remote(self.base_model,dict(zip(self.tunables.keys(),p)),exp,ode_method=ode_method) for exp in experiments for p in parameters])
        return np.array(costs)
    
@ray.remote(num_cpus=NUM_CORES)
def _single_cost_ray(base_model:adm.Model,
                 parameters:dict,
                 experiment:core.Experiment,
                 ode_method:str="BDF")->float:
    """
    This function is called by ray to calculate the cost of a configuration.
    """
    base_model=base_model.copy()

    ic=experiment.initial_concentrations.copy()
    ic.update({k:experiment.data[0,idx] for idx,k in enumerate(experiment.variables)})
    
    for k,v in base_model._ic.items():
        if k in parameters:
            ic[k]=parameters[k]
    
    for k,v in base_model.model_parameters.items():
        if k in parameters:
            base_model.model_parameters[k]=parameters[k]
            
    
    base_model.update_parameters(initial_conditions=ic)
    
    for k,v in base_model.base_parameters.items():
        if k in parameters:
            base_model.base_parameters[k]=parameters[k]
            
    base_model.control_state={k:experiment.initial_concentrations[k] for k in experiment.constants}
    base_model.feed=experiment.feed
    solution=base_model.solve_model(np.array(experiment.time),method=ode_method).y[[base_model.species.index(i) for i in experiment.variables],:]
    return np.sum(np.square(solution.T-experiment.data))

    
def _single_cost(base_model:adm.Model,
                    parameters:dict,
                    experiment:core.Experiment,
                    var_type:str="model_parameters",
                    ode_method:str="LSODA")->float:
    
        ic=experiment.initial_concentrations.copy()
        ic.update({k:experiment.data[0,idx] for idx,k in enumerate(experiment.variables)})
        for k,v in base_model.initial_conditions.items():
            ic[k]=parameters.get(k,v)
        base_model.update_parameters(**{var_type:parameters})
        base_model.update_parameters(initial_conditions=ic)
        base_model.base_parameters=experiment.base_parameters
        base_model.control_state={k:experiment.initial_concentrations[k] for k in experiment.constants}
        base_model.feed=experiment.feed
        solution=base_model.solve_model(np.array(experiment.time),method=ode_method).y[[base_model.species.index(i) for i in experiment.variables],:]
        return np.sum(np.square(solution.T-experiment.data))
    

            
    

            
        

    

def validate_model(model:adm.Model,data:core.Experiment|Iterable[core.Experiment],plot:bool=False,show_extra_states:Iterable[str]|None=None,ode_solver:str="Radau")->tuple[dict[str,pd.DataFrame],plotly.graph_objs.Figure|None]:
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
    fig=None

    
    if isinstance(data,Iterable):
        

        fig=go.Figure()
        ic=pd.concat([pd.DataFrame(i.initial_concentrations,index=[0]) for  i in data ]).mean().to_dict()
        
        ic.update({k:data[0].data[0,idx] for idx,k in enumerate(data[0].variables) })
        model.control_state={k:data[0].initial_concentrations[k] for k in data[0].constants}
        model.update_parameters(base_parameters=data[0].base_parameters)
        model.update_parameters(initial_conditions=ic)
        all_time_points=np.array(sorted(list(set(sum([exp.time for exp in data],start=[])))))
        solution=model.solve_model(all_time_points,method=ode_solver)
        out={"model":pd.DataFrame(solution.y[[model.species.index(i) for i in data[0].variables],:].T,index=np.array(all_time_points).tolist(),columns=data[0].variables)}
        for idx,variable in enumerate(data[0].variables):
            fig.add_trace(go.Scatter(x=out["model"].index,y=out["model"][variable],name=variable,mode="lines",line=dict(
                    color=pallet[idx])))
        df=[]
        for data_ in data:
            conc=pd.DataFrame(data_.data,columns=data_.variables,index=data_.time)
            conc["time"]=conc.index
            conc=conc.melt(id_vars=["time"])
            df.append(conc)
        df=pd.concat(df)
        comps=df["variable"].unique()
        times=df["time"].unique()

        df_grouped=df.groupby(["time","variable"])
        for idx,comp in enumerate(comps):
            t=[]
            y=[]
            e=[]
            for time in times:
                t.append(time)
                y.append(df_grouped.get_group((time,comp)).mean(numeric_only=True)["value"])
                e.append(df_grouped.get_group((time,comp)).std(numeric_only=True)["value"]/np.sqrt(df_grouped.get_group((time,comp)).shape[0]))
            fig.add_trace(go.Scatter(x=t,
                                    y=y,
                                    error_y=dict(
                                    type='data', # value of error bar given in data coordinates
                                    array=e,
                                    color=pallet[idx],
                                    visible=True),
                                    mode="markers",
                                    marker=dict(
                                        color=pallet[idx]
                                    ),
                                    name=comp+" Observed"
                                    ))
            
        if show_extra_states:
            for idx,extra in enumerate(show_extra_states):
                fig.add_trace(go.Scatter(x=out["model"].index,y=solution.y[model.species.index(extra)],name=extra,mode="lines",line=dict(
                    color=pallet[idx+len(data[0].variables)])))
        if plot:
            fig.show(renderer="svg")


    else:
  
        ic=data.initial_concentrations.copy()
        ic.update({k:data.data[0,idx] for idx,k in enumerate(data.variables) })
        model.control_state={k:data.initial_concentrations[k] for k in data.constants}
        model.update_parameters(base_parameters=data.base_parameters)
        model.update_parameters(initial_conditions=ic)
        solution=model.solve_model(np.array(data.time),method=ode_solver)
        out={"model":pd.DataFrame(solution.y[[model.species.index(i) for i in data.variables],:].T,index=np.array(data.time).tolist(),columns=data.variables),
             "data":pd.DataFrame(data.data,index=data.time,columns=data.variables)}
        if plot:
            fig=go.Figure()
            for idx,variable in enumerate(data.variables):
                fig.add_trace(go.Scatter(x=out["model"].index,y=out["model"][variable],name=variable,mode="lines",line=dict(
                    color=pallet[idx])))
                fig.add_trace(go.Scatter(x=out["data"].index,y=out["data"][variable],name=variable+" observed",mode="markers",marker=dict(
                    color=pallet[idx])))
            if show_extra_states:
                for idx,extra in enumerate(show_extra_states):
                    fig.add_trace(go.Scatter(x=out["model"].index,y=solution.y[model.species.index(extra)],name=extra,mode="lines",line=dict(
                        color=pallet[idx+len(data.variables)])))
            fig.show(renderer="svg")
        
    return  out,fig

def calculate_fit_stats(model:adm.Model,data:Iterable[core.Experiment])->Validation:
    """This function calculates RMSE and R-squared metrics on a set of experiment objects
    #AIC
    """
    x=[]
    y=[]
    for study in data:
        formatted_data=validate_model(model,study)[0]
        model_,data_=formatted_data["model"],formatted_data["data"]
        for column in model_.columns:
            x.extend(model_[column])
            y.extend(data_[column])
    x=np.array(x)
    y=np.array(y)
    return Validation(r_squared=1-(np.sum(np.square(y-x))/np.sum(np.square(y-np.mean(y)))),rmse=np.sqrt(np.sum(np.square(y-x))))
      
            
        




    
if __name__ == "__main__":
    import utils
    params=utils.load_multiple_json_files(configs.E_ADM_2_LOCAL)
    model=adm.Model(
    initial_conditions=params.initial_conditions,
    inlet_conditions=params.inlet_conditions,
    model_parameters=params.model_parameters,
    reactions=params.reactions,
    species=params.species,
    feed=adm.DEFAULT_FEED,
    base_parameters=params.base_parameters,
    control_state={},
    build_stoichiometric_matrix=adm.build_e_adm_2_stoichiometric_matrix,
    ode_system=adm.e_adm_2_ode_sys,
    )
    db=core.Database(configs.Database())
    exp=db.get_experiment_from_experiments_db("name","")[:3]
    tuner=NNSurrogateTuner(
        base_model=model,
        train_data=exp[:3],
        tuneables={"k_m_bu":(1,100),"k_m_ac":(1,100)},
         history_file_path=pathlib.Path("./test_parallel.pkl"),
            n_steps=100,
    )
    tuner.optimize(parallel=True,parallel_framework="native",ode_method="BDF")

    
    
    
   
