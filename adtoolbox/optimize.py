from openbox import Optimizer,ParallelOptimizer
import core,configs
from openbox import space as sp
from adtoolbox.core import Experiment
import pygad
import adm
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
             exp_std:float=5
             )->None:

        self.base_model = base_model
        if set(tuneables.keys())-set(base_model.__getattribute__(var_type).keys()):
            raise ValueError("Tuneable parameters not in model parameters.")
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
                ic.update({self.base_model.species[k]:experiment.data[0,idx] for idx,k in enumerate(experiment.variables)})
                _model= self.base_model.copy()
                _model.update_parameters(**{self.var_type:parameters})
                _model.update_parameters(initial_conditions=ic)
                solution=_model.solve_model(np.array(experiment.time),method=ode_method).y[experiment.variables,:]
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

        
    def optimize(self, perturbation_method:str="random",ode_method="LSODA", **kwargs)->dict:
        costs=[]
        if not self._initialized:
            for pop in self._popuplation:
                costs.append(self._cost(dict(zip(self.tunables.keys(),pop)),ode_method=ode_method))
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
            new_cost=self._cost(dict(zip(self.tunables.keys(),new_params)),ode_method=ode_method)
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

            
    

            
        

    

def validate_model(model:adm.Model,data:core.Experiment|Iterable[core.Experiment],plot:bool=False,show_extra_states:Iterable[str]|None=None)->tuple[dict[str,pd.DataFrame],plotly.graph_objs.Figure|None]:
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
    if isinstance(data,Experiment):
  
        ic=data.initial_concentrations.copy()
        ic.update({model.species[k]:data.data[0,idx] for idx,k in enumerate(data.variables) })
        model.update_parameters(initial_conditions=ic)
        solution=model.solve_model(np.array(data.time))
        out={"model":pd.DataFrame(solution.y[data.variables,:].T,index=np.array(data.time).tolist(),columns=data.variables),
             "data":pd.DataFrame(data.data,index=data.time,columns=data.variables)}
        if plot:
            fig=go.Figure()
            for idx,variable in enumerate(data.variables):
                fig.add_trace(go.Scatter(x=out["model"].index,y=out["model"][variable],name=model.species[variable],mode="lines",line=dict(
                    color=pallet[idx])))
                fig.add_trace(go.Scatter(x=out["data"].index,y=out["data"][variable],name=model.species[variable]+" observed",mode="markers",marker=dict(
                    color=pallet[idx])))
            if show_extra_states:
                for idx,extra in enumerate(show_extra_states):
                    fig.add_trace(go.Scatter(x=out["model"].index,y=solution.y[model.species.index(extra)],name=model.species[model.species.index(extra)],mode="lines",line=dict(
                        color=pallet[idx+len(data.variables)])))
            fig.show(renderer="svg")
    
    elif isinstance(data,Iterable):
        
        if plot:
            fig=go.Figure()
        ic=pd.concat([pd.DataFrame(i.initial_concentrations,index=[0]) for  i in data ]).mean().to_dict()
        ic.update({model.species[k]:data[0].data[0,idx] for idx,k in enumerate(data[0].variables) })
        model.update_parameters(initial_conditions=ic)
        solution=model.solve_model(np.array(data[0].time))
        out={"model":pd.DataFrame(solution.y[data[0].variables,:].T,index=np.array(data[0].time).tolist(),columns=data[0].variables)}
        for idx,variable in enumerate(data[0].variables):
            fig.add_trace(go.Scatter(x=out["model"].index,y=out["model"][variable],name=model.species[variable],mode="lines",line=dict(
                    color=pallet[idx])))
        df=[]
        for data_ in data:
            conc=pd.DataFrame(data_.data,columns=[model.species[i] for i in data_.variables],index=data_.time)
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
                y.append(df_grouped.get_group((time,comp)).mean()["value"])
                e.append(df_grouped.get_group((time,comp)).std()["value"]/np.sqrt(df_grouped.get_group((time,comp)).shape[0]))
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
                fig.add_trace(go.Scatter(x=out["model"].index,y=solution.y[model.species.index(extra)],name=model.species[model.species.index(extra)],mode="lines",line=dict(
                    color=pallet[idx+len(data[0].variables)])))
        fig.show(renderer="svg")
    else:
        raise TypeError("Data argument should be either a core. Experiment object or an iterable of core.Experiment objects.")
        
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
            x.append(model_[column])
            y.append(data_[column])
    x=np.array(x)
    y=np.array(y)
    return Validation(r_squared=1-(np.sum(np.square(y-x))/np.sum(np.square(y-np.mean(y)))),rmse=np.sqrt(np.sum(np.square(y-x))))
      
            
        




    
if __name__ == "__main__":
    pass
    # with open(configs.Database().initial_conditions) as file:
    #     initial_conditions=json.load(file)
    # study=core.Experiment(
    #     "test_study",
    #     time=[0,1,2,3,4,5],
    #     variables=[10,12],
    #     data=[[1,2,3,4,5,6],[2,3,4,5,6,7]],
    #     reference="test_reference"
    # )
    # with open(configs.Database().initial_conditions) as file:
    #     initial_conditions=json.load(file)
    
    # with open(configs.Database().base_parameters) as file:
    #     base_parameters=json.load(file)
    
    # with open(configs.Database().model_parameters) as file:
    #     model_parameters=json.load(file)
    
    # with open(configs.Database().inlet_conditions) as file:
    #     inlet_conditions=json.load(file)
    
    # with open(configs.Database().reactions) as file:
    #     reactions=json.load(file)
    
    # with open(configs.Database().species) as file:
    #     species=json.load(file)

    # tuner=Tuner(
    #     base_model=adm.Model(initial_conditions=initial_conditions,
    #                          base_parameters=base_parameters,
    #                          model_parameters=model_parameters,
    #                         inlet_conditions=inlet_conditions,
    #                         feed=adm.DEFAULT_FEED,
    #                         reactions=reactions,
    #                         species=species,
    #                         ode_system=adm.modified_adm_ode_sys,
    #                         build_stoichiometric_matrix=adm.build_modified_adm_stoichiometric_matrix,
    #                         name="test_model"),
                        
    #                 train_data=[study],
    #                 tuneables={"Y_Me_h2":(0,1),},
    #                 fitness_mode="equalized",
    #                 var_type="model_parameters"
    #          )
    # hist=tuner.optimize(max_runs=2)
   
