import adm
import configs
import torch
from typing import Iterable, List, Tuple, Dict, Any
import torch.nn as nn


class AD_LSTM:
    def __init__(
        initial_conditions:dict,
        time_res:float,
        species:Iterable,
        reactions:Iterable,

    ):
        pass

class AD_Dense:
    
    def __init__(
        initial_conditions:dict,
        time_res:float,
        duration:float,
        species:Iterable,
        reactions:Iterable,
        ):
        pass
    

    def step(self)->torch.FloatTensor:
        next_state=self.predict_next()
        self.state=next_state.copy().detach()
        return next_state   

    def simulate(self)->torch.FloatTensor:
        self.reset()
        simulation=[self.state]
        t=0
        while t<self.duration:
            simulation.append(self.step())
            t+=self.time_res

        return torch.FloatTensor(simulation)

