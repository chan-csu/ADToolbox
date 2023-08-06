from adtoolbox import __version__
from adtoolbox import core,configs,optimize,adm
import numpy as np
import json
from pytest import fixture
import pytest

def test_version():
    assert __version__ == '0.1.0'
    
###FIXTURES###

@fixture
def experimental_data():
    with open(configs.Database().initial_conditions) as file:
        initial_conditions=json.load(file)
    study=core.Experiment(
        "test_study",
        initial_conditions=initial_conditions,
        time=[0,1,2,3,4,5],
        variables=[10,12],
        data=[[1,2,3,4,5,6],[2,3,4,5,6,7]],
        reference="test_reference"
    )
    return study


    
###--optimize module--###
@pytest.mark.long
def test_build_optimizer_object(experimental_data):
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
    
    

    tuner=optimize.Tuner(
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

                        
                    train_data=[experimental_data],
                    tuneables={"Y_Me_h2":(0,1),
                               "Y_Me_h2":(0,1),
                               },
                    fitness_mode="equalized",
                    var_type="model_parameters"
             )
    hist=tuner.optimize(max_runs=1)
