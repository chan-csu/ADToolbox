# Parameter tuning

This is a how-to tutorial on how to train ADM parametrs using experimental data.

First, you need to import the packages that you need.

```

from adtoolbox import ADM,Main_Dir
import pygad
import numpy as np
import os
import json 
import pandas as pd

```
Now you are ready to start:

## 1- Load the default parameters

This script is build for tuning the parameters of Modified-ADM using a genetic algorithm. 
First you need to load the current default parameters.

```

Parameters_Dir=os.path.join(Main_Dir,"Database","Modified_ADM")
Base_Parameters=os.path.join(Parameters_Dir,"Modified_ADM_Base_Parameters.json")
IC=os.path.join(Parameters_Dir,"Modified_ADM_Initial_Conditions.json")
Inlet_C=os.path.join(Parameters_Dir,"Modified_ADM_Inlet_Conditions.json")
Model_Parameters=os.path.join(Parameters_Dir,"Modified_ADM_Model_Parameters.json")
Reactions=os.path.join(Parameters_Dir,"Modified_ADM_Reactions.json")
Species=os.path.join(Parameters_Dir,"Modified_ADM_Species.json")

with open(Base_Parameters) as f:
    Base_Parameters=json.load(f)
with open(IC) as f:
    IC=json.load(f)
with open(Inlet_C) as f:
    Inlet_C=json.load(f)
with open(Model_Parameters) as f:
    Model_Parameters=json.load(f)
with open(Reactions) as f:
    Reactions=json.load(f)
with open(Species) as f:
    Species=json.load(f)



```
## 2- Load the experimental data

In the next step you need to load your experimental data that you need to train parameters.
The data must be in a CSV file following this format(All in ADM Units):

|time|ADM_State_var1|ADM_State_var2|....|ADM_State_varM|           
|----|--------------|--------------|----|--------------|
|0.  | 0.1          | 0.15         | ...| 0.0|

```

Exp_Data=pd.read_csv(os.path.join("/Users/parsaghadermarzi/Desktop/Academics/Projects/Anaerobic_Digestion_Modeling/ADToolBox/Optimization","Experimental_Data.csv"),delimiter=",")

Inds=[]
for i in Exp_Data.columns[1:]:
    Inds.append(Species.index(i))
for state in Exp_Data.columns[1:]:
    IC[state]=Exp_Data[state][0]


```
## 3- Determine the parameters to be trained

Now you have to determine what parameters you want to tune using the experimental data that you have:

```

Params_list=["n_cap","f_ch_TSS","f_et_bu"]

gene_space=[[0,5],[0,1],[0,1]]

```

## 4- Start Optimizing

```

model=ADM.Model(Model_Parameters,Base_Parameters,IC,Inlet_C,Reactions,Species,ODE_System=ADM.Modified_ADM1_ODE_Sys,Build_Stoiciometric_Matrix=ADM.Build_Modified_ADM1_Stoiciometric_Matrix)

def fitness_func(solution, solution_idx):
    for i in range(len(Params_list)):
        Model_Parameters[Params_list[i]]=solution[i]
    model.Model_Parameters=Model_Parameters
    Sol = model.Solve_Model(
        (0, Exp_Data["time"].iloc[-1]+1), model.Initial_Conditions[:, 0], Exp_Data["time"].to_list())

    # Calculate the fitness
    fitness = 1.0 /(np.square(Sol.y[Inds,:] - Exp_Data.to_numpy()[:,1:].T)).mean()
    return fitness

fitness_function = fitness_func

num_generations = 50
num_parents_mating = 4

sol_per_pop = 8
num_genes = len(Params_list)

init_range_low = 0.0000000000000001
init_range_high =5

parent_selection_type = "sss"
keep_parents = 1

crossover_type = "single_point"

mutation_type = "random"
mutation_percent_genes = 10

ga_instance = pygad.GA(num_generations=num_generations,
                       num_parents_mating=num_parents_mating,
                       fitness_func=fitness_function,
                       sol_per_pop=sol_per_pop,
                       num_genes=num_genes,
                       init_range_low=init_range_low,
                       init_range_high=init_range_high,
                       parent_selection_type=parent_selection_type,
                       keep_parents=keep_parents,
                       crossover_type=crossover_type,
                       mutation_type=mutation_type,
                       mutation_percent_genes=mutation_percent_genes
                       ,gene_space=gene_space)

ga_instance.run()

```

## 5- Display the Results

You can plot the results:

```

solution, solution_fitness, solution_idx = ga_instance.best_solution()
print("Parameters of the best solution : {solution}".format(solution=solution))
print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))
ga_instance.plot_fitness()

```