import configs
import subprocess
import pathlib
import os
import core
from typing import Iterable, Any,Union
def wrap_for_slurm(command:str,run:bool,save:bool,config:configs.Utils())->str:
    """
    This is a function that wraps a bash script in a slurm script.
    All resource allocation configuration is obtained from config argument
    Args:
        command: bash command to run in python string format
        run: if True, the slurm script is executed
        save: if True, the slurm script is saved to the path specified in config
        config: Configs.Utils object that detemines the template form and slurm options
    
    Returns:
        The slurm script as a string
    """
    with open(config.slurm_template,'r') as f:
        slurm_template = f.read()
    
    command = command.replace("#!/bin/bash","") 
    slurm_script = slurm_template.replace("<command>",command)
    slurm_script = slurm_script.replace("<sample_name>",config.slurm_job_name)
    slurm_script = slurm_script.replace("<wall_time>",config.slurm_wall_time)
    slurm_script = slurm_script.replace("<executer>",config.slurm_executer)
    slurm_script = slurm_script.replace("<job_name>",config.slurm_job_name)
    slurm_script = slurm_script.replace("<sample_outlog>",config.slurm_outlog)
    slurm_script = slurm_script.replace("<memory>",config.slurm_memory)
    slurm_script = slurm_script.replace("<cpus>",config.slurm_cpus)
    
    if save:
        with open(config.slurm_save_dir,'w') as f:
            f.write(slurm_script)

    if run:
        workingdir=pathlib.Path(config.slurm_save_dir).parents[1]
        slurmfile=pathlib.Path(config.slurm_save_dir)
        subprocess.run([f"sbatch {str(slurmfile.relative_to(workingdir))}"],cwd=workingdir,shell=True)        
    return slurm_script

def fasta_to_dict(fasta:str)->dict:
    """
    This function converts a fasta file to a dictionary
    Args:
        fasta: the affress fasta file as a string
    
    Returns:
        A dictionary with the fasta labels as keys and the fasta sequences as values
    """
    Dictionary={}
    
    with open(fasta,'r') as f:
        for line in f:
            if line.startswith('>'):
                label = line.lstrip(">").strip()
                Dictionary[label]=""
            else:
                Dictionary[label] += line.strip()
    return Dictionary

def dict_to_fasta(dictionary:dict,fasta:str)->None:
    """
    This function converts a dictionary to a fasta file
    Args:
        dictionary: the dictionary to convert
        fasta: the address of the fasta file to save
    
    Returns:
        None
    """
    with open(fasta,'w') as f:
        for label,seq in dictionary.items():
            f.write(f">{label}\n{seq}\n")



def extract_zipped_file(input_file:str,container:str="None",configs=configs.Utils())->str:
    """
    This function extracts a zipped file to a directory
    Args:
        container: The container to run the script in. If None, the script is run locally    
    Returns:
        None
    """
    if container == "None":
        script=f"gzip -d {input_file}"
    
    elif container =="singularity":
        script=f"singularity exec -B {input_file}:{input_file} {configs.adtoolbox_singularity} gzip -d {input_file}"
        
    elif container =="docker":
        script=f"docker run -v {input_file}:{input_file} {configs.adtoolbox_docker} gzip -d {input_file}"
        
    return script,

def generate_batch_script(
    generator_function:callable,
    number_of_batches:int,
    input_series:list[list],
    input_var:list[str],
    container:str="None",
    save:Union[str,None]=None,
    run:bool=False,
    header:str="#!/bin/bash\n",
    **kwargs)->tuple:
    """
    This is a general function that generates an iterable of bash scripts for running a function
    that creates a bash script on an iterable of inputs.
    
    """
    batch_size = len(list(zip(*input_series)))//number_of_batches+1
    batches=[]
    for i in range(len(input_series)):
        batches.append([input_series[i][j:j+batch_size] for j in range(0,len(input_series[i]),batch_size)])
    scripts=[]
    
    for inputs in zip(*batches):
        script=header
        for i,item in enumerate(zip(*inputs)):
            for ind,inp in enumerate(input_var):
                kwargs.update(dict([(inp,item[ind])]))
            script_=generator_function(container=container,**kwargs)[0]
            script=script+script_+"\n"
        if save:
            if not os.path.exists(save):
                os.makedirs(save)
            with open(os.path.join(save,f"{item[0]}.sh"),'w') as f:
                f.write(script)
        scripts.append(script)
    if run:
        for script in scripts:
            subprocess.run([script],shell=True)
            
    return tuple(scripts)
    

    
if __name__ == "__main__":
    # accessions=["AAA"+str(i) for i in range(100)]
    # names=["BBB"+str(i) for i in range(100)]
    # answers=generate_batch_script(
    #     generator_function=core.Metagenomics(configs.Metagenomics()).align_genome_to_protein_db,
    #     number_of_batches=5,
    #     input_series=[accessions,names],
    #     input_var=["address",'name'],
    #     container="singularity",
    #     save=os.path.join("/Users/parsaghadermarzi","Desktop","test_bash_generator",))
    pass