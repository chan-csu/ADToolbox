import configs
import subprocess
import pathlib
import os
from typing import Iterable
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

def run_batch(
    series:Iterable,
    number_of_batches:int,
    function:callable,
    ):
    pass

# def file_finder(base_dir:str,pattern:str,save:bool=True)->dict[str:Iterable[str]]:
#     """
#     This function finds files in subdirectories of a base directory, and generates a dictionary
#     with subdirectory names as keys and file addresses as values.
    
#     """
#     base_dir=pathlib.Path(base_dir)
#     files_dict={}
#     for sub_dir in base_dir.iterdir():
#         if sub_dir.is_dir():
#             files_dict[sub_dir.name]=list(sub_dir.rglob(pattern))
def extract_zipped_file(zipped_file:str,container:str="None")->str:
    """
    This function extracts a zipped file to a directory
    Args:
        zipped_file: the address of the zipped file as a string
        container: The container to run the script in. If None, the script is run locally    
    Returns:
        None
    """
    if container == "None":
        script=f"gzip -d {zipped_file}"
    
    elif container =="singularity":
        script=f"singularity exec -B {zipped_file}:{zipped_file} {configs.singularity_container} gzip -d {zipped_file}"
        
    elif container =="docker":
        script=f"docker run -v {zipped_file}:{zipped_file} {configs.docker_container} gzip -d {zipped_file}"
        
    return script