import configs
import subprocess
import pathlib
import os

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

                
                


