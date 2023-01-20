import Configs

def wrap_for_slurm(command:str,run:bool,save:bool,config:Configs.Utils())->str:
    """
    This is a function that wraps a bash script in a slurm script.
    All resource allocation configuration is obtained from config argument
    Args:
        command: bash command to run in python string format
        run: if True, the slurm script is executed
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
    slurm_script = slurm_script.replace("<cpus>",config.slurm_cpus)
    
    if save:
        with open(config.slurm_outlog,'w') as f:
            f.write(slurm_script)

    if run:
        import subprocess
        subprocess.run(slurm_script)
    return slurm_script


    
    