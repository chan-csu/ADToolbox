import configs
import subprocess
import pathlib
import json
import pandas as pd
import os
from typing import Iterable, Any,Union
from warnings import warn
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
    
def get_sample_metadata_from_accession(accession:str,save:Union[None,str]=None)->dict:
    """This function returns a dictionary that contains the sample metadata from the SRA database.
    The function uses the bio command line tool to get the sample metadata. The function also saves the sample metadata if save=True.
    if the download somehow fails, the function will return a dictionary with default values like -1 and "Unknown".
    
    Requires:
        <R>project_accession</R>
        <R>Configs.Metagenomics</R>
        
    Satisfies:
        <S>sample_metadata</S>
    
    Args:
        accession (str): The accession number of the SRA project or run
        save (bool): If True, the sample metadata will be saved in the SRA work directory
        
    Returns:
        sample_metadata (dict): A dictionary that contains the sample metadata
    """
    
    metadata={}
    if os.name=="nt":
        res=subprocess.run(["bio","search",accession,"--all"],shell=True,capture_output=True)
    else:
        res=subprocess.run([f"""bio search {accession} --all"""],shell=True,capture_output=True)
    try:
        res=json.loads(res.stdout.decode("utf-8"))[0]
        metadata["host"]=res.setdefault("host","Unknown")
        metadata["library_strategy"]=res["library_strategy"].lower() if res["library_strategy"] else "Unknown"
        metadata["library_layout"]=res.setdefault("library_layout","Unknown").lower()
        metadata["base_count"]=float(res.setdefault("base_count",-1))
        metadata["read_count"]=float(res.setdefault("read_count",-1))
        avg_read_len=int(metadata["base_count"]/metadata["read_count"])
        metadata["read_length"]=avg_read_len if metadata["library_layout"]=="single" else avg_read_len/2
        metadata["sample_accession"]=accession
        metadata["study_accession"]=res["study_accession"]
    except:
        metadata["host"]="Unknown"
        metadata["library_strategy"]="Unknown"
        metadata["library_layout"]="Unknown"
        metadata["base_count"]=-1
        metadata["read_count"]=-1
        metadata["read_length"]=-1
        metadata["sample_accession"]=accession
        metadata["study_accession"]="Unknown"
    if save:
        with open(save,"w") as f:
            json.dump(metadata,f)
    return metadata

def create_mmseqs_database(fasta_db:str,mmseqs_db:str,container:str="None",save:Union[str,None]=None,run:bool=True,config=configs.Config())->str:
    """This function converts any database in fasta format for mmseqs2.
    Specifically, this function can be used for the protein database and is highly recomended for
    shotgun metagenomics samples."""
    fasta_db=os.path.abspath(fasta_db)
    mmseqs_db=os.path.abspath(mmseqs_db)
    bashscript = f"mmseqs createdb {fasta_db} {mmseqs_db}"
    db_name_path=pathlib.Path(mmseqs_db)
    if container == "None":
        pass
    
    elif container == "singularity":
        bashscript = f"singularity exec --bind {fasta_db}:{fasta_db},{str(db_name_path.parent)}:{str(db_name_path.parent)} {config.adtoolbox_singularity} {bashscript}"
    
    elif container == "docker":
        bashscript = f"docker run -v {fasta_db}:{fasta_db} -v {db_name_path.parent}:{db_name_path.parent} {config.adtoolbox_docker} {bashscript}"
    
    else:
        raise ValueError("Invalid container type. Please choose between None, singularity and docker.")
    if save:
        with open(save,'w') as f:
            f.write(bashscript)
    if run:
        subprocess.run(bashscript,shell=True)
    
    return bashscript

def index_mmseqs_db(mmseqs_db:str,container:str="None",save:Union[str,None]=None,run:bool=True,config=configs.Config)->str:
    """This function indexes any database in fasta format for mmseqs2.
    Specifically, this function can be used for the protein database and is highly recomended for
    shotgun metagenomics samples."""
    mmseqs_db=os.path.abspath(mmseqs_db)
    bashscript = f"mmseqs createindex {mmseqs_db} {mmseqs_db}"
    db_name_path=pathlib.Path(mmseqs_db)
    if container != "None":
        pass
    
    elif container == "singularity":
        bashscript = f"singularity exec --bind {str(pathlib.Path(mmseqs_db).parent)}:{str(pathlib.Path(mmseqs_db).parent)},{db_name_path.parent}:{db_name_path.parent} {config.adtoolbox_singularity} {bashscript}"
    
    elif container == "docker":
        bashscript = f"docker run -v {mmseqs_db}:{mmseqs_db} -v {db_name_path.parent}:{db_name_path.parent} {config.adtoolbox_docker} {bashscript}"
    
    else:
        raise ValueError("Invalid container type. Please choose between None, singularity and docker.")
    if save:
        with open(save,'w') as f:
            f.write(bashscript)
    if run:
        subprocess.run(bashscript,shell=True)
    
    return bashscript

def mmseqs_search(
    query_db:str,
    target_db:str,
    results_db:str,
    container:str="None",
    save:Union[str,None]=None,
    run:bool=True,
    config=configs.Config()
                        )->str:
    """This function searches a query database against a target database using mmseqs2.
    Specifically, this function can be used for the protein database and is highly recomended for
    shotgun metagenomics samples."""
    query_db=os.path.abspath(query_db)
    target_db=os.path.abspath(target_db)
    results_db=os.path.abspath(results_db)
    bashscript = f"mmseqs search {query_db} {target_db} {results_db} tmp"
    if container == "None":
        pass 
    
    elif container == "singularity":
        bashscript = f"singularity exec --bind {query_db}:{query_db},{str(pathlib.Path(target_db).parent)}:{str(pathlib.Path(target_db).parent)},{results_db}:{results_db} {config.adtoolbox_singularity} {bashscript}"
    
    elif container == "docker":
        bashscript = f"docker run -v {query_db}:{query_db} -v {str(pathlib.Path(target_db).parent)}:{str(pathlib.Path(target_db).parent)} -v {results_db}:{results_db} {config.adtoolbox_docker} {bashscript}"

    else:
        raise ValueError("Invalid container type. Please choose between None, singularity and docker.")

    if save:
        with open(save,'w') as f:
            f.write(bashscript) 
    if run:
        subprocess.run(bashscript,shell=True)
    
    return bashscript


def mmseqs_result_db_to_tsv(query_db:str,target_db:str,results_db:str,tsv_file:str,container:str="None",save:Union[str,None]=None,run:bool=True,config=configs.Config())->str:
    """This function converts the results of mmseqs search to a tsv file.
    Specifically, this function can be used for the protein database and is highly recomended for
    shotgun metagenomics samples."""
    query_db=os.path.abspath(query_db)
    target_db=os.path.abspath(target_db)
    results_db=os.path.abspath(results_db)
    tsv_file=os.path.abspath(tsv_file)
    bashscript = f"mmseqs convertalis {query_db} {target_db} {results_db} {tsv_file}"
    
    if container == "None":
        pass
    
    elif container == "singularity":
        bashscript = f"singularity exec --bind {query_db}:{query_db},{target_db}:{target_db},{results_db}:{results_db} {config.adtoolbox_singularity} {bashscript}"
    
    elif container == "docker":
        bashscript = f"docker run -v {query_db}:{query_db} -v {str(pathlib.Path(target_db).parent)}:{str(pathlib.Path(target_db).parent)} -v {results_db}:{results_db} {config.adtoolbox_docker} {bashscript}"
    
    else:
        raise ValueError("Invalid container type. Please choose between None, singularity and docker.")
    if save:
        with open(save,'w') as f:
            f.write(bashscript)
    if run:
        subprocess.run(bashscript,shell=True)
    
    return bashscript

def make_json_from_genomes(input_dir:str,output_dir:str)->dict:
    """
    This function takes a csv file directory. The CSV file must include three columns:
    1- genome_id: Identifier for the genome, preferrably; NCBI ID
    2- NCBI_Name: NCBI taxonomy name for the genome: Does not need to be in a specific format
    3- Genome_Dir: Absolute path to the fasta files: NOT .gz
    """
    genomes_json={}
    genomes_table=pd.read_table(input_dir,delimiter=",")
    genomes_table.set_index("genome_id",inplace=True)
    for[gi] in genomes_table.index:
        genomes_json[gi]={}
        genomes_json[gi]["NCBI_Name"]= genomes_table.loc[gi,"NCBI_Name"]
        genomes_json[gi]["Genome_Dir"]= genomes_table.loc[gi,"Genome_Dir"]
    with open(output_dir,"w") as fp:
        json.dump(genomes_json,fp)
    return

def needs_repair(func):
    def to_be_repaired(*args,**kwargs):
        warn("This function is not optimized yet or have issues in running. Please use it with caution.")
        return func(*args,**kwargs)
    return to_be_repaired
        
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