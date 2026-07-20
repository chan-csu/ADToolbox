from distutils.log import warn
import subprocess
import os
from collections import UserDict
import time
import json
import logging
import numpy as np
import re
import requests
import shlex
from requests.adapters import HTTPAdapter
import utils
import configs
from requests.packages.urllib3.util.retry import Retry
from requests.exceptions import Timeout
from datetime import datetime
from pathlib import Path
from collections import Counter
from collections import namedtuple
import pathlib
import asyncio
import gzip
import tomllib
import configs
from rich.progress import track,Progress
import rich
from typing import Iterable
from typing import Union
from dataclasses import dataclass
import dataclasses
from stats import scaler
from utils import (wrap_for_slurm,
                   fasta_to_dict,
                   extract_zipped_file,
                   needs_repair,
                   create_mmseqs_database,
                   index_mmseqs_db,
                   mmseqs_search,
                   mmseqs_result_db_to_tsv) 
import polars as pl
# import doctest
# doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)


def _read_json_records(path: str | os.PathLike) -> list[dict]:
    with open(path) as handle:
        payload = json.load(handle)
    if isinstance(payload, list):
        return payload
    if isinstance(payload, dict):
        return list(payload.values())
    raise ValueError(f"Expected JSON records in {path}")


def _empty_csv(path: str | os.PathLike, columns: list[str], separator: str = "\t") -> None:
    pathlib.Path(path).parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(schema={column: pl.Utf8 for column in columns}).write_csv(path, separator=separator)


def _read_table(path: str | os.PathLike, separator: str = "\t") -> pl.DataFrame:
    return pl.read_csv(path, separator=separator, infer_schema_length=None)


@dataclass
class PipelineTask:
    sample_name: str
    step_name: str
    backend: str
    command: str
    status: str
    dependencies: list[str] = dataclasses.field(default_factory=list)
    sbatch: str | None = None
    job_id: str | None = None
    submission: str | None = None


class PipelineTaskManager:
    """Small execution manager for local and Slurm pipeline tasks."""

    def __init__(
        self,
        *,
        sample_name: str,
        output_dir: str | os.PathLike,
        execution_profile: dict,
        logger: logging.Logger,
    ):
        self.sample_name = sample_name
        self.output_dir = pathlib.Path(output_dir)
        self.execution_profile = execution_profile
        self.logger = logger
        self.events_path = self.output_dir / "task_events.jsonl"
        self.events_path.parent.mkdir(parents=True, exist_ok=True)

    def record(self, event: str, task: PipelineTask, **payload) -> None:
        record = {
            "time": datetime.now().isoformat(timespec="seconds"),
            "event": event,
            "sample": task.sample_name,
            "step": task.step_name,
            "backend": task.backend,
            "status": task.status,
            "command": task.command,
            "sbatch": task.sbatch,
            "job_id": task.job_id,
            "dependencies": task.dependencies,
            **payload,
        }
        with open(self.events_path, "a") as f:
            f.write(json.dumps(record, default=str, sort_keys=True) + "\n")

    @staticmethod
    def parse_slurm_job_id(submission: str) -> str | None:
        match = re.search(r"\b(\d+)(?:\.\d+)?\b", submission or "")
        return match.group(1) if match else None

    @staticmethod
    def dependency_job_ids(dependencies: Iterable[dict] | None) -> list[str]:
        job_ids = []
        for dependency in dependencies or []:
            job_id = dependency.get("job_id") if isinstance(dependency, dict) else None
            if job_id:
                job_ids.append(str(job_id))
        return job_ids

    @staticmethod
    def _truthy(value) -> bool:
        if isinstance(value, bool):
            return value
        if value is None:
            return False
        return str(value).strip().lower() in {"1", "true", "yes", "on"}

    def wait_for_slurm_capacity(self, job_name_prefix: str) -> None:
        global_slurm = self.execution_profile.get("slurm", {})
        max_jobs = global_slurm.get("max_concurrent_jobs")
        if not max_jobs:
            return
        max_jobs = int(max_jobs)
        user = global_slurm.get("user") or os.environ.get("USER")
        if not user:
            return
        while True:
            completed = subprocess.run(
                ["squeue", "-h", "-u", str(user), "-o", "%j"],
                capture_output=True,
                text=True,
            )
            if completed.returncode:
                self.logger.warning("Could not query Slurm capacity with squeue; submitting without throttling")
                return
            active = [
                name for name in completed.stdout.splitlines()
                if name.startswith(job_name_prefix)
            ]
            if len(active) < max_jobs:
                return
            self.logger.info(
                "Slurm capacity reached for %s: %s/%s active jobs; waiting",
                job_name_prefix,
                len(active),
                max_jobs,
            )
            time.sleep(int(global_slurm.get("poll_seconds", 30)))

    def submit_slurm_job(
        self,
        *,
        sbatch_path: str | os.PathLike,
        task: PipelineTask,
        job_name_prefix: str,
        attempt: int,
    ) -> tuple[str, str | None]:
        self.wait_for_slurm_capacity(job_name_prefix)
        completed = subprocess.run(["sbatch", str(sbatch_path)], capture_output=True, text=True)
        submission = completed.stdout.strip()
        if completed.returncode:
            message = completed.stderr.strip() or submission or f"sbatch exited with status {completed.returncode}"
            task.status = "failed"
            self.record("submission_failed", task, attempt=attempt, returncode=completed.returncode, message=message)
            raise RuntimeError(f"Could not submit Slurm step {task.step_name}: {message}")
        job_id = self.parse_slurm_job_id(submission)
        task.status = "submitted"
        task.submission = submission
        task.job_id = job_id
        self.record("submitted", task, attempt=attempt, submission=submission)
        return submission, job_id

    def slurm_job_state(self, job_id: str) -> str | None:
        try:
            completed = subprocess.run(
                ["sacct", "-j", str(job_id), "--format=State", "--noheader", "--parsable2"],
                capture_output=True,
                text=True,
            )
        except FileNotFoundError:
            completed = None
        if completed is not None and completed.returncode == 0:
            states = [line.split("|", 1)[0].strip().split()[0] for line in completed.stdout.splitlines() if line.strip()]
            if states:
                for state in states:
                    if state in {"FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE", "REVOKED", "SPECIAL_EXIT"}:
                        return state
                if all(state == "COMPLETED" for state in states):
                    return "COMPLETED"
                return states[0]

        try:
            queued = subprocess.run(
                ["squeue", "-h", "-j", str(job_id), "-o", "%T"],
                capture_output=True,
                text=True,
            )
        except FileNotFoundError:
            return None
        if queued.returncode == 0 and queued.stdout.strip():
            return queued.stdout.splitlines()[0].strip().split()[0]

        return None

    def wait_for_slurm_terminal_state(self, job_id: str, poll_seconds: int) -> str | None:
        terminal_states = {
            "COMPLETED",
            "FAILED",
            "CANCELLED",
            "TIMEOUT",
            "OUT_OF_MEMORY",
            "NODE_FAIL",
            "PREEMPTED",
            "BOOT_FAIL",
            "DEADLINE",
            "REVOKED",
            "SPECIAL_EXIT",
        }
        while True:
            state = self.slurm_job_state(job_id)
            if state is None:
                self.logger.warning("Could not determine Slurm state for job %s; stopping retry monitor", job_id)
                return None
            if state in terminal_states:
                return state
            self.logger.info("Slurm job %s is %s; waiting %s seconds", job_id, state, poll_seconds)
            time.sleep(poll_seconds)


class MetagenomicsWorkflowState:
    """Persistent state and event log for resumable metagenomics batches."""

    TERMINAL_SLURM_FAILURES = {
        "FAILED",
        "CANCELLED",
        "TIMEOUT",
        "OUT_OF_MEMORY",
        "NODE_FAIL",
        "PREEMPTED",
        "BOOT_FAIL",
        "DEADLINE",
        "REVOKED",
        "SPECIAL_EXIT",
    }
    ACTIVE_SLURM_STATES = {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "SUSPENDED", "REQUEUED"}

    def __init__(self, output_dir: str | os.PathLike):
        self.output_dir = pathlib.Path(output_dir)
        self.state_path = self.output_dir / "workflow_state.json"
        self.events_path = self.output_dir / "workflow_events.jsonl"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.state = self._load()

    def _load(self) -> dict:
        if self.state_path.exists():
            with open(self.state_path) as f:
                state = json.load(f)
            state.setdefault("version", 1)
            state.setdefault("samples", {})
            return state
        return {"version": 1, "samples": {}}

    def save(self) -> None:
        with open(self.state_path, "w") as f:
            json.dump(self.state, f, indent=2, sort_keys=True, default=str)

    def sample(self, sample_name: str) -> dict:
        sample = self.state["samples"].setdefault(sample_name, {"stages": {}})
        sample.setdefault("stages", {})
        return sample

    def stage(self, sample_name: str, stage_name: str) -> dict:
        return self.sample(sample_name)["stages"].setdefault(stage_name, {})

    def record(
        self,
        sample_name: str,
        stage_name: str,
        status: str,
        *,
        artifact: dict | None = None,
        message: str | None = None,
        paths: dict | None = None,
    ) -> dict:
        previous = self.stage(sample_name, stage_name)
        attempts = int(previous.get("attempts", 0))
        if status == "submitted" and previous.get("job_id") != (artifact or {}).get("job_id"):
            attempts += 1
        entry = {
            **previous,
            "status": status,
            "updated_at": datetime.now().isoformat(timespec="seconds"),
            "attempts": attempts,
        }
        if artifact:
            entry["artifact"] = artifact
            if artifact.get("job_id"):
                entry["job_id"] = artifact["job_id"]
            if artifact.get("sbatch"):
                entry["sbatch"] = artifact["sbatch"]
        if message:
            entry["message"] = message
        if paths:
            entry["paths"] = paths
        self.sample(sample_name)["stages"][stage_name] = entry
        self._event(sample_name, stage_name, status, entry, message=message)
        self.save()
        return entry

    def _event(self, sample_name: str, stage_name: str, status: str, entry: dict, *, message: str | None = None) -> None:
        record = {
            "time": datetime.now().isoformat(timespec="seconds"),
            "sample": sample_name,
            "stage": stage_name,
            "status": status,
            "attempts": entry.get("attempts", 0),
            "job_id": entry.get("job_id"),
            "message": message,
        }
        with open(self.events_path, "a") as f:
            f.write(json.dumps(record, default=str, sort_keys=True) + "\n")

    def active_submission(
        self,
        sample_name: str,
        stage_name: str,
        *,
        slurm_checker: PipelineTaskManager | None = None,
    ) -> bool:
        entry = self.stage(sample_name, stage_name)
        if entry.get("status") not in {"submitted", "running", "monitoring"}:
            return False
        job_id = entry.get("job_id")
        if not job_id or slurm_checker is None:
            return True
        state = slurm_checker.slurm_job_state(str(job_id))
        if state is None:
            return True
        if state == "COMPLETED":
            self.record(sample_name, stage_name, "completed", message=f"Slurm job {job_id} completed")
            return False
        if state in self.TERMINAL_SLURM_FAILURES:
            self.record(sample_name, stage_name, "failed", message=f"Slurm job {job_id} ended as {state}")
            return False
        return state in self.ACTIVE_SLURM_STATES

@dataclass
class Feed:

    """
    The Feed class is used to store the feed information, and later use it in the e_adm model.
    all the entered numbers must in percentages. Carbohudrates, lipids, and proteins and si must sum up to 100, 
    and they form the total dissolved solids. Carbohydrates, lipids, proteins, and xi must sum up to 100, and they form the total suspended solids.
    
    IMPORTANT: It is assumed that lipid, proteins and carbohydrates have the same fraction in soluble and insoluble fractions.
    
    Args:
        name (str): A unique name for the feed.
        carbohydrates (float): percentage of carbohydrates in the feed.
        lipids (float): percentage of lipids in the feed.
        proteins (float): percentage of proteins in the feed.
        tss (float): percentage of total COD in the form of suspended solids.
        si (float): percentage of percentage of soluble inorganics in the TDS.
        xi (float): percentage of percentage of insoluble inorganics in the TSS.
        reference (str, optional): A reference for the feed data. Defaults to ''.    
    
    Examples:
        >>> feed=Feed(name="Test",carbohydrates=25,lipids=25,proteins=25,si=25,xi=25,tss=70)
        >>> assert feed.ch_tss==feed.lip_tss==feed.prot_tss==feed.xi_tss==0.25
        
    """
    # total_cod:float Transfer to base parameters
    name:str            # A unique name for the feed
    carbohydrates:float # percentage of carbohydrates in the feed
    lipids:float        # percentage of lipids in the feed
    proteins:float      # percentage of proteins in the feed
    tss:float           # percentage of total COD in the form of suspended solids
    si:float            # percentage of percentage of soluble inorganics in the TDS
    xi:float            # percentage of percentage of insoluble inorganics in the TSS
    reference:str=''    # A reference for the feed data

    def __post_init__(self):
        if self.carbohydrates+self.lipids+self.proteins>100:
            raise ValueError("The sum of the percentages must less than 100")
        if self.carbohydrates+self.lipids+self.proteins+self.si<1:
            warn("The sum of lipids, carbohydrates, proteins are suspiciously low! Make sure youhave input the numbers in percentages!")
        li_ch_pr=self.carbohydrates+self.lipids+self.proteins
        without_xi=100-self.xi
        self.ch_tss=self.carbohydrates/li_ch_pr*without_xi/100
        self.lip_tss=self.lipids/li_ch_pr*without_xi/100
        self.prot_tss=self.proteins/li_ch_pr*without_xi/100
        self.xi_tss=self.xi/100
        without_si=100-self.si
        self.ch_tds=self.carbohydrates/li_ch_pr*without_si/100
        self.lip_tds=self.lipids/li_ch_pr*without_si/100
        self.prot_tds=self.proteins/li_ch_pr*without_si/100
        self.si_tds=self.si/100
    
    def to_dict(self)->dict:
        return {"name":self.name,
                "carbohydrates":self.carbohydrates,
                "lipids":self.lipids,
                "proteins":self.proteins,
                "tss":self.tss,
                "si":self.si,
                "xi":self.xi,
                "reference":self.reference}

          
@dataclass
class Experiment:
    """
    This class creates an interface for the experimental data to be used in different places in ADToolbox.
    First you should give each experiment a name. Time must be a list of time points in days, and there must be a time 0 point assinged to each experiment.
    variables must be a list of integers that represent the variables that are the index of the ADM species that we have concentration data for.
    data must be a list of lists. Each list in the list must be a list of concentrations for each species at each time point.
    IMPORTANT: The order of the species in the data list must match the order of the species in the variables list.
    if there are specific initial concentrations for the ADM species, they can be passed as a dictionary to the initial_concentrations argument.
    reference is an optional argument that can be used to provide a reference for the experimental data. If using the database module 
    to query for Experiment objects you can query by name or reference or model_type. So, having a descriptive reference can be useful for querying as well.
    default model name is "e_adm". This can be changed by passing a different model name to the model_name argument. This also helps with querying.
    
    Args:
        name (str): A unique name for the experiment.
        time (list): A list of time points in days.
        variables (list): A list of strings that represent the species that are in the used model (most commonly ADM) that we have concentration data for.
        data (list): A list of lists. Each list in the list must be a list of concentrations for each species at each time point.
        feed (Feed): An instance of feed class
        initial_concentrations (dict, optional): A dictionary of initial concentrations for the ADM species. Defaults to {}.
        base_parameters (dict, optional): A dictionary of base parameters for the model. Defaults to {}.
        constants (list, optional): A list of strings that represent the species that are in the used model (most commonly ADM) that are held constant during the simulations. Their value will come from the initial_concentrations. Defaults to [].
        reference (str, optional): A reference for the experimental data. Defaults to ''.
        model_name (str, optional): The name of the model that the experimental data is for. Defaults to "e_adm".
    
    Examples:
        >>> import json
        >>> with open(configs.Database().adm_parameters["species"],"r") as f:
        ...     species=json.load(f)
        >>> feed=Feed(name="Test",carbohydrates=25,lipids=25,proteins=25,si=25,xi=25,tss=70)
        >>> exp=Experiment(name="Test",time=[0,1,2],variables=["S_su","S_aa"],data=[[1,2,3],[4,5,6]],feed=Feed,reference="Test reference")
        
    """
    name:str
    time: list[float]
    variables: list[str]
    data: list[list[float]]
    feed: Feed
    initial_concentrations: dict[str,float] = dataclasses.field(default_factory=dict)
    base_parameters: dict[str,float] = dataclasses.field(default_factory=dict)
    constants: list[str] = dataclasses.field(default_factory=tuple)
    reference: str = ""
    model_name: str = "e_adm"
    
    
    def __post_init__(self):
        self.data=np.array(self.data).T
        self.validate()
    
    def validate(self):
        assert len(self.time)==self.data.shape[0], "Number of time points must match number of rows in data."
        assert len(self.variables)==self.data.shape[1] , "Number of variables must match number of columns in data."
        assert self.time[0]==0, "Time must start at 0."
        return "successful"
    
    def to_dict(self):
        return {"name":self.name,
                "time":self.time,
                "variables":self.variables,
                "data":self.data.T.tolist(),
                "feed":self.feed.to_dict(),
                "initial_concentrations":self.initial_concentrations,
                "base_parameters":self.base_parameters,
                "constants":self.constants,
                "reference":self.reference,
                "model_name":self.model_name}
    

    
    
    

@dataclass
class MetagenomicsStudy:
    """
    This class is used to communicate between the metagenomics studies database and the ADM model.
    
    Args:
        name (str): The name of the metagenomics study. Its okay if it is not unique.
        study_type (str): The type of the metagenomics study. It can be "amplicon" or "WGS".
        microbiome (str): The microbiome that the metagenomics study is about.
        sample_accession (str): The SRA sample accession number of the metagenomics study. This must be unique.
        comments (str): Any comments that you want to add to the metagenomics study.
        study_accession (str): The SRA study accession number of the metagenomics study.   
    
    Examples:
        >>> study=MetagenomicsStudy(name="Test",study_type="WGS",microbiome="test_microbiome",sample_accession="test_accession",comments="test_comments",study_accession="test_study_accession")
        >>> assert study.name=="Test"

    """
    name:str
    study_type:str
    microbiome:str
    sample_accession:str
    comments:str
    study_accession:str
    
    def to_dict(self)->dict:
        return {"name":self.name,
                "study_type":self.study_type,
                "microbiome":self.microbiome,
                "sample_accession":self.sample_accession,
                "comments":self.comments,
                "study_accession":self.study_accession}

class Reaction:
    """
    This class provides a simple interface between information about biochemical reactions and multiple functionalities of ADToolbox.
    In order to instantiate a reaction object, you need to pass a dictionary of the reaction information.
    This dictionary must include 'name','stoichiometry' keys. This follows the format of the seed database.
    stoichiometry must be formatted like seed database. The seed database format is as follows:
    stoichiometry: '-1:cpd00079:0:0:\"D-glucose-6-phosphate\";1:cpd00072:0:0:\"D-fructose-6-phosphate\"'

    Args:
        data (dict): A dictionary containing the reaction information. This follows the format of the seed database.


    Examples:
        >>> A={"name":'D-glucose-6-phosphate aldose-ketose-isomerase',"stoichiometry":'-1:cpd00079:0:0:\"D-glucose-6-phosphate\";1:cpd00072:0:0:\"D-fructose-6-phosphate\"'}
        >>> a=Reaction(A)
        >>> print(a)
        D-glucose-6-phosphate aldose-ketose-isomerase

    """
    def __init__(self, data:dict)->None:
        self.data = data

    def __str__(self)->str:
        return self.data['name']

    @property
    def stoichiometry(self)->dict:
        """
        Returns the stoichiometry of the reaction by the seed id of the compounds as key and the
        stoichiometric coefficient as value.
        Examples:
            >>> A={"name":'D-glucose-6-phosphate aldose-ketose-isomerase',"stoichiometry":'-1:cpd00079:0:0:\"D-glucose-6-phosphate\";1:cpd00072:0:0:\"D-fructose-6-phosphate\"'}
            >>> a=Reaction(A)
            >>> a.stoichiometry=={'cpd00079': -1, 'cpd00072': 1}
            True
        
        Args:
            self (Reaction): An instance of the Reaction.

        Returns:
            dict: The stoichiometry of the reaction 
        """
        return {compound.split(':')[1]:float(compound.split(':')[0]) for compound in self.data['stoichiometry'].split(';') }


class Metabolite:
    """
    This class provides a simple interface between information about metabolites and multiple functionalities of ADToolbox.
    In order to instantiate a metabolite object, you need to pass a dictionary of the metabolite information.
    This dictionary must include 'name','mass','formula' keys. This follows the format of the seed database.
    formula must be formatted like seed database. The seed database format is as follows:
    formula: 'C6H12O6'
    Possibly the main advantage of instantiating a metabolite object is that it provides a COD attribute that can be used to convert
    the concentration of the metabolite from g/l to gCOD/l. This is useful for comparing the experimental data with the model outputs.

    Args:
        data (dict): A dictionary containing the metabolite information. This follows the format of the seed database.


    Examples:
        >>> A={"name":"methane","mass":16,"formula":"CH4"}
        >>> a=Metabolite(A)
        >>> print(a)
        methane

    """

    def __init__(self, data):
        self.data = data
        self.cod = self.cod_calc()
        self.mw= self.data.get('mass',None)

    def __str__(self) -> str:
        return self.data['name']

    def cod_calc(self,add_h:float=0,add_c:float=0,add_o:float=0)->float:
        """
        Calculates the conversion rates for g/l -> gCOD/l
        In some cases we would like to add extra atoms for COD calculations
        For example, model seed biochemistry database only uses acetate instead of acetic acid.
        The 1 hydrogen difference changes the COD conversion rate. For this reason we can add extra atoms to the formula
        to calculate the COD conversion rate without changing anything else.
        
        Args:
            add_h (float): The number of extra hydrogen atoms to add to the formula for COD calculation.
            add_c (float): The number of extra carbon atoms to add to the formula for COD calculation.
            add_o (float): The number of extra oxygen atoms to add to the formula for COD calculation.

        Examples:
            >>> A={"name":"methane","mass":16,"formula":"CH4"}
            >>> a=Metabolite(A)
            >>> a.cod
            4.0

        Args:
            self (Metabolite): An instance of the Metabolite class: Note

        Returns:
            float: COD conversion from g/l to gCOD/l

        """
        if self.data['formula'] and self.data['mass']:
            contents = {}
            atoms = ["H", "C", "O"]
            mw = self.data['mass']+add_h*1+add_c*12+add_o*16
            for atom in atoms:
                pattern = atom + r'\d*'
                if re.search(pattern, self.data['formula']):
                    if len(re.search(pattern, self.data['formula']).group()[1:]) == 0:
                        contents[atom] = 1
                    else:
                        contents[atom] = int(
                            re.search(pattern, self.data['formula']).group()[1:])
                else:
                    contents[atom] = 0
            contents['H']+=add_h
            contents['C']+=add_c
            contents['O']+=add_o
            cod_conv=1/mw*(contents['H']+4*contents['C']-2*contents['O'])/4*32
            return cod_conv

        else:
            return 'None'

        


class SeedDB:

    """
    This class is designed to interact with seed database. The main advantage of using this class is that it can be used to instantiate
    a reaction and metabolite object, and it provides extra functionalities that rely on information in the seed database. For example, 
    If there is a chemical formula assigned to a metabolite in the seed database, then the informattion about the COD of that metabolite
    can be computed using the chemical formula. 
    
    Args:
        config (configs.SeedDB): An instance of the SeedDB class in the configs module. This class contains the information about the seed database.
    
    Examples:
        >>> seed_db=SeedDB(configs.Database())
        >>> assert seed_db.compound_db==configs.Database().compound_db
        >>> assert seed_db.reaction_db==configs.Database().reaction_db

    """

    def __init__(self, config:configs.Database) -> None:
        
        self.reaction_db = config.reaction_db
        self.compound_db = config.compound_db

    def instantiate_rxns(self, seed_id:str)->Reaction:
        """
        This method is used to instantiate reaction objects from the seed database.
        in order to instantiate a reaction object, you need to pass the seed identifier for that reaction.
        
        Args:
            seed_id (str): The seed identifier for the reaction.
    
        Returns:
            Reaction: An instance of the Reaction class.
            
        Required Configs:
            - config.reaction_db
        
        Examples:
            >>> db_conf=configs.Database()
            >>> seed_db=SeedDB(db_conf)
            >>> rxn=seed_db.instantiate_rxns("rxn00558")
            >>> assert rxn.data["name"]=="D-glucose-6-phosphate aldose-ketose-isomerase"
        """
        records = _read_json_records(self.reaction_db)
        return Reaction(data=next(record for record in records if record.get("id") == seed_id))

    def instantiate_metabs(self, seed_id:str)->Metabolite:
        """
        This method is used to instantiate metabolite objects from the seed database.
        In order to instantiate a metabolite object, you need to pass the seed identifier for that metabolite.

        Args:
            seed_id (str): The seed identifier for the metabolite.
        
        Returns:
            Metabolite: An instance of the Metabolite class. 
        
        Required Configs:
            - config.compound_db
        
        Examples:
            >>> seed_db=SeedDB(configs.Database())
            >>> metab=seed_db.instantiate_metabs("cpd01024")
            >>> assert metab.cod==4.0
        """
        records = _read_json_records(self.compound_db)
        return Metabolite(data=next(record for record in records if record.get("id") == seed_id))

    def get_seed_rxn_from_ec(self, ec_number:str)->list:
        """
        This method is used to get the seed reaction identifiers for a given EC number.

        Args:
            ec_number (str): The EC number.
        
        Returns:
            list: A list of seed reaction identifiers.
        
        Required Configs:
            - config.reaction_db
        
        Examples:
            >>> seed_db=SeedDB(config=configs.Database())
            >>> seed_rxn_list=seed_db.get_seed_rxn_from_ec("1.1.1.1")
            >>> assert len(seed_rxn_list)>0
        
        """
        seen = set()
        matches = []
        for record in _read_json_records(self.reaction_db):
            if ec_number not in (record.get("ec_numbers") or []):
                continue
            if record.get("id") in seen:
                continue
            seen.add(record.get("id"))
            matches.append(record)
        return matches
            

class Database:

    '''
    This class is designed to supply any data requirement for ADToolbox. All functionalisties for saving, loading, and querying data are implemented here.
    ADToolbox in general contains the following databases:
    
    - The seed reaction database
    
    - The seed compound database
    
    - ADToolbox's Feed database
    
    - ADToolbox's Metagenomics studies database
    
    - ADToolbox's Experimental data database
    
    - ADToolbox's Protein database
    
    - ADToolbox's Reaction database
    
    - GTDB-tk database for bacterial and archaeal 16s rRNA sequences
    
    - ADM and e_adm model parameters
    
    This class is instantiated with a configs.Database object. This object contains the paths to all the databases that ADToolbox uses.
    Please refer to the documentation of each method for more information on the required configurations.
    
    Args:
        config (configs.Database, optional): A configs.Database object. Defaults to configs.Database().
    
    Examples:
        >>> db=Database(config=configs.Database())
        >>> assert type(db)==Database and type(db.config)==configs.Database

    '''
    def __init__(self, config:configs.Database|None=None)->None:
        self.config = config or configs.Database()


    def initialize_protein_db(self)->None:
        """This function intializes ADToolbox's protein database by creating an empty fasta file.
        Be careful, this will overwrite any existing file with the same name.
        Logically, this needs method needs config.protein_db to be defined.
        
        Required Configs:
            - config.protein_db
            --------
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==False # This is just to make sure that the following lines create the file
            >>> db=Database(config=configs.Database(protein_db=os.path.join(Main_Dir,"protein_test_db.fasta"))) # point to a test non-existing file
            >>> db.initialize_protein_db() # initialize the protein database
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==True # check if the file is created
            >>> os.remove(os.path.join(Main_Dir,"protein_test_db.fasta")) # remove the file to clean up
        """

        if not (pathlib.Path(self.config.protein_db).parent).exists():
            pathlib.Path(self.config.protein_db).parent.mkdir(parents=True)
        with open(self.config.protein_db, 'w') as f:
            pass
    
    def initialize_reaction_db(self)->None:
        r"""This function intializes ADToolbox's reaction database by creating an empty tsv file.
        Be careful, this will overwrite any existing file with the same name.
        
        Required Configs:
            - config.reaction_db
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"reaction_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(reaction_db=os.path.join(Main_Dir,"reaction_test_db.tsv")))
            >>> db.initialize_reaction_db()
            >>> assert pd.read_table(os.path.join(Main_Dir,"reaction_test_db.tsv"),delimiter="\t").shape[0]==0
            >>> assert set(pd.read_csv(os.path.join(Main_Dir,"reaction_test_db.tsv"),delimiter="\t").columns)==set(["ec_numbers","seed_ids","reaction_names","adm1_reaction","e_adm_reactions","pathways"])
            >>> os.remove(os.path.join(Main_Dir,"reaction_test_db.tsv"))
        
        """
        _empty_csv(self.config.reaction_db, ["ec_numbers","seed_ids","reaction_names","adm1_reaction","e_adm_reactions","pathways"])
        
    def initialize_feed_db(self)->None:
        r"""This function intializes ADToolbox's Feed database by creating an empty tsv file.
        Be careful, this will overwrite any existing file with the same name.
        
        Required Configs:
            - config.feed_db
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(feed_db=os.path.join(Main_Dir,"feed_test_db.tsv")))
            >>> db.initialize_feed_db()
            >>> assert pd.read_table(os.path.join(Main_Dir,"feed_test_db.tsv"),delimiter='\t').shape[0]==0
            >>> assert set(pd.read_table(os.path.join(Main_Dir,"feed_test_db.tsv"),delimiter='\t').columns)==set(["name","carbohydrates","lipids","proteins","tss","si","xi","reference"])
            >>> os.remove(os.path.join(Main_Dir,"feed_test_db.tsv"))
        
        """
        _empty_csv(self.config.feed_db, ["name","carbohydrates","lipids","proteins","tss","si","xi","reference"])
    
    def initialize_metagenomics_studies_db(self)->None:
        r"""This function intializes ADToolbox's Metagenomics studies database by creating an empty tsv file.
        Be careful, this will overwrite any existing file with the same name.
        
        Required Configs:
            - config.metagenomics_studies_db
        
        Examples:
            >>> import os
            >>> local_dir={'metagenomics_studies':os.path.join(Main_Dir,"test","metagenomics_test_db.tsv")}
            >>> db=Database(config=configs.Database(studies_local=local_dir))
            >>> db.initialize_metagenomics_studies_db()
            >>> assert pd.read_table(local_dir['metagenomics_studies'],delimiter="\t").shape[0]==0
            >>> assert set(pd.read_table(local_dir['metagenomics_studies'],delimiter="\t").columns)==set(["name","study_type","microbiome","sample_accession","comments","study_accession"])
            >>> os.remove(local_dir['metagenomics_studies'])
         
        """
        _empty_csv(self.config.studies_local["metagenomics_studies"], ["name","study_type","microbiome","sample_accession","comments","study_accession"])
        
    def initialize_experimental_data_db(self)->None:
        """This function intializes ADToolbox's experimental data database by creating an empty json file.
        Be careful, this will overwrite any existing file with the same name.

        Required Configs:
            - config.experimental_data_db
        
        Examples:
            >>> import os,json
            >>> local_dir={'experimental_data_db':os.path.join(Main_Dir,"test","experiments_test_db.json")}
            >>> db=Database(config=configs.Database(studies_local=local_dir))
            >>> db.initialize_experimental_data_db()
            >>> assert pd.read_json(local_dir['experimental_data_db']).shape[0]==0
            >>> with open(local_dir['experimental_data_db'],"r") as f:
            ...     assert json.load(f)==[]
            >>> os.remove(local_dir['experimental_data_db'])
        """
        if not (pathlib.Path(self.config.studies_local["experimental_data_db"]).parent).exists():
            pathlib.Path(self.config.studies_local["experimental_data_db"]).parent.mkdir(parents=True)
        with open(self.config.studies_local["experimental_data_db"], "w") as handle:
            json.dump([], handle)
        
    
    def filter_seed_from_ec(self, 
                            ec_list:list[str],
                            save:bool=False) -> tuple:
        """
        This function takes a list of EC numbers and filters the seed database to find the seed reactions that have the EC numbers in their EC number list.
        This will help to trim the large seed database to a smaller one that only contains the reactions that are relevant to the AD process.

        Args:
            ec_list (list[str]): A list of EC numbers.
            save (bool, optional): Whether to save the filtered seed database or not. Defaults to False.
        
        Returns:
            tuple: A tuple containing the filtered seed reaction database and the seed compound database, respectively.
        
        Required Configs:
        
            - config.reaction_db
            --------
            - config.compound_db
            --------
            - config.local_reaction_db
            --------
            - config.local_compound_db
            --------
            
            
        Examples:
            >>> db=Database()
            >>> seed_rxn_db,seed_compound_db=db.filter_seed_from_ec(["1.1.1.1","1.1.1.2"])
            >>> assert len(seed_rxn_db)>0 and len(seed_compound_db)>0
            >>> assert pd.read_json(configs.Database().reaction_db).shape[0]>pd.DataFrame(seed_rxn_db).shape[0]
        """
        seed_rxn_db = [
            record for record in _read_json_records(self.config.reaction_db)
            if any(ec in (record.get("ec_numbers") or []) for ec in ec_list)
        ]
        stoichiometry_ids = {
            metabolite
            for record in seed_rxn_db
            for metabolite in (record.get("stoichiometry") or [])
        }
        seed_compound_db = [
            record for record in _read_json_records(self.config.compound_db)
            if record.get("id") in stoichiometry_ids
        ]
        if save:
            with open(self.config.local_reaction_db, "w") as handle:
                json.dump(seed_rxn_db, handle)
            with open(self.config.local_compound_db, "w") as handle:
                json.dump(seed_compound_db, handle)
        return seed_rxn_db, seed_compound_db
        
            

    def get_protein_seqs_from_uniprot(self, uniprot_id:str) -> str:
        """
        This function takes a uniprot id and fetches the protein sequence from Uniprot.

        Args:
            uniprot_id (str): The uniprot id of the protein.
        
            
        Returns:
            str: The protein sequence.
        
        Examples:
            >>> db=Database()
            >>> seq=db.get_protein_seqs_from_uniprot("P0A9P0")
            >>> assert type(seq)==str and len(seq)>0
        """
        Base_URL = "https://rest.uniprot.org/uniprotkb/"
        session = requests.Session()
        retry = Retry(connect=3, backoff_factor=0.5)
        adapter = HTTPAdapter(max_retries=retry)
        session.mount('http://', adapter)
        try:
            file = session.get(
                f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta", timeout=10)
        except:
            print("Could not fetch the sequence! Trying again ...")
            while True:
                time.sleep(5)
                file = session.get(Base_URL+uniprot_id+".fasta", timeout=10)
                if file.ok:
                    break
           
        return ''.join(file.text.split('\n')[1:-1])
   
    def proteins_from_ec(self,ec_number:str) -> dict:
        """
        This function returns a dictionary of protein sequences for a given EC number.
        The keys are the uniprot ids and ec number compatible with ADToolbox protein database
        and the values are the protein sequences. Since ADToolbox deals with microbial process,
        only bacterial and archaeal proteins are considered.

        Args:
            ec_number (str): The EC number.
        
        Returns:
            dict: A dictionary of protein sequences.
            
        Examples:
            >>> db=Database()
            >>> protein_seqs=db.proteins_from_ec("1.1.1.1")
            >>> assert len(protein_seqs)>0
            >>> assert list(protein_seqs.keys())[0].split("|")[1]=="1.1.1.1"
        """
        session = requests.Session()
        retry = Retry(connect=3, backoff_factor=0.5)
        adapter = HTTPAdapter(max_retries=retry)
        session.mount('http://', adapter)
        protein_seqs={}
        try:
            file = session.get(
                f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28ec%3A{ec_number}%29%20AND%20%28reviewed%3Atrue%29%20NOT%20%28taxonomy_id%3A2759%29%29", timeout=30)
        except requests.exceptions.HTTPError or requests.exceptions.ConnectionError:
            print("Request Error! Trying again ...")
            time.sleep(30)
            file = session.get(
                f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28ec%3A{ec_number}%29%20AND%20%28reviewed%3Atrue%29%20NOT%20%28taxonomy_id%3A2759%29%29", timeout=30)
        # This alsp does a sanity chec
        except Exception:
            print('Something went wrong!')
        text = file.text
        if text:
            text=text.split('>')
            text.remove("")
            for seq in text:
                protein_seqs.update([(seq.split("\n")[0].split("|")[1]+"|"+ec_number, "".join(seq.split("\n")[1:]))])
                
        
        return protein_seqs


    def build_protein_db_from_reactions_db(self):
        r"""
        This function builds the protein database from the reaction database.
        It takes the reaction database and finds the protein sequences for each EC number in the reaction database.
        Then it saves the protein sequences in a fasta file.

        Required Configs:
            - config.reaction_db
            --------
            - config.protein_db
            --------
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==False
            >>> assert os.path.exists(os.path.join(Main_Dir,"reaction_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(protein_db=os.path.join(Main_Dir,"protein_test_db.fasta"),reaction_db=os.path.join(Main_Dir,"reaction_test_db.tsv")))
            >>> reaction_db=pd.DataFrame(columns=["EC_Numbers","Seed Ids","Reaction Names","ADM1_Reaction","e_adm_Reactions","Pathways"])
            >>> reaction_db.loc[0,"EC_Numbers"]="1.1.1.1"
            >>> reaction_db.to_csv(os.path.join(Main_Dir,"reaction_test_db.tsv"),index=False,sep="\t")
            >>> db.build_protein_db_from_reactions_db()
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==True
            >>> assert os.path.exists(os.path.join(Main_Dir,"reaction_test_db.tsv"))==True
            >>> assert os.path.getsize(os.path.join(Main_Dir,"protein_test_db.fasta"))>0
            >>> os.remove(os.path.join(Main_Dir,"protein_test_db.fasta"))
            >>> os.remove(os.path.join(Main_Dir,"reaction_test_db.tsv"))
        """
        rxn_db=_read_table(self.config.reaction_db)
        ec_numbers=list(set(rxn_db["EC_Numbers"].to_list()))
        protein_seqs={}
        for ec in ec_numbers:
            protein_seqs.update(self.proteins_from_ec(ec))
        with open(self.config.protein_db,"w") as f:
            for key,value in protein_seqs.items():
                f.write(">"+key+"\n")
                f.write(value+"\n")

    def cazy_ec(self)->list:
        pass
          
    def add_protein_to_protein_db(self, protein_id:str, header_tail:str)->None:
        """
        This funciton adds a protein sequence to the protein database. It takes a uniprot id and an EC number it is assigned to 
        and adds the corresponding protein sequence to the protein database.
        
        Required Configs:
            - config.protein_db

        Args:
            protein_id (str): The uniprot id of the protein.
            header_tail (str): A text to append to the header of the entry in the database.
                In ADToolbox it is better to use ec number for compatibility with downstream functions.
        
    
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==False
            >>> db=Database(config=configs.Database(protein_db=os.path.join(Main_Dir,"protein_test_db.fasta")))
            >>> db.add_protein_to_protein_db("P0A9P0","1.2.3.4")
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==True
            >>> assert os.path.getsize(os.path.join(Main_Dir,"protein_test_db.fasta"))>0
            >>> import utils
            >>> assert len(utils.fasta_to_dict(os.path.join(Main_Dir,"protein_test_db.fasta")))>0
            >>> os.remove(os.path.join(Main_Dir,"protein_test_db.fasta"))
        """
        if not os.path.exists(self.config.protein_db):
            self.initialize_protein_db()
        with open(self.config.protein_db,"a") as f:
            f.write(">"+protein_id+"|"+header_tail+"\n")
            f.write(self.get_protein_seqs_from_uniprot(protein_id)+"\n")
            
    def add_proteins_from_ecnumbers_to_protein_db(self, ec_numbers:list)->None:
        """
        This function adds protein sequences to the protein database from a list of EC numbers.
        It takes a list of EC numbers and finds the protein sequences for each EC number in the list.
        Then it saves the protein sequences in a fasta file.
        
        Required Configs:
            - config.protein_db
        
        Args:
            ec_numbers (list): A list of EC numbers.
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==False
            >>> db=Database(config=configs.Database(protein_db=os.path.join(Main_Dir,"protein_test_db.fasta")))
            >>> db.add_proteins_from_ecnumbers_to_protein_db(["1.1.1.1","1.1.1.2"])
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==True
            >>> import utils
            >>> assert len(utils.fasta_to_dict(os.path.join(Main_Dir,"protein_test_db.fasta")))>0
            >>> os.remove(os.path.join(Main_Dir,"protein_test_db.fasta"))
        """
        if not os.path.exists(self.config.protein_db):
            self.initialize_protein_db()
        
        protein_seqs={}
        for ec in ec_numbers:
            protein_seqs.update(self.proteins_from_ec(ec))
        
        with open(self.config.protein_db,"a") as f:
            for key,value in protein_seqs.items():
                f.write(">"+key+"\n")
                f.write(value+"\n")
        
    def add_feed_to_feed_db(self,feed:Feed)->None:
        r"""
        This function adds a feed to the feed database. It takes the feed name and the feed composition and adds them to the feed database.

        Required Configs:
            - config.feed_db

        Args:
            feed (Feed): An instance of the Feed class.
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(feed_db=os.path.join(Main_Dir,"feed_test_db.tsv")))
            >>> feed=Feed(name="test_feed",carbohydrates=10,lipids=20,proteins=30,tss=80,si=10,xi=30,reference="test")
            >>> db.add_feed_to_feed_db(feed)
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==True
            >>> assert pd.read_table(os.path.join(Main_Dir,"feed_test_db.tsv"),delimiter="\t").shape[0]>0
            >>> os.remove(os.path.join(Main_Dir,"feed_test_db.tsv"))
        
        """
        if not os.path.exists(self.config.feed_db):
            self.initialize_feed_db()
            
        feed_db = _read_table(self.config.feed_db)
        if feed.name in feed_db["name"].to_list():
            raise ValueError("Feed already exists in the database.")
        pl.concat([feed_db, pl.DataFrame([feed.to_dict()])], how="diagonal_relaxed").write_csv(self.config.feed_db, separator="\t")
    
    def remove_feed_from_feed_db(self,field_name:str,query:str)->None:
        r"""
        This function removes studyes that contain the query in the given column, field name, from the feed database.

        Required Configs:
            - config.feed_db
        
        Args:
            field_name (str): The name of the column to query.
            query (str): The query string.
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(feed_db=os.path.join(Main_Dir,"feed_test_db.tsv")))
            >>> feed=Feed(name="test_feed",carbohydrates=10,lipids=20,proteins=30,tss=80,si=10,xi=30,reference="test")
            >>> db.add_feed_to_feed_db(feed)
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==True
            >>> assert pd.read_table(os.path.join(Main_Dir,"feed_test_db.tsv"),delimiter="\t").shape[0]>0
            >>> db.remove_feed_from_feed_db("name","test_feed")
            >>> assert pd.read_table(os.path.join(Main_Dir,"feed_test_db.tsv"),delimiter="\t").shape[0]==0
            >>> os.remove(os.path.join(Main_Dir,"feed_test_db.tsv"))
        
        """
        if not os.path.exists(self.config.feed_db):
            raise FileNotFoundError("Feed database does not exist!")
        
        
        feed_db=_read_table(self.config.feed_db)
        feed_db.filter(~pl.col(field_name).str.contains(query, literal=True)).write_csv(self.config.feed_db, separator="\t")
        
    def get_feed_from_feed_db(self,field_name:str,query:str)->list[Feed]:
        r"""
        This function returns a feed from the feed database. It takes the query string and the column name to query and returns the feed that contains the query string in the given column.

        Required Configs:
            - config.feed_db
        
        Args:
            field_name (str): The name of the column to query.
            query (str): The query string.
        
        Returns:
            Feed: An instance of the Feed class.
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(feed_db=os.path.join(Main_Dir,"feed_test_db.tsv")))
            >>> feed=Feed(name="test_feed",carbohydrates=10,lipids=20,proteins=30,tss=80,si=10,xi=30,reference="test")
            >>> db.add_feed_to_feed_db(feed)
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==True
            >>> assert pd.read_table(os.path.join(Main_Dir,"feed_test_db.tsv"),delimiter="\t").shape[0]>0
            >>> feed=db.get_feed_from_feed_db("name","test_feed")
            >>> assert feed[0].name=="test_feed"
            >>> os.remove(os.path.join(Main_Dir,"feed_test_db.tsv"))
        
        """
        if not os.path.exists(self.config.feed_db):
            raise FileNotFoundError("Feed database does not exist!")
        
        feed_db=_read_table(self.config.feed_db)
        feed_db=feed_db.filter(pl.col(field_name).str.contains(query, literal=True))
        return [Feed(**feed) for feed in feed_db.to_dicts()]
    
    def add_metagenomics_study_to_metagenomics_studies_db(self,metagenomics_study:MetagenomicsStudy)->None:
        r"""
        This function adds a metagenomics study to the metagenomics studies database. It takes a metagenomics study and adds it to the metagenomics studies database.
        
        Required Configs:
            - config.metagenomics_studies_db
        
        Args:
            metagenomics_study (MetagenomicsStudy): An instance of the MetagenomicsStudy class.

        Examples:
            >>> import os
            >>> local_dir={'metagenomics_studies':os.path.join(Main_Dir,"test","metagenomics_test_db.tsv")}
            >>> db=Database(config=configs.Database(studies_local=local_dir))
            >>> metagenomics_study=MetagenomicsStudy(name="test_study",study_type="metagenomics",microbiome="anaerobic digester",sample_accession="test",comments="test",study_accession="test")
            >>> db.add_metagenomics_study_to_metagenomics_studies_db(metagenomics_study)
            >>> assert os.path.exists(local_dir['metagenomics_studies'])==True
            >>> assert pd.read_table(local_dir['metagenomics_studies'],delimiter="\t").shape[0]>0
            >>> os.remove(local_dir['metagenomics_studies'])
        """
        if not os.path.exists(self.config.studies_local["metagenomics_studies"]):
            self.initialize_metagenomics_studies_db()
        metagenomics_studies_db=_read_table(self.config.studies_local["metagenomics_studies"])
        pl.concat([metagenomics_studies_db, pl.DataFrame([metagenomics_study.to_dict()])], how="diagonal_relaxed").write_csv(
            self.config.studies_local["metagenomics_studies"], separator="\t"
        )
    
    def remove_metagenomics_study_from_metagenomics_studies_db(self,field_name:str,query:str)->None:
        r"""
        This function removes studies that contain the query in the given column, field name, from the metagenomics studies database.

        Required Configs:
            - config.metagenomics_studies_db

        Args:
            field_name (str): The name of the column to query.
            query (str): The query string.

        Examples:
            >>> import os,json
            >>> local_dir={'metagenomics_studies':os.path.join(Main_Dir,"test","metagenomics_test_db.json")}
            >>> db=Database(config=configs.Database(studies_local=local_dir))
            >>> metagenomics_study=MetagenomicsStudy(name="test_study",study_type="metagenomics",microbiome="anaerobic digester",sample_accession="test",comments="test",study_accession="test")
            >>> db.add_metagenomics_study_to_metagenomics_studies_db(metagenomics_study)
            >>> assert os.path.exists(local_dir['metagenomics_studies'])==True
            >>> assert pd.read_table(local_dir['metagenomics_studies'],delimiter="\t").shape[0]>0
            >>> db.remove_metagenomics_study_from_metagenomics_studies_db("name","test_study")
            >>> assert pd.read_table(local_dir['metagenomics_studies'],delimiter="\t").shape[0]==0
            >>> os.remove(local_dir['metagenomics_studies'])
        """
        if not os.path.exists(self.config.studies_local["metagenomics_studies"]):
            raise FileNotFoundError("Metagenomics studies database does not exist!")

        metagenomics_studies_db=_read_table(self.config.studies_local["metagenomics_studies"])
        metagenomics_studies_db.filter(~pl.col(field_name).str.contains(query, literal=True)).write_csv(
            self.config.studies_local["metagenomics_studies"], separator="\t"
        )
    
    def get_metagenomics_study_from_metagenomics_studies_db(self,field_name:str,query:str)->list[MetagenomicsStudy]:
        r"""
        This function returns a metagenomics study from the metagenomics studies database. It takes the query string and the column name to query and returns the metagenomics study that contains the query string in the given column.

        Required Configs:
            - config.metagenomics_studies_db
        
        Args:
            field_name (str): The name of the column to query.
            query (str): The query string.
        
        Returns:
            MetagenomicsStudy: An instance of the MetagenomicsStudy class.
        
        Examples:
            >>> import os
            >>> local_dir={'metagenomics_studies':os.path.join(Main_Dir,"test","metagenomics_test_db.tsv")}
            >>> db=Database(config=configs.Database(studies_local=local_dir))
            >>> metagenomics_study=MetagenomicsStudy(name="test_study",study_type="metagenomics",microbiome="anaerobic digester",sample_accession="test",comments="test",study_accession="test")
            >>> db.add_metagenomics_study_to_metagenomics_studies_db(metagenomics_study)
            >>> assert os.path.exists(local_dir['metagenomics_studies'])==True
            >>> assert pd.read_table(local_dir['metagenomics_studies'],delimiter="\t").shape[0]>0
            >>> metagenomics_study=db.get_metagenomics_study_from_metagenomics_studies_db("name","test_study")
            >>> assert metagenomics_study[0].name=="test_study"
            >>> os.remove(local_dir['metagenomics_studies'])
        """
        if not os.path.exists(self.config.studies_local["metagenomics_studies"]):
            raise FileNotFoundError("Metagenomics studies database does not exist!")

        metagenomics_studies_db=_read_table(self.config.studies_local["metagenomics_studies"])
        metagenomics_studies_db=metagenomics_studies_db.filter(pl.col(field_name).str.contains(query, literal=True))
        return [MetagenomicsStudy(**metagenomics_study) for metagenomics_study in metagenomics_studies_db.to_dicts()]
    
    def add_experiment_to_experiments_db(self,experiment:Experiment,force:bool=False)->None:
        r"""
        This function adds an experiment to the experiments database. It takes an experiment and adds it to the experiments database.
        
        Required Configs:
            - config.experimental_data_db
        
        Args:
            experiment (Experiment): An instance of the Experiment class.
        
        Examples:
            >>> import os,json
            >>> local_dir={'experimental_data_db':os.path.join(Main_Dir,"test","experiments_test_db.json")}
            >>> if os.path.exists(local_dir['experimental_data_db']):
            ...     os.remove(local_dir['experimental_data_db'])
            >>> db=Database(config=configs.Database(studies_local=local_dir))
            >>> feed=Feed(name="test_feed",carbohydrates=10,lipids=20,proteins=30,tss=80,si=10,xi=30,reference="test")
            >>> experiment = Experiment(name="test2_study",initial_concentrations={}, time=[0, 1, 2], variables=["S_bu","S_ac"], data=[[1, 2, 3], [4, 5, 6]],feed=feed,reference="test")
            >>> db.add_experiment_to_experiments_db(experiment,force=True)
            >>> assert os.path.exists(local_dir['experimental_data_db'])==True
            >>> assert os.path.getsize(local_dir['experimental_data_db'])>0
            >>> os.remove(local_dir['experimental_data_db'])

        """
        if not os.path.exists(self.config.studies_local["experimental_data_db"]):
            self.initialize_experimental_data_db()
        if force:
            for exp in self.get_experiment_from_experiments_db("name",experiment.name):
                self.remove_experiment_from_experiments_db("name",exp.name)
            
        else:
            if self.get_experiment_from_experiments_db("name",experiment.name):
                raise ValueError("Experiment already exists in the database!")
        
        with open(self.config.studies_local["experimental_data_db"],"r") as f:
            experiments_db=json.load(f)
        experiments_db.append(experiment.to_dict())
        with open(self.config.studies_local["experimental_data_db"],"w") as f:
            json.dump(experiments_db,f)
        
    def remove_experiment_from_experiments_db(self,field_name:str,query:str)->None:
        r"""
        This function removes experiments that contain the query in the given column, field name, from the experiments database.

        Required Configs:
            - config.experimental_data_db

        Args:
            field_name (str): The name of the column to query.
            query (str): The query string.

        Examples:
            >>> import os,json
            >>> local_dir={'experimental_data_db':os.path.join(Main_Dir,"test","experiments_test_db.json")}
            >>> db=Database(config=configs.Database(studies_local=local_dir))
            >>> feed=Feed(name="test_feed",carbohydrates=10,lipids=20,proteins=30,tss=80,si=10,xi=30,reference="test")
            >>> experiment = Experiment(name="test2_study",initial_concentrations={}, time=[0, 1, 2], variables=["S_bu","S_ac"], data=[[1, 2, 3], [4, 5, 6]],feed=feed,reference="test")
            >>> db.add_experiment_to_experiments_db(experiment,force=True)
            >>> assert os.path.exists(local_dir['experimental_data_db'])==True
            >>> assert os.path.getsize(local_dir['experimental_data_db'])>0
            >>> db.remove_experiment_from_experiments_db("name","test2_study")
            >>> assert pd.read_json(local_dir['experimental_data_db']).shape[0]==0
            >>> os.remove(local_dir['experimental_data_db'])

        """
        if not os.path.exists(self.config.studies_local["experimental_data_db"]):
            raise FileNotFoundError("Experimental data database does not exist!")

        with open(self.config.studies_local["experimental_data_db"],"r") as f:
            experiments_db=json.load(f)
        experiments_db=[experiment for experiment in experiments_db if query not in experiment[field_name]]
        with open(self.config.studies_local["experimental_data_db"],"w") as f:
            json.dump(experiments_db,f)

    def get_experiment_from_experiments_db(self,field_name:str,query:str)->list[Experiment]:
        r"""
        This function returns an experiment from the experiments database. It takes the query string and the column name to query and returns the experiment that contains the query string in the given column.

        Required Configs:
            - config.experimental_data_db
        
        Args:
            field_name (str): The name of the column to query.
            query (str): The query string.
        
        Returns:
            Experiment: An instance of the Experiment class.
        
        Examples:
            >>> import os,json
            >>> local_dir={'experimental_data_db':os.path.join(Main_Dir,"test","experiments_test_db.json")}
            >>> db=Database(config=configs.Database(studies_local=local_dir))
            >>> feed=Feed(name="test_feed",carbohydrates=10,lipids=20,proteins=30,tss=80,si=10,xi=30,reference="test")
            >>> experiment = Experiment(name="test2_study",initial_concentrations={}, time=[0, 1, 2], variables=["S_bu","S_ac"], data=[[1, 2, 3], [4, 5, 6]],feed=feed,reference="test")
            >>> db.add_experiment_to_experiments_db(experiment,force=True)
            >>> assert os.path.exists(local_dir['experimental_data_db'])==True
            >>> assert os.path.getsize(local_dir['experimental_data_db'])>0
            >>> experiment=db.get_experiment_from_experiments_db("name","test2_study")
            >>> assert experiment[0].name=="test2_study"
            >>> os.remove(local_dir['experimental_data_db'])
        """


        
        if not os.path.exists(self.config.studies_local["experimental_data_db"]):
            raise FileNotFoundError("Experimental data database does not exist!")

        with open(self.config.studies_local["experimental_data_db"],"r") as f:
            experiments_db=json.load(f)
        experiments_db=[experiment for experiment in experiments_db if query in experiment[field_name]]
        for experiment in experiments_db:
            experiment["feed"]=Feed(**experiment["feed"])
        return [Experiment(**experiment) for experiment in experiments_db]
        
    def build_mmseqs_database(self,container:str="None")->str:
        """Builds an indexed mmseqs database from the ADToolbox's fasta protein database.
        
        Required Configs:
            - config.protein_db
            - config.adtoolbox_singularity
            - config.adtoolbox_docker
            
        Args:
            container (str, optional): The container to run the script with. Defaults to "None".
        Returns:
            str: The script to build the mmseqs database.
            
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==False
            >>> db=Database(config=configs.Database(protein_db=os.path.join(Main_Dir,"protein_test_db.fasta")))
            >>> db.add_protein_to_protein_db("P0A9P0","x,x,x,x")
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==True
            >>> assert os.path.getsize(os.path.join(Main_Dir,"protein_test_db.fasta"))>0
            >>> script=db.build_mmseqs_database()
            >>> assert script=="mmseqs createdb "+str(os.path.join(Main_Dir,"protein_test_db.fasta"))+" "+str(db.config.protein_db_mmseqs)
            >>> os.remove(os.path.join(Main_Dir,"protein_test_db.fasta"))
        
        """
        script=create_mmseqs_database(self.config.protein_db,
                                      self.config.protein_db_mmseqs,
                                      container=container,
                                      run=False,
                                      config=self.config)

    
        return script


    def download_adm_parameters(self,verbose:bool=True)->None:
        """
        Downloads the parameters needed for running ADM models in ADToolbox.
        
        Required Configs:
            - config.adm_parameters_base_dir
            - config.adm_parameters_urls
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"adm_parameters_test"))==False
            >>> E_ADM_LOCAL={"model_parameters":os.path.join(Main_Dir, "test","ADM_Parameters","e_adm_model_parameters.json"),\
	        "base_parameters":os.path.join(Main_Dir, "test","ADM_Parameters","e_adm_base_parameters.json"),\
	        "initial_conditions":os.path.join(Main_Dir, "test","ADM_Parameters","e_adm_initial_conditions.json"),\
	        "inlet_conditions":os.path.join(Main_Dir, "test","ADM_Parameters","e_adm_inlet_conditions.json"),\
	        "reactions":os.path.join(Main_Dir, "test","ADM_Parameters","e_adm_reactions.json"),\
	        "species":os.path.join(Main_Dir, "test","ADM_Parameters","e_adm_species.json")}
            >>> db=Database(config=configs.Database(adm_parameters=E_ADM_LOCAL))
            >>> db.download_adm_parameters(verbose=False) 
            >>> assert os.path.exists(os.path.join(Main_Dir,"test"))==True
            >>> assert len(os.listdir(os.path.join(Main_Dir,"test")))>0
            >>> os.system("rm -r "+os.path.join(Main_Dir,"test"))
            0
        
        Args:
        
            verbose (bool, optional): Whether to print the progress or not. Defaults to True.
        
        
        """
        for param in self.config.adm_parameters.keys():
            if not pathlib.Path(self.config.adm_parameters[param]).parent.exists():
                os.makedirs(pathlib.Path(self.config.adm_parameters[param]).parent)
            r = requests.get(self.config.adm_parameters_urls[param], allow_redirects=True)
            with open(self.config.adm_parameters[param], 'wb') as f:
                f.write(r.content)
            if verbose:
                rich.print(f"[green]{param} downloaded to {self.config.adm_parameters[param]}")
        
    def download_seed_databases(self,verbose:bool=True) -> None :
        """This function will download the seed databases, both compound and reaction databases.

        Required Configs:
            - config.seed_rxn_url
            - config.seed_compound_url
            - config.reaction_db
            - config.compound_db

        Args:
            verbose (bool, optional): Whether to print the progress or not. Defaults to True.
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"seed_compound.json"))==False
            >>> db=Database(config=configs.Database(reaction_db=os.path.join(Main_Dir,"seed_rxn.json"),compound_db=os.path.join(Main_Dir,"seed_compound.json")))
            >>> db.download_seed_databases(verbose=False)
            >>> assert os.path.exists(os.path.join(Main_Dir,"seed_rxn.json"))==True
            >>> assert os.path.exists(os.path.join(Main_Dir,"seed_compound.json"))==True
            >>> os.remove(os.path.join(Main_Dir,"seed_rxn.json"))
            >>> os.remove(os.path.join(Main_Dir,"seed_compound.json"))
        """
        r = requests.get(self.config.seed_rxn_url, allow_redirects=True,stream=True)
        if not os.path.exists(Path(self.config.reaction_db).parent):
            os.makedirs(Path(self.config.reaction_db).parent)
        with open(self.config.reaction_db, 'wb') as f:
            f.write(r.content)
        if verbose:
            rich.print(f"[green]Reaction database downloaded to {self.config.reaction_db}")
        r=requests.get(self.config.seed_compound_url,allow_redirects=True,stream=True)
        with open(self.config.compound_db, 'wb') as f:
            f.write(r.content)
        if verbose:
            rich.print(f"[green]Compound database downloaded to {self.config.compound_db}")

    def download_protein_database(self, verbose:bool=True) -> None:
        """
        Downloads the prebuilt protein database from the remote repository.
        
        Required Configs:
            - config.protein_db_url
            - config.protein_db
        
        Args:
            verbose (bool, optional): Whether to print the progress or not. Defaults to True.
            
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==False
            >>> db=Database(config=configs.Database(protein_db=os.path.join(Main_Dir,"protein_test_db.fasta")))
            >>> db.download_protein_database(verbose=False)
            >>> assert os.path.exists(os.path.join(Main_Dir,"protein_test_db.fasta"))==True
            >>> assert os.path.getsize(os.path.join(Main_Dir,"protein_test_db.fasta"))>0
            >>> os.remove(os.path.join(Main_Dir,"protein_test_db.fasta"))
        """
        r = requests.get(self.config.protein_db_url, allow_redirects=True)
        
        if not os.path.exists(Path(self.config.protein_db).parent):
            os.makedirs(Path(self.config.protein_db).parent)
        
        with open(self.config.protein_db, 'wb') as f:
            f.write(r.content)
        if verbose:
            rich.print(f"[green]Protein database downloaded to {self.config.protein_db}")
        
    def download_reaction_database(self,verbose:bool=True)->None:
        """
        This function will download the reaction database from the remote repository.
        
        Required Configs:
            - config.adtoolbox_rxn_db_url
            - config.csv_reaction_db
        
        Args:
            verbose (bool, optional): Whether to print the progress or not. Defaults to True.

        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"reaction_test_db.csv"))==False
            >>> db=Database(config=configs.Database(csv_reaction_db=os.path.join(Main_Dir,"reaction_test_db.csv")))
            >>> db.download_reaction_database(verbose=False)
            >>> assert os.path.exists(os.path.join(Main_Dir,"reaction_test_db.csv"))==True
            >>> assert os.path.getsize(os.path.join(Main_Dir,"reaction_test_db.csv"))>0
            >>> os.remove(os.path.join(Main_Dir,"reaction_test_db.csv"))
        """
    
        r = requests.get(self.config.adtoolbox_rxn_db_url, allow_redirects=True)
        
        if not os.path.exists(Path(self.config.csv_reaction_db).parent):
            os.makedirs(Path(self.config.csv_reaction_db).parent)

        with open(self.config.csv_reaction_db, 'wb') as f:
            f.write(r.content)
        if verbose:
            rich.print(f"[green]Reaction database downloaded to {self.config.csv_reaction_db}")

    
    def download_feed_database(self,verbose:bool=True)-> None:
        """
        This function will download the feed database from the remote repository.

        Required Configs:
            - config.feed_db_url
            - config.feed_db
        
        Args:
            verbose (bool, optional): Whether to print the progress or not. Defaults to True.

        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==False
            >>> db=Database(config=configs.Database(feed_db=os.path.join(Main_Dir,"feed_test_db.tsv")))
            >>> db.download_feed_database(verbose=False)
            >>> assert os.path.exists(os.path.join(Main_Dir,"feed_test_db.tsv"))==True
            >>> assert os.path.getsize(os.path.join(Main_Dir,"feed_test_db.tsv"))>0
            >>> os.remove(os.path.join(Main_Dir,"feed_test_db.tsv"))
        """
        r = requests.get(self.config.feed_db_url, allow_redirects=True)
        
        if not os.path.exists(Path(self.config.feed_db).parent):
            os.makedirs(Path(self.config.feed_db).parent)
        
        with open(self.config.feed_db, 'wb') as f:
            f.write(r.content)
        if verbose:
            rich.print(f"[green]Feed database downloaded to {self.config.feed_db}")
    
    def download_studies_database(self,verbose:bool=True)->None:
        """
        This function will download the required files for studies functionality.

        Args:
            verbose (bool, optional): Whether to print the progress or not. Defaults to True.
        
        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"studies_test_db.tsv"))==False
            >>> STUDIES_LOCAL={"metagenomics_studies":os.path.join(Main_Dir,"Database","test","Studies","metagenomics_studies.tsv"),\
	                            "experimental_data_db":os.path.join(Main_Dir,"Database","test","Studies","experimental_data_references.json")}
            >>> db=Database(config=configs.Database(studies_local=STUDIES_LOCAL))
            >>> db.download_studies_database(verbose=False)
            >>> assert os.path.exists(STUDIES_LOCAL['metagenomics_studies'])==True
            >>> assert os.path.getsize(STUDIES_LOCAL['metagenomics_studies'])>0
            >>> os.remove(STUDIES_LOCAL['metagenomics_studies'])
        """
        for i in ["metagenomics_studies","experimental_data_db"]:
            r = requests.get(self.config.studies_remote[i], allow_redirects=True)
            if not os.path.exists(Path(self.config.studies_local[i]).parent):
                os.makedirs(Path(self.config.studies_local[i]).parent)
            with open(self.config.studies_local[i], 'wb') as f:
                f.write(r.content)
            
            if verbose:
                rich.print(f"[bold green]Downloaded {self.config.studies_remote[i]}[/bold green]")
    
    def download_amplicon_to_genome_db(self,verbose:bool=True):
        """
        This function will automatically download the GTDB-tk database for genome assignment.
        
        Required Configs:
            - config.amplicon_to_genome_db
            - config.amplicon_to_genome_urls
        
        Args:
            verbose (bool, optional): Whether to print the progress or not. Defaults to True.

        Examples:
            >>> import os
            >>> assert os.path.exists(os.path.join(Main_Dir,"amplicon_to_genome_test_db"))==False
            >>> db=Database(config=configs.Database(amplicon_to_genome_db=os.path.join(Main_Dir,"amplicon_to_genome_test_db")))
            >>> db.download_amplicon_to_genome_db(verbose=False)
            >>> assert os.path.exists(os.path.join(Main_Dir,"amplicon_to_genome_test_db"))==True
            >>> assert len(os.listdir(os.path.join(Main_Dir,"amplicon_to_genome_test_db")))>0
            >>> os.system("rm -r "+os.path.join(Main_Dir,"amplicon_to_genome_test_db"))
            0
        """
        if not os.path.exists(self.config.amplicon_to_genome_db):
            os.mkdir(self.config.amplicon_to_genome_db)

        url = self.config.amplicon_to_genome_urls
        if verbose:
            for keys in ['Version', 'MD5SUM', 'FILE_DESCRIPTIONS']:
                with requests.get(url[keys], allow_redirects=True, stream=True) as r:
                    total_size = int(r.headers.get('content-length', 0))
                    block_size = 1024
                    with Progress() as progress:
                        task1 = progress.add_task("Downloading " + keys, total=total_size)
                        with open(os.path.join(self.config.amplicon_to_genome_db, keys), 'wb') as f:
                            for data in r.iter_content(block_size):
                                progress.update(task1, advance=len(data))
                                f.write(data)
            with requests.get(url['metadata_field_desc'], allow_redirects=True, stream=True) as r:
                total_size = int(r.headers.get('content-length', 0))
                block_size = 1024
                with Progress() as progress:
                    task1 = progress.add_task("Downloading metadata_field_desc.tsv", total=total_size)
                    with open(os.path.join(self.config.amplicon_to_genome_db, 'metadata_field_desc.tsv'), 'wb') as f:
                        for data in r.iter_content(block_size):
                            progress.update(task1, advance=len(data))
                            f.write(data)

            for keys in ['bac120_ssu']:
                with requests.get(url[keys], allow_redirects=True, stream=True) as r:
                    total_size = int(r.headers.get('content-length', 0))
                    block_size = 1024
                    with Progress() as progress:
                        task1 = progress.add_task("Downloading " + keys, total=total_size)
                        with open(os.path.join(self.config.amplicon_to_genome_db, url[keys].split("/")[-1]), 'wb') as f:
                            for data in r.iter_content(block_size):
                                progress.update(task1, advance=len(data))
                                f.write(data)
                with gzip.open(os.path.join(self.config.amplicon_to_genome_db, url[keys].split("/")[-1]),"r") as f_in:
                    with open(os.path.join(self.config.amplicon_to_genome_db, url[keys].split("/")[-1].replace(".gz","")),"wb") as f_out:
                        f_out.write(f_in.read())

                
                os.remove(os.path.join(self.config.amplicon_to_genome_db, url[keys].split("/")[-1]))
        else:
            for keys in ['Version', 'MD5SUM', 'FILE_DESCRIPTIONS']:
                with requests.get(url[keys], allow_redirects=True, stream=False) as r:
                    with open(os.path.join(self.config.amplicon_to_genome_db, keys), 'wb') as f:
                        f.write(r.content)
            with requests.get(url['metadata_field_desc'], allow_redirects=True, stream=False) as r:
                with open(os.path.join(self.config.amplicon_to_genome_db, 'metadata_field_desc.tsv'), 'wb') as f:
                    f.write(r.content)
            for keys in [ 'bac120_ssu']:
                with requests.get(url[keys], allow_redirects=True, stream=False) as r:
                    with open(os.path.join(self.config.amplicon_to_genome_db, url[keys].split("/")[-1]), 'wb') as f:
                        f.write(r.content)
                with gzip.open(os.path.join(self.config.amplicon_to_genome_db, url[keys].split("/")[-1]),"r") as f_in:
                    with open(os.path.join(self.config.amplicon_to_genome_db, url[keys].split("/")[-1].replace(".gz","")),"wb") as f_out:
                        f_out.write(f_in.read())
        if verbose:
            rich.print("[bold green]Downloaded all the required files for Amplicon to Genome functionality.[/bold green]")
                    
                        

            
    def download_all_databases(self,verbose:bool=True)->None:
        """
        This function will download all the required databases for all the functionalities of ADToolbox.
        NOTE: each method that this function calls is individually tested so it is skipped from testing!

        Args:
            verbose (bool, optional): Whether to print the progress or not. Defaults to True.

        Required Configs:
            - config.adm_parameters_base_dir
            - config.adm_parameters_urls
            - config.seed_rxn_url
            - config.seed_compound_url
            - config.reaction_db
            - config.compound_db
            - config.protein_db_url
            - config.protein_db
            - config.adtoolbox_rxn_db_url
            - config.csv_reaction_db
            - config.feed_db_url
            - config.feed_db
            - config.amplicon_to_genome_db
            - config.amplicon_to_genome_urls
            - config.studies_db
            - config.studies_urls
            
        Examples:
            >>> import os # doctest: +SKIP
            >>> db=Database(config=configs.Database()) # doctest: +SKIP
            >>> db.download_all_databases(verbose=False) # doctest: +SKIP

        """

        self.download_seed_databases(verbose=verbose)
        self.download_adm_parameters(verbose=verbose)
        self.download_protein_database(verbose=verbose)
        self.download_reaction_database(verbose=verbose)
        self.download_feed_database(verbose=verbose)
        self.download_studies_database(verbose=verbose)
        self.download_amplicon_to_genome_db(verbose=verbose)
        

class Metagenomics:

    """
    This is the main class for Metagenomics functionality of ADToolbox. This class contains all the methods required for metagenomics analysis 
    that ADToolbox offers.
    """
    def __init__(self,config:configs.Metagenomics)->None:
        """In order to instntiate an object from this class, you need to provide a metagenomics configs object from the configs module : configs.Metagenomics.
        Information for inputs and of each method is then obtained from the corresponding configs object. The following example shows how to instantiate an object from this class
        using the default configs object:
        
        Examples:
            >>> from adtoolbox import core, configs
            >>> config=configs.Metagenomics() ### This uses default arguments. Refer to configs module for more information.
            >>> metagenomics=core.Metagenomics(config)
            >>> assert type(metagenomics)==core.Metagenomics
        
        Args:
            config (configs.Metagenomics): A metagenomics configs object from configs module.
        
        Returns:
            None
        """
        self.config=config

    #### NEEDS to BE FIXED        
    def find_top_taxa(
        self,
        sample_name:str,
        treshold:Union[int,float],
        mode:str='top_k',
        )->dict:
        """
        This function needs three amplicon outputs:
        1. feature table: This is the abundance of each feature in each sample (TSV).
        2. taxonomy table: This is the taxonomy of each feature (TSV). 
        3. rep seqs: This is the representative sequence of each feature (fasta).
        It then finds the top k features or features that form specific percentile of the community of the sample.
        
        Required Configs:
        
            config.feature_table_dir: The path to the feature table tsv file.
            ---------
            config.taxonomy_table_dir: The path to the taxonomy table tsv file.
            ---------
            config.rep_seq_fasta: The path to the representative sequence fasta file.
            ---------
        
        Args:
            sample_name (str): The name of the sample.
            treshold (int, float): The threshold for the top k or the percentile.
            mode (str, optional): Whether to find the top k features or features that form specific percentile of the community of the sample. Defaults to 'top_k'. Options: 'top_k', 'percentile'.
        
        Returns:
            dict: A dictionary of the top k features and their taxonomy.
        """
        ### Load all the required files
        feature_table = pl.read_csv(self.config.feature_table_dir, separator="\t", skip_rows=1, infer_schema_length=0)
        taxonomy_table = pl.read_csv(self.config.taxonomy_table_dir, separator="\t", infer_schema_length=0)
        repseqs=fasta_to_dict(self.config.rep_seq_fasta)
        ### End Loading
        if mode == 'top_k':
            sorted_df=feature_table.with_columns(pl.col(sample_name).cast(pl.Float64, strict=False).fill_null(0.0)).sort(sample_name, descending=True)
            top_featureids=sorted_df['#OTU ID'].head(treshold).to_list()
            top_taxa=[
                taxonomy_table.filter(pl.col('Feature ID') == featureid)['Taxon'][0]
                for featureid in top_featureids
            ]
            top_repseqs=[repseqs[featureid] for featureid in top_featureids]
            total = sorted_df[sample_name].sum()
            top_abundances=[float(value) / float(total) for value in sorted_df[sample_name].head(treshold).to_list()]
            
        elif mode == 'percentile':
            total = feature_table.select(pl.col(sample_name).cast(pl.Float64, strict=False).fill_null(0.0).sum()).item()
            sorted_df=(
                feature_table
                .with_columns((pl.col(sample_name).cast(pl.Float64, strict=False).fill_null(0.0) / total).alias(sample_name))
                .sort(sample_name, descending=True)
                .with_columns((pl.col(sample_name).cum_sum() * 100).alias("cumsum"))
            )
            sorted_df_filtered=sorted_df.filter(pl.col("cumsum") <= treshold)
            top_featureids=sorted_df_filtered['#OTU ID'].to_list()
            top_taxa=[
                taxonomy_table.filter(pl.col('Feature ID') == featureid)['Taxon'][0]
                for featureid in top_featureids
            ]
            top_repseqs=[repseqs[featureid] for featureid in top_featureids]
            top_abundances=sorted_df_filtered[sample_name].to_list()
        else:
            raise ValueError("mode must be either 'top_k' or 'percentile'")
        
        return {'top_featureids':top_featureids,'top_taxa':top_taxa,'top_repseqs':top_repseqs,'top_abundances':top_abundances}    
        
    
    def align_to_gtdb(self,
                      query_dir:str,
                      output_dir:str,
                      container:str="None",
                      image: str | None = None)->tuple[str]:
        r"""This function takes the representative sequences of the top k features and generates the script to
        align these feature sequences to gtdb using VSEARCH. If you intend to run this you either
        need to have VSEARCH installed or run it with a container option. You can use either the docker or singularity
        as container options. Otherwise you can use None and run it with the assumption that VSEARCH is installed.
        If you only want the script and not to run it, set run to False.

        Required Configs:
        
            ---------
            config.gtdb_dir_fasta: The path to the gtdb fasta database.
            ---------
            config.vsearch_similarity: The similarity threshold for the alignment to be used by VSEARCH.
            ---------
            config.vsearch_threads: The number of threads to be used by VSEARCH.
            ---------
            config.adtoolbox_docker: The name of the docker image to be used by ADToolbox (Only if using Docker as container).
            ---------
            config.adtoolbox_singularity: The name of the singularity image to be used by ADToolbox (Only if using Singularity as container).
            ---------
        Examples:
            >>> import os
            >>> query_dir=os.path.join(Main_Dir,"test","query.fa")
            >>> output_dir=os.path.join(Main_Dir,"test")
            >>> conf=configs.Metagenomics(vsearch_similarity=0.8)
            >>> conf.gtdb_dir_fasta=(os.path.join(Main_Dir,"db.fa"))
            >>> obj = Metagenomics(conf) 
            >>> assert obj.align_to_gtdb(query_dir, output_dir, container="docker")[0] == f'docker run -v {output_dir}:{output_dir} -v {conf.gtdb_dir_fasta}:{conf.gtdb_dir_fasta} -v {query_dir}:{query_dir} parsaghadermazi/adtoolbox:x86 vsearch --top_hits_only --blast6out {output_dir}/matches.blast --usearch_global {query_dir} --db {conf.gtdb_dir_fasta} --id {conf.vsearch_similarity} --threads 4 --alnout {output_dir}/Alignments --top_hits_only\n'

        Args:
            container (str, optional): The container to use. Defaults to "None".
        
        Returns:
            str: The script that is supposed to be running later.
        """
        ### Load all the required files
        alignment_dir = str(pathlib.Path(os.path.join(output_dir,'Alignments')).absolute())
        match_table=str(pathlib.Path(os.path.join(output_dir,'matches.blast')))
        if self.config.gtdb_dir_fasta is None:
            raise FileNotFoundError(
                "No GTDB/amplicon-to-genome FASTA was found. "
                f"Looked under {self.config.amplicon2genome_db!r} for pattern {self.config.gtdb_dir!r}. "
                "Pass --amplicon-to-genome-db to a directory containing an SSU FASTA, "
                "or run `adtoolbox Database download-amplicon-to-genome-dbs` first."
            )
        gtdb_dir_fasta=str(pathlib.Path(self.config.gtdb_dir_fasta))
        ### End Loading
        query=query_dir
        dirs=[output_dir,
            gtdb_dir_fasta,
            query
            ]
        for dir in dirs:
            if not pathlib.Path(dir).exists():
                os.mkdir(dir)
        if container=="None":
            bash_script=('vsearch --top_hits_only --blast6out '+
                        match_table+
                        ' --usearch_global '+ query +
                        ' --db '+ gtdb_dir_fasta +
                        ' --id ' +str(self.config.vsearch_similarity) +
                        ' --threads '+str(self.config.vsearch_threads)+
                        ' --alnout '+ alignment_dir +
                        ' --top_hits_only'+'\n')
        
        if container=="docker":
            bash_script='docker run '
            for dir in dirs:
                bash_script+=('-v '+dir+':'+dir+' ')
            
            bash_script += (self._container_image(container, image)+' vsearch --top_hits_only --blast6out '+
                        match_table+
                        ' --usearch_global '+ query +
                        ' --db '+ gtdb_dir_fasta +
                        ' --id ' +str(self.config.vsearch_similarity) +
                        ' --threads '+str(self.config.vsearch_threads)+
                        ' --alnout '+ alignment_dir +
                        ' --top_hits_only'+'\n')
        
        if container in {"singularity", "apptainer"}:
            runtime = "apptainer" if container == "apptainer" else "singularity"
            bash_script=f'{runtime} exec '
            for dir in dirs:
                bash_script+=('-B '+str(dir)+':'+str(dir)+' ')
            
            bash_script += (self._container_image(container, image)+' vsearch --top_hits_only --blast6out '+
                        match_table+
                        ' --usearch_global '+ str(query) +
                        ' --db '+ gtdb_dir_fasta +
                        ' --id ' +str(self.config.vsearch_similarity) +
                        ' --threads '+str(self.config.vsearch_threads)+
                        ' --alnout '+ alignment_dir +
                        ' --top_hits_only'+'\n')
        return bash_script,
    
    
    
    def get_genomes_from_gtdb_alignment(self,alignment_dir:str)->dict:
        r"""This function takes the alignment file generated from the align_to_gtdb function and generates the the genome information
        using the GTDB-Tk. In the outputted dictionary, the keys are feature ids and the values are the representative genomes.

        Examples:
            >>> import os
            >>> alignments=["0034e3d0368ec41aeb7f346b434f8d46","GB_GCA_937889405.1~CALAPC010000121.1","100.0","253","0","0"	,"1","253",	"1","1009",	"-1","0"]
            >>> hit="\t".join(alignments)
            >>> output= os.path.join(Main_Dir,"test","matches.blast")
            >>> with open(output,"w") as f:
            ...    f.write(hit)
            101
            >>> obj=Metagenomics(configs.Metagenomics())
            >>> obj.get_genomes_from_gtdb_alignment(output)
            {'0034e3d0368ec41aeb7f346b434f8d46': 'GCA_937889405.1'}
            
    
        Args:
            alignment_dir (str): The path to the alignment file generated by align_to_gtdb.
        """
        aligned = pl.read_csv(
            alignment_dir,
            separator="\t",
            has_header=False,
            infer_schema_length=0,
        )
        if aligned.is_empty():
            return {}
        query_col, target_col = aligned.columns[:2]
        aligned = (
            aligned
            .unique(subset=[query_col], keep="first")
            .with_columns(
                pl.col(target_col)
                .map_elements(lambda value: ("_".join(str(value).split("_")[1:])).split("~")[0], return_dtype=pl.Utf8)
                .alias("genome_id")
            )
        )
        return dict(zip(aligned[query_col].to_list(), aligned["genome_id"].to_list()))

    @staticmethod
    def _ncbi_genome_path(identifier: str) -> tuple[str, str]:
        accession = identifier.strip()
        match = re.match(r"^(GC[AF])_?(\d+)(?:\.(\d+))?$", accession)
        if not match:
            raise ValueError(f"Genome accession must look like GCF_000146505.1 or GCA_937889405.1: {identifier}")
        prefix, digits, version = match.groups()
        normalized = f"{prefix}_{digits}.{version}" if version else f"{prefix}_{digits}"
        chunks = "/".join(digits[i:i + 3] for i in range(0, len(digits), 3))
        return normalized, f"{prefix}/{chunks}"
    
    
    def download_genome(self,identifier:str,output_dir:str,container:str="None")-> str:
        r"""This function downloads the genomes from NCBI using the refseq/genbank identifiers.
        Note that this function uses rsync to download the genomes. 

        Required Configs:
            config.genomes_base_dir: The path to the base directory where the genomes will be saved.
            ---------
            config.adtoolbox_docker: The name of the docker image to be used by ADToolbox (Only if using Docker as container).
            ---------
            config.adtoolbox_singularity: The name of the singularity image to be used by ADToolbox (Only if using Singularity as container).
            ---------
        Examples:
            >>> import os 
            >>> genome_identifier="GCA_937889405.1"
            >>> output=os.path.join(Main_Dir,"test")
            >>> obj = Metagenomics(configs.Metagenomics()) 
            >>> assert obj.download_genome(identifier=genome_identifier,output_dir= output)[0] == 'rsync -avz --progress rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/937/889/405 '+output

        Args:
            identifier (str): The identifier for the genome. It can be either refseq or genbank.
            output_dir (str): The directory where the genome should be downloaded.
            container (str, optional): The container to use. Defaults to "None". You may select from "None", "docker", "singularity".
        
        Returns:
            str: The bash script that is used to download the genomes or to be used to download the genomes.

        """
        accession, accession_path = self._ncbi_genome_path(identifier)
        genome_dir = pathlib.Path(output_dir)
        base_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{accession_path}"
        script_body = f"""set -euo pipefail
mkdir -p {shlex.quote(str(genome_dir))}
base_url={shlex.quote(base_url)}
accession={shlex.quote(accession)}
out_dir={shlex.quote(str(genome_dir))}
assembly_dir=$(curl -fsSL "$base_url/" | sed -n 's/.*href="\\([^"]*'${{accession}}'[^"]*\\/\\)".*/\\1/p' | head -n 1)
if [ -z "$assembly_dir" ]; then
  echo "Could not find assembly directory for $accession under $base_url" >&2
  exit 1
fi
assembly_name=${{assembly_dir%/}}
mkdir -p "$out_dir/$assembly_name"
genome_file=$(curl -fsSL "$base_url/$assembly_dir" | sed -n 's/.*href="\\([^"]*_genomic\\.fna\\.gz\\)".*/\\1/p' | head -n 1)
if [ -z "$genome_file" ]; then
  echo "Could not find genomic FASTA for $accession under $base_url/$assembly_dir" >&2
  exit 1
fi
curl -fL "$base_url/$assembly_dir$genome_file" -o "$out_dir/$assembly_name/$genome_file"
"""

        if container=="None":
            bash_script = script_body
        elif container=="docker":
            bash_script = (
                f"docker run -v {shlex.quote(str(genome_dir))}:{shlex.quote(str(genome_dir))} "
                f"{self.config.adtoolbox_docker} bash -lc {shlex.quote(script_body)}"
            )
        elif container=="singularity":
            bash_script = (
                f"singularity exec -B {shlex.quote(str(genome_dir))}:{shlex.quote(str(genome_dir))} "
                f"{self.config.adtoolbox_singularity} bash -lc {shlex.quote(script_body)}"
            )
        else:
            raise ValueError("container must be one of: None, docker, singularity")
        
        return bash_script,
    
    def async_genome_downloader(self,identifiers:Iterable[str],batch_size:float=10,container:str="None"):
        sem=asyncio.Semaphore(batch_size)
        asyncio.run(self._collect_coros(identifiers=identifiers,semaphore=sem,container=container))
        
    async def _collect_coros(self,identifiers:Iterable[str],semaphore:asyncio.Semaphore,container:str="None"):
        await asyncio.gather(*[self._genome_dl_coro(identifier=i,semaphore=semaphore,container=container) for i in identifiers])
        
    async def _genome_dl_coro(self,identifier:str,semaphore:asyncio.Semaphore,container:str="None")->None:
        async with semaphore:
            await asyncio.create_subprocess_exec(*self.download_genome(identifier=identifier,container=container).split(" "))
    
    def extract_genome_info(self,
                            base_dir:str,
                            endpattern:str="genomic.fna.gz",
                            filters:dict={
                                          "INCLUDE":[],
                                          "EXCLUDE":["cds","rna"],
                                            })->dict[str,str]:
        """This function extracts the genome information from the genomes base directory. The output
        is a dictionary where the keys are the genome IDs and the values are the paths to the genome files.
        
        Required Configs:
            None
        ---------
        Args:
            base_dir (str): The path to the base directory where the genomes are saved.
            endpattern (str, optional): The end pattern of the genome files. Defaults to "genomic.fna.gz".
            filters (dict, optional): The filters to be applied to the genome files. This filter must be a 
            dictionary with two keys: INCLUDE and EXCLUDE. The values of these keys must be lists of strings.
            Defaults to {"INCLUDE":[],"EXCLUDE":["cds","rna"]}. This defult is compatible with the genomes downloaded
            from NCBI i.e. only change this if you are providing your own genomes with different file name conventions.
        Returns:
            dict[str,str]: A dictionary containing the address of the genomes that are downloaded or to be downloaded.
        """
        genome_df = self.extract_genome_info_df(base_dir, endpattern=endpattern, filters=filters)
        return dict(zip(genome_df["genome_id"].to_list(), genome_df["path"].to_list()))

    def extract_genome_info_df(
            self,
            base_dir: str,
            endpattern: str = "genomic.fna.gz",
            filters: dict | None = None,
            ) -> pl.DataFrame:
        filters = filters or {"INCLUDE": [], "EXCLUDE": ["cds", "rna", "protein"]}
        base_path = pathlib.Path(base_dir)
        records = []
        if not base_path.exists():
            return pl.DataFrame(schema={
                "genome_id": pl.Utf8,
                "assembly_accession": pl.Utf8,
                "assembly_name": pl.Utf8,
                "path": pl.Utf8,
            })
        for candidate in sorted(base_path.rglob(f"*{endpattern}")):
            name = candidate.name
            if not all(text in name for text in filters.get("INCLUDE", [])):
                continue
            if any(text in name for text in filters.get("EXCLUDE", [])):
                continue
            genome_id = name.replace("_genomic.fna.gz", "").replace("_genomic.fna", "")
            accession_match = re.match(r"^(GC[AF]_\d+\.\d+)", genome_id)
            assembly_accession = accession_match.group(1) if accession_match else genome_id
            records.append({
                "genome_id": genome_id,
                "assembly_accession": assembly_accession,
                "assembly_name": candidate.parent.name,
                "path": str(candidate.resolve()),
            })
        return pl.DataFrame(records, schema={
            "genome_id": pl.Utf8,
            "assembly_accession": pl.Utf8,
            "assembly_name": pl.Utf8,
            "path": pl.Utf8,
        })
     
    def align_genome_to_protein_db(
            self,
            address:str,
            outdir:str,
            name:str,
            container:str="None",
            image: str | None = None,
            )->tuple[str,str]:
        r"""
        This is a function that will align a genome to the Protein Database of the ADToolbox using mmseqs2.
        If you want to save the scripts, set save to True. Note that the alignment tables will be saved in any case.
        Note that this function uses mmseqs2 to align the genomes to the protein database. So, to run this function without
        any container you need to have mmseqs2 installed on your system. However, if you want to run this function with a container,
        you need to have the container installed on your system. You may select from "None", "docker", "singularity".

        Requires:
            config.protein_db: The address of the protein database of the ADToolbox.
            ---------
            config.adtoolbox_docker: The name of the docker image to be used by ADToolbox (Only if using Docker as container).
            ---------
            config.adtoolbox_singularity: The name of the singularity image to be used by ADToolbox (Only if using Singularity as container).
            ---------
        Example: 
            >>> import os
            >>> output=os.path.join(Main_Dir,"test")
            >>> address=os.path.join(Main_Dir,"test","genome_sequence.fa")
            >>> obj = Metagenomics(configs.Metagenomics()) 
            >>> name="test_align"
            >>> assert obj.align_genome_to_protein_db(address=address,outdir=output,name=name,container="docker")[0].split(" ")==['docker', 'run', '', '-v', address+':'+address, '-v', configs.Database().protein_db+':'+configs.Database().protein_db, '-v', output+':'+output, 'parsaghadermazi/adtoolbox:x86', '', 'mmseqs', 'easy-search', address, configs.Database().protein_db, output+'/Alignment_Results_mmseq_test_align.tsv', 'tmpfiles', '--format-mode', '4', '\n\n']

        Args:
            address (str): The address of the genome to be aligned.
            outdir (str): The output directory where the alignment files will be saved.
            name (str): The name of the genome.
            container (str, optional): The container to use. Defaults to "None". You may select from "None", "docker", "singularity".

        Returns:
            str: A dictionary containing the alignment files.
            str: The bash script that is used to align the genomes or to be used to align the genomes.
        """
 
        if container=="None":
            bash_script = ""
            alignment_file=os.path.join(outdir,"Alignment_Results_mmseq_"+name+".tsv")
            bash_script += "mmseqs easy-search " + \
                address + " " + \
                self.config.protein_db + " " + \
                alignment_file+ ' tmp --format-mode 4 '+"\n\n"
        
        if container=="docker":
            bash_script = ""
            alignment_file=os.path.join(outdir,"Alignment_Results_mmseq_"+name+".tsv")
            bash_script +="docker run "+ \
            " -v "+address+":"+address+ \
            " -v "+self.config.protein_db+":"+self.config.protein_db+ \
            " -v "+outdir+":"+outdir+ \
            f" {self._container_image(container, image)}  mmseqs easy-search " + \
                address + " " + \
                self.config.protein_db + " " + \
                alignment_file+' tmpfiles --format-mode 4 '+"\n\n"

        if container in {"singularity", "apptainer"}:
            runtime = "apptainer" if container == "apptainer" else "singularity"
            bash_script = ""
            alignment_file=os.path.join(outdir,"Alignment_Results_mmseq_"+name+".tsv")
            bash_script +=f"{runtime} exec "+ \
            " -B "+address+":"+address+ \
            " -B "+self.config.protein_db+":"+self.config.protein_db+ \
            " -B "+outdir+":"+outdir+ \
            f" {self._container_image(container, image)}  mmseqs easy-search " + \
                address + " " + \
                self.config.protein_db + " " + \
                alignment_file+' tmpfiles --format-mode 4 '+"\n\n"
        
        return  bash_script,alignment_file

    def align_short_reads_to_protein_db(self,
                                        query_seq:str,
                                        alignment_file_name:str,
                                        container:str="None",
                                        image: str | None = None,
                                        )->tuple[str,str]:
        r"""This function aligns shotgun short reads to the protein database of the ADToolbox using mmseqs2.
        mmseqs wrappers in utils are used to perform this task. The result of this task is an alignment table.
        
        Required Configs:
        
            protein_db_mmseqs (str): The address of the existing/to be created protein database of the ADToolbox for mmseqs.
            --------
        Example:
            >>> import os
            >>> query_seq= os.path.join(Main_Dir,"test","query.seq")
            >>> alignment_file= "output_alignment.txt"
            >>> protein_db_mmseqs=os.path.join(Main_Dir,"test","protein_db_mmseqs")
            >>> with open(protein_db_mmseqs,'w') as f:
            ...     f.write("")
            0
            >>> obj = Metagenomics(configs.Metagenomics(protein_db_mmseqs=protein_db_mmseqs)) 
            >>> assert obj.align_short_reads_to_protein_db(query_seq=query_seq, alignment_file_name= protein_db_mmseqs, container="docker")[0] == 'docker run -v /home/parsa/ADresearch/test/query.seq:/home/parsa/ADresearch/test/query.seq -v /home/parsa/ADresearch/test:/home/parsa/ADresearch/test parsaghadermazi/adtoolbox:x86 mmseqs createdb /home/parsa/ADresearch/test/query.seq /home/parsa/ADresearch/test/query\ndocker run -v /home/parsa/ADresearch/test/query:/home/parsa/ADresearch/test/query -v /home/parsa/ADresearch/test:/home/parsa/ADresearch/test -v /home/parsa/ADresearch/test/protein_db_mmseqs:/home/parsa/ADresearch/test/protein_db_mmseqs parsaghadermazi/adtoolbox:x86 mmseqs search /home/parsa/ADresearch/test/query /home/parsa/ADresearch/test/protein_db_mmseqs /home/parsa/ADresearch/test/protein_db_mmseqs /home/parsa/ADresearch/test/tmp\ndocker run -v /home/parsa/ADresearch/test/query:/home/parsa/ADresearch/test/query -v /home/parsa/ADresearch/test:/home/parsa/ADresearch/test -v /home/parsa/ADresearch/test/protein_db_mmseqs:/home/parsa/ADresearch/test/protein_db_mmseqs parsaghadermazi/adtoolbox:x86 mmseqs convertalis /home/parsa/ADresearch/test/query /home/parsa/ADresearch/test/protein_db_mmseqs /home/parsa/ADresearch/test/protein_db_mmseqs /home/parsa/ADresearch/test/protein_db_mmseqs.tsv --format-mode 4\n'
            >>> os.remove(protein_db_mmseqs)

        Args:
            query_seq (str): The address of the query sequence.
            alignment_file_name (str): The name of the alignment file.
            container (str, optional): The container to use. Defaults to "None". You may select from "None", "docker", "singularity".


        Returns:
            str: The bash script that is used to align the genomes or to be used to align the genomes.
            str: The address of the alignment file.
        """
        if not pathlib.Path(self.config.protein_db_mmseqs).exists():
            raise FileNotFoundError("""The protein database of the ADToolbox for mmseqs is not found. Please build it first
                                    using Database.build_mmseqs_database method.""")
        path_query=pathlib.Path(query_seq)
        old_image_values = {}
        if image:
            old_image_values = {
                "adtoolbox_docker": self.config.adtoolbox_docker,
                "adtoolbox_singularity": self.config.adtoolbox_singularity,
            }
            if container == "docker":
                self.config.adtoolbox_docker = image
            elif container in {"singularity", "apptainer"}:
                self.config.adtoolbox_singularity = image
        script = ""
        try:
            script += create_mmseqs_database(query_seq,str(path_query.parent/path_query.name.split(".")[0]),container=container,save=None,run=False,config=self.config)+"\n"
            script += mmseqs_search(
                query_db=str(path_query.parent/path_query.name.split(".")[0]),
                target_db=self.config.protein_db_mmseqs,
                results_db=path_query.parent/alignment_file_name,
                run=False,
                save=None,
                container=container,
                config=self.config,
            )+"\n"
            script += mmseqs_result_db_to_tsv(
                query_db=str(path_query.parent/path_query.name.split(".")[0]),
                target_db=self.config.protein_db_mmseqs,
                results_db=path_query.parent/alignment_file_name,
                tsv_file=path_query.parent/(alignment_file_name+".tsv"),
                container=container,
                save=None,
                run=False,
                config=self.config,
            )+"\n"
        finally:
            for key, value in old_image_values.items():
                setattr(self.config, key, value)
        return script,path_query.parent/(alignment_file_name+".tsv")
    
    def extract_ec_from_alignment(self,alignment_file:str)->dict[str,int]:
        r"""
        This function extracts the number of times an EC number is found in the alignment file when aligned to ADToolbox protein database.
        
        Required Configs:
            config.e_value: The e-value threshold for the filtering the alignment table.
            ---------
            config.bit_score: The bit score threshold for the filtering the alignment table.
            ---------
            config.ec_counts_from_alignment: The address of the json file that the results will be saved in.
            ---------

        Example:
            >>> import os
            >>> output = os.path.join(Main_Dir, "test", "mmseqs_alignments.tsv")
            >>> alignments = ["CP001673.1", "Q8YNF9|1.4.4.2", "0.566", "2859", "414", "0", "1021521", "1024379", "27", "982", "0.000E+00", "1109"]
            >>> headers = ["query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend ", "tstart", "tend", "evalue", "bits"]
            >>> hit = "\t".join(alignments)
            >>> headers_tab = "\t".join(headers)
            >>> combine = headers_tab + "\n" + hit
            >>> with open(output, "w") as f:
            ...     f.write(combine)
            161
            >>> obj = Metagenomics(configs.Metagenomics())
            >>> obj.extract_ec_from_alignment(output) 
            {'1.4.4.2': 1}
    
        Args:
            alignment_file (str): The address of the alignment file.
        
        Returns:
            dict: A dictionary of EC numbers and their counts.

        """
        df=(pl.scan_csv(alignment_file, separator='\t',) 
        .filter((pl.col("evalue") < self.config.e_value)&(pl.col("bits") >self.config.bit_score))
        .with_columns((pl.col("target").str.split_exact("|",1).struct[1].alias("EC")))
        .unique(["query","EC"],keep="first")
        .group_by("EC")
        .len()
        ).collect()
        records = df.to_dicts()
        if not records:
            return {}
        return {str(row["EC"]): int(row["len"]) for row in records}
    
    def get_cod_from_ec_counts(self,ec_counts:dict)->dict:
        r"""This function takes a json file that comtains ec counts and converts it to ADM microbial agents counts.
        Required Configs:
            config.adm_mapping : A dictionary that maps ADM reactions to ADM microbial agents.
            ---------
            config.csv_reaction_db : The address of the reaction database of ADToolbox.
            ---------
            config.adm_cod_from_ec  : The address of the json file that the results will be saved in.
            ---------
        Example:
            >>> import os
            >>> output = os.path.join(Main_Dir, "test", "cod_from_ec_counts")
            >>> alignments = ["CP001673.1", "Q8YNF9|1.4.4.2", "0.566", "2859", "414", "0", "1021521", "1024379", "27", "982", "0.000E+00", "1109"]
            >>> headers = ["query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend ", "tstart", "tend", "evalue", "bits"]
            >>> hit = "\t".join(alignments)
            >>> headers_tab = "\t".join(headers)
            >>> combine = headers_tab + "\n" + hit
            >>> with open(output, "w") as f:
            ...     f.write(combine)
            161
            >>> obj = Metagenomics(configs.Metagenomics())
            >>> obj.get_cod_from_ec_counts(output) 
            [('ADM_microbe_1', 161), ('ADM_microbe_2', 0), ('ADM_microbe_3', 0)]

        Args:
            ec_counts (dict): A dictionary containing the counts for each ec number.  
        Returns:
            dict: A dictionary containing the ADM microbial agents counts.
        """
        reaction_rows = {
            row["EC_Numbers"]: row
            for row in pl.read_csv(self.config.csv_reaction_db, separator=",", infer_schema_length=0)
            .unique(subset=["EC_Numbers"], keep="first")
            .to_dicts()
        }
        adm_reactions_agents = {k:0 for k in self.config.adm_mapping.keys()}
        for ec in ec_counts.keys():
            if ec not in reaction_rows:
                continue
            l=reaction_rows[ec]["e_adm_Reactions"].split("|")
            for adm_rxn in l: 
                adm_reactions_agents[adm_rxn]+=ec_counts[ec]
        adm_microbial_agents={}
        for k,v in self.config.adm_mapping.items():
            adm_microbial_agents[v]=adm_reactions_agents[k]
        return adm_microbial_agents

    def _sample_pipeline_logger(self, sample_name: str, output_dir: str | os.PathLike, verbose: bool = True) -> logging.Logger:
        output_path = pathlib.Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        logger = logging.getLogger(f"adtoolbox.metagenomics.{sample_name}")
        logger.setLevel(logging.INFO)
        logger.propagate = False
        logger.handlers.clear()

        formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
        file_handler = logging.FileHandler(output_path / "pipeline.log")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if verbose:
            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(formatter)
            logger.addHandler(stream_handler)
        return logger

    @staticmethod
    def _load_execution_profile(profile: str | os.PathLike | dict | None) -> dict:
        default_profile = {
            "backend": "local",
            "container": "None",
            "image": None,
            "slurm": {},
            "steps": {},
        }
        if profile is None:
            return default_profile
        if isinstance(profile, (str, os.PathLike)):
            with open(profile, "rb") as f:
                loaded_profile = tomllib.load(f)
        else:
            loaded_profile = dict(profile)

        merged = dict(default_profile)
        merged.update({key: value for key, value in loaded_profile.items() if key not in {"slurm", "steps"}})
        merged["slurm"] = {**default_profile["slurm"], **loaded_profile.get("slurm", {})}
        merged["steps"] = loaded_profile.get("steps", {})
        return merged

    @staticmethod
    def _step_settings(execution_profile: dict, step_name: str) -> dict:
        return execution_profile.get("steps", {}).get(step_name, {})

    def _step_container(self, execution_profile: dict, step_name: str, fallback: str) -> str:
        step_settings = self._step_settings(execution_profile, step_name)
        return str(step_settings.get("container", execution_profile.get("container", fallback)))

    def _step_image(self, execution_profile: dict, step_name: str, container: str) -> str | None:
        step_settings = self._step_settings(execution_profile, step_name)
        image = step_settings.get("image", execution_profile.get("image"))
        if image:
            return str(image)
        if container == "docker":
            return self.config.adtoolbox_docker
        if container in {"singularity", "apptainer"}:
            return self.config.adtoolbox_singularity
        return None

    @staticmethod
    def _slurm_script(
        command: str,
        *,
        job_name: str,
        log_file: str | os.PathLike,
        global_slurm: dict,
        step_settings: dict,
        dependency_job_ids: Iterable[str] | None = None,
    ) -> str:
        cpus = step_settings.get("cpus", global_slurm.get("cpus", 1))
        memory = step_settings.get("memory", global_slurm.get("memory", "8G"))
        wall_time = step_settings.get("time", global_slurm.get("time", "01:00:00"))
        partition = step_settings.get("partition", global_slurm.get("partition"))
        account = step_settings.get("account", global_slurm.get("account"))
        qos = step_settings.get("qos", global_slurm.get("qos"))
        max_retries = int(step_settings.get("retries", global_slurm.get("retries", 0)))
        retry_delay = int(step_settings.get("retry_delay_seconds", global_slurm.get("retry_delay_seconds", 60)))
        requeue = PipelineTaskManager._truthy(step_settings.get("requeue", global_slurm.get("requeue", False)))

        lines = [
            "#!/bin/bash",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --cpus-per-task={cpus}",
            f"#SBATCH --mem={memory}",
            f"#SBATCH --time={wall_time}",
            f"#SBATCH --output={log_file}",
        ]
        if partition:
            lines.append(f"#SBATCH --partition={partition}")
        if account:
            lines.append(f"#SBATCH --account={account}")
        if qos:
            lines.append(f"#SBATCH --qos={qos}")
        if requeue:
            lines.append("#SBATCH --requeue")
        dependency_job_ids = list(dependency_job_ids or [])
        if dependency_job_ids:
            lines.append(f"#SBATCH --dependency=afterok:{':'.join(dependency_job_ids)}")
        for option in step_settings.get("extra_sbatch", global_slurm.get("extra_sbatch", [])):
            lines.append(f"#SBATCH {option}")

        command_lines = command.strip().splitlines()
        if max_retries > 0:
            lines.extend(
                [
                    "",
                    "set -uo pipefail",
                    "attempt=0",
                    f"max_retries={max_retries}",
                    f"retry_delay_seconds={retry_delay}",
                    "while true; do",
                    '  echo "ADToolbox step attempt $((attempt + 1))/$((max_retries + 1))"',
                    "  set +e",
                    "  (",
                ]
            )
            lines.extend([f"    {line}" for line in command_lines])
            lines.extend(
                [
                    "  )",
                    "  status=$?",
                    "  set -e",
                    '  if [ "$status" -eq 0 ]; then',
                    "    exit 0",
                    "  fi",
                    '  if [ "$attempt" -ge "$max_retries" ]; then',
                    '    echo "ADToolbox step failed after $((max_retries + 1)) attempt(s)" >&2',
                    '    exit "$status"',
                    "  fi",
                    '  echo "ADToolbox step failed with status $status; retrying in ${retry_delay_seconds}s" >&2',
                    "  attempt=$((attempt + 1))",
                    '  sleep "$retry_delay_seconds"',
                    "done",
                    "",
                ]
            )
        else:
            lines.extend(["", "set -euo pipefail", command.strip(), ""])
        return "\n".join(lines)

    def _execute_step(
        self,
        script: str,
        *,
        step_name: str,
        sample_name: str,
        output_dir: str | os.PathLike,
        logger: logging.Logger,
        execute: bool,
        execution_profile: dict,
        dependencies: Iterable[dict] | None = None,
    ) -> dict:
        output_path = pathlib.Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        step_settings = self._step_settings(execution_profile, step_name)
        backend = str(step_settings.get("backend", execution_profile.get("backend", "local"))).lower()
        command_path = output_path / f"{step_name}.sh"
        command_path.write_text("#!/bin/bash\nset -euo pipefail\n" + script.strip() + "\n")
        manager = PipelineTaskManager(
            sample_name=sample_name,
            output_dir=output_path,
            execution_profile=execution_profile,
            logger=logger,
        )
        dependency_job_ids = manager.dependency_job_ids(dependencies)
        dependency_steps = [
            str(dependency.get("step", dependency.get("step_name", "")))
            for dependency in dependencies or []
            if isinstance(dependency, dict)
        ]
        task = PipelineTask(
            sample_name=sample_name,
            step_name=step_name,
            backend=backend,
            command=str(command_path),
            status="prepared",
            dependencies=dependency_job_ids,
        )

        artifact = {
            "step": step_name,
            "backend": backend,
            "command": str(command_path),
            "executed": execute,
            "status": "prepared",
            "dependencies": dependency_steps,
            "dependency_job_ids": dependency_job_ids,
            "task_events": str(manager.events_path),
        }
        manager.record("prepared", task)
        if backend == "local":
            if execute:
                logger.info("Running local step %s", step_name)
                task.status = "running"
                artifact["status"] = "running"
                manager.record("running", task)
                completed = subprocess.run(script, shell=True, capture_output=True, text=True)
                if completed.stdout:
                    logger.info(completed.stdout.strip())
                if completed.stderr:
                    logger.error(completed.stderr.strip())
                if completed.returncode:
                    message = completed.stderr.strip() or completed.stdout.strip() or f"Step exited with status {completed.returncode}"
                    if len(message) > 1200:
                        message = message[-1200:]
                    task.status = "failed"
                    artifact["status"] = "failed"
                    manager.record("failed", task, returncode=completed.returncode, message=message)
                    raise RuntimeError(f"Step {step_name} failed. See {command_path}. Last output: {message}")
                task.status = "completed"
                artifact["status"] = "completed"
                manager.record("completed", task, returncode=completed.returncode)
            else:
                logger.info("Prepared local step %s at %s", step_name, command_path)
            return artifact

        if backend == "slurm":
            slurm_dir = output_path / "slurm"
            slurm_dir.mkdir(parents=True, exist_ok=True)
            job_name = str(step_settings.get("job_name", f"adtoolbox_{sample_name}_{step_name}"))
            sbatch_path = slurm_dir / f"{step_name}.sbatch"
            slurm_log = slurm_dir / f"{step_name}.%j.out"
            task.sbatch = str(sbatch_path)
            sbatch_path.write_text(
                self._slurm_script(
                    script,
                    job_name=job_name,
                    log_file=slurm_log,
                    global_slurm=execution_profile.get("slurm", {}),
                    step_settings=step_settings,
                    dependency_job_ids=dependency_job_ids,
                )
            )
            artifact["sbatch"] = str(sbatch_path)
            artifact["job_name"] = job_name
            if execute:
                job_prefix = str(execution_profile.get("slurm", {}).get("job_name_prefix", f"adtoolbox_{sample_name}"))
                slurm_settings = execution_profile.get("slurm", {})
                max_retries = int(step_settings.get("retries", slurm_settings.get("retries", 0)))
                wait_for_completion = manager._truthy(
                    step_settings.get("wait_for_completion", slurm_settings.get("wait_for_completion", False))
                )
                poll_seconds = int(step_settings.get("poll_seconds", slurm_settings.get("poll_seconds", 30)))
                logger.info("Submitting Slurm step %s with sbatch %s", step_name, sbatch_path)
                submission, job_id = manager.submit_slurm_job(
                    sbatch_path=sbatch_path,
                    task=task,
                    job_name_prefix=job_prefix,
                    attempt=0,
                )
                artifact["status"] = "submitted"
                artifact["submission"] = submission
                artifact["job_id"] = job_id
                artifact["retries"] = max_retries
                artifact["retry_mode"] = "slurm_job_wrapper" if max_retries > 0 else "none"

                if wait_for_completion and job_id:
                    task.status = "monitoring"
                    artifact["status"] = "monitoring"
                    manager.record("monitoring", task, attempt=0)
                    state = manager.wait_for_slurm_terminal_state(job_id, poll_seconds)
                    artifact["slurm_state"] = state
                    if state == "COMPLETED":
                        task.status = "completed"
                        artifact["status"] = "completed"
                        manager.record("completed", task, attempt=0, slurm_state=state)
                    elif state is None:
                        task.status = "submitted"
                        artifact["status"] = "submitted"
                        manager.record("monitor_unknown", task, attempt=0)
                    else:
                        task.status = "failed"
                        artifact["status"] = "failed"
                        manager.record("failed", task, attempt=0, slurm_state=state)
                        raise RuntimeError(f"Slurm step {step_name} failed as {state}. See {sbatch_path}.")
                elif wait_for_completion and not job_id:
                    manager.record("monitor_skipped", task, attempt=0, reason="missing_job_id")
                    logger.warning("Could not parse Slurm job id for %s; wait_for_completion is disabled for this step", step_name)
            else:
                logger.info("Prepared Slurm step %s at %s", step_name, sbatch_path)
                manager.record("prepared_slurm", task)
            return artifact

        raise ValueError("Execution backend must be local or slurm")

    def _apply_step_config_settings(self, execution_profile: dict, step_name: str) -> dict:
        old_values = {}
        for key, value in self._step_settings(execution_profile, step_name).get("settings", {}).items():
            if hasattr(self.config, key):
                old_values[key] = getattr(self.config, key)
                setattr(self.config, key, value)
        return old_values

    def _restore_config_settings(self, old_values: dict) -> None:
        for key, value in old_values.items():
            setattr(self.config, key, value)

    @staticmethod
    def _quote(value: str | os.PathLike) -> str:
        return shlex.quote(str(value))

    def _container_image(self, container: str, image: str | None = None) -> str:
        if image:
            return image
        if container == "docker":
            return self.config.adtoolbox_docker
        if container in {"singularity", "apptainer"}:
            return self.config.adtoolbox_singularity
        return ""

    def _wrap_external_command(
        self,
        command: str,
        *,
        container: str,
        mounts: Iterable[str | os.PathLike],
        image: str | None = None,
    ) -> str:
        container = str(container)
        if container == "None":
            return command

        mount_dirs = []
        for mount in mounts:
            if mount is None:
                continue
            path = pathlib.Path(mount)
            mount_dir = path if path.suffix == "" else path.parent
            mount_dirs.append(str(mount_dir.absolute()))
        mount_dirs = sorted(set(mount_dirs))

        if container == "docker":
            mount_args = " ".join(f"-v {self._quote(path)}:{self._quote(path)}" for path in mount_dirs)
            return f"docker run {mount_args} {self._container_image(container, image)} sh -lc {self._quote(command)}"
        if container in {"singularity", "apptainer"}:
            runtime = "apptainer" if container == "apptainer" else "singularity"
            bind_args = " ".join(f"-B {self._quote(path)}:{self._quote(path)}" for path in mount_dirs)
            return f"{runtime} exec {bind_args} {self._container_image(container, image)} sh -lc {self._quote(command)}"
        raise ValueError("container must be one of: None, docker, singularity, apptainer")

    @staticmethod
    def _write_json(path: str | os.PathLike, payload: dict) -> str:
        path = pathlib.Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(payload, f, indent=2, sort_keys=True, default=str)
        return str(path)

    @staticmethod
    def _write_tall_mapping(
        path: str | os.PathLike,
        mapping: dict,
        *,
        sample_name: str,
        key_name: str,
        value_name: str,
        value_dtype: pl.DataType = pl.Float64,
    ) -> str:
        rows = [
            {"sample": sample_name, key_name: str(key), value_name: value}
            for key, value in mapping.items()
        ]
        frame = pl.DataFrame(
            rows,
            schema={"sample": pl.Utf8, key_name: pl.Utf8, value_name: value_dtype},
        )
        path = pathlib.Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.write_csv(path)
        return str(path)

    @staticmethod
    def _write_tall_nested_profile(
        path: str | os.PathLike,
        nested: dict[str, dict[str, float]],
        *,
        sample_name: str,
        entity_name: str,
    ) -> str:
        rows = [
            {"sample": sample_name, entity_name: str(entity), "group": str(group), "value": float(value)}
            for entity, profile in nested.items()
            for group, value in profile.items()
        ]
        frame = pl.DataFrame(
            rows,
            schema={"sample": pl.Utf8, entity_name: pl.Utf8, "group": pl.Utf8, "value": pl.Float64},
        )
        path = pathlib.Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.write_csv(path)
        return str(path)

    @staticmethod
    def _read_json_or_table(path: str | os.PathLike) -> dict:
        path = pathlib.Path(path)
        if path.suffix.lower() == ".json":
            with open(path) as f:
                return json.load(f)

        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        table = pl.read_csv(path, separator=sep, infer_schema_length=0)
        if table.width < 2:
            raise ValueError(f"{path} must contain at least two columns")
        key_candidates = ["genome_id", "feature_id", "ec", "group", table.columns[0]]
        value_candidates = ["abundance", "count", "value", table.columns[-1]]
        key_col = next(column for column in key_candidates if column in table.columns)
        value_col = next(column for column in value_candidates if column in table.columns)
        return {
            str(row[key_col]): float(row[value_col])
            for row in table.select([key_col, value_col]).to_dicts()
        }

    @staticmethod
    def _read_sample_manifest(path: str | os.PathLike) -> list[dict]:
        path = pathlib.Path(path)
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        table = pl.read_csv(path, separator=sep, infer_schema_length=0).fill_null("")
        return table.to_dicts()

    @staticmethod
    def _row_value(row: dict, *keys: str) -> str | None:
        normalized = {str(key).lower(): value for key, value in row.items()}
        for key in keys:
            value = normalized.get(key.lower())
            if value not in (None, ""):
                return str(value)
        return None

    @staticmethod
    def _row_bool(row: dict, key: str, default: bool = True) -> bool:
        value = Metagenomics._row_value(row, key)
        if value is None:
            return default
        return value.strip().lower() not in {"0", "false", "no", "single", "single-end"}

    @staticmethod
    def _normalize_profile(profile: dict[str, float], keys: Iterable[str]) -> dict[str, float]:
        normalized = {key: float(profile.get(key, 0) or 0) for key in keys}
        total = sum(value for value in normalized.values() if value > 0)
        if total > 0:
            normalized = {key: value / total for key, value in normalized.items()}
        return normalized

    def _reaction_column(self, reaction_db) -> str:
        for column in ("e_adm_Reactions", "e_adm_reactions", "Modified_ADM_Reactions"):
            if column in reaction_db.columns:
                return column
        raise ValueError("Reaction database must contain an e-ADM reaction mapping column")

    def cod_from_ec_counts(self, ec_counts: dict[str, int | float], normalize: bool = True) -> dict[str, float]:
        reaction_db = pl.read_csv(self.config.csv_reaction_db).drop_nulls("EC_Numbers")
        reaction_column = self._reaction_column(reaction_db)
        reaction_lookup = {
            str(row["EC_Numbers"]): row[reaction_column]
            for row in reaction_db.unique(subset=["EC_Numbers"], keep="first").select(["EC_Numbers", reaction_column]).to_dicts()
        }
        reaction_scores = {reaction: 0.0 for reaction in self.config.adm_mapping}

        for ec, count in ec_counts.items():
            mapped_reactions = reaction_lookup.get(str(ec))
            if mapped_reactions in (None, ""):
                continue
            for reaction in str(mapped_reactions).split("|"):
                reaction = reaction.strip()
                if reaction in reaction_scores:
                    reaction_scores[reaction] += float(count)

        microbe_scores = {group: 0.0 for group in self.config.adm_mapping.values()}
        for reaction, group in self.config.adm_mapping.items():
            microbe_scores[group] += reaction_scores.get(reaction, 0.0)
        return self._normalize_profile(microbe_scores, microbe_scores) if normalize else microbe_scores

    def cod_from_alignment(self, alignment_file: str | os.PathLike, normalize: bool = True) -> dict[str, float]:
        return self.cod_from_ec_counts(self.extract_ec_from_alignment(str(alignment_file)), normalize=normalize)

    def _alignment_files_from_path(self, path: str | os.PathLike) -> dict[str, str]:
        path = pathlib.Path(path)
        if path.is_file() and path.suffix.lower() == ".json":
            with open(path) as f:
                return {str(key): str(value) for key, value in json.load(f).items()}
        if path.is_file():
            return {path.stem.replace("Alignment_Results_mmseq_", ""): str(path)}
        if not path.exists():
            raise FileNotFoundError(path)
        alignments = {}
        for alignment in sorted(path.glob("Alignment_Results_mmseq_*.tsv")):
            name = alignment.stem.replace("Alignment_Results_mmseq_", "")
            alignments[name] = str(alignment)
            alignments.setdefault(name.split("~", 1)[0], str(alignment))
        if not alignments:
            raise FileNotFoundError(f"No Alignment_Results_mmseq_*.tsv files found in {path}")
        return alignments

    def _genome_files_from_dir(self, genomes_dir: str | os.PathLike) -> dict[str, str]:
        genome_info = self.extract_genome_info(str(genomes_dir))
        normalized = dict(genome_info)
        for key, value in genome_info.items():
            normalized.setdefault(key.replace("_genomic", ""), value)
            normalized.setdefault(key.split("_genomic")[0], value)
            accession_match = re.match(r"^(GC[AF]_\d+\.\d+)", key)
            if accession_match:
                normalized.setdefault(accession_match.group(1), value)
        return normalized

    def _run_shell_script(self, script: str, logger: logging.Logger, execute: bool) -> None:
        if not execute:
            logger.info("Prepared command: %s", script.strip())
            return
        logger.info("Running command: %s", script.strip())
        subprocess.run(script, shell=True, check=True)

    def _write_sample_repseqs(
        self,
        rep_seqs: str | os.PathLike,
        feature_abundances: dict[str, float],
        output_fasta: str | os.PathLike,
    ) -> str:
        sequences = {}
        for header, sequence in fasta_to_dict(str(rep_seqs)).items():
            feature_id = header.split(";", 1)[0].split(None, 1)[0]
            sequences.setdefault(feature_id, sequence)
        selected = {feature: sequences[feature] for feature in feature_abundances if feature in sequences}
        if not selected:
            raise ValueError("None of the selected feature IDs were found in the representative sequence FASTA")
        utils.dict_to_fasta(selected, str(output_fasta))
        return str(output_fasta)

    def trim_amplicon_reads(
        self,
        *,
        read_1: str | os.PathLike,
        read_2: str | os.PathLike | None,
        output_dir: str | os.PathLike,
        sample_name: str,
        forward_primer: str | None = None,
        reverse_primer: str | None = None,
        adapter_1: str | None = None,
        adapter_2: str | None = None,
        minimum_length: int = 100,
        quality_cutoff: str | int | None = None,
        threads: int = 1,
        container: str = "None",
        image: str | None = None,
    ) -> tuple[str, dict[str, str | None]]:
        output_path = pathlib.Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        trimmed_1 = output_path / f"{sample_name}_trimmed_R1.fastq.gz"
        trimmed_2 = output_path / f"{sample_name}_trimmed_R2.fastq.gz" if read_2 else None

        adapter_1 = adapter_1 or forward_primer
        adapter_2 = adapter_2 or reverse_primer

        report_prefix = output_path / f"{sample_name}_fastp"
        args = [
            "fastp",
            "-w",
            str(threads),
            "--length_required",
            str(minimum_length),
            "-i",
            str(read_1),
            "-o",
            str(trimmed_1),
            "--json",
            str(report_prefix.with_suffix(".json")),
            "--html",
            str(report_prefix.with_suffix(".html")),
        ]
        if read_2:
            args.extend(
                [
                    "-I",
                    str(read_2),
                    "-O",
                    str(trimmed_2),
                ]
            )
            if not adapter_1 and not adapter_2:
                args.append("--detect_adapter_for_pe")
        if quality_cutoff is not None:
            args.extend(["--qualified_quality_phred", str(quality_cutoff).split(",", 1)[0]])
        if adapter_1:
            args.extend(["--adapter_sequence", adapter_1])
        if adapter_2 and read_2:
            args.extend(["--adapter_sequence_r2", adapter_2])

        command = " ".join(self._quote(arg) for arg in args)
        script = self._wrap_external_command(
            command,
            container=container,
            mounts=[read_1, read_2, output_path],
            image=image,
        )
        return script + "\n", {"read_1": str(trimmed_1), "read_2": str(trimmed_2) if trimmed_2 else None}

    def build_amplicon_features(
        self,
        *,
        read_1: str | os.PathLike,
        read_2: str | os.PathLike | None,
        output_dir: str | os.PathLike,
        sample_name: str,
        identity: float = 0.97,
        maxee: float = 1.0,
        minimum_length: int = 100,
        min_unique_size: int = 2,
        chimera_filter: bool = True,
        threads: int = 1,
        container: str = "None",
        image: str | None = None,
    ) -> tuple[str, dict[str, str]]:
        output_path = pathlib.Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        merged = output_path / f"{sample_name}_merged.fastq"
        filtered = output_path / f"{sample_name}_filtered.fasta"
        uniques = output_path / f"{sample_name}_uniques.fasta"
        denoised = output_path / f"{sample_name}_denoised.fasta"
        rep_seqs = output_path / "rep-seqs.fasta"
        feature_table = output_path / "feature-table.tsv"
        raw_feature_table = output_path / f"{sample_name}_vsearch_otutab.tsv"
        derep_uc = output_path / f"{sample_name}_derep.uc"

        commands = []
        if read_2:
            commands.append(
                " ".join(
                    self._quote(arg)
                    for arg in [
                        "vsearch",
                        "--fastq_mergepairs",
                        str(read_1),
                        "--reverse",
                        str(read_2),
                        "--fastqout",
                        str(merged),
                        "--threads",
                        str(threads),
                    ]
                )
            )
            fastq_input = merged
        else:
            fastq_input = pathlib.Path(read_1)

        commands.append(
            " ".join(
                self._quote(arg)
                for arg in [
                    "vsearch",
                    "--fastq_filter",
                    str(fastq_input),
                    "--fastq_maxee",
                    str(maxee),
                    "--fastq_minlen",
                    str(minimum_length),
                    "--fastaout",
                    str(filtered),
                ]
            )
        )
        commands.append(
            " ".join(
                self._quote(arg)
                for arg in [
                    "vsearch",
                    "--derep_fulllength",
                    str(filtered),
                    "--output",
                    str(uniques),
                    "--sizeout",
                    "--minuniquesize",
                    str(min_unique_size),
                    "--uc",
                    str(derep_uc),
                ]
            )
        )
        commands.append(
            " ".join(
                self._quote(arg)
                for arg in [
                    "vsearch",
                    "--cluster_unoise",
                    str(uniques),
                    "--centroids",
                    str(denoised),
                    "--minsize",
                    str(min_unique_size),
                ]
            )
        )
        if chimera_filter:
            commands.append(
                " ".join(
                    self._quote(arg)
                    for arg in [
                        "vsearch",
                        "--uchime3_denovo",
                        str(denoised),
                        "--nonchimeras",
                        str(rep_seqs),
                    ]
                )
            )
        else:
            commands.append(f"cp {self._quote(denoised)} {self._quote(rep_seqs)}")
        commands.append(
            " ".join(
                self._quote(arg)
                for arg in [
                    "vsearch",
                    "--usearch_global",
                    str(filtered),
                    "--db",
                    str(rep_seqs),
                    "--id",
                    str(identity),
                    "--otutabout",
                    str(raw_feature_table),
                    "--threads",
                    str(threads),
                ]
            )
        )
        commands.append(
            f"awk -v sample={self._quote(sample_name)} "
            f"'BEGIN{{FS=OFS=\"\\t\"}} NR==1{{if (NF >= 2) $2=sample; print; next}} {{print}}' "
            f"{self._quote(raw_feature_table)} > {self._quote(feature_table)}"
        )

        command = "\n".join(commands)
        script = self._wrap_external_command(
            command,
            container=container,
            mounts=[read_1, read_2, output_path],
            image=image,
        )
        return script + "\n", {"feature_table": str(feature_table), "rep_seqs": str(rep_seqs)}

    def run_trim_reads_step(
        self,
        *,
        sample_name: str,
        output_dir: str | os.PathLike,
        read_1: str | os.PathLike,
        read_2: str | os.PathLike | None = None,
        step_output_dir: str | os.PathLike | None = None,
        forward_primer: str | None = None,
        reverse_primer: str | None = None,
        adapter_1: str | None = None,
        adapter_2: str | None = None,
        minimum_length: int = 100,
        quality_cutoff: str | int | None = None,
        container: str = "None",
        execute: bool = False,
        verbose: bool = True,
        execution_profile: str | os.PathLike | dict | None = None,
        dependencies: Iterable[dict] | None = None,
    ) -> dict:
        """Prepare or run the amplicon read trimming step for one sample."""
        sample_dir = pathlib.Path(output_dir) / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)
        step_output = pathlib.Path(step_output_dir) if step_output_dir else sample_dir / "amplicon_preprocess"
        logger = self._sample_pipeline_logger(sample_name, sample_dir, verbose=verbose)
        execution_profile = self._load_execution_profile(execution_profile)

        step_name = "trim_reads"
        settings = self._step_settings(execution_profile, step_name).get("settings", {})
        step_container = self._step_container(execution_profile, step_name, container)
        script, trimmed_reads = self.trim_amplicon_reads(
            read_1=read_1,
            read_2=read_2,
            output_dir=step_output,
            sample_name=sample_name,
            forward_primer=settings.get("forward_primer", forward_primer),
            reverse_primer=settings.get("reverse_primer", reverse_primer),
            adapter_1=settings.get("adapter_1", adapter_1),
            adapter_2=settings.get("adapter_2", adapter_2),
            minimum_length=int(settings.get("minimum_length", minimum_length)),
            quality_cutoff=settings.get("quality_cutoff", quality_cutoff),
            threads=int(settings.get("threads", self._step_settings(execution_profile, step_name).get("cpus", 1))),
            container=step_container,
            image=self._step_image(execution_profile, step_name, step_container),
        )
        artifact = self._execute_step(
            script,
            step_name=step_name,
            sample_name=sample_name,
            output_dir=step_output.parent if step_output_dir else sample_dir,
            logger=logger,
            execute=execute,
            execution_profile=execution_profile,
            dependencies=dependencies,
        )
        return {
            "sample_name": sample_name,
            "step": step_name,
            "output_dir": str(sample_dir),
            "execute": execute,
            "artifacts": {
                step_name: artifact,
                "trimmed_reads": trimmed_reads,
            },
        }

    def _sra_download_script(
        self,
        *,
        accession: str,
        target_dir: str | os.PathLike,
        container: str = "None",
        image: str | None = None,
    ) -> tuple[str, dict[str, str]]:
        target_path = pathlib.Path(target_dir)
        accession_dir = target_path / accession
        sra_file = accession_dir / f"{accession}.sra"
        read_1 = accession_dir / f"{accession}_1.fastq"
        read_2 = accession_dir / f"{accession}_2.fastq"
        command = f"""set -e
mkdir -p {self._quote(accession_dir)}
if command -v prefetch >/dev/null 2>&1 && command -v fasterq-dump >/dev/null 2>&1; then
  prefetch {self._quote(accession)} -O {self._quote(target_path)} --max-size 100000000
  fasterq-dump {self._quote(sra_file)} -O {self._quote(accession_dir)} --split-3 --temp {self._quote(accession_dir)}
  rm -f {self._quote(sra_file)}
elif command -v curl >/dev/null 2>&1; then
  ena_url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession={self._quote(accession)}&result=read_run&fields=fastq_ftp&format=tsv&download=false"
  fastq_ftp=$(curl -fsSL "$ena_url" | awk -F '\\t' 'NR==2 {{print $NF}}')
  if [ -z "$fastq_ftp" ]; then
    echo "Could not find ENA FASTQ URLs for {self._quote(accession)}" >&2
    exit 1
  fi
  old_ifs=$IFS
  IFS=';'
  for ftp_path in $fastq_ftp; do
    IFS=$old_ifs
    file_name=$(basename "$ftp_path")
    case "$ftp_path" in
      ftp://*|https://*) download_url="$ftp_path" ;;
      *) download_url="ftp://$ftp_path" ;;
    esac
    https_url=$(printf '%s' "$download_url" | sed 's#^ftp://#https://#')
    echo "Downloading $file_name from ENA"
    curl -fsSL --retry 3 --connect-timeout 30 "$https_url" -o {self._quote(str(accession_dir))}/"$file_name" || \
      curl -fsSL --retry 3 --connect-timeout 30 "$download_url" -o {self._quote(str(accession_dir))}/"$file_name"
    IFS=';'
  done
  IFS=$old_ifs
else
  echo "SRA download requires either prefetch/fasterq-dump from SRA Toolkit or curl for ENA FASTQ download." >&2
  exit 127
fi"""
        script = self._wrap_external_command(
            command,
            container=container,
            mounts=[target_path],
            image=image,
        )
        return script + "\n", {"read_1": str(read_1), "read_2": str(read_2)}

    @staticmethod
    def _resolved_sra_reads(accession: str, target_dir: str | os.PathLike, *, paired: bool) -> dict[str, str | None]:
        accession_dir = pathlib.Path(target_dir) / accession
        read_1_candidates = [accession_dir / f"{accession}_1.fastq", accession_dir / f"{accession}_1.fastq.gz"]
        read_2_candidates = [accession_dir / f"{accession}_2.fastq", accession_dir / f"{accession}_2.fastq.gz"]
        single_candidates = [accession_dir / f"{accession}.fastq", accession_dir / f"{accession}.fastq.gz"]
        read_1 = next((path for path in read_1_candidates if path.exists()), None)
        read_2 = next((path for path in read_2_candidates if path.exists()), None)
        single = next((path for path in single_candidates if path.exists()), None)
        if read_1 is not None and read_2 is not None:
            return {"read_1": str(read_1), "read_2": str(read_2)}
        if single is not None:
            return {"read_1": str(single), "read_2": None}
        if read_1 is not None and not paired:
            return {"read_1": str(read_1), "read_2": None}
        expected = ", ".join(str(path) for path in read_1_candidates + read_2_candidates + single_candidates)
        raise FileNotFoundError(f"SRA download finished, but no FASTQ output was found. Expected {expected}")

    def run_sra_download_step(
        self,
        *,
        sample_name: str,
        output_dir: str | os.PathLike,
        accession: str,
        sra_dir: str | os.PathLike | None = None,
        paired: bool = True,
        container: str = "None",
        execute: bool = False,
        verbose: bool = True,
        execution_profile: str | os.PathLike | dict | None = None,
        dependencies: Iterable[dict] | None = None,
    ) -> dict:
        """Prepare or run an SRA download step for one sample."""
        sample_dir = pathlib.Path(output_dir) / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)
        target_dir = pathlib.Path(sra_dir) if sra_dir else pathlib.Path(output_dir) / "sra"
        target_dir.mkdir(parents=True, exist_ok=True)
        logger = self._sample_pipeline_logger(sample_name, sample_dir, verbose=verbose)
        execution_profile = self._load_execution_profile(execution_profile)
        step_name = "download_sra"
        step_container = self._step_container(execution_profile, step_name, container)
        script, reads = self._sra_download_script(
            accession=accession,
            target_dir=target_dir,
            container=step_container,
            image=self._step_image(execution_profile, step_name, step_container),
        )
        if not paired:
            reads["read_2"] = None
        artifact = self._execute_step(
            script,
            step_name=step_name,
            sample_name=sample_name,
            output_dir=sample_dir,
            logger=logger,
            execute=execute,
            execution_profile=execution_profile,
            dependencies=dependencies,
        )
        if execute and artifact["status"] in {"completed", "running"}:
            reads = self._resolved_sra_reads(accession, target_dir, paired=paired)
        elif not paired:
            reads["read_2"] = None
        return {
            "sample_name": sample_name,
            "step": step_name,
            "output_dir": str(sample_dir),
            "execute": execute,
            "artifacts": {
                step_name: artifact,
                "accession": accession,
                "reads": reads,
            },
        }

    def run_build_amplicon_features_step(
        self,
        *,
        sample_name: str,
        output_dir: str | os.PathLike,
        read_1: str | os.PathLike,
        read_2: str | os.PathLike | None = None,
        step_output_dir: str | os.PathLike | None = None,
        identity: float = 0.97,
        maxee: float = 1.0,
        minimum_length: int = 100,
        min_unique_size: int = 2,
        chimera_filter: bool = True,
        container: str = "None",
        execute: bool = False,
        verbose: bool = True,
        execution_profile: str | os.PathLike | dict | None = None,
        dependencies: Iterable[dict] | None = None,
    ) -> dict:
        """Prepare or run the VSEARCH feature-table step for one sample."""
        sample_dir = pathlib.Path(output_dir) / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)
        step_output = pathlib.Path(step_output_dir) if step_output_dir else sample_dir / "amplicon_preprocess"
        logger = self._sample_pipeline_logger(sample_name, sample_dir, verbose=verbose)
        execution_profile = self._load_execution_profile(execution_profile)

        step_name = "build_amplicon_features"
        settings = self._step_settings(execution_profile, step_name).get("settings", {})
        step_container = self._step_container(execution_profile, step_name, container)
        script, feature_artifacts = self.build_amplicon_features(
            read_1=read_1,
            read_2=read_2,
            output_dir=step_output,
            sample_name=sample_name,
            identity=float(settings.get("identity", identity)),
            maxee=float(settings.get("maxee", maxee)),
            minimum_length=int(settings.get("minimum_length", minimum_length)),
            min_unique_size=int(settings.get("min_unique_size", min_unique_size)),
            chimera_filter=bool(settings.get("chimera_filter", chimera_filter)),
            threads=int(settings.get("threads", self._step_settings(execution_profile, step_name).get("cpus", 1))),
            container=step_container,
            image=self._step_image(execution_profile, step_name, step_container),
        )
        artifact = self._execute_step(
            script,
            step_name=step_name,
            sample_name=sample_name,
            output_dir=step_output.parent if step_output_dir else sample_dir,
            logger=logger,
            execute=execute,
            execution_profile=execution_profile,
            dependencies=dependencies,
        )
        return {
            "sample_name": sample_name,
            "step": step_name,
            "output_dir": str(sample_dir),
            "execute": execute,
            "artifacts": {
                step_name: artifact,
                **feature_artifacts,
            },
        }

    def run_short_read_alignment_step(
        self,
        *,
        sample_name: str,
        output_dir: str | os.PathLike,
        reads: str | os.PathLike,
        container: str = "None",
        execute: bool = False,
        verbose: bool = True,
        execution_profile: str | os.PathLike | dict | None = None,
        dependencies: Iterable[dict] | None = None,
    ) -> dict:
        """Prepare or run the MMseqs shotgun-read alignment step."""
        sample_dir = pathlib.Path(output_dir) / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)
        logger = self._sample_pipeline_logger(sample_name, sample_dir, verbose=verbose)
        execution_profile = self._load_execution_profile(execution_profile)
        step_name = "align_short_reads"
        step_container = self._step_container(execution_profile, step_name, container)
        old_settings = self._apply_step_config_settings(execution_profile, step_name)
        try:
            script, alignment_path = self.align_short_reads_to_protein_db(
                str(reads),
                f"{sample_name}_mmseq",
                container=step_container,
                image=self._step_image(execution_profile, step_name, step_container),
            )
        finally:
            self._restore_config_settings(old_settings)
        artifact = self._execute_step(
            script,
            step_name=step_name,
            sample_name=sample_name,
            output_dir=sample_dir,
            logger=logger,
            execute=execute,
            execution_profile=execution_profile,
        )
        return {
            "sample_name": sample_name,
            "step": step_name,
            "output_dir": str(sample_dir),
            "execute": execute,
            "artifacts": {
                step_name: artifact,
                "alignment_file": str(alignment_path),
            },
        }

    def run_gtdb_alignment_step(
        self,
        *,
        sample_name: str,
        output_dir: str | os.PathLike,
        query_fasta: str | os.PathLike,
        container: str = "None",
        execute: bool = False,
        verbose: bool = True,
        execution_profile: str | os.PathLike | dict | None = None,
    ) -> dict:
        """Prepare or run the VSEARCH representative-sequence to GTDB step."""
        sample_dir = pathlib.Path(output_dir) / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)
        logger = self._sample_pipeline_logger(sample_name, sample_dir, verbose=verbose)
        execution_profile = self._load_execution_profile(execution_profile)
        step_name = "align_to_gtdb"
        step_container = self._step_container(execution_profile, step_name, container)
        old_settings = self._apply_step_config_settings(execution_profile, step_name)
        try:
            script = self.align_to_gtdb(
                str(query_fasta),
                str(sample_dir),
                container=step_container,
                image=self._step_image(execution_profile, step_name, step_container),
            )[0]
        finally:
            self._restore_config_settings(old_settings)
        artifact = self._execute_step(
            script,
            step_name=step_name,
            sample_name=sample_name,
            output_dir=sample_dir,
            logger=logger,
            execute=execute,
            execution_profile=execution_profile,
        )
        return {
            "sample_name": sample_name,
            "step": step_name,
            "output_dir": str(sample_dir),
            "execute": execute,
            "artifacts": {
                step_name: artifact,
                "matches": str(sample_dir / "matches.blast"),
            },
        }

    def run_genome_alignment_step(
        self,
        *,
        sample_name: str,
        output_dir: str | os.PathLike,
        genome_name: str,
        genome_file: str | os.PathLike,
        alignment_dir: str | os.PathLike | None = None,
        container: str = "None",
        execute: bool = False,
        verbose: bool = True,
        execution_profile: str | os.PathLike | dict | None = None,
    ) -> dict:
        """Prepare or run one genome-to-protein-database MMseqs alignment step."""
        sample_dir = pathlib.Path(output_dir) / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)
        alignment_path = pathlib.Path(alignment_dir) if alignment_dir else sample_dir / "genome_alignments"
        alignment_path.mkdir(parents=True, exist_ok=True)
        logger = self._sample_pipeline_logger(sample_name, sample_dir, verbose=verbose)
        execution_profile = self._load_execution_profile(execution_profile)
        profile_step_name = "align_genome"
        step_name = f"align_genome_{genome_name}"
        step_container = self._step_container(execution_profile, profile_step_name, container)
        old_settings = self._apply_step_config_settings(execution_profile, profile_step_name)
        try:
            script, alignment = self.align_genome_to_protein_db(
                str(genome_file),
                str(alignment_path),
                genome_name,
                container=step_container,
                image=self._step_image(execution_profile, profile_step_name, step_container),
            )
        finally:
            self._restore_config_settings(old_settings)
        artifact = self._execute_step(
            script,
            step_name=step_name,
            sample_name=sample_name,
            output_dir=sample_dir,
            logger=logger,
            execute=execute,
            execution_profile={
                **execution_profile,
                "steps": {
                    **execution_profile.get("steps", {}),
                    step_name: self._step_settings(execution_profile, profile_step_name),
                },
            },
        )
        return {
            "sample_name": sample_name,
            "step": step_name,
            "output_dir": str(sample_dir),
            "execute": execute,
            "artifacts": {
                step_name: artifact,
                "alignment_file": str(alignment),
            },
        }

    def preprocess_amplicon_sample(
        self,
        *,
        sample_name: str,
        output_dir: str | os.PathLike,
        read_1: str | os.PathLike,
        read_2: str | os.PathLike | None = None,
        forward_primer: str | None = None,
        reverse_primer: str | None = None,
        adapter_1: str | None = None,
        adapter_2: str | None = None,
        minimum_length: int = 100,
        quality_cutoff: str | int | None = None,
        quality_maxee: float = 1.0,
        identity: float = 0.97,
        min_unique_size: int = 2,
        chimera_filter: bool = True,
        container: str = "None",
        execute: bool = False,
        verbose: bool = True,
        execution_profile: str | os.PathLike | dict | None = None,
        dependencies: Iterable[dict] | None = None,
    ) -> dict:
        output_path = pathlib.Path(output_dir) / sample_name
        scratch_path = output_path / "scratch"
        preprocess_path = scratch_path / "amplicon_preprocess"
        preprocess_path.mkdir(parents=True, exist_ok=True)
        logger = self._sample_pipeline_logger(sample_name, output_path, verbose=verbose)
        execution_profile = self._load_execution_profile(execution_profile)
        result = {
            "sample_name": sample_name,
            "output_dir": str(output_path),
            "execute": execute,
            "artifacts": {},
        }

        trim_result = self.run_trim_reads_step(
            sample_name=sample_name,
            output_dir=output_dir,
            read_1=read_1,
            read_2=read_2,
            step_output_dir=preprocess_path,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            adapter_1=adapter_1,
            adapter_2=adapter_2,
            minimum_length=minimum_length,
            quality_cutoff=quality_cutoff,
            container=container,
            execute=execute,
            verbose=verbose,
            execution_profile=execution_profile,
            dependencies=dependencies,
        )
        result["artifacts"].update(trim_result["artifacts"])
        trimmed_reads = trim_result["artifacts"]["trimmed_reads"]

        feature_result = self.run_build_amplicon_features_step(
            sample_name=sample_name,
            output_dir=output_dir,
            read_1=trimmed_reads["read_1"],
            read_2=trimmed_reads["read_2"],
            step_output_dir=preprocess_path,
            identity=identity,
            maxee=quality_maxee,
            minimum_length=minimum_length,
            min_unique_size=min_unique_size,
            chimera_filter=chimera_filter,
            container=container,
            execute=execute,
            verbose=verbose,
            execution_profile=execution_profile,
            dependencies=[trim_result["artifacts"]["trim_reads"]],
        )
        result["artifacts"].update(feature_result["artifacts"])
        result["artifacts"]["preprocess_manifest"] = self._write_json(scratch_path / "preprocess_artifacts.json", result["artifacts"])
        logger.info("Finished amplicon preprocessing: %s", sample_name)
        return result

    def aggregate_genome_cod(
        self,
        genome_cods: dict[str, dict[str, float]],
        genome_abundances: dict[str, float],
        normalize: bool = True,
    ) -> dict[str, float]:
        groups = set(self.config.adm_mapping.values())
        aggregated = {group: 0.0 for group in groups}
        abundance_total = sum(float(value) for value in genome_abundances.values() if float(value) > 0)
        if abundance_total <= 0:
            raise ValueError("Genome abundances must contain at least one positive value")

        for genome, abundance in genome_abundances.items():
            if genome not in genome_cods:
                continue
            weight = float(abundance) / abundance_total
            for group, value in genome_cods[genome].items():
                aggregated[group] = aggregated.get(group, 0.0) + float(value) * weight
        return self._normalize_profile(aggregated, groups) if normalize else aggregated

    def sample_to_cod(
        self,
        sample_name: str,
        output_dir: str | os.PathLike,
        *,
        mode: str,
        alignment_file: str | os.PathLike | None = None,
        reads: str | os.PathLike | None = None,
        read_1: str | os.PathLike | None = None,
        read_2: str | os.PathLike | None = None,
        genome_abundances: str | os.PathLike | dict[str, float] | None = None,
        genome_alignments: str | os.PathLike | dict[str, str] | None = None,
        genomes_dir: str | os.PathLike | None = None,
        feature_table: str | os.PathLike | None = None,
        rep_seqs: str | os.PathLike | None = None,
        gtdb_matches: str | os.PathLike | None = None,
        forward_primer: str | None = None,
        reverse_primer: str | None = None,
        adapter_1: str | None = None,
        adapter_2: str | None = None,
        minimum_length: int = 100,
        quality_cutoff: str | int | None = None,
        quality_maxee: float = 1.0,
        identity: float = 0.97,
        min_unique_size: int = 2,
        chimera_filter: bool = True,
        top_k: int = -1,
        container: str = "None",
        execute: bool = False,
        normalize: bool = True,
        verbose: bool = True,
        execution_profile: str | os.PathLike | dict | None = None,
        dependencies: Iterable[dict] | None = None,
    ) -> dict:
        """Run one metagenomics-to-eADM-COD pipeline for a single sample.

        The method writes tall CSV tables for table-shaped artifacts,
        plus ``provenance.json`` and ``pipeline.log`` under ``output_dir``.
        """
        output_path = pathlib.Path(output_dir) / sample_name
        output_path.mkdir(parents=True, exist_ok=True)
        scratch_path = output_path / "scratch"
        scratch_path.mkdir(parents=True, exist_ok=True)
        logger = self._sample_pipeline_logger(sample_name, output_path, verbose=verbose)
        execution_profile = self._load_execution_profile(execution_profile)
        logger.info("Starting sample COD pipeline: sample=%s mode=%s", sample_name, mode)

        result = {
            "sample_name": sample_name,
            "mode": mode,
            "output_dir": str(output_path),
            "execute": execute,
            "artifacts": {},
        }

        if mode == "shotgun-alignment":
            if alignment_file is None:
                raise ValueError("alignment_file is required for shotgun-alignment mode")
            logger.info("Converting shotgun alignment to EC counts and COD profile")
            ec_counts = self.extract_ec_from_alignment(str(alignment_file))
            cod_profile = self.cod_from_ec_counts(ec_counts, normalize=normalize)
            result["artifacts"]["ec_counts"] = self._write_tall_mapping(output_path / "ec_counts.csv", ec_counts, sample_name=sample_name, key_name="ec", value_name="count", value_dtype=pl.Int64)
            result["artifacts"]["cod_profile"] = self._write_tall_mapping(output_path / "cod_profile.csv", cod_profile, sample_name=sample_name, key_name="group", value_name="value")

        elif mode == "shotgun-reads":
            if reads is None:
                raise ValueError("reads is required for shotgun-reads mode")
            logger.info("Aligning shotgun reads to protein database")
            step_name = "align_short_reads"
            old_settings = self._apply_step_config_settings(execution_profile, step_name)
            try:
                script, alignment_path = self.align_short_reads_to_protein_db(
                    str(reads),
                    f"{sample_name}_mmseq",
                    container=self._step_container(execution_profile, step_name, container),
                )
            finally:
                self._restore_config_settings(old_settings)
            result["artifacts"][step_name] = self._execute_step(
                script,
                step_name=step_name,
                sample_name=sample_name,
                output_dir=scratch_path,
                logger=logger,
                execute=execute,
                execution_profile=execution_profile,
            )
            if not pathlib.Path(alignment_path).exists():
                logger.info("Alignment output is not available yet: %s", alignment_path)
                cod_profile = {}
                result["status"] = "waiting_for_alignment"
            else:
                ec_counts = self.extract_ec_from_alignment(str(alignment_path))
                cod_profile = self.cod_from_ec_counts(ec_counts, normalize=normalize)
                result["artifacts"]["ec_counts"] = self._write_tall_mapping(output_path / "ec_counts.csv", ec_counts, sample_name=sample_name, key_name="ec", value_name="count", value_dtype=pl.Int64)
                result["artifacts"]["cod_profile"] = self._write_tall_mapping(output_path / "cod_profile.csv", cod_profile, sample_name=sample_name, key_name="group", value_name="value")

        elif mode == "genome-alignments":
            if genome_abundances is None or genome_alignments is None:
                raise ValueError("genome_abundances and genome_alignments are required for genome-alignments mode")
            abundances = genome_abundances if isinstance(genome_abundances, dict) else self._read_json_or_table(genome_abundances)
            alignments = genome_alignments if isinstance(genome_alignments, dict) else self._alignment_files_from_path(genome_alignments)
            logger.info("Converting %s genome alignments to genome-level COD profiles", len(alignments))
            genome_cods = {
                genome: self.cod_from_alignment(alignment, normalize=normalize)
                for genome, alignment in alignments.items()
                if genome in abundances
            }
            cod_profile = self.aggregate_genome_cod(genome_cods, abundances, normalize=normalize)
            result["artifacts"]["genome_cods"] = self._write_tall_nested_profile(output_path / "genome_cods.csv", genome_cods, sample_name=sample_name, entity_name="genome_id")
            result["artifacts"]["cod_profile"] = self._write_tall_mapping(output_path / "cod_profile.csv", cod_profile, sample_name=sample_name, key_name="group", value_name="value")

        elif mode == "amplicon-reads":
            if read_1 is None:
                raise ValueError("read_1 is required for amplicon-reads mode")
            logger.info("Preprocessing amplicon reads before COD conversion")
            preprocess_result = self.preprocess_amplicon_sample(
                sample_name=sample_name,
                output_dir=output_dir,
                read_1=read_1,
                read_2=read_2,
                forward_primer=forward_primer,
                reverse_primer=reverse_primer,
                adapter_1=adapter_1,
                adapter_2=adapter_2,
                minimum_length=minimum_length,
                quality_cutoff=quality_cutoff,
                quality_maxee=quality_maxee,
                identity=identity,
                min_unique_size=min_unique_size,
                chimera_filter=chimera_filter,
                container=container,
                execute=execute,
                verbose=verbose,
                execution_profile=execution_profile,
                dependencies=dependencies,
            )
            result["artifacts"]["preprocess"] = preprocess_result["artifacts"]
            feature_table = preprocess_result["artifacts"]["feature_table"]
            rep_seqs = preprocess_result["artifacts"]["rep_seqs"]

            if not pathlib.Path(feature_table).exists() or not pathlib.Path(rep_seqs).exists():
                logger.info("Amplicon feature outputs are not available yet")
                cod_profile = {}
                result["status"] = "waiting_for_preprocess"
            else:
                downstream_dependencies = [preprocess_result["artifacts"]["build_amplicon_features"]]
                downstream_result = self.sample_to_cod(
                    sample_name=sample_name,
                    output_dir=output_dir,
                    mode="amplicon",
                    genome_alignments=genome_alignments,
                    genomes_dir=genomes_dir,
                    feature_table=feature_table,
                    rep_seqs=rep_seqs,
                    gtdb_matches=gtdb_matches,
                    top_k=top_k,
                    container=container,
                    execute=execute,
                    normalize=normalize,
                    verbose=verbose,
                    execution_profile=execution_profile,
                    dependencies=downstream_dependencies,
                )
                cod_profile = downstream_result["cod_profile"]
                result["artifacts"].update(downstream_result["artifacts"])

        elif mode == "amplicon":
            if feature_table is None or rep_seqs is None:
                raise ValueError("feature_table and rep_seqs are required for amplicon mode")
            feature_abundances = self.extract_relative_abundances(str(feature_table), sample_names=[sample_name], top_k=top_k)[sample_name]
            result["artifacts"]["feature_abundances"] = self._write_tall_mapping(output_path / "feature_abundances.csv", feature_abundances, sample_name=sample_name, key_name="feature_id", value_name="abundance")
            sample_repseqs = self._write_sample_repseqs(rep_seqs, feature_abundances, scratch_path / "sample_repseqs.fasta")
            result["artifacts"]["sample_repseqs"] = sample_repseqs

            matches_path = pathlib.Path(gtdb_matches) if gtdb_matches else scratch_path / "matches.blast"
            if gtdb_matches is None:
                logger.info("Aligning sample representative sequences to GTDB")
                step_name = "align_to_gtdb"
                step_container = self._step_container(execution_profile, step_name, container)
                old_settings = self._apply_step_config_settings(execution_profile, step_name)
                try:
                    script = self.align_to_gtdb(
                        sample_repseqs,
                        str(scratch_path),
                        container=step_container,
                        image=self._step_image(execution_profile, step_name, step_container),
                    )[0]
                finally:
                    self._restore_config_settings(old_settings)
                result["artifacts"][step_name] = self._execute_step(
                    script,
                    step_name=step_name,
                    sample_name=sample_name,
                    output_dir=scratch_path,
                    logger=logger,
                    execute=execute,
                    execution_profile=execution_profile,
                    dependencies=dependencies,
                )
            if not matches_path.exists():
                logger.info("GTDB matches are not available yet: %s", matches_path)
                cod_profile = {}
                result["status"] = "waiting_for_gtdb_alignment"
            else:
                representative_genomes = self.get_genomes_from_gtdb_alignment(str(matches_path))
                genome_abund = {}
                for feature, abundance in feature_abundances.items():
                    genome = representative_genomes.get(feature)
                    if genome:
                        genome_abund[genome] = genome_abund.get(genome, 0.0) + float(abundance)
                result["artifacts"]["representative_genomes"] = self._write_tall_mapping(output_path / "representative_genomes.csv", representative_genomes, sample_name=sample_name, key_name="feature_id", value_name="genome_id", value_dtype=pl.Utf8)
                result["artifacts"]["genome_abundances"] = self._write_tall_mapping(output_path / "genome_abundances.csv", genome_abund, sample_name=sample_name, key_name="genome_id", value_name="abundance")

                if genome_alignments:
                    alignments = genome_alignments if isinstance(genome_alignments, dict) else self._alignment_files_from_path(genome_alignments)
                elif genomes_dir:
                    genome_files = self._genome_files_from_dir(genomes_dir)
                    alignments = {}
                    commands = {}
                    alignment_dir = scratch_path / "genome_alignments"
                    alignment_dir.mkdir(parents=True, exist_ok=True)
                    missing_genome_fastas = []
                    genome_alignment_dependencies = (
                        [result["artifacts"]["align_to_gtdb"]]
                        if result["artifacts"].get("align_to_gtdb")
                        else None
                    )
                    for genome in genome_abund:
                        if genome not in genome_files:
                            logger.info("No genome FASTA found for %s", genome)
                            missing_genome_fastas.append(genome)
                            continue
                        step_name = f"align_genome_{genome}"
                        profile_step_name = "align_genome"
                        step_container = self._step_container(execution_profile, profile_step_name, container)
                        old_settings = self._apply_step_config_settings(execution_profile, profile_step_name)
                        try:
                            script, alignment = self.align_genome_to_protein_db(
                                genome_files[genome],
                                str(alignment_dir),
                                genome,
                                container=step_container,
                                image=self._step_image(execution_profile, profile_step_name, step_container),
                            )
                        finally:
                            self._restore_config_settings(old_settings)
                        commands[genome] = script
                        result["artifacts"][step_name] = self._execute_step(
                            script,
                            step_name=step_name,
                            sample_name=sample_name,
                            output_dir=scratch_path,
                            logger=logger,
                            execute=execute,
                            execution_profile={
                                **execution_profile,
                                "steps": {
                                    **execution_profile.get("steps", {}),
                                    step_name: self._step_settings(execution_profile, profile_step_name),
                                },
                            },
                            dependencies=genome_alignment_dependencies,
                        )
                        alignments[genome] = alignment
                    result["artifacts"]["genome_alignment_scripts"] = self._write_json(scratch_path / "genome_alignment_commands.json", commands)
                    if missing_genome_fastas:
                        result["artifacts"]["missing_genome_fastas"] = missing_genome_fastas
                else:
                    raise ValueError("amplicon mode requires genome_alignments or genomes_dir after GTDB matching")

                genome_cods = {
                    genome: self.cod_from_alignment(alignment, normalize=normalize)
                    for genome, alignment in alignments.items()
                    if genome in genome_abund and pathlib.Path(alignment).exists()
                }
                missing_alignments = [
                    genome
                    for genome, alignment in alignments.items()
                    if genome in genome_abund and not pathlib.Path(alignment).exists()
                ]
                if missing_alignments and not genome_cods:
                    logger.info(
                        "Genome alignment outputs are not available yet for %s genome(s)",
                        len(missing_alignments),
                    )
                    cod_profile = {}
                    result["status"] = "waiting_for_genome_alignment"
                    result["artifacts"]["missing_genome_alignments"] = missing_alignments
                elif result["artifacts"].get("missing_genome_fastas") and not genome_cods:
                    logger.info(
                        "Genome FASTA files are missing for %s genome(s)",
                        len(result["artifacts"]["missing_genome_fastas"]),
                    )
                    cod_profile = {}
                    result["status"] = "waiting_for_genome_fasta"
                else:
                    cod_profile = self.aggregate_genome_cod(genome_cods, genome_abund, normalize=normalize) if genome_cods else {}
                    result["artifacts"]["genome_cods"] = self._write_tall_nested_profile(output_path / "genome_cods.csv", genome_cods, sample_name=sample_name, entity_name="genome_id")
                    result["artifacts"]["cod_profile"] = self._write_tall_mapping(output_path / "cod_profile.csv", cod_profile, sample_name=sample_name, key_name="group", value_name="value")

        else:
            raise ValueError("mode must be one of: shotgun-alignment, shotgun-reads, genome-alignments, amplicon, amplicon-reads")

        result["cod_profile"] = cod_profile
        result["artifacts"]["provenance"] = self._write_json(
            output_path / "provenance.json",
            {
                "sample_name": sample_name,
                "mode": mode,
                "execute": execute,
                "container": container,
                "normalize": normalize,
                "reaction_db": self.config.csv_reaction_db,
                "protein_db": self.config.protein_db,
                "bit_score": self.config.bit_score,
                "e_value": self.config.e_value,
                "artifacts": result["artifacts"],
            },
        )
        if result.get("status", "").startswith("waiting_for_"):
            logger.info("Sample COD pipeline is waiting for submitted outputs: %s (%s)", sample_name, result["status"])
        else:
            logger.info("Finished sample COD pipeline: %s", sample_name)
        return result

    def batch_sample_to_cod(
        self,
        *,
        manifest: str | os.PathLike,
        output_dir: str | os.PathLike,
        input_type: str | None = None,
        sra_dir: str | os.PathLike | None = None,
        stage: str = "all",
        amplicon_to_genome_db: str | os.PathLike | None = None,
        genome_alignments: str | os.PathLike | dict[str, str] | None = None,
        genomes_dir: str | os.PathLike | None = None,
        gtdb_matches_dir: str | os.PathLike | None = None,
        forward_primer: str | None = None,
        reverse_primer: str | None = None,
        adapter_1: str | None = None,
        adapter_2: str | None = None,
        minimum_length: int = 100,
        quality_cutoff: str | int | None = None,
        quality_maxee: float = 1.0,
        identity: float = 0.97,
        min_unique_size: int = 2,
        chimera_filter: bool = True,
        top_k: int = -1,
        container: str = "None",
        execute: bool = False,
        normalize: bool = True,
        verbose: bool = True,
        execution_profile: str | os.PathLike | dict | None = None,
    ) -> dict:
        """Run a manifest-driven batch of amplicon samples to COD artifacts.

        Manifest rows may provide either an SRA ``accession`` or FASTQ paths
        in ``read_1``/``read_2``. The ``stage`` argument can be ``download``,
        ``preprocess``, ``cod``, or ``all``.
        """
        if stage not in {"download", "preprocess", "cod", "all"}:
            raise ValueError("stage must be one of: download, preprocess, cod, all")
        if input_type not in {None, "sra", "reads"}:
            raise ValueError("input_type must be 'sra' or 'reads'")
        if amplicon_to_genome_db is not None:
            self.config.amplicon2genome_db = str(amplicon_to_genome_db)
            matches = list(pathlib.Path(amplicon_to_genome_db).rglob(self.config.gtdb_dir))
            self.config.gtdb_dir_fasta = str(matches[0]) if matches else self.config.gtdb_dir_fasta

        output_path = pathlib.Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        rows = self._read_sample_manifest(manifest)
        workflow = MetagenomicsWorkflowState(output_path)
        execution_profile_data = self._load_execution_profile(execution_profile)
        results = {
            "manifest": str(manifest),
            "output_dir": str(output_path),
            "stage": stage,
            "execute": execute,
            "workflow_state": str(workflow.state_path),
            "workflow_events": str(workflow.events_path),
            "samples": {},
        }

        def _csv_has_rows(path: str | os.PathLike) -> bool:
            path = pathlib.Path(path)
            if not path.exists():
                return False
            try:
                return pl.read_csv(path, infer_schema_length=0).height > 0
            except Exception:
                return False

        def _artifact_status(artifacts: Iterable[dict]) -> str:
            statuses = {
                str(artifact.get("status", "prepared"))
                for artifact in artifacts
                if isinstance(artifact, dict)
            }
            if "failed" in statuses:
                return "failed"
            if statuses & {"submitted", "monitoring", "running"}:
                return "submitted"
            if statuses == {"completed"}:
                return "completed"
            return "prepared"

        for row in rows:
            accession = self._row_value(row, "accession", "sra", "run")
            read_1 = self._row_value(row, "read_1", "read1", "forward", "fastq_1", "fastq1")
            read_2 = self._row_value(row, "read_2", "read2", "reverse", "fastq_2", "fastq2")
            if input_type == "sra" and not accession:
                raise ValueError(f"Sample row {row} needs an accession when input_type='sra'")
            if input_type == "reads" and not read_1:
                raise ValueError(f"Sample row {row} needs read_1 when input_type='reads'")
            sample_name = self._row_value(row, "sample_name", "sample", "name")
            if sample_name is None:
                if accession:
                    sample_name = accession
                elif read_1:
                    sample_name = pathlib.Path(read_1).name.split(".")[0]
                else:
                    raise ValueError("Each manifest row must include sample/sample_name, accession, or read_1")
            paired = self._row_bool(row, "paired", default=read_2 is not None or accession is not None)
            sample_dir = output_path / sample_name
            scratch_dir = sample_dir / "scratch"
            preprocess_dir = scratch_dir / "amplicon_preprocess"
            feature_table = self._row_value(row, "feature_table", "feature-table") or str(preprocess_dir / "feature-table.tsv")
            rep_seqs = self._row_value(row, "rep_seqs", "rep-seqs", "representative_sequences") or str(preprocess_dir / "rep-seqs.fasta")
            cod_profile_path = sample_dir / "cod_profile.csv"
            sample_result = {"input": dict(row), "artifacts": {}, "stages": {}}
            sample_dependencies: list[dict] = []
            logger = self._sample_pipeline_logger(sample_name, sample_dir, verbose=verbose)
            slurm_checker = PipelineTaskManager(
                sample_name=sample_name,
                output_dir=sample_dir,
                execution_profile=execution_profile_data,
                logger=logger,
            )

            def record_stage(stage_name: str, status_value: str, *, artifact: dict | None = None, message: str | None = None, paths: dict | None = None) -> dict:
                entry = workflow.record(
                    sample_name,
                    stage_name,
                    status_value,
                    artifact=artifact,
                    message=message,
                    paths=paths,
                )
                sample_result["stages"][stage_name] = entry
                return entry

            def active(stage_name: str) -> bool:
                return workflow.active_submission(sample_name, stage_name, slurm_checker=slurm_checker)

            if accession and not read_1:
                target_sra_dir = pathlib.Path(sra_dir) if sra_dir else output_path / "sra"
                try:
                    resolved_reads = self._resolved_sra_reads(accession, target_sra_dir, paired=paired)
                    read_1 = resolved_reads["read_1"]
                    read_2 = resolved_reads["read_2"]
                    record_stage(
                        "download_sra",
                        "completed",
                        message="cached FASTQ files found",
                        paths=resolved_reads,
                    )
                except FileNotFoundError:
                    pass

            if accession and not read_1 and stage in {"download", "preprocess", "all"}:
                if active("download_sra"):
                    record_stage("download_sra", "submitted", message="download is already active")
                    sample_result["status"] = "waiting_for_download"
                    results["samples"][sample_name] = sample_result
                    continue
                else:
                    download_result = self.run_sra_download_step(
                        sample_name=sample_name,
                        output_dir=output_path,
                        accession=accession,
                        sra_dir=sra_dir,
                        paired=paired,
                        container=container,
                        execute=execute,
                        verbose=verbose,
                        execution_profile=execution_profile_data,
                    )
                    sample_result["artifacts"]["download"] = download_result["artifacts"]
                    download_artifact = download_result["artifacts"]["download_sra"]
                    read_1 = download_result["artifacts"]["reads"]["read_1"]
                    read_2 = download_result["artifacts"]["reads"]["read_2"]
                    sample_dependencies = [download_artifact]
                    record_stage(
                        "download_sra",
                        download_artifact.get("status", "prepared"),
                        artifact=download_artifact,
                        paths=download_result["artifacts"]["reads"],
                    )
                    if download_artifact.get("status") in {"submitted", "prepared", "running", "monitoring"} and not pathlib.Path(read_1).exists():
                        sample_result["status"] = "waiting_for_download"
                        results["samples"][sample_name] = sample_result
                        continue

            if stage == "download":
                sample_result.setdefault("status", "completed" if read_1 else "submitted")
                results["samples"][sample_name] = sample_result
                continue

            preprocess_ready = pathlib.Path(feature_table).exists() and pathlib.Path(rep_seqs).exists()
            if stage in {"preprocess", "all"}:
                if preprocess_ready:
                    record_stage(
                        "preprocess",
                        "completed",
                        message="cached feature table and representative sequences found",
                        paths={"feature_table": feature_table, "rep_seqs": rep_seqs},
                    )
                elif active("preprocess"):
                    record_stage("preprocess", "submitted", message="preprocessing is already active")
                    sample_result["status"] = "waiting_for_preprocess"
                    results["samples"][sample_name] = sample_result
                    continue
                elif read_1 and (not execute or pathlib.Path(read_1).exists()):
                    preprocess_result = self.preprocess_amplicon_sample(
                        sample_name=sample_name,
                        output_dir=output_path,
                        read_1=read_1,
                        read_2=read_2,
                        forward_primer=forward_primer,
                        reverse_primer=reverse_primer,
                        adapter_1=adapter_1,
                        adapter_2=adapter_2,
                        minimum_length=minimum_length,
                        quality_cutoff=quality_cutoff,
                        quality_maxee=quality_maxee,
                        identity=identity,
                        min_unique_size=min_unique_size,
                        chimera_filter=chimera_filter,
                        container=container,
                        execute=execute,
                        verbose=verbose,
                        execution_profile=execution_profile_data,
                        dependencies=sample_dependencies,
                    )
                    sample_result["artifacts"]["preprocess"] = preprocess_result["artifacts"]
                    preprocess_artifacts = [
                        preprocess_result["artifacts"].get("trim_reads"),
                        preprocess_result["artifacts"].get("build_amplicon_features"),
                    ]
                    preprocess_status = _artifact_status(preprocess_artifacts)
                    record_stage(
                        "preprocess",
                        preprocess_status,
                        artifact={"steps": preprocess_artifacts},
                        paths={"feature_table": feature_table, "rep_seqs": rep_seqs},
                    )
                    if not pathlib.Path(feature_table).exists() or not pathlib.Path(rep_seqs).exists():
                        sample_result["status"] = "waiting_for_preprocess"
                        results["samples"][sample_name] = sample_result
                        continue
                    preprocess_ready = True
                else:
                    if input_type == "reads":
                        raise ValueError(f"Sample {sample_name} needs an existing read_1 file before preprocessing")
                    record_stage("download_sra", "waiting", message="FASTQ files are not available yet")
                    sample_result["status"] = "waiting_for_download"
                    results["samples"][sample_name] = sample_result
                    continue

            if stage == "preprocess":
                sample_result.setdefault("status", "completed" if preprocess_ready else "waiting_for_preprocess")
                results["samples"][sample_name] = sample_result
                continue

            row_gtdb_matches = self._row_value(row, "gtdb_matches", "matches")
            if row_gtdb_matches is None and gtdb_matches_dir is not None:
                candidate = pathlib.Path(gtdb_matches_dir) / sample_name / "matches.blast"
                if candidate.exists():
                    row_gtdb_matches = str(candidate)
            if row_gtdb_matches is None:
                candidate = scratch_dir / "matches.blast"
                if candidate.exists():
                    row_gtdb_matches = str(candidate)

            if not pathlib.Path(feature_table).exists() or not pathlib.Path(rep_seqs).exists():
                record_stage(
                    "preprocess",
                    "waiting",
                    message="feature table or representative sequences are not available yet",
                    paths={"feature_table": feature_table, "rep_seqs": rep_seqs},
                )
                sample_result["status"] = "waiting_for_preprocess"
                results["samples"][sample_name] = sample_result
                continue

            if _csv_has_rows(cod_profile_path):
                record_stage("cod", "completed", message="cached COD profile found", paths={"cod_profile": str(cod_profile_path)})
                sample_result["status"] = "completed"
                sample_result["artifacts"]["cod"] = {"cod_profile": str(cod_profile_path)}
                results["samples"][sample_name] = sample_result
                continue

            if row_gtdb_matches is None and active("align_to_gtdb"):
                record_stage("align_to_gtdb", "submitted", message="GTDB alignment is already active")
                sample_result["status"] = "waiting_for_gtdb_alignment"
                results["samples"][sample_name] = sample_result
                continue

            if row_gtdb_matches is not None and active("align_genome"):
                record_stage("align_genome", "submitted", message="genome alignments are already active")
                sample_result["status"] = "waiting_for_genome_alignment"
                results["samples"][sample_name] = sample_result
                continue

            if stage == "cod":
                cod_result = self.sample_to_cod(
                    sample_name=sample_name,
                    output_dir=output_path,
                    mode="amplicon",
                    genome_alignments=genome_alignments,
                    genomes_dir=genomes_dir,
                    feature_table=feature_table,
                    rep_seqs=rep_seqs,
                    gtdb_matches=row_gtdb_matches,
                    top_k=top_k,
                    container=container,
                    execute=execute,
                    normalize=normalize,
                    verbose=verbose,
                    execution_profile=execution_profile_data,
                    dependencies=sample_dependencies,
                )
            else:
                cod_result = self.sample_to_cod(
                    sample_name=sample_name,
                    output_dir=output_path,
                    mode="amplicon",
                    genome_alignments=genome_alignments,
                    genomes_dir=genomes_dir,
                    feature_table=feature_table,
                    rep_seqs=rep_seqs,
                    gtdb_matches=row_gtdb_matches,
                    top_k=top_k,
                    container=container,
                    execute=execute,
                    normalize=normalize,
                    verbose=verbose,
                    execution_profile=execution_profile_data,
                    dependencies=sample_dependencies,
                )
            sample_result["artifacts"]["cod"] = cod_result["artifacts"]
            sample_result["cod_profile"] = cod_result.get("cod_profile", {})
            if cod_result["artifacts"].get("align_to_gtdb"):
                gtdb_artifact = cod_result["artifacts"]["align_to_gtdb"]
                record_stage(
                    "align_to_gtdb",
                    gtdb_artifact.get("status", "prepared"),
                    artifact=gtdb_artifact,
                    paths={"matches": str(scratch_dir / "matches.blast")},
                )
            genome_artifacts = [
                artifact
                for name, artifact in cod_result["artifacts"].items()
                if str(name).startswith("align_genome_") and isinstance(artifact, dict)
            ]
            if genome_artifacts:
                record_stage("align_genome", _artifact_status(genome_artifacts), artifact={"steps": genome_artifacts})
            if cod_result.get("status"):
                sample_result["status"] = cod_result["status"]
                waiting_stage = {
                    "waiting_for_gtdb_alignment": "align_to_gtdb",
                    "waiting_for_genome_alignment": "align_genome",
                    "waiting_for_genome_fasta": "genome_fasta",
                }.get(cod_result["status"], "cod")
                current_stage = workflow.stage(sample_name, waiting_stage)
                if current_stage.get("status") in {"submitted", "running", "monitoring"}:
                    sample_result["stages"][waiting_stage] = current_stage
                else:
                    record_stage(waiting_stage, "waiting", message=cod_result["status"])
            elif _csv_has_rows(cod_profile_path):
                sample_result["status"] = "completed"
                record_stage("cod", "completed", paths={"cod_profile": str(cod_profile_path)})
            else:
                sample_result["status"] = "waiting_for_cod"
                record_stage("cod", "waiting", message="COD profile was not produced")
            results["samples"][sample_name] = sample_result

        results["summary"] = self._write_json(output_path / "batch_summary.json", results)
        return results
    
    def calculate_group_abundances(self,elements_feature_abundances:dict[str,dict],rel_abund:dict[str,dict])->dict[str,dict[str,float]]:
        """
        This method is defined to calculate the features for each sample given:
        1) The relative abundances of the genomes in each sample:
            - In this dictionary the keys are the sample names and the values are dictionaries where the keys are the genome names and the values are the relative abundances of the genomes in the sample.
        2) The relative abundances of the elements in each genome.
            - In this dictionary the keys are the genome names and the values are dictionaries where the keys are the element names and the values are the relative abundances of the elements in the genome.

        Required Configs:
            None
        
        Args:
            elements_feature_abundances (dict[str,dict]): A dictionary containing the relative abundances of the elements in each genome.
            rel_abund (dict[str,dict]): A dictionary containing the relative abundances of the genomes in each sample.
        
        Returns:
            dict[str,dict[str,float]]: A dictionary containing the relative abundances of the elements in each sample.
        """
        out={}
        features = sorted({feature for abundances in elements_feature_abundances.values() for feature in abundances})
        for sample,abunds in rel_abund.items():
            weighted = {feature: 0.0 for feature in features}
            for element, abundance in abunds.items():
                for feature, value in elements_feature_abundances.get(element, {}).items():
                    weighted[feature] = weighted.get(feature, 0.0) + float(value) * float(abundance)
            out[sample]=scaler(pl.DataFrame([weighted])).to_dicts()[0]
        return out
    
    def extract_relative_abundances(self,feature_table_dir:str,sample_names:Union[list[str],None]=None,top_k:int=-1)->dict:
        
        r"""
        This method extracts the relative abundances of the features in each sample from a TSV feature table.
        VSEARCH OTU tables and exported BIOM-style TSV tables are both supported.
        NOTE: The final feature abundances sum to 1 for each sample.
        Required Configs:
            None
        Args:
            feature_table_dir (str): The path to the feature table.
            sample_names (Union[list[str],None], optional): The list of sample names. to be considered. If None, all the samples will be considered. Defaults to None.
            top_k (int, optional): The number of top features to be used. If -1, all the features will be used. Defaults to -1.

        Returns:
            dict: A dictionary containing the relative abundances of the features in each sample.
        """
        with open(feature_table_dir) as f:
            first_line = f.readline()
        skiprows = 1 if first_line.startswith("#") and not first_line.startswith("#OTU ID") else 0
        feature_table = pl.read_csv(
            feature_table_dir,
            separator="\t",
            skip_rows=skiprows,
            infer_schema_length=0,
        )
        if "#OTU ID" not in feature_table.columns:
            raise ValueError("Feature table must contain a '#OTU ID' column")
        if sample_names is None:
            sample_names = [column for column in feature_table.columns if column != "#OTU ID"]
        relative_abundances={sample:[] for sample in sample_names}
        if top_k == -1:
            top_k = feature_table.height
        for sample in sample_names:
            if sample not in feature_table.columns:
                raise ValueError(f"Sample {sample} not found in feature table")
            top_features = (
                feature_table
                .select([
                    pl.col("#OTU ID").cast(pl.Utf8).alias("feature_id"),
                    pl.col(sample).cast(pl.Float64, strict=False).fill_null(0.0).alias("abundance"),
                ])
                .sort("abundance", descending=True)
                .head(top_k)
            )
            total = top_features["abundance"].sum()
            if not total:
                relative_abundances[sample] = {}
                continue
            relative_abundances[sample] = {
                row["feature_id"]: float(row["abundance"]) / float(total)
                for row in top_features.to_dicts()
            }
        return relative_abundances
    
    def assign_ec_to_genome(self,alignment_file:str)->dict:
        r"""
        This function takes an alignment file and assigns the EC numbers to the genomes based on the alignment file,
        and the e-adm groupings of the EC numbers. The output is a dictionary where the keys e-adm reactions and the values are the EC numbers,
        that are found in the genome and are grouped under the e-adm reaction.
        
        Example: 
            >>> import os
            >>> output = os.path.join(Main_Dir, "test", "ec_to_genome")
            >>> alignments = ["CP001673.1", "Q8YNF9|1.4.4.2", "0.566", "2859", "414", "0", "1021521", "1024379", "27", "982", "0.000E+00", "1109"]
            >>> headers = ["query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend ", "tstart", "tend", "evalue", "bits"]
            >>> hit = "\t".join(alignments)
            >>> headers_tab = "\t".join(headers)
            >>> combine = headers_tab + "\n" + hit
            >>> with open(output, "w") as f:
            ...     f.write(combine)
            161
            >>> obj = Metagenomics(configs.Metagenomics())
            >>> obj.assign_ec_to_genome(output) 

        Args:
            alignment_file (str): The address of the alignment file.
            
        Returns:
            dict: A dictionary containing the e-adm reactions and the EC numbers that are found in the genome and are grouped under the e-adm reaction.
        """

        aligntable = pl.read_csv(alignment_file, separator="\t", infer_schema_length=0).with_columns([
            pl.col("bits").cast(pl.Float64, strict=False).fill_null(0.0),
            pl.col("evalue").cast(pl.Float64, strict=False).fill_null(float("inf")),
        ])
        aligntable = aligntable.filter((pl.col("bits") > self.config.bit_score) & (pl.col("evalue") < self.config.e_value))

        ec_align_list = (
            aligntable
            .select(pl.col("target").str.split_exact("|", 1).struct.field("field_1").alias("ec"))
            .drop_nulls()
            .unique()
            ["ec"]
            .to_list()
        )

        metadatatable = (
            pl.read_csv(self.config.csv_reaction_db, separator=",", infer_schema_length=0)
            .unique(subset=["EC_Numbers"], keep="first")
            .select(["EC_Numbers","Modified_ADM_Reactions"])
            .drop_nulls()
            .filter(pl.col("EC_Numbers").is_in(ec_align_list))
        )
        adm_reactions=sorted({
            reaction
            for row in metadatatable.to_dicts()
            for reaction in row["Modified_ADM_Reactions"].split("|")
        })
        adm_to_ecs={}
        for reaction in adm_reactions:
            adm_to_ecs[reaction]=[
                row["EC_Numbers"]
                for row in metadatatable.to_dicts()
                if reaction in row["Modified_ADM_Reactions"].split("|")
            ]
            
        return adm_to_ecs

    


    def seqs_from_sra(self,accession:str,target_dir:str,container:str="None",**kwargs)-> tuple[str,dict]:
        """ 
        This method downloads the fastq files from the SRA database using the accession number (ONLY SAMPLE ACCESSION AND NOT PROJECT ACCESSION) of the project or run.
        The method uses the fasterq-dump tool to download the fastq files. This method also extracts the sample metadata from the SRA database for future use.
        #NOTE In order for this method to work without any container, you need to have the SRA toolkit installed on your system or
        at least have prefetch and fasterq-dump installed on your system. For more information on how to install the SRA toolkit, please refer to the following link:
        https://github.com/ncbi/sra-tools.
        
        Any additional keyword
    
        Required Configs:
            None
        
        
        Args:
            accession (str): The accession number of the SRA project or run
            target_dir (str): The directory where the fastq files will be downloaded
            container (str, optional): The containerization tool that will be used to run the bash scripts. Defaults to "None". Options are "None","docker","singularity"
        
        Returns:
            prefetch_script (str): The bash script that will be used to download the SRA files in python string format
            sample_metadata (dict): A dictionary that contains the sample metadata
    
        """   
        if container=="None":
            prefetch_script=f"""prefetch {accession} -O {target_dir} --max-size 100000000\n"""
            acc_folder=pathlib.Path(target_dir)/accession
            fasterq_dump_script=""
            sra_file=acc_folder/(accession+".sra")
            fasterq_dump_script+=f"fasterq-dump {sra_file} -O {acc_folder} --split-files --temp {acc_folder}\n"
            fasterq_dump_script+=f"rm {sra_file}"
            prefetch_script+=fasterq_dump_script
 
        
        elif container=="docker":
            prefetch_script=""""""
            prefetch_script+=f"docker run -v {target_dir}:{target_dir} {self.config.adtoolbox_docker} prefetch {accession} -O {target_dir} --max-size 100000000\n"
            acc_folder=pathlib.Path(target_dir)/accession
            fasterq_dump_script=""
            sra_file=acc_folder/(accession+".sra")
            fasterq_dump_script+=f"docker run -v {target_dir}:{target_dir} {self.config.adtoolbox_docker} fasterq-dump {sra_file} -O {acc_folder} --split-files --temp {acc_folder}\n"
            fasterq_dump_script+=f"docker run -v {target_dir}:{target_dir} {self.config.adtoolbox_docker} rm {sra_file}"
            prefetch_script+=fasterq_dump_script
        
        elif container=="singularity":
            prefetch_script=""""""
            prefetch_script+=f"singularity exec {self.config.adtoolbox_singularity} prefetch {accession} -O {target_dir} --max-size 100000000\n"
            acc_folder=pathlib.Path(target_dir)/accession
            fasterq_dump_script=""
            sra_file=acc_folder/(accession+".sra")
            fasterq_dump_script+=f"singularity exec {self.config.adtoolbox_singularity} fasterq-dump {sra_file} -O {acc_folder} --split-files --temp {acc_folder}\n"
            fasterq_dump_script+=f"singularity exec {self.config.adtoolbox_singularity} rm {sra_file}"
            prefetch_script+=fasterq_dump_script
       
        sample_metadata=utils.get_sample_metadata_from_accession(accession)      
            
        
        return prefetch_script,sample_metadata     
            
    def merge_paired_sequences(self,read_1:str,read_2:str,outputfile:str,container:str="None",**kwargs)->tuple[str]:
        """ This method merges the paired end reads using the fastp tool. Note that any additional keyword arguments will be passed to the fastp tool.   
        Required Configs:
            None
        Args:
            read_1 (str): The directory of the forward reads file.
            read_2 (str): The directory of the reverse reads file.
            outputfile (str): The directory of the output file.
            container (str, optional): The containerization tool that will be used to run the bash scripts. Defaults to "None". Options are "None","docker","singularity".
        
        Returns:
            str: The bash script that will be used to merge the paired end reads in python string format
        """
        if container=="None":
            bash_script=f"""fastp -i {read_1} -I {read_2}  -m --merged_out {outputfile}"""
            for key,value in kwargs.items():
                key_=key.replace("_","-")
                bash_script+=f" --{key_} {value} "
                
        elif container=="docker":
            bash_script=f"""docker run -v {read_1}:{read_1} -v {read_2}:{read_2} -v {outputfile}:{outputfile} {self.config.adtoolbox_docker} fastp -i {read_1} -I {read_2}  -m --merged_out {outputfile}"""
            for key,value in kwargs.items():
                key_=key.replace("_","-")
                bash_script+=f" --{key_} {value} "
        
        elif container=="singularity":
            bash_script=f"""singularity exec {self.config.adtoolbox_singularity} fastp -i {read_1} -I {read_2}  -m --merged_out {outputfile}"""
            for key,value in kwargs.items():
                key_=key.replace("_","-")
                bash_script+=f" --{key_} {value} "
                    
        return bash_script,
    
class Annotation:
    
    def __init__(self,config:configs.Annotation):
        self.config=config
        
    def annotate_with_metacyc(self,
                                alignment_file:str,
                                )->dict[str,dict[str,str]]:
        """"
        This function annotates the genomes with the metabolic pathways using the metabolic pathways database.
        Required Configs:
            None
        Args:
            alignment_file (str): The path to the alignment file.
        Returns:
            dict: A dictionary containing the metabolic pathways  and their coverage in based on
            the input alignment file.
        """
        ### Creating the metacyc full dictionry from the protein database
        metacyc_full_dict={}
        annotation_dict={}
        with open(self.config.metacyc_protein_db, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    header = line[1:].strip()
                    pathway, reaction, _ = header.split('|')
                    if pathway not in metacyc_full_dict:
                        metacyc_full_dict[pathway] = set()
                    metacyc_full_dict[pathway].add(reaction)
        alignment_table=pl.read_csv(alignment_file, separator="\t", infer_schema_length=0)
        grouped: dict[str, dict[str, set[str]]] = {}
        for target in alignment_table["target"].to_list():
            pathway, reaction, *_ = target.split("|")
            grouped.setdefault(pathway, {"reactions": set(), "extras": set()})
            grouped[pathway]["reactions"].add(reaction)
            if _:
                grouped[pathway]["extras"].add(_[0])
        for pathway, values in grouped.items():
            annotation_dict.setdefault(pathway,{})["reactions"]=values["reactions"]
            annotation_dict.setdefault(pathway,{})["coverage"]=len(values["reactions"])/len(metacyc_full_dict[pathway])
            annotation_dict.setdefault(pathway,{})["all_reactions"]=metacyc_full_dict[pathway]
            
        return annotation_dict
            

        


if __name__ == "__main__":
    annot_conf=configs.Annotation()
    annot=Annotation(annot_conf)
    annot.annotate_with_metacyc("/Users/parsaghadermarzi/Downloads/sach.m8")

    
    
