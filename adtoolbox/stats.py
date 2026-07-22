"""
This module provides statistical functions for ADToolbox.
"""
from threading import Thread,Semaphore
from scipy.spatial import distance
from typing import Iterable
import polars as pl


def _row_names(df: pl.DataFrame) -> list[str]:
    return df.get_column(df.columns[0]).cast(pl.Utf8).to_list() if df.columns and df.columns[0] == "name" else [str(i) for i in range(df.height)]


def _numeric_rows(df: pl.DataFrame) -> tuple[list[str], list[list[float]]]:
    frame = df
    names = _row_names(frame)
    if frame.columns and frame.columns[0] == "name":
        frame = frame.drop("name")
    numeric = frame.select(pl.all().cast(pl.Float64, strict=False).fill_null(0.0))
    return names, numeric.rows()


class Pair(Thread):
    """
    This class is used to calculate the distance between two vectors. it is a subclass of Thread.
    This design is used to do multi-threading given that multiple pairs of vectors can be calculated simultaneously.
    """
    def __init__(self,u:Iterable,v:Iterable,semaphore:Semaphore):
        """
        Given two vectors, this function calculates the bray-curtis distance between them.
        Args:
            u: the first vector
            v: the second vector
            semaphore: the semaphore used to control the number of running threads (This is automatically handled by the calculate_dist function)
        
        """
        super().__init__()
        self.u=u
        self.v=v
        self.semaphore=semaphore
    def run(self):
        """
        This method is called to start the thread.
        """
        with self.semaphore:
            self.res=distance.braycurtis(self.u,self.v)

def calculate_dist(df:pl.DataFrame,feature_axis:int=1,threads:int=8)->pl.DataFrame:
    """
    This function calculates the bray-curtis distance between each pair of vectors in a dataframe.
    Args:
        df: the dataframe containing the vectors
        feature_axis: the axis of the dataframe that contains the vectors
        threads: the number of threads used to calculate the distances
    Returns:
        a dataframe containing the distances between each pair of vectors
    """
    if feature_axis==0:
        df=df.transpose(include_header=True, header_name="name")
    sem=Semaphore(threads)
    names, rows = _numeric_rows(df)
    res = [[0.0 for _ in names] for _ in names]
    pairs=[]
    tasks=[]
    for i in range(len(rows)):
        for j in range(i+1,len(rows)):
            pairs.append((i,j))
            tasks.append(Pair(rows[i],rows[j],sem))
    for task in tasks:
        task.start()
    for idx,task in enumerate(tasks):    
        task.join()
        i1,j1=pairs[idx]
        res[i1][j1]=task.res
        res[j1][i1]=task.res
    
    return pl.DataFrame([{"name": name, **{names[j]: row[j] for j in range(len(names))}} for name, row in zip(names, res)])

def scaler(df:pl.DataFrame,feature_axis:int=1,inplace:bool=False)->pl.DataFrame|None:
    """
    This function scales the vectors in a dataframe to have unit norm.
    Args:
        df (pd.DataFrame): the dataframe containing the vectors
        feature_axis (int): the axis of the dataframe that contains the vectors
        inplace (bool): whether to scale the vectors in place
    Returns:
        (pd.DataFrame|None): the scaled dataframe; if inplace is True, returns None
    """
    if feature_axis==0:
        df=df.transpose(include_header=True, header_name="name")
    columns = df.columns
    scaled_rows = []
    for row in df.select(pl.all().cast(pl.Float64, strict=False).fill_null(0.0)).to_dicts():
        total = sum(row.values())
        scaled_rows.append({key: (value / total if total else 0.0) for key, value in row.items()})
    if inplace:
        return None
    return pl.DataFrame(scaled_rows, schema={column: pl.Float64 for column in columns})
    
