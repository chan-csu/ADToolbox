"""
This module provides statistical functions for ADToolbox.
"""
import pandas as pd
from threading import Thread,Semaphore
from scipy.spatial import distance
from typing import Iterable
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

def calculate_dist(df:pd.DataFrame,feature_axis:int=1,threads:int=8)->pd.DataFrame:
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
        df=df.T
    sem=Semaphore(threads)
    res=pd.DataFrame(index=df.index,columns=df.index)
    pairs=[]
    tasks=[]
    for i in range(df.shape[0]):
        res.iloc[i,i]=0
        for j in range(i+1,df.shape[0]):
            pairs.append((i,j))
            tasks.append(Pair(df.iloc[i,:],df.iloc[j,:],sem))
    for task in tasks:
        task.start()
    for idx,task in enumerate(tasks):    
        task.join()
        i1,j1=pairs[idx]
        res.iloc[i1,j1]=task.res
        res.iloc[j1,i1]=task.res
    
    return res

def scaler(df:pd.DataFrame,feature_axis:int=1,inplace:bool=False)->pd.DataFrame|None:
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
        df=df.T
    if inplace:
        df=df.apply(lambda x:x/x.sum(),axis=1)
        return None
    else:
        return df.apply(lambda x:x/x.sum(),axis=1)
    

