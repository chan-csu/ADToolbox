class Pair(Thread):
    def __init__(self,u,v):
        super().__init__()
        self.u=u
        self.v=v
    def run(self):
        self.res=distance.braycurtis(self.u,self.v)
def calculate_dist(df:pd.DataFrame,feature_axis:int=1,threads:int=8)->pd.DataFrame:
    res=pd.DataFrame(index=df.index,columns=df.index)
    pairs=[]
    tasks=[]
    for i in range(df.shape[0]):
        res.iloc[i,i]=0
        for j in range(i+1,df.shape[0]):
            pairs.append((i,j))
            tasks.append(Pair(df.iloc[i,:],df.iloc[j,:]))
    for task in tasks:
        task.start()
    for idx,task in enumerate(tasks):    
        task.join()
        i1,j1=pairs[idx]
        res.iloc[i1,j1]=task.res
        res.iloc[j1,i1]=task.res
    
    return res
        