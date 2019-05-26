import numpy as np
from scipy.spatial.distance import pdist, squareform
import pandas as pd

def dotsum(X,Y):
    return sum(np.dot(x,y) for x, y in zip(X,Y))

def nei(X, Y): ##for speedup: memorize dotsum(X,X)
    return -np.log(dotsum(X, Y)/np.sqrt(dotsum(X, X) * dotsum(Y,Y)))

def nei2(X, Y, grouping):
    """expects 2 vectors and a list for nr columns per locus (except the last)"""
    return nei(np.split(X, np.cumsum(grouping)), np.split(Y, np.cumsum(grouping)))

def  neiDF(df, grouping, format=squareform):
    """handles a pandas dataframe with rows as vectors"""
    return format(pdist(df.values, lambda X, Y: nei2(X,Y,grouping)))

if __name__ == "__main__":
    ## 3 loci
    X = [.25, .6, .15,   .3, .7,   .5, .5]
    Y = [.24, .6, .16,   .2, .8,   .4, .6]
    Z = [.62, .0, .38,   .5, .5,   .1, .9]

    df = pd.DataFrame([X, Y, Z], index = list("XYZ"))
    print (neiDF(df, [3,2]))

    #print nei(X,Y)
    #print nei(X,Z)
    #print nei(Y,Z)



