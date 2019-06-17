from alleleFreqAPI import PopulationTree
from importlib import reload; import polyBootstrap; reload(polyBootstrap)
from polyBootstrap import Gene, Poly, Population, PopulationSet
import pandas as pd
import pickle

ptree = PopulationTree()
ptree.fix()

with open('../Data/pops_A-B-DQB1_min50.pcl', 'rb') as f:
    populations = pickle.load(f)
ps = PopulationSet(populations)
ps.bootstrap(useAllLoci=True, basename='dmtest', bootstraps=1)
ps.visualizeDM()

