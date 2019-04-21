import numpy as np
import pandas as pd
import pickle
from polyBootstrap import Population, PopulationSet
from neiGeneticDistance import *
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
import pdb
import sys

def bootstrap(ps, jobID, basename='majorityTree', treebuilder='nj', bootstraps=1, outgroup=None):
    """treebuilder could be nj/upgma, outgroup: a population name or 'midpoint'"""
    allLoci = set()
    for pop in ps.populations:
        allLoci = allLoci.union(pop.allpolySites)
    allLoci = list(allLoci)  ## sort it? Reduce to independent sites (www.pnas.org/content/93/23/13429, run LD?)

    sites = len(allLoci)
    trees = []
    for bootstrap in range(bootstraps):
        selectedLoci0 = np.random.choice(range(len(allLoci)), sites, replace=True)
        selectedLoci = [allLoci[l] for l in selectedLoci0]
        df = pd.DataFrame([pop.bootstrap(selectedLoci) for pop in ps.populations], index=ps.popnames)
        dmNei = neiDF(df, [5] * (sites - 1))
        ## annoying conversion, BioPython couldnt be just more compatible with scipy/pdist?
        dmTriangular = [list(dmNei[i, :(i + 1)]) for i in range(len(dmNei))]
        try:
            m = _DistanceMatrix(ps.popnames, dmTriangular)
        except ValueError:
            pdb.set_trace()
        constructor = DistanceTreeConstructor()  # could've passed treebuilder here too
        tree = getattr(constructor, treebuilder)(m)
        if outgroup == 'midpoint':
            tree.root_at_midpoint()
        elif not outgroup is None:
            tree.root_with_outgroup({'name': outgroup})
        trees.append(tree)
    filename = '../Data/%s_%s_%s_%s_%04d.pcl' % (basename, bootstraps, treebuilder, len(ps.populations), jobID)
    with open(filename, "wb") as w:
        pickle.dump(trees, w)


jobID = int(sys.argv[-1])
basename = sys.argv[-2]
popfile = sys.argv[-3] # popsAB.pcl

with open(popfile, "rb") as f:
    ps = PopulationSet(pickle.load(f))

bootstrap(ps, jobID, basename=basename, treebuilder='nj', outgroup='midpoint')