"""
What is the best tree we can derive from AlleleFrequencies.net data?

Input: a selection of MHC genes, gets allele frequencies from the scraped spreadsheet (see scrape*.py)
Output:
This script simply uses the functionality leveraged by polyBootstrap.py, where the heart of the calculation lies.
Here we only produce Population objects from the AlleleFrequencies spreadsheet,
compatible with polyBootstrap

"""

from polyBootstrap import Gene, Poly, Population, PopulationSet
from alleleFreqAPI import PopulationTree
import pandas as pd
import pickle
import numpy as np
allGenes = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
## TODO: read this from command line args
with open('../Data/hlaAllPoly.pcl', 'rb') as f: genes = pickle.load(f) 
for gene in allGenes:
    if not gene in genes:
        genes[gene] = Gene(gene=gene)

ptree = PopulationTree()
ptree.fix()

stripsum = lambda x: [e[4:] for e in x]
#selectedGenes = ['A', 'B', 'DRB1', 'DQB1']
minPopSize = 50 
for nrGenes in []:#range(2,7): CHAMGE!!!!
    print(nrGenes)
    ps = ptree.exploreGoodCombinations(nrGenes=nrGenes)[-3:]
    geneCombinations = [stripsum(com[1:]) for com in ps] 

    for selectedGenes in geneCombinations:
        ptree.findPopsWithData(selectedGenes, minPopSize=minPopSize)
        pops = [Population(pdict, genes, pname) for (pname, pdict) in ptree.makeAFdicts(selectedGenes)] ## CHANGE!!!

        ps = PopulationSet(pops)
        filename = "../Data/pops_%s_min%s.pcl" % ("-".join(selectedGenes), minPopSize)
        with open(filename, "wb") as w:
            pickle.dump( ps.populations,w)
        print("Wrote %s" % filename)

#Often to slow to run on a single machine. dump population data as pickle
#see parallelBootstrap.py, which is started by the bsub script runAllBootstraps.sh
#TODO: gather all trees using
#ps.bootstrap(bootstraps=100, basename='maj_A', treebuilder='nj', outgroup='midpoint') ## run parallel


'''
#Deprecated code, intended for single run
for popIdx in range(len(ptree.df)):
    row = ptree.df.iloc[popIdx]
    if type(row.PopName) != str or len(row.PopName)<2: continue

    try:
        alleleFreqs = {gene: {} for gene in selectedGenes}
        for c in starNames:
            gene = c.split('*')[0]
            if gene in selectedGenes:
                if c in genes[gene].polymorphSites:
                    alleleFreqs[gene][c] = row[c]# todo: maybe normalize percentages, AF not consistent
                elif row[c]>0.01:
                    print("Warning: No info on %s (%s)" % (c, row[c]))
        pops.append(Population(alleleFreqs, genes, name=row.PopName))
        print ("Successfully added stats for %s" % row.PopName)
    except:
        print ("something went wrong with %s" % row.PopName)

ps = PopulationSet(pops)
ps.bootstrap(bootstraps=100, basename='maj_A', treebuilder='nj', outgroup='midpoint') ## run parallel
'''