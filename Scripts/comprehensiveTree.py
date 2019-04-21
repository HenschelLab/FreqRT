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

## TODO: read this from command line args
selectedGenes = ['A', 'B', 'DQB1']
minPopSize = 95 # so to include UAE

with open('../Data/hlaABPoly.pcl', 'rb') as f: genes = pickle.load(f) ## TODO: precalc other genes!
for gene in selectedGenes:
    if not gene in genes:
        genes[gene] = Gene(gene=gene)

genes['DQB1'] = Gene(gene='DQB1')

ptree = PopulationTree()
ptree.fix()
ptree.findPopsWithData(selectedGenes, minPopSize=minPopSize)

pops = [Population(pdict, genes, pname) for (pname, pdict) in ptree.makeAFdicts(selectedGenes)] ## CHANGE!!!

ps = PopulationSet(pops)
filename = "../Data/pops_%s_min%s.pcl" % ("-".join(selectedGenes), minPopSize)
with open(filename, "wb") as w:
    pickle.dump(ps.populations,w)
print("Wrote %s" % filename)

#Often to slow to run on a single machine. dump population data as pickle
#see parallelBootstrap.py, which is started by the bsub script runAllBootstraps.sh
#TODO: gather all trees using
#ps.bootstrap(bootstraps=100, basename='maj_A', treebuilder='nj', outgroup='midpoint') ## run parallel


'''
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