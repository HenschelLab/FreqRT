import pickle
from Bio.Phylo import Consensus
from Bio import Phylo
import glob
import sys


"""
This tool creates a consensus tree from individual trees. 
It simply collects all pickled trees that were calculated by individual bootstrap runs of
the script parallelBootstrap.py, which dumps a list of trees that were created with BioPython's
DistanceTreeConstructor. 

Alternative consensus construction that BioPython provides:
 * strict_consensus
 * adam_consensus
They don't provide bootstrap support values, though.

#run on HPC
#source activate bio3
"""

datadir = "/research/gutsybugs/HLA/Data/Trees"

treebase = "maj_ABC_1_nj_92"
if len(sys.path)>2:
    treebase = sys.argv[-1]
#maj_AB_min95_1_nj_246_0001.pcl

treefiles =  glob.glob("%s/%s*.pcl" % (datadir, treebase))
trees = []
for treefile in treefiles:
    with open(treefile, "rb") as tf:
        trees += pickle.load(tf)

majorityTree = Consensus.majority_consensus(trees)

Phylo.write(majorityTree, '../Data/%s.nwk' % treebase, format='newick')

Phylo.draw_ascii(majorityTree)

