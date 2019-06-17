"Generating Violin plots of confidence distributions from branch support values in bootstrapped trees"

import glob
from Bio import Phylo
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

#    tf = '/Users/ahenschel/Github/FreqRT/Results/Trees2/pops_A-B-C-DPB1-DQB1-DRB1_min50_1_upgma_14_cut0.nwk'

def getConfidence(clade):
    """little recursive function traversing a tree (starting from tree.root, 
    which is of type Clade, and so are all children)"""
    confidence = []
    if not clade.confidence is None:
        confidence.append(clade.confidence)
    for child in clade.clades:
        confidence += getConfidence(child)
    return confidence

class PopInfo:
    def __init__(self, filename):
        tree = next(Phylo.parse(filename,'newick'))   
        self.data = np.array(getConfidence(tree.root))
        _, genes, _,_, self.method, pops = os.path.basename(filename).split('_')[:6]
        self.pops = int(pops)
        self.genes = genes.split('-')
        #self.repr = len(self.genes), self.method, self.pops 
        self.repr = f'{genes}/{self.pops}' 
        
dfdata = []
pis = [PopInfo(tf) for tf in glob.glob(f'../Results/Trees2/*_cut0.nwk')]
pis.sort(key=lambda p: (p.pops, p.method))
datas = [pi.data for pi in pis]
labels = [pi.repr for pi in pis]
for pi in pis:
    dfdata += [[confidence, pi.repr, pi.method] for confidence in pi.data]
df1 = pd.DataFrame(dfdata, columns=['confidence', 'genes', 'method'])
ax = sns.violinplot(x='genes', y='confidence', data=df1, cut=0, width=1., hue='method', split=True) 
ax.set_xticklabels(labels[::2], rotation=90)
""" parts = plt.violinplot(datas, range(len(datas)))
for idx, pc in enumerate(parts['bodies']):
    color = '#D43F3F' if idx%2 else '#3F3FD4'
    pc.set_edgecolor(color)
    pc.set_color(color)
    pc.set_alpha(0.7)

plt.xticks(range(len(labels)), labels, rotation='vertical')
plt.show() """