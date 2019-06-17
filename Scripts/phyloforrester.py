from Bio import Phylo
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os

texTableHead="""\\begin{table}[H]
\\caption{Average bootstrap confidence for the best gene combination, 
thresholded by confidence cutoff %s\\%%.}
\\centering
\\label{tab:genecombinations} 
"""
#plt.figure(num=None, figsize=(18, 12), dpi=300, facecolor='w', edgecolor='k')
def infoFromName(fn):
    f = os.path.splitext(os.path.basename(fn))[0].split('_')
    return (f[1], f[1].count('-')+1, int(f[6][3:]), f[4], f[5])
def parse(clade):
    confidences = []
    for c in clade.clades:
        #if not c.confidence is None:
        confidences.append(c.confidence)
        confidences += parse(c)
    return confidences

texDir = '/Users/ahenschel/Github/FreqRT/ManuscriptBootstrapping/Tables'
nwkFiles = [n for n in sys.argv[1:] if n.endswith('.nwk')]
if not nwkFiles:
    print("Provide newick files through command line, e.g.: Data/*.nwk")
else:
    data = []
    for nwkFile in nwkFiles: 
        tree = Phylo.read(nwkFile, format='newick')
        rec = infoFromName(nwkFile)
        confidences = np.array(parse(tree.root), dtype=np.float)
        nrNans = np.sum(np.isnan(confidences))
        nrNansPct = nrNans/len(confidences)
        rec += (np.nanmean(confidences), np.nanmedian(confidences), nrNans, nrNansPct)
        data.append(rec)
        if True:
            plt.clf()
            Phylo.draw(tree)
            plt.savefig(f'{nwkFile[:-4]}.svg', figsize=(18, 12), dpi=300)
            plt.close()
            #input("Press enter")
        #break
    columns = 'Loci NrLoci Thresh Build Pops AvgConf MedConf Miss MissPct'.split()
    df = pd.DataFrame(data, columns=columns)
    
if False:
    with open(f'{texDir}/treeConfidenceStatsMissSort.tex', 'w') as tex: 
        for threshold in [0,25,50,75]:
            df0 = df[df.Thresh==threshold].sort_values('MissPct', ascending=True)
            df0 = df0[columns[:2]+columns[3:]]
            tex.write(texTableHead % threshold)
            df0.to_latex(tex, index=False, na_rep='-', float_format='%.2f') 
            tex.write('\\end{table}\n\n\n')
#pops_A-B-DQA1-DQB1-DRB1_min50_1_upgma_26_cut75.nwk