import pandas as pd
from collections import defaultdict, Counter
import os
import pdb
import numpy as np
import pickle
from neiGeneticDistance import *
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
from Bio.Phylo import Consensus
from Bio import Phylo
import io

#imgtDir = '/research/gutsybugs/HLA/Data/IMGT/alignments' ## adjust this!
imgtDir = '/Users/ahenschel/Github/FreqRT/Data/IMGT/alignments'
vcfDir = '../Data/'
resultDir = '../Results/Trees'

allGenes = ['A', 'B', 'C', 'DMA', 'DMB', 'DOA', 'DOB', 'DPA1', 'DPA2', 'DPB1',
            'DPB2', 'DQA1', 'DQB1', 'DRA', 'DRB1', 'DRB3', 'DRB4', 'E', 'F', 'G', 'H',
            'HFE', 'J', 'K', 'L', 'MICA', 'MICB', 'P', 'T', 'TAP1', 'TAP2', 'V', 'W', 'Y']
bases = list('ACGT.')
vcfdtypes = {'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 
             'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str}
def rowMgmt(key, seq):
    return tuple([':'.join(key), key[0]] + list(seq.replace('|', '')))

class Gene:
    def __init__(self, gene='A', seqtype='gen', parse=True):
        self.alifile = os.path.join(imgtDir, '%s_%s.txt' % (gene, seqtype))
        if parse:
            self.parseAlignment()
            self.getPolymorphSites()
        
    def parseAlignment(self):
        ## first creates a (default) dict with sequence id
        ## Interpretation:
        lc = 0
        seqs = defaultdict(str)
        for line in open(self.alifile):
            lc += 1
            # cutting the crap
            if lc < 6: continue
            if line.strip() == '|': continue
            if line.strip().startswith('AA codon'): continue
            if not line.strip(): continue
            fields = line.strip().split()
            if fields[0].lower() in ['gdna', 'prot', 'cdna', 'please']: continue

            id = tuple(fields[0].split(':'))
            if len(seqs) == 0:
                self.reference = id
            seqs[id] += ''.join(fields[1:])

        try:
            assert len(set(map(len, seqs.values()))) == 1 ## all seqs have to have equal length!!!
        except AssertionError:
            print ("Warning! Sequences have different lengths")
        ## Beware: the reference sequence still contains |
        self.spliceSites = np.cumsum(list(map(len, seqs[self.reference].split('|'))))#[idx for idx, pos in enumerate(seqs[self.reference]) if pos=='|']
        self.refseq0 = seqs[self.reference] ## hopefully this is always a copy, not a pointer
        self.refseq = self.refseq0.replace('|','')
        # probably not efficient, but only happens once
        
        ## taking care of the ref seq, ie. formatting it like the rest (. means consensus)
        seqs[self.reference] = ''.join(['.' if nt=='.' else '-' for nt in self.refseq])
        self.df = pd.DataFrame.from_records([rowMgmt(key, seqs[key]) for key in sorted(seqs.keys())])
        self.tableOffset = 2 ## 2 columns for meta data
        self.alleletypes = Counter(self.df[1])
        if self.df.shape[1] - self.tableOffset != len(self.refseq):
            pdb.set_trace()

    def getPolymorphSites(self):
        """
        A.polymorphSites {'A*01': <Poly>, ... 'A*80': <Poly>}"""

        self.polymorphSites = {atype: Poly(atype, self.df, self.tableOffset, self.refseq) for atype in self.alleletypes.keys()}

    def totalPolymorphSites(self):
        return sum([len(allele.loci) for allele in self.polymorphSites.values()])

class Poly:
    """For an allele (eg: A*01) maintains the record of polymorphism for each polymorph locus, eg.:
    {310: Counter({'G': 1, 'C': 70}),
    320: Counter({'C': 70}),
    449: Counter({'A': 1, 'C': 70}),
    461: Counter({'G': 70}),
    """
    def __init__(self, atype, df, offset, refseq, countIndels=False):
        subali = df[df[1] == atype]
        gene = atype.split('*')[0]
        polySites = set('ACGT.') if countIndels else set('ACGT')
        self.loci = {}
        for col in range(offset, df.shape[1]):
            stats = Counter(subali.iloc[:, col])
            if set(stats.keys()).intersection(polySites): ## polymorphic
                #pdb.set_trace()
                refNt = refseq[col - offset]
                ## acc. to convention cant have an explicit nt that does not differ from ref
                if refNt !='.' and refNt in stats.keys():
                    pdb.set_trace()
                if '-' in stats:
                    stats[refNt] = stats.pop('-')
                self.loci[(gene, col - offset)] = stats
    #def entropy(self):

class Population:
    def __init__(self, alleleFreqs, genes, name='Pop', run=True, keepData=True): ## options: 'core exons', 'nonsilent', 'all'
        def normalized(afdict, gene):
            ## first, filter out those, we don't have data for, shouldn't happen (that much)
            d = {k:v for k,v in afdict.items() if k in genes[gene].polymorphSites} ## TODO: give warning
            dropped = [(k,v) for (k,v) in afdict.items() if (not k in d) and v>0.01]
            if dropped: print("Warning: dropped %s (Pop %s, gene %s)" % (dropped, name, gene))
            return {k:v/sum(d.values()) for k,v in d.items()}
        self.alleleFreqs = {gene: normalized(af, gene) for gene, af in alleleFreqs.items()}
        self.name = name
        if run:
            self.getAllPolySites(genes)
            self.makePWM(genes)
        if keepData: ## this could be large if population is pickled. Needed in bootstrap
            self.genes = genes

    def getAllPolySites(self, genes):
        allSites = set()
        for gene, geneAFs in self.alleleFreqs.items():
            for atype in geneAFs.keys():
                allSites = allSites.union(genes[gene].polymorphSites[atype].loci.keys())
        self.allpolySites = sorted(list(allSites))
    def makePWM(self, genes): ## admixing probabilities for all encountered allele types according to observed allele freqs.
        pwm = np.zeros((5, len(self.allpolySites))) ## A C G T .
        for idx, (gene, site) in enumerate(self.allpolySites):
            for atype, freq in self.alleleFreqs[gene].items():
                refNt = genes[gene].refseq[site]
                refNtpos = bases.index(refNt)
                if (gene, site) in genes[gene].polymorphSites[atype].loci:
                    stats = genes[gene].polymorphSites[atype].loci[(gene, site)]
                    total = float(np.sum(list(stats.values())) - stats['*'])

                    probDist = np.array([stats[nt]/total for nt in bases])
                    try:
                        assert np.isclose(probDist.sum(), 1)
                    except AssertionError:
                        pdb.set_trace()
                    pwm[:,idx] += freq * probDist

                else:
                    a = np.zeros(5)
                    a[refNtpos] = 1
                    pwm[:, idx] += freq * a
        ## Make a nice dataframe out of it, good for merging, viewing ...
        self.pwm = pd.DataFrame(pwm.T, index=self.allpolySites, columns=bases)
    def bootstrap(self, selectedLoci, afbased=True):
        pwm = self.pwm if afbased else self.pwmVCF
        if not "genes" in vars(self):
            print ("Warning: needs gene data for variable loci, will probably fail!!!")
        def getProbDist(locus):
            if locus in pwm.index:
                return pwm.loc[[locus]].iloc[0] ## this required rewriting since locus became a tuple
            ## slightly redundant code, see makePWM
            ## locus not variable in this population
            refNt = self.genes[locus[0]].refseq[locus[1]]
            refNtpos = bases.index(refNt)
            a = np.zeros(5)
            a[refNtpos] = 1.
            return a
        return np.hstack([getProbDist(locus) for locus in selectedLoci]) ##

    def readVCF(self, genesShort):
        ## TODO: fix for all genes
        pwm = []
        for gene in genesShort: ## FIX:
            with open(f'{vcfDir}chr6_1000g_HLA_{gene}.vcf', 'r') as f:
                lines = [l for l in f if not l.startswith('##')] ## everythin in memory :-(
            vcf = pd.read_csv(io.StringIO(''.join(lines)), dtype=vcfdtypes, \
                            sep='\t').rename(columns={'#CHROM': 'CHROM'})
            vcf = vcf[['POS','REF','ALT']+self.sampleIDs] ## just select pop samples # like vcf subset  
            snvMask = (vcf['REF'].str.len() == 1) & (vcf['ALT'].str.len() == 1) # only biallelic SNVs
            self.vcf = vcf.loc[snvMask]

            for i, r in self.vcf.iterrows():
                row = [(gene, r.POS), 0,0,0,0,0]
                alleleCount = Counter(''.join(list(r)[3:]))
                total = alleleCount['0'] + alleleCount['1']
                refNtpos = bases.index(r.REF)
                row[refNtpos+1] = alleleCount['0']/total
                refNtpos = bases.index(r.ALT)
                row[refNtpos+1] = alleleCount['1']/total
                pwm.append(row)
        self.pwmVCF = pd.DataFrame(pwm, columns=['index']+bases)
        self.pwmVCF.set_index('index', inplace=True)
        self.allpolySitesVCF = self.pwmVCF.index.to_list() ## no overwriting, compatible with phylo construction
    
class PopulationSet:
    def __init__(self, populations):
        self.populations = populations
        self.popnames = [pop.name for pop in self.populations]
    def bootstrap(self, afbased=True, basename='majorityTree', treebuilder='nj', bootstraps=1000, outgroup=None):
        """treebuilder could be nj/upgma, outgroup: a population name or 'midpoint'"""
        ## allLoci: all loci that are variable in at least one population
        allpolySites, pwm = {True:  ['allpolySites', 'pwm'],
                             False: ['allpolySitesVCF', 'pwmVCF'],}[afbased]
        allLoci = set()
        for pop in self.populations:
            allLoci = allLoci.union(getattr(pop, allpolySites))
        allLoci = list(allLoci) ## sort it? Reduce to independent sites (www.pnas.org/content/93/23/13429, run LD?)
        
        sites = len(allLoci)
        trees = []
        print ("Bootstrapping, rounds:", end=' ')
        for bootstrap in range(bootstraps): ## see also parallelized version
            print(bootstrap, end=' ')
            selectedLoci0 = np.random.choice(range(len(allLoci)), sites, replace=True)
            selectedLoci = [allLoci[l] for l in selectedLoci0]
            df = pd.DataFrame([pop.bootstrap(selectedLoci, afbased) for pop in self.populations], index=self.popnames)
            #import pdb; pdb.set_trace()
            dmNei = neiDF(df, [5]*(sites-1))
            ## annoying conversion, BioPython couldnt be just more compatible with scipy/pdist?
            dmTriangular = [list(dmNei[i, :(i + 1)]) for i in range(len(dmNei))]

            m = _DistanceMatrix(self.popnames, dmTriangular)
            constructor = DistanceTreeConstructor() # could've passed treebuilder here too
            tree = getattr(constructor, treebuilder)(m)
            if outgroup == 'midpoint':
                tree.root_at_midpoint()
            elif not outgroup is None:
                tree.root_with_outgroup({'name': outgroup})
            trees.append(tree) ## use nj!
        ## debug info:
        print(f'selectedLoci: {selectedLoci[:30]}')
        ## see https://biopython.org/wiki/Phylo, turned out to be more suitable than dendropy/sumtrees
        self.majorityTree = Consensus.majority_consensus(trees) ## also consider strict_consensus and adam_consensus (but they don't have bootstrap support values)
        treefile = '%s/%s_%s_%s_%s.nwk' %(resultDir, basename, bootstraps,treebuilder, len(self.populations))
        Phylo.write(self.majorityTree, treefile, format='newick')
        print(f'wrote {treefile}')
        Phylo.draw_ascii(self.majorityTree)

        #Phylo.draw(self.majorityTree)

if __name__ == "__main__":
    #genes = {g: Gene(gene=g) for g in ['A', 'B']}
    #Count number of
    with open('../Data/hlaABPoly.pcl', 'rb') as f: genes = pickle.load(f)
    pop1 = Population({'A': {'A*36': 0.2, 'A*80': 0.8}, 'B': {'B*07': .9, 'B*08': .1}}, genes, 'Pop1')
    pop2 = Population({'A': {'A*36': 0.1, 'A*80': 0.9}, 'B': {'B*07': .93, 'B*08': .07}}, genes, 'Pop2')

    pop3 = Population({'A': {'A*36': 0.9, 'A*80': 0.1}, 'B': {'B*07': .5, 'B*08': .5}}, genes, 'Pop3')
    pop4 = Population({'A': {'A*36': 0.9, 'A*80': 0.1}, 'B': {'B*07': .9, 'B*08': .1}}, genes, 'Pop4')
    pop5 = Population({'A': {'A*36': 0.8, 'A*80': 0.2}, 'B': {'B*07': .4, 'B*08': .6}}, genes, 'Pop5')

    ps = PopulationSet([pop1, pop2, pop3, pop4, pop5])
    ps.bootstrap(bootstraps=10, treebuilder='upgma', outgroup='Pop4')
