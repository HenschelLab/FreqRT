import pandas as pd
from collections import Counter
from importlib import reload; import polyBootstrap; reload(polyBootstrap) ## when debugging, 2b removed
from polyBootstrap import Population, PopulationSet, Gene, Poly
import pickle
import subprocess
import alleleFreqAPI; reload(alleleFreqAPI)
from alleleFreqAPI import PopulationTree

## assuming that this script is run from the Scripts dir
## table derived from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt, which includes 4 digits
df = pd.read_csv('../Data/20181129_HLA_types_full_1000_Genomes_Project_panelA.txt', sep='\t').replace('None', pd.np.nan)  
genes = ['HLA-A 1', 'HLA-A 2', 'HLA-B 1', 'HLA-B 2', 'HLA-C 1', 'HLA-C 2', 'HLA-DQB1 1', 'HLA-DQB1 2', 'HLA-DRB1 1', 'HLA-DRB1 2']

cols = ['Region', 'Population', 'Sample ID'] + genes[:6]
genesShort = [g.split()[0].split('-')[-1] for g in genes[:6:2]]  
df = df[cols].dropna(how='any')
## reducing the rows to those samples for which we actually have vcf data
vcfdir = '../Data'
availableSamples = subprocess.check_output(f'grep "^#CHROM" {vcfdir}/chr6_1000g_HLA_A.vcf', shell=True).decode().split('\t')[9:]
mask = df['Sample ID'].isin(availableSamples)
df = df.loc[mask]  ## power of pandas, loving it!!
"""
http://www.internationalgenome.org/faq/which-populations-are-part-your-study
         Counter({('AFR', 'ACB'): 96,
         ('AFR', 'ASW'): 59,
         ('AFR', 'ESN'): 96,
         ('AFR', 'GWD'): 112,
         ('AFR', 'LWK'): 98,
         ('AFR', 'MSL'): 85,
         ('AFR', 'YRI'): 108,
         ('AMR', 'CLM'): 94,
         ('AMR', 'MXL'): 64,
         ('AMR', 'PEL'): 85,
         ('AMR', 'PUR'): 104,
         ('EAS', 'CDX'): 93,
         ('EAS', 'CHB'): 103,
         ('EAS', 'CHS'): 105,
         ('EAS', 'JPT'): 104,
         ('EAS', 'KHV'): 99,
         ('EUR', 'CEU'): 99,
         ('EUR', 'FIN'): 99,
         ('EUR', 'GBR'): 91,
         ('EUR', 'IBS'): 107,
         ('EUR', 'TSI'): 107,
         ('SAS', 'BEB'): 86,
         ('SAS', 'GIH'): 102,
         ('SAS', 'ITU'): 102,
         ('SAS', 'PJL'): 96,
         ('SAS', 'STU'): 102})
         """

def makePop(pop, df, afbased=True):
    print (f"Dealing with Population {pop}")
    df1 = df[df.Population == pop]
    ## jackknifing: compose populations randomly, equal size
    df2 = df1.sample(n=popSize, random_state=1)
    ##memorize df2['Sample ID']
    pdict = {}
    for gene in genesShort:
        tmpdict = Counter(df2[f'HLA-{gene} 1']) + Counter(df2[f'HLA-{gene} 2']) 
        total = sum(tmpdict.values())
        pdict[gene] = {f'{gene}*{k}': v/total for k,v in tmpdict.items()}
    population = Population(pdict, genesData, pop)
    population.sampleIDs = list(df2['Sample ID'])
    if not afbased:
        population.readVCF(genesShort)
    return population

def runPhylip(ps, basename):
    ptree = PopulationTree(precalc=False)
    ptree.df = ps.makeAFtable() ## emulate the previously used spreadsheet format for af's
    ptree.fix(renameCols=False, lastCol=ptree.df.columns.values[-3])
    ptree.runPhylipGenDist(basename, selectedLoci = genesShort, selectedPops=False)
    ptree.runPhylipNeighbor() 
    ptree.renamingPopulationsInTree()

popSelection = 'YRI GBR IBS CEU PJL PUR ACB CLM GIH PEL JPT'.split()
popCounts = Counter(df.Population)
popSize = min([popCounts[pop] for pop in popSelection])
#popSize = 90 

with open('../Data/hlaAllPoly.pcl', 'rb') as f: 
    genesData = pickle.load(f) 

ps = PopulationSet([makePop(pop, df, afbased=False) for pop in popSelection])
basenamePhylip  = f'phylip_{"-".join(genesShort)}_ps{popSize}_pops{len(popSelection)}' 
basenameAF  = f'afbased_{"-".join(genesShort)}_ps{popSize}_pops{len(popSelection)}'
basenameVCF = f'vcfbased_{"-".join(genesShort)}_ps{popSize}_pops{len(popSelection)}'

runPhylip(ps, basenamePhylip) ## files end up in ../GenDist, needs to be created on first use!
#ps.bootstrap(bootstraps=100, basename=basenameAF, treebuilder='nj', outgroup='YRI')
#ps.bootstrap(bootstraps=100, basename=basenameVCF, treebuilder='nj', outgroup='YRI', afbased=False)
