import pandas as pd
from collections import Counter
from importlib import reload; import polyBootstrap; reload(polyBootstrap) ## when debugging, 2b removed
from polyBootstrap import Population, PopulationSet, Gene, Poly
import pickle
import subprocess

## assuming that this script is run from the Scripts dir
## table derived from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt, which includes 4 digits
df = pd.read_csv('../Data/20181129_HLA_types_full_1000_Genomes_Project_panelA.txt', sep='\t').replace('None', pd.np.nan)  
genes = ['HLA-A 1', 'HLA-A 2', 'HLA-B 1', 'HLA-B 2', 'HLA-C 1', 'HLA-C 2', 'HLA-DQB1 1', 'HLA-DQB1 2', 'HLA-DRB1 1', 'HLA-DRB1 2']
#genes = 'A B C DQB1 DRB1'.split()[:3] 
cols = ['Region', 'Population', 'Sample ID'] + genes[:4]
genesShort = [g.split()[0].split('-')[-1] for g in genes[:4:2]]  
df = df[cols].dropna(how='any')
## reducing the rows to those samples for which we actually have vcf data
vcfdir = '../Data'
availableSamples = subprocess.check_output(f'grep "^#CHROM" {vcfdir}/chr6_1000g_HLA_A.vcf', shell=True).decode().split('\t')[9:]
mask = df['Sample ID'].isin(availableSamples)
df = df.loc[mask]  ## power of pandas, loving it!!
"""
http://www.internationalgenome.org/faq/which-populations-are-part-your-study
Counter({('AFR', 'ACB'): 98,
         ('AFR', 'ASW'): 68,
         ('AFR', 'ESN'): 111,
         ('AFR', 'GWD'): 120,
         ('AFR', 'LWK'): 106,
         ('AFR', 'MSL'): 98,
         ('AFR', 'YRI'): 111,
         ('AMR', 'CLM'): 105,
         ('AMR', 'MXL'): 70,
         ('AMR', 'PEL'): 91,
         ('AMR', 'PUR'): 107,
         ('EAS', 'CDX'): 108,
         ('EAS', 'CHB'): 108,
         ('EAS', 'CHS'): 112,
         ('EAS', 'JPT'): 105,
         ('EAS', 'KHV'): 103,
         ('EUR', 'CEU'): 102,
         ('EUR', 'FIN'): 105,
         ('EUR', 'GBR'): 102,
         ('EUR', 'IBS'): 108,
         ('EUR', 'TSI'): 112,
         ('SAS', 'BEB'): 103,
         ('SAS', 'GIH'): 109,
         ('SAS', 'ITU'): 112,
         ('SAS', 'PJL'): 108,
         ('SAS', 'STU'): 111})"""

def makePop(pop, df, vcfbased=False):
    df1 = df[df.Population == pop]
    df2 = df1.sample(n=popSize, random_state=1)
    ##memorize df2['Sample ID']
    pdict = {}
    for gene in genesShort:
        tmpdict = Counter(df2[f'HLA-{gene} 1']) + Counter(df2[f'HLA-{gene} 2']) 
        total = sum(tmpdict.values())
        pdict[gene] = {f'{gene}*{k}': v/total for k,v in tmpdict.items()}
    population = Population(pdict, genesData, pop)
    population.sampleIDs = list(df2['Sample ID'])
    if vcfbased:
        population.readVCF(genesShort)
    return population

highConfidenceTreeSelection = 'YRI FIN IBS PJL CLM JPT'.split()
popSize = 60 
## jackknifing: compose populations randomly, equal size
with open('../Data/hlaABPoly.pcl', 'rb') as f: 
    genesData = pickle.load(f)
#p = makePop('YRI', df)
#p.readVCF(['A', 'C'])

ps = PopulationSet([makePop(pop, df, vcfbased=True) for pop in highConfidenceTreeSelection])
ps.bootstrap(bootstraps=100, basename='vcfbased', treebuilder='upgma', outgroup='YRI')
