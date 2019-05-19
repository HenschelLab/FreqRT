import pandas as pd

## assuming that this script is run from the Scripts dir
## table derived from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt, which includes 4 digits
df = pd.read_csv('../Data/20181129_HLA_types_full_1000_Genomes_Project_panelA.txt', sep='\t') 

