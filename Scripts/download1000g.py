import os

"""gene positions collected from GeneCards
link to exons, maybe also useful
https://genecards.weizmann.ac.il/geneloc-bin/exon_struct.pl?disp_name=HLA-DMB&chr_nr=6"""
chromPositions = '''chr6 29909037 29913661 A +
chr6 31321649 31324989 B +
chr6 31236526 31239907 C -
chr6 32916390 32936871 DMA -
chr6 33032346 33048555 DMB -
chr6 33043703 33057473 DPB1 +
chr6 32595956 32614839 DQA1 +
chr6 32407619 32412823 DRA +
chr6 32546546 32557625 DRB1 -
chr6 30457183 30461982 E +
chr6 29690552 29706305 F +
chr6 29794744 29798902 G +'''.split('\n')

#download 1000G data for selected genes:

for chromPos in chromPositions:
    chrom, pos1, pos2, gene, strand = chromPos.split()
        #| perl vcf-subset -c HG00098 | bgzip -c > HG00098.vcf.gz') ##selecting individuals
    #os.system(f'tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr17.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz 17:1471000-1472000')
    os.system(f'tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 6:{pos1}-{pos2} > chr6_1000g_HLA_{gene}.vcf')
