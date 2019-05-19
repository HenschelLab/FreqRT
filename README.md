# FreqRT

(pronounce fr&#0230;kkert, danish for cheeky person, also affectionate)
The goal is to provide a method that allows to draw phylogenetic trees from few but hyperpolymorphic genes, while doing reliable bootstrapping.
For example genes belonging to the MHC have many variable loci, leading to many different alleles.
Conventional methods for however are not suitable for low numbers of marker
genes (Phylip, DISPAN, check GenPop), since they perform bootstrapping with gene loci as units.

E.g., in
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0169929
the authors claim that the provided phylogeny (in Figure 2) has 100%
bootstrap values. This tree however is just based on a single gene
(HLA-DRB1), and the used tool (DISPAN) "resamples" exactly from that
one gene only, thus inevitably producing identical bootstraps. In
other words, if only allele frequencies from a single gene are used,
conventional bootstrapping performs resampling with replacement, draws
repeatedly from only a single locus. Trivially, this single locus is
in agreement with itself, and the tools report 100\% branch support
for all clades in the tree, which leads to false confidence.

Sample usage:

look at comprehensiveTree.py
