# Repermulize

This is an accompanying package for phylogenize to run the repermulize module.

Install with bioconda using ```conda install -b bioconda repermulize``` and run with ```library(repermulize)```

Repermulize is an R package that uses a phenotype vector across species, a phylogenetic tree, and a gene matrix such as gene presence/absence or pangenome prevalence across those same species. First, the function aligns the species shared by the tree, the phenotype, and the gene matrix, then prunes the tree and reorders the phenotype so everything matches tip-for-tip. Then, the phenotype is converted into phylogenetically independent contrasts (PICs), so the downstream regression is done in phylogenetically corrected space rather than treating species as independent samples. The same PIC conversion is also applied gene-by-gene to the predictor matrix unless the genes are already supplied as PICs.
