### popgen-filters-snakemake
## A pipeline that filters large bcfs by several popgen metrics. 
## Runs filtering steps in parallel so as not to take your entire life


## Currently supports filtering by the following:
- filter by missingness, depth, allele freq, quality, etc  (bcftools view)
- remove SNPs with mismapping reads (ngsParalog calcR)
- prune SNPs in LD, accounting for pop structure (PCAone ancestry-adjusted LD statistic)


