#!/usr/bin/env Rscript
## usage: Rscript --vanilla mismapped-LR_significance_test.R ngsParalog_LR_output  

#args = commandArgs(trailingOnly=TRUE)
#infile = args[1]
#outfile = paste0(infile,".nonmismapped_sites.txt")
#
#lr <- read.table(args[1]) # read in ngsParalog calcLR output
#lr$pval <- 0.5*pchisq(lr$V5,df=1,lower.tail=FALSE) # append column of p-values
#lr$pval.adj <- p.adjust(lr$pval, method="bonferroni") # p-values adjusted for number of tested sites

# The 7th column of the lr data.frame is the adjusted p-value for rejecting the null hypothesis that reads
# covering the site derive from a single locus. Of course you can use any p-value adjustment of your
# choosing, e.g. "fdr".

# generate list of sites that don't show evidence of mismapping at 0.05 significance level:
#qc.sites <- lr[-which(lr$pval.adj < 0.001),1:2]
#write.table(qc.sites,file=outfile,row.names=F)


library(data.table)

args = commandArgs(trailingOnly = TRUE)
infile = args[1]
outfile = paste0(infile, ".nonmismapped_sites.txt")

# Read in the data using fread (faster than read.table)
lr <- fread(infile)

# Append columns for p-value and adjusted p-value
lr$pval <- 0.5 * pchisq(lr$V5, df = 1, lower.tail = FALSE)
lr$pval.adj <- p.adjust(lr$pval, method = "bonferroni")

# Subset the columns
lr.pvals <- lr[,.(V1,V2,pval,pval.adj)]
# Write to a file
fwrite(lr.pvals, file = paste0(infile,".pvals.txt"),sep='\t')

# Filter out sites that show evidence of mismapping
# null hypothesis: reads covering site map uniquely / derive from one locus
# alternative hypothesis: reads covering site mismap
# so keep all sites >= p-value threshold
# using a p-value threshold of 0.001
qc.sites <- lr[pval.adj >= 0.001, .(V1, V2)]  # Efficient subsetting with data.table
# Write the filtered sites to a file
fwrite(qc.sites, file = outfile,sep='\t',col.names=F)

