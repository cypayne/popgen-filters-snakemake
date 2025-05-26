## Visualize variant and individual quality statistics
## produced from vcftools output
## adapted from:
## https://speciationgenomics.github.io/filtering_vcfs/
## cyp X-2024

library(readr)
library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
stat_dir=args[1]
file_prefix=args[2]

cat(sprintf("VCFTOOLS SUMMARY STATISTICS\n"))
cat(sprintf("The following have been calculated on a random subset of ~100k variants\n\n"))

### variant statistics

## variant quality
lqual_file <- paste0(stat_dir,"/",file_prefix,".lqual")
if(file.exists(lqual_file)) {
    var_qual <- read_delim(lqual_file, delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
    a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() 

    ggsave(a, filename=paste0(stat_dir,"/",file_prefix,".lqual.pdf"), width=8, height=8, bg="transparent")

    cat(sprintf("SITE QUALITY\n"))
    cat(sprintf("# sites before QUAL filter: %s\n", dim(var_qual)[1]))
    cat(sprintf("# sites after QUAL>30 filter: %s\n", dim(subset(var_qual,qual>30))[1]))
    cat(sprintf("Distribution of per site quality:\n"))
    summary(var_qual$qual)
}

## variant mean depth
## e.g. mean of read depth across individuals per site
ldepth_mean_file <- paste0(stat_dir,"/",file_prefix,".ldepth.mean")
if(file.exists(ldepth_mean_file)) {
    var_depth <- read_delim(ldepth_mean_file, delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
    a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + xlim(0,50) 

    ggsave(a, filename=paste0(stat_dir,"/",file_prefix,".ldepth-mean.pdf"), width=8, height=8, bg="transparent")

    cat(sprintf("\nSITE MEAN DEPTH\n"))
    cat(sprintf("Distribution of per site mean depths:\n"))
    summary(var_depth$mean_depth)
}

## variant missingness
lmiss_file <- paste0(stat_dir,"/",file_prefix,".lmiss")
if(file.exists(lmiss_file)) {
    var_miss <- read_delim(lmiss_file, delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
    a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() + xlim(0,1)

    ggsave(a, filename=paste0(stat_dir,"/",file_prefix,".lmiss.pdf"), width=8, height=8, bg="transparent")

    cat(sprintf("\nSITE MISSINGNESS\n"))
    cat(sprintf("Distribution of per site proportion of missing data:\n"))
    summary(var_miss$fmiss)
}

## minor allele frequency
frq_file <- paste0(stat_dir,"/",file_prefix,".frq")
if(file.exists(frq_file)) {
    var_freq <- read_delim(frq_file, delim = "\t", col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
    var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
    a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()

    ggsave(a, filename=paste0(stat_dir,"/",file_prefix,".frq.pdf"), width=8, height=8, bg="transparent")

    cat(sprintf("\nSITE MINOR ALLELE FREQUENCY\n"))
    cat(sprintf("Distribution of per site minor allele frequencies:\n"))
    summary(var_freq$maf)
}

### individual statistics

## depth per individual
idepth_file <- paste0(stat_dir,"/",file_prefix,".idepth")
if(file.exists(idepth_file)) {
    ind_depth <- read_delim(idepth_file, delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
    a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light() 

    ggsave(a, filename=paste0(stat_dir,"/",file_prefix,".idepth.pdf"), width=8, height=8, bg="transparent")

    cat(sprintf("\nINDIVIDUAL READ DEPTH\n"))
    cat(sprintf("Distribution of per individual read depth:\n"))
    summary(ind_depth$depth)
}

## proportion of missing data per individual
imiss_file <- paste0(stat_dir,"/",file_prefix,".imiss")
if(file.exists(imiss_file)) {
    ind_miss  <- read_delim(imiss_file, delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
    a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()

    ggsave(a, filename=paste0(stat_dir,"/",file_prefix,".imiss.pdf"), width=8, height=8, bg="transparent")

    cat(sprintf("\nINDIVIDUAL MISSINGNESS\n"))
    cat(sprintf("Distribution of per individual proportion of missing data:\n"))
    summary(ind_miss$fmiss)
}

## heterozygosity and inbreeding coefficient per individual
het_file <- paste0(stat_dir,"/",file_prefix,".het")
if(file.exists(het_file)) {
    ind_het <- read_delim( het_file, delim = "\t", col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)

    ind_het$heterozygosity <- (as.numeric(ind_het$nsites) - as.numeric(ind_het$ho)) / as.numeric(ind_het$nsites)
    write_tsv( ind_het, paste0(stat_dir,"/",file_prefix,".heterozygosity.txt")) 

    a <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + theme_light()

    ggsave(a, filename=paste0(stat_dir,"/",file_prefix,".inbreeding-coef.pdf"), width=8, height=8, bg="transparent")

    cat(sprintf("\nINDIVIDUAL HETEROZYGOSITY\n"))
    cat(sprintf("Distribution of per individual heterozygosity:\n"))
    print(summary(ind_het$heterozygosity))

    cat(sprintf("\nINDIVIDUAL INBREEDING COEFFICIENT\n"))
    cat(sprintf("Distribution of per individual inbreeding coefficient:\n"))
    summary(ind_het$f)
}
