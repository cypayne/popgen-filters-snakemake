# popgen-filters-snakemake
A pipeline that filters variants by several popgen metrics and is appropriate for low-coverage datasets. Runs filtering steps in parallel fairly efficiently so as to not ruin your life. Takes the common BCF as an input format so as to not ruin your life.

I developed this pipeline so that you can take unfiltered variants (called from low or high coverage data) and pop them through this so 
they are ready for downstream analysis! It includes cool tools that I've folded into my filtering routine and speeds up computational
time by running them in parrallel on little chunks of the genome and then stitching results together at the end. 

For NOAA folks, here's how you might use this pipeline:
Eric's awesome [mega-non-model-wgs-snakeflow](https://github.com/eriqande/mega-non-model-wgs-snakeflow) BCF > 
popgen-filters-snakemake BCF > [mega-post-bcf-exploratory-snakeflows](https://github.com/eriqande/mega-post-bcf-exploratory-snakeflows)

## Currently supports filtering by the following
- filter by missingness, depth, minor allele frequency, quality, etc  (bcftools view)
- remove SNPs with mismapping reads (ngsParalog calcR)
- prune SNPs in LD, accounting for population structure (PCAone ancestry-adjusted LD statistic)

## Inputs

Required:
- BCF [binary variant file]

Optional, depending on filter:
- Mismapping reads filter with ngsParalog
  - bam_list.txt [list of bams] 
  - chrom_lengths.sorted.txt [list of chromosomes and lengths] 
- Ancestry-adjusted LD pruning filter with PCAone
  - no additional files needed 


## How To Use

1. This is a snakemake pipeline that uses Mamba/Conda environments. So, install Mamba if you haven't already, and 
create a conda environment with your snakemake specs (see the [Snakemake Installation guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)). 
I developed and tested this pipeline with snakemake v7.32.4 and python v3.11, so I'd recommend creating an environment with 
those versions! Proceed with a later version at your own risk!
```
mamba create -c conda-forge -c bioconda -n popgen-filters-snakemake snakemake-7.32.4
mamba install python=3.11
``` 

2. Git clone this repository into your working directory
```
git clone https://github.com/cypayne/popgen-filters-snakemake.git
```

3. Edit the provided config.yaml file to reflect your data and filtering needs.

4. If you are using a cluster with a job scheduler like Slurm, and would like to schedule the snakemake
jobs to take advantage of compute nodes, you can create a profile specifying the nodes, 
resource limits, etc specific to your cluster 
(you can use the provided hpcc-profiles/slurm/sedna profile as a guide - you may not need to change much).

5. Activate your snakemake Mamba environment 
```
mamba activate popgen-filters-snakemake
```

6. Run your pipeline!
On cluster with job scheduler, take advantage of compute nodes with:
```
snakemake -p --profile hpcc-profiles/slurm/cluster --configfile config.yaml 
```
On local machine:
```
snakemake -p --configfile config.yaml
```


## Configuration - config.yaml 
You can edit the provided config.yaml file to reflect your data and filtering needs.

### Standard population genetics filters

If you just want to filter by standard popgen filters with bcftools (i.e. all other filters are set to False),
then you only need to pass in a bcf file with your unfiltered variants. You'll also want to set the parameters
you want bcftools to filter by. You should set the parameters for what you want to include in the filtered 
variant set (the filtering step uses bcftools view -i). Only biallelic SNPs meeting these criteria will be
retained. Here's what the example filters by:
- FILTER="PASS" = keep variants that pass filter
- MAF>0.05 = keep variants with a minor allele frequency above 5%
- F_MISSING<0.1 = keep variants with genotypes for at least 90% of samples (<10% missing data) 
- FMT/DP>3 & FMT/DP<20 = keep variants with more than 3 reads but less than 20 reads for every sample
- AVG(FMT/DP)>3 & AVG(FMT/DP)<20 = keep variants with on average more than 3 reads but less than 20 reads across samples 
- QUAL>30 = keep variants with a quality score greater than 30

```
## in config.yaml
bcf: inputs/all.bcf
bcftools_filters: 'FILTER="PASS" & MAF>0.05 & F_MISSING<0.1 & FMT/DP>3 & FMT/DP<20 & AVG(FMT/DP)>3 & AVG(FMT/DP)<20 & QUAL>30' 
```

### Ancestry-adjusted LD pruning with PCAone 
If you want to filter by linkage disequilibrium for a sample set that may or may not have population structure with
[PCAone](https://github.com/Zilong-Li/PCAone), you don't need to pass in any new input files, but you will need to 
set some parameters. You will need to decide the window size to calculate the LD statistic for and the r^2 threshold
for pruning SNPs (i.e. in this example, within a 1Mbp window, any SNP with an LD r^2 greater than 0.2 with an upstream 
SNP will be dropped). 
```
## in config.yaml
PCAone: True
# define window-size within which to calculate LD
ld_window_bp: 100000
# define r^2 ld threshold for filtering
ld_threshold_r2: 0.2
```

### Mismapping reads filter with ngsParalog
If you want to filter by mismapped reads with [ngsParalog](https://github.com/tplinderoth/ngsParalog), you 
will need to download ngsParalog separately and provide the path to the executable. You can download 
with the following steps:
```
git clone https://github.com/tplinderoth/ngsParalog
cd ngsParalog; make
```
You will need to pass in a list of the paths to the bam files and the list of chromosomes with their lengths 
sorted in descending order. The pipeline will automate everything else.
```
## in config.yaml
ngsParalog: True
# path to ngsParalog executable
ngsParalog_path: "/home/cpayne/programs/ngsParalog"
# list of bam files
bam_list: inputs/bam_list.txt
# list of chromosomes and lengths in descending order
chrom_lengths_list: inputs/chrom_lengths.sorted.txt  
```
The bam list has the path to each sample's bam file, one per line, and no header. The path can be relative 
to the pipeline directory or absolute.
```
inputs/bams/sample1.bam
inputs/bams/sample2.bam
inputs/bams/sample3.bam
```
The chromosome lengths list has two TAB-delimited columns, one with chromosome name and one with the length.
The chromosomes are sorted in descending order by length. There is no header. 
```
chrom1  4000000
chrom2  3500000
chrom3  2000000
```
