


## BEFORE filtering, let's run the per site stats including
## - quality per site
## - mean depth per site
## - proportion of missingness per site
## can use these to set appropriate filters

## XXX to randomly sample ~100k sites from the summary stat output 
# awk 'BEGIN {srand()} !/^$/ { if (rand() <= .0001 || FNR==1) print $0}' pre_filter.ldepth.mean > out

## randomly subsample 1% of sites
rule pre_random_subsample_bcf:
    input:
        config['bcf'],
    output:
        "results/vcftools_summary_stats/pre_filter_stats/pre_filter.subsample.vcf",
    conda:
        "../envs/vcftools.yaml",
    shell:
        "bcftools view -v snps -m 2 -M 2 {input} | vcfrandomsample -r 0.01 > {output} "  

## output various summary statistics
rule pre_filter_stats:
    input:
        "results/vcftools_summary_stats/pre_filter_stats/pre_filter.subsample.vcf", 
    output:
        "results/vcftools_summary_stats/pre_filter_stats/pre_filter.subsample.ldepth.mean",
        "results/vcftools_summary_stats/pre_filter_stats/pre_filter.subsample.lqual",
        "results/vcftools_summary_stats/pre_filter_stats/pre_filter.subsample.lmiss",
        "results/vcftools_summary_stats/pre_filter_stats/pre_filter.subsample.idepth",
        "results/vcftools_summary_stats/pre_filter_stats/pre_filter.subsample.imiss",
        "results/vcftools_summary_stats/pre_filter_stats/pre_filter.subsample.frq",
        "results/vcftools_summary_stats/pre_filter_stats/pre_filter.subsample.het",
    params:
        out_prefix="results/vcftools_summary_stats/pre_filter_stats/pre_filter.subsample",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --depth --out {params.out_prefix} ;"
        "vcftools --gzvcf {input} --site-mean-depth --out {params.out_prefix} ;"
        "vcftools --gzvcf {input} --site-quality --out {params.out_prefix} ;"
        "vcftools --gzvcf {input} --missing-indv --out {params.out_prefix} ;"
        "vcftools --gzvcf {input} --missing-site --out {params.out_prefix} ;"
        "vcftools --gzvcf {input} --freq2 --out {params.out_prefix} --max-alleles 2 ;"
        "vcftools --gzvcf {input} --het --out {params.out_prefix} "


## visualize and summarize stats
rule pre_summarize_stats:
    input:
        "results/vcftools_summary_stats/pre_filter_stats/pre_filter.subsample.vcf",
        "results/vcftools_summary_stats/pre_filter_stats/pre_filter.subsample.ldepth.mean",
    output:
        "results/vcftools_summary_stats/pre_filter_stats/pre_filter.subsample.summary.txt",
    params:
        stat_dir="results/vcftools_summary_stats/pre_filter_stats",
        in_prefix="pre_filter.subsample",
    conda:
        "../envs/R.yaml",
    shell:
        "Rscript --vanilla workflow/scripts/summarize_stats.R {params.stat_dir} {params.in_prefix} > {output}"

## generate vcftools summary stats for post snp filtered bcf 
use rule pre_random_subsample_bcf as post_snp_random_subsample_bcf with:
    input: 
        "results/snp_filter/bi-snp.filtered.bcf", 
    output:
        "results/vcftools_summary_stats/post_snp_filter_stats/bi-snp.filtered.subsample.vcf",

use rule pre_filter_stats as post_snp_filter_stats with:
    input:
        "results/vcftools_summary_stats/post_snp_filter_stats/bi-snp.filtered.subsample.vcf",
    output: 
        "results/vcftools_summary_stats/post_snp_filter_stats/bi-snp.filtered.subsample.het",
    params: 
        out_prefix="results/vcftools_summary_stats/post_snp_filter_stats/bi-snp.filtered.subsample"

## post summary visuals and stats
use rule pre_summarize_stats as post_summarize_stats with:
    input:
        "results/vcftools_summary_stats/post_snp_filter_stats/bi-snp.filtered.subsample.het",
    output:
        "results/vcftools_summary_stats/post_snp_filter_stats/bi-snp.filtered.subsample.summary.txt",
    params:
        stat_dir="results/vcftools_summary_stats/post_snp_filter_stats",
        in_prefix="bi-snp.filtered.subsample"
        

## follow this:
## https://speciationgenomics.github.io/filtering_vcfs/
## output all the stats
