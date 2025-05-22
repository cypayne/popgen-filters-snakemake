
## filter by linkage disequilibrium 
## using a population structure aware adjustment 
## implemented in PCAone

## make plink binary inputs
rule generate_plink_binaries:
    input: 
        "results/snp_filter/bi-snp.filtered.bcf",
    output:
        "results/plink_binaries/{chrom}.bed",
        "results/plink_binaries/{chrom}.bim",
        "results/plink_binaries/{chrom}.fam",
    params:
        out_prefix="results/plink_binaries/{chrom}",
    conda:
        "../envs/plink.yaml",
    shell:
        "plink2 --bcf {input} --make-bed --chr {wildcards.chrom} --maf 0.01 "
        "--allow-extra-chr --set-all-var-ids @:# --out {params.out_prefix} || "
        "touch {output} "

## calculate ancestry-adjusted linkage disequilibrium matrix 
## assumes that first 3 PCs in PCA capture population structure
## runs through different SVD options, since some can fail depending
## on the size of the scatter data 
rule adjusted_ld_matrix:
    input: 
        bed="results/plink_binaries/{chrom}.bed",
        bim="results/plink_binaries/{chrom}.bim",
        fam="results/plink_binaries/{chrom}.fam",
    output:
        "results/ld_filter/ld_matrix/{chrom}.residuals",
        "results/ld_filter/ld_matrix/{chrom}.mbim",
    params:
        in_prefix="results/plink_binaries/{chrom}",
        out_prefix="results/ld_filter/ld_matrix/{chrom}",
    conda:
        "../envs/pcaone.yaml",
    shell:
        "if [ -s {input.bed} ]; then " 
        "    PCAone -b {params.in_prefix} -k 3 --ld -o {params.out_prefix} --svd 2 || "
        "    PCAone -b {params.in_prefix} -k 3 --ld -o {params.out_prefix} --svd 0 || "
        "    PCAone -b {params.in_prefix} -k 3 --ld -o {params.out_prefix} --svd 1 ; "
        "else "
        "touch {output} ; "
        "fi " 


## report ld statistics
rule output_ld_stats:
    input:
        resid="results/ld_filter/ld_matrix/{chrom}.residuals",
        mbim="results/ld_filter/ld_matrix/{chrom}.mbim",
    output:
        "results/ld_filter/ld_stats/{chrom}.ld.gz",
    params:
        out_prefix="results/ld_filter/ld_stats/{chrom}",
        ld_bp=config['ld_window_bp'],
    conda:
        "../envs/pcaone.yaml",
    shell:
        "if [ -s {input.resid} ] ; then "
        "    PCAone -B {input.resid} --match-bim {input.mbim} --ld-bp {params.ld_bp} "
        "    --print-r2 -o {params.out_prefix} --svd 2 ; "
        "else "
        "    touch {output} ; "
        "fi "


## combine ld stats
rule combine_ld_stats:
    input:
        expand("results/ld_filter/ld_stats/{chrom}.ld.gz", chrom=unique_chroms),
    output: 
        "results/ld_filter/bi-snp.filtered.ld.gz",
    shell:
        "cat {input} > {output} "


## calculate LD r^2 to plot the r^2 decay curve
rule calc_r2_decay:
    input:
        "results/ld_filter/bi-snp.filtered.ld.gz",
    output:
        "results/ld_filter/ld_r2_decay_curve.txt",
    params:
        ld_bp=config['ld_window_bp'],
    conda:
        #"../envs/R.yaml",
        "../envs/python.yaml",
    shell:
        "python workflow/scripts/calc_decay_bin.py {input} {output} "
        "--window {params.ld_bp} --bins 100 "
        #"Rscript --vanilla workflow/scripts/calc_decay_bin.R {input} {output} "
        #"--window {params.ld_bp} --bins 100 "

        
## prune by ancestry-adjusted LD 
## keep all variants in ld.prune.in 
rule prune_ld:
    input:
        resid="results/ld_filter/ld_matrix/{chrom}.residuals",
        mbim="results/ld_filter/ld_matrix/{chrom}.mbim",
    output:
        "results/ld_filter/ld_prune/{chrom}.ld.prune.in",
    params:
        out_prefix="results/ld_filter/ld_prune/{chrom}",
        ld_bp=config['ld_window_bp'],
        ld_r2=config['ld_threshold_r2'],
    conda:
        "../envs/pcaone.yaml",
    shell:
        "if [ -s {input.resid} ]; then "
        "    PCAone -B {input.resid} --match-bim {input.mbim} --ld-bp {params.ld_bp} "
        "    --ld-r2 {params.ld_r2} -o {params.out_prefix} --svd 2 ; "
        "else "
        "    touch {output} ; "
        "fi "


## combine the ld.prune.in files
rule combine_pruned:
    input:
        expand("results/ld_filter/ld_prune/{chrom}.ld.prune.in",chrom=unique_chroms),
    output:
        "results/ld_filter/bi-snp.filtered.ld.prune.in",
    shell:
        "cat {input} > {output} "


## create a positions file to subset bcf
rule make_ld_pos_file:
    input:
        "results/ld_filter/bi-snp.filtered.ld.prune.in",
    output: 
        "results/ld_filter/bi-snp.filtered.ld.prune.pos",
    shell:
        "cut -d$'\t' -f1,4 {input} > {output}.temp ; "  
        "sort -k 1,1 -k 2,2 {output}.temp > {output} ; "
        "rm {output}.temp "


## filter bcf by ancestry-adjusted ld
rule filter_bcf_by_ld:
    input:
        bcf="results/snp_filter/bi-snp.filtered.bcf",
        pos="results/ld_filter/bi-snp.filtered.ld.prune.pos",
    output: 
        bcf="results/filtered_bcfs/ld.bi-snp.filtered.bcf",
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bcftools view -R {input.pos} {input.bcf} -Ob -o {output.bcf} ; "
        "bcftools index -c {output.bcf} "
     

## get stats for filtered bcf
rule ld_bcf_stats:
    input:
        bcf="results/filtered_bcfs/ld.bi-snp.filtered.bcf",
    output:
        stats="results/filtered_bcfs/ld.bi-snp.filtered.bcf.stats",
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bcftools stats {input.bcf} > {output.stats} "
