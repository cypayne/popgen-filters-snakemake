        
## make scatters file
#checkpoint make_scatters_file:
#    input:
#        chr_lengths=config["chr_lengths_file"],
#    output:
#        scatters="results/scatters/scatters_list.txt",
#    params:
#        chrom_num=config['chrom_num'],
#        scatter_length=config["scatter_length"]
#    shell:
#        "python3 ../scripts/make_scatters_file_lengths.py {input.chr_lengths} "
#        "{params.chr_num} {params.scatter_length} " 

## make scatter bed files
rule make_scatter_files:
    input:
        chr_file=config['chrom_lengths_list'],
    output:
        bed="results/scatter_beds/{chrom}_{scat}.bed",
    run:
        create_scatter_bed_files(input.chr_file)
 
## ngsParalog takes samtools mpileup
## create snp positions list 
rule make_pos_list:
    input:
        bcf="results/snp_filter/bi-snp.filtered.bcf",
        csi="results/snp_filter/bi-snp.filtered.bcf.csi",
        bed="results/scatter_beds/{chrom}_{scat}.bed",
    output:
        pos_list="results/positions/{chrom}_{scat}.pos",
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bcftools query -f '%CHROM\t%POS\n' -R {input.bed} {input.bcf} > {output.pos_list} "


## generate the pileup files for each chromosome 
rule prep_mpileup:
    input:
        bam_list=config['bam_list'],
        pos_list="results/positions/{chrom}_{scat}.pos",
    output: 
        pileup="results/mpileup/{chrom}_{scat}.pileup",
    conda:
        "../envs/samtools.yaml"
    log:
        "results/logs/mpileup/{chrom}_{scat}.log"
    benchmark:
        "results/benchmarks/mpileup/{chrom}_{scat}.bmk"
    shell:
        "samtools mpileup -b {input.bam_list} -l {input.pos_list} "
        "-r {wildcards.chrom} -q 0 -Q 0 --ff UNMAP,DUP -o {output.pileup} > {log} "


## calculate likelihood ratios of mismapping reads at each site
## run ngsParalog calcLR
rule calculate_LR: 
    input:
        pileup="results/mpileup/{chrom}_{scat}.pileup",
    output:
        lr="results/ngsParalog/{chrom}_{scat}.lr", 
    params:
        path=config['ngsParalog_path'],
    resources:
        time="20:00:00",
    shell:
        "if [ -s {input.pileup} ] ; then "
        "   {params.path}/ngsParalog calcLR -infile {input.pileup} "
        "   -outfile {output.lr} -minQ 30 -minind 1 -mincov 1 ; "
        "else "
        "   touch {output.lr} ; "
        "fi "


## combine per chromosome likelihood ratios
rule combine_LRs:
    input:
        expand("results/ngsParalog/{chrom}_{scat}.lr",zip,chrom=chrom_scatter_table.chrom,scat=chrom_scatter_table.scatter_name),
    output:
        "results/ngsParalog/bi-snp.filtered.lr",
    shell:
        "cat {input} > {output} "

## make list of sites that don't show evidence of mismapping 
## at 0.05 significance level
rule list_passing_sites:
    input:
        lr="results/ngsParalog/bi-snp.filtered.lr",
    output:
        pass_list="results/ngsParalog/bi-snp.filtered.lr.nonmismapped_sites.txt",
    conda:
        "../envs/R.yaml",
    shell:
        "Rscript --vanilla workflow/scripts/mismapped-LR_significance_test.R {input.lr} "


## average LRs in windows
rule avg_lr_in_windows:
    input:
        lr="results/ngsParalog/bi-snp.filtered.lr",
    output:
        windowed="results/ngsParalog/windows.bi-snp.filtered.lr",
    conda:
        "../envs/R.yaml",
    shell:
        "python3 workflow/scripts/avg_by_window.LR.py 10000 {input.lr} {output.windowed}"


## filter by nonmismaped sites
rule filter_bcf_by_nomismap: 
    input:
        bcf="results/snp_filter/bi-snp.filtered.bcf",
        pos="results/ngsParalog/bi-snp.filtered.lr.nonmismapped_sites.txt",
    output:
        bcf="results/filtered_bcfs/nomismap.bi-snp.filtered.bcf",
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bcftools view -R {input.pos} {input.bcf} -Ob -o {output.bcf} ; "
        "bcftools index -c {output.bcf} "


## get stats for filtered bcf
rule nomismap_bcf_stats:
    input:
        bcf="results/filtered_bcfs/nomismap.bi-snp.filtered.bcf",
    output:
        stats="results/filtered_bcfs/nomismap.bi-snp.filtered.bcf.stats",
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bcftools stats {input.bcf} > {output.stats} "


## calculate per sample stats
## includes heterozyosity, transition, transversion, etc info
rule nommismap_bcf_stats_per_sample:
    input:
        bcf="results/filtered_bcfs/nomismap.bi-snp.filtered.bcf",
    output:
        stats="results/filtered_bcfs/nomismap.bi-snp.filtered.bcf.sample-stats",
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bcftools stats -s - {input.bcf} | grep '^PSC' -B 1 > {output.stats} "  
