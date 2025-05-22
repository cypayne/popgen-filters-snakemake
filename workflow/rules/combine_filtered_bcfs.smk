
## if filtering by both mismapping reads and LD, 
## intersect the two bcfs
## i.e. keep only the SNPs in both bcfs 
rule combine_filtered_bcfs: 
    input:
        ld="results/filtered_bcfs/ld.bi-snp.filtered.bcf",
        mismap="results/filtered_bcfs/nomismap.bi-snp.filtered.bcf",
    output:
        "results/filtered_bcfs/ld.nomismap.bi-snp.filtered.bcf",
    conda:
        "../envs/bcftools.yaml",
    params:
        tempdir="results/filtered_bcfs/bcftools_intersect",
    shell:
        "bcftools isec {input.ld} {input.mismap} -Ob -p {params.tempdir} ; "
        "mv {params.tempdir}/0002.bcf {output} ; "
        "mv {params.tempdir}/0002.bcf.csi {output}.csi ; "
        "rm -r {params.tempdir} ; "


## get stats for combined bcf
rule combined_bcf_stats:
    input:
        bcf="results/filtered_bcfs/ld.nomismap.bi-snp.filtered.bcf",
    output:
        stats="results/filtered_bcfs/ld.nomismap.bi-snp.filtered.bcf.stats",
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bcftools stats {input.bcf} > {output.stats} "


## calculate per sample stats
## includes heterozyosity, transition, transversion, etc info
## BUT keep in mind this reports hard calls, i.e. does not use
## genotype likelihoods and therefore isn't ideal for low cov data
rule combined_bcf_stats_per_sample:
    input:
        bcf="results/filtered_bcfs/ld.nomismap.bi-snp.filtered.bcf",
    output:
        stats="results/filtered_bcfs/ld.nomismap.bi-snp.filtered.bcf.sample-stats",
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bcftools stats -s - {input.bcf} | grep '^PSC' -B 1 > {output.stats} "  
