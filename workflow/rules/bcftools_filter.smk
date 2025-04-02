
## filter bcf by biallelic snps
rule snp_filter:
    input:
        bcf=config['bcf'],
    output:
        bcf="results/snp_filter/bi-snp.filtered.bcf",
        csi="results/snp_filter/bi-snp.filtered.bcf.csi",
    params:
        filters=config['bcftools_filters'],
    conda:
        "../envs/bcftools.yaml",
    log:
        "results/logs/snp_filter/snp_filter.log"
    shell:
        "bcftools view -v snps -m 2 -M 2 "
        "-i '{params.filters}' {input.bcf} -Ob -o {output.bcf} > {log}; " 
        "bcftools index -c {output.bcf} " 

## get stats for filtered bcf
rule bcf_stats:
    input:
        bcf="results/snp_filter/bi-snp.filtered.bcf",
        csi="results/snp_filter/bi-snp.filtered.bcf.csi",
    output:
        stats="results/snp_filter/bi-snp.filtered.bcf.stats",
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bcftools stats {input.bcf} > {output.stats} "
