
## thin filtered bcf by physical distance between snps 
rule thin_filter:
    input:
        bcf="results/snp_filter/bi-snp.filtered.bcf",
    output:
        bcf="results/thin_filter/thinned.bi-snp.filtered.bcf",
    params:
        thin_bp=config['thin_distance_bp'],
    conda:
        "../envs/python.yaml",
    log:
        "results/logs/thin_filter/thin_filter.log"
    shell:
        "python3 workflow/scripts/thin_snps_by_distance_bp.py {input.bcf} {output.bcf} -N {params.thin_bp} "

## get stats for thinned bcf
rule thin_bcf_stats:
    input:
        bcf="results/thin_filter/thinned.bi-snp.filtered.bcf",
    output:
        csi="results/thin_filter/thinned.bi-snp.filtered.bcf.csi",
        stats="results/thin_filter/thinned.bi-snp.filtered.bcf.stats",
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bcftools index {input.bcf} ; "
        "bcftools stats {input.bcf} > {output.stats} "


## thin ngsParalog filtered bcf, if applicable
use rule thin_filter as ngsparalog_thin_filter with:
    input:
        bcf="results/filtered_bcfs/nomismap.bi-snp.filtered.bcf",
    output:
        bcf="results/filtered_bcfs/thin.nomismap.bi-snp.filtered.bcf",

## get stats for thinned ngsParalog filtered bcf, if applicable
use rule thin_bcf_stats as ngsparalog_thin_bcf_stats with:
    input:
        bcf="results/filtered_bcfs/thin.nomismap.bi-snp.filtered.bcf",
    output:
        csi="results/filtered_bcfs/thin.nomismap.bi-snp.filtered.bcf.csi",
        stats="results/filtered_bcfs/thin.nomismap.bi-snp.filtered.bcf.stats",

## calculate per sample stats
## includes heterozyosity, transition, transversion, etc info
rule thin_nommismap_bcf_stats_per_sample:
    input:
        bcf="results/filtered_bcfs/thin.nomismap.bi-snp.filtered.bcf",
    output:
        stats="results/filtered_bcfs/thin.nomismap.bi-snp.filtered.bcf.sample-stats",
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bcftools stats -s - {input.bcf} | grep '^PSC' -B 1 > {output.stats} "  


