process BCFTOOLS_QUERY_STATS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.15.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.15.1--h0ea216a_0' :
        'quay.io/biocontainers/bcftools:1.15.1--h0ea216a_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*_GT.txt"), emit: gt
    tuple val(meta), path("*_DP.txt"), emit: dp
    tuple val(meta), path("*_VD.txt"), emit: vd
    tuple val(meta), path("*_RD.txt"), emit: rd
    tuple val(meta), path("*_DP_mean.txt"), emit: dp_mean
    tuple val(meta), path("*_VD_mean.txt"), emit: vd_mean
    tuple val(meta), path("*_RD_mean.txt"), emit: rd_mean
    tuple val(meta), path("*_AD_mean.txt"), emit: ad_mean
    tuple val(meta), path("*_VAF.txt"), emit: vaf
    tuple val(meta), path("*_header.txt"), emit: programs
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract sample names from the VCF and split to get program names
    programs=\$(bcftools query -l ${vcf} | sed 's/.*\\.//' | tr '\n' '\t' | sed 's/\\t\$//')

    echo -e "\$programs" > ${prefix}_header.txt

    # Get genotype, depth and variant allele depth
    bcftools query -f '[\\t%GT]\\n' ${vcf} | sed 's/\\t//1' | sed 's/\\.\\/\\.//g' | sed 's/1\\/0/0\\/1/g' > ${prefix}_GT_body.txt
    cat ${prefix}_header.txt ${prefix}_GT_body.txt > ${prefix}_GT.txt

    bcftools query -f '[\\t%DP]\\n' ${vcf} | sed 's/\\t//1' | sed 's/\\.//g' > ${prefix}_DP.txt
    bcftools query -f '[\\t%AD{1}]\\n' ${vcf} | sed 's/\\t//1' | sed 's/\\.//g' > ${prefix}_VD.txt
    bcftools query -f '[\\t%AD{0}]\\n' ${vcf} | sed 's/\\t//1' | sed 's/\\.//g' > ${prefix}_RD.txt

    # Calculate the average depth and variant allele depth
    awk -v OFMT=%.0f '{sum = 0; for (i = 1; i <= NF; i++) sum += \$i; sum /= NF; print sum}' ${prefix}_DP.txt > ${prefix}_DP_mean.txt
    awk -v OFMT=%.0f '{sum = 0; for (i = 1; i <= NF; i++) sum += \$i; sum /= NF; print sum}' ${prefix}_VD.txt > ${prefix}_VD_mean.txt
    awk -v OFMT=%.0f '{sum = 0; for (i = 1; i <= NF; i++) sum += \$i; sum /= NF; print sum}' ${prefix}_RD.txt > ${prefix}_RD_mean.txt
    paste -d, ${prefix}_RD_mean.txt ${prefix}_VD_mean.txt > ${prefix}_AD_mean.txt

    # Calculate variant allele depth (VAD)
    paste ${prefix}_VD_mean.txt ${prefix}_DP_mean.txt | awk -v OFMT=%.2f '{if (\$2 == 0) print "-nan"; else print(\$1/\$2)}' > ${prefix}_VAF.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}