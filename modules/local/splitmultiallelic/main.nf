process SPLITMULTIALLELIC {
    tag "${meta.id}_${meta.program}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(fasta)
    val(program)

    output:
    tuple val(meta), path("${prefix}.biallelic.${program}.vcf.gz"), path("${prefix}.biallelic.${program}.vcf.gz.tbi"), emit: biallelic_renamed_vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    //program = task.ext.prefix ?: "${meta.program}"

    """
    bcftools query -l ${vcf} > samples.${program}.txt
    sed -e 's/\$/.${program}/' -i samples.${program}.txt

    bcftools reheader --samples samples.${program}.txt -o ${prefix}.renamed.${program}.vcf.gz ${vcf}
    tabix -p vcf ${prefix}.renamed.${program}.vcf.gz

    bcftools norm -m - -c s -f ${fasta} ${prefix}.renamed.${program}.vcf.gz | bcftools view --min-ac=1 -O z -o ${prefix}.biallelic.${program}.vcf.gz 
    tabix -p vcf ${prefix}.biallelic.${program}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}