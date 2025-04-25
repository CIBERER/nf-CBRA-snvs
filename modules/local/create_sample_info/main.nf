process CREATE_SAMPLE_INFO {
    tag "${meta.id}_${meta.program}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(GT_consensus), path(GT_discordances), path(AD_mean), path(DP_mean), path(RD_mean), path(VD_mean), path(VAF), path(SF), path(programs)
    val (assembly)
    
    output:
    tuple val(meta), path("${prefix}_FORMAT_SAMPLE.txt"), emit: final_vcf
    //tuple val(meta), path("*.final.vcf.gz"), path("*.final.vcf.gz.tbi"), emit: final_vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools query -f '[:%GT:%DP:%AD{1}]\\n' ${vcf} > ${prefix}_FORMAT_SF.txt
    paste -d ":" ${GT_consensus} ${AD_mean} ${DP_mean} ${VAF} ${SF} ${GT_discordances} > ${prefix}_FORMAT_JOIN.txt
    paste ${prefix}_FORMAT_JOIN.txt ${prefix}_FORMAT_SF.txt | tr -d '\\t' > ${prefix}_FORMAT_SAMPLE.txt


    # Read the programs from the file and format the output
    PROGRAMS=\$(head -n 1 ${programs})
    OUTPUT=()
    for PROGRAM in \$PROGRAMS; do
        OUTPUT+=("\${PROGRAM}_GT:\${PROGRAM}_DP:\${PROGRAM}_VD")
    done
    format_string=\$(echo "GT:AD:DP:VAF:SF:GD:\${OUTPUT[*]}" | tr ' ' ':')

    bcftools view -H ${vcf} | cut -f 1-8 | sed "s/\$/\t\$format_string/" > ${prefix}_VCF_CONTENT.txt

    bcftools view -h ${vcf} | grep "##" > ${prefix}.${assembly}.final.vcf

    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}