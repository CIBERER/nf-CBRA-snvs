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
    tuple val(meta), path("${prefix}.${assembly}.final.vcf.gz"), emit: final_vcf
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

    echo "##FORMAT=<ID=SF,Number=1,Type=String,Description=\\"Software\\">" >> ${prefix}.${assembly}.final.vcf
    echo "##FORMAT=<ID=GD,Number=1,Type=String,Description=\\"Genotype discordances. 0 same genotype and 1 different genotype\\">" >> ${prefix}.${assembly}.final.vcf

    # Loop through each program in the list
    for PROGRAM in \$PROGRAMS; do
      # Append the desired FORMAT lines to the temporary file
      echo "##FORMAT=<ID=\${PROGRAM}_GT,Number=1,Type=String,Description=\"\${PROGRAM} genotype\">" >> ${prefix}.${assembly}.final.vcf
      echo "##FORMAT=<ID=\${PROGRAM}_DP,Number=1,Type=String,Description=\"\${PROGRAM} depth\">" >> ${prefix}.${assembly}.final.vcf
      echo "##FORMAT=<ID=\${PROGRAM}_VD,Number=1,Type=String,Description=\"\${PROGRAM} Variant frequency\">" >> ${prefix}.${assembly}.final.vcf
    done

    bcftools view -h ${vcf} | tail -n 1 | cut -f 1-10 | cut -d"." -f1 >> ${prefix}.${assembly}.final.vcf ## vcf header

    paste -d "\\t" ${prefix}_VCF_CONTENT.txt ${prefix}_FORMAT_SAMPLE.txt >> ${prefix}.${assembly}.final.vcf ## add the content of the vcf to the final vcf
    bgzip -c ${prefix}.${assembly}.final.vcf > ${prefix}.${assembly}.final.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}