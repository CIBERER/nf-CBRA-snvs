process ADD_VAF_TRIO {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.15.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.15.1--h0ea216a_0' :
        'quay.io/biocontainers/bcftools:1.15.1--h0ea216a_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:

    tuple val(meta), path("*.final.vcf.gz"), path("*.final.vcf.gz.tbi"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    bcftools query -f '%FORMAT\n' ${vcf} | cut -f 1 | awk '{print \$0 ":VAF"}'  > ${prefix}_FORMAT.txt

    for sample in \$(bcftools query -l ${vcf})
	do
        bcftools view -s \${sample} ${vcf} | bcftools query -f '%FORMAT\\n' | cut -f 2 > \${sample}_FORMAT.txt

		bcftools view -s \${sample} ${vcf} | bcftools query -f '[\\t%DP]\\n' | sed 's/\\t//1' | sed 's/\\.//g' > \${sample}_DP.txt
        bcftools view -s \${sample} ${vcf} | bcftools query -f '[\\t%AD{1}]\\n'  | sed 's/\\t//1' | sed 's/\\.//g' > \${sample}_VD.txt
        bcftools view -s \${sample} ${vcf} | bcftools query -f '[\\t%AD{0}]\\n'  | sed 's/\\t//1' | sed 's/\\.//g' > \${sample}_RD.txt

        # Calculate variant allele depth (VAD)
        paste \${sample}_VD.txt \${sample}_DP.txt | \
            awk -v OFMT=%.2f '{
                vd = (\$1 == "" ? 0 : \$1)
                dp = (\$2 == "" ? 0 : \$2)
                if (dp == 0) print "-nan"; else print(vd/dp)
            }' > \${sample}_VAF.txt

        paste -d ":" \${sample}_FORMAT.txt \${sample}_VAF.txt | tr -d '\\t' > \${sample}_FORMAT_SAMPLE.txt

        paste ${prefix}_FORMAT.txt \${sample}_FORMAT_SAMPLE.txt > ${prefix}_FORMAT_tmp.txt
        cat ${prefix}_FORMAT_tmp.txt > ${prefix}_FORMAT.txt

	done

    bcftools view -H ${vcf} | cut -f 1-8 > ${prefix}_VCF_CONTENT.txt

    bcftools view -h ${vcf} | grep "##" > ${prefix}.final.vcf
    echo "##FORMAT=<ID=VAF,Number=1,Type=String,Description=\\"Variant Allele Frequency\\">" >> ${prefix}.final.vcf
    bcftools view -h  ${vcf} | tail -1 >> ${prefix}.final.vcf

    paste -d "\\t" ${prefix}_VCF_CONTENT.txt ${prefix}_FORMAT.txt >> ${prefix}.final.vcf ## add the content of the vcf to the final vcf
    bgzip -c ${prefix}.final.vcf > ${prefix}.final.vcf.gz
    bcftools index -t ${prefix}.final.vcf.gz

    """
}