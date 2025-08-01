process FORMAT2INFO {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}.vcf_to_annotate.vcf.gz"), emit: vcf_to_annotate
    tuple val(meta), path("${prefix}.fields.txt"), emit: fields
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    FORMAT=\$(gunzip -c ${vcf} | awk '/^#[^#]/{getline; gsub(":", " ", \$9); print \$9}') 

	bcftools view -h ${vcf} | grep "##" > ${prefix}.vcf_to_annotate.vcf
	echo "##INFO=<ID=variant_id,Number=.,Type=String,Description=\\"variant identification\\">" >> ${prefix}.vcf_to_annotate.vcf
	echo "##INFO=<ID=Original_pos,Number=.,Type=String,Description=\\"original position\\">" >> ${prefix}.vcf_to_annotate.vcf
	
    for sample in \$(bcftools query -l ${vcf})
	do
		for field in \${FORMAT[@]}
		do
			echo "##INFO=<ID=\${sample}_\${field},Number=.,Type=String,Description=\\"\${sample} \${field}\\">" >> ${prefix}.vcf_to_annotate.vcf
		done
	done
	
    bcftools view -h ${vcf} | grep "#CHROM" | cut -f1-8 >> ${prefix}.vcf_to_annotate.vcf

	newfields=""
	for field in \${FORMAT[@]}; do newfields="\$(echo "\${newfields}[;%SAMPLE\\_\${field}=%\${field}]")"; done
	bcftools query -f "variant_id=%CHROM\\_%POS\\_%REF\\_%ALT;Original_pos=%POS;\${newfields}\\n" -u ${vcf} | sed 's/;//2' | sed 's/,/_/g' > new_info.txt
	bcftools view -H ${vcf} | cut -f1-8 > old_info.txt
	paste -d ';' old_info.txt new_info.txt >> ${prefix}.vcf_to_annotate.vcf
	
	if bcftools view -h ${vcf} | grep -Fq hiConfDeNovo; then
		fields=",hiConfDeNovo,loConfDeNovo,variant_id,Original_pos"
	else
		fields=",variant_id,Original_pos"
	fi
		 
	for sample in \$(bcftools query -l ${vcf}); do for field in \${FORMAT[@]}; do fields="\$(echo "\${fields},\${sample}_\${field}")"; done; done
	echo \${fields} > ${prefix}.fields.txt

    bgzip -c ${prefix}.vcf_to_annotate.vcf > ${prefix}.vcf_to_annotate.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}