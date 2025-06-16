process AUTOMAP {
    tag "$meta.id"
    label 'process_single'
    errorStrategy 'retry'

    // Conda is not supported
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/yolandabq/automap:1' :
        'docker.io/yolandabq/automap:1' }"



    input:
    tuple val(meta), path(vcf)
    val automap_assembly
    path projectDir

    output:
    tuple val(meta), path("*HomRegions*.t*"), emit: roh_automap_file
    tuple val(meta), path("*HomRegions*.pdf"), emit: roh_automap_pdf, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    if ( task.attempt == 1 )
        """
        if [[ \$(bcftools query -l ${vcf} | wc -l) -gt 1 ]]; then
		    for sample in \$(bcftools query -l ${vcf}); do

			    bcftools view -s \${sample} -O v -o \${sample}.indv.vcf ${vcf}
					
			    bash /app/AutoMap_v1.3.sh \\
			    --vcf \${sample}.indv.vcf \\
			    --out . \\
			    --genome ${automap_assembly}

			    mv \${sample}/*HomRegions* .

		    done
			
	    else
					
		    bcftools view -O v -o ${vcf}.vcf ${vcf}

		    bash /app/AutoMap_v1.3.sh \\
		    --vcf ${vcf}.vcf \\
		    --out . \\
		    --genome ${automap_assembly}

		    mv \$(bcftools query -l ${vcf})/*HomRegions* .

	    fi
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
	    """
    
    else
        """
        echo "Less than 10k variants (with quality)" > ${prefix}_no_automap_HomRegions.txt
		
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
        """
}

