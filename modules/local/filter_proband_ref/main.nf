process FILTER_PROBAND_REF {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.15.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.15.1--h0ea216a_0' :
        'quay.io/biocontainers/bcftools:1.15.1--h0ea216a_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(ped)

    output:
    tuple val(meta), path("*.filtered.vcf.gz"),  path("*.filtered.vcf.gz.tbi"),  emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract proband ID from PED file (assuming it's in the second column)
    proband_id=\$(head -1 ${ped} | cut -f2)
    
    # Debug: Check what samples are in the VCF
    echo "Samples in VCF:"
    bcftools query -l ${vcf}
    echo "Proband ID from PED: \${proband_id}"
    
    # Get the sample index (0-based) for the proband
    sample_index=\$(bcftools query -l ${vcf} | grep -n "^\${proband_id}\$" | cut -d: -f1)
    sample_index=\$((sample_index - 1))  # Convert to 0-based index
    echo "Sample index: \${sample_index}"
    
    # Validate that we found the sample
    if [ "\${sample_index}" -lt 0 ]; then
        echo "ERROR: Sample \${proband_id} not found in VCF"
        exit 1
    fi
    
    # Filter variants where the proband is not homozygous reference (0/0)
    # Using single quotes to avoid escaping issues
    bcftools view -i 'GT['\${sample_index}']!="0/0" && GT['\${sample_index}']!="./."' ${vcf} -Oz -o ${prefix}.filtered.vcf.gz
    
    # Index the filtered VCF
    bcftools index -t ${prefix}.filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}