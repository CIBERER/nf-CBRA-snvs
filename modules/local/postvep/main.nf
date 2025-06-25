process POSTVEP {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/yolandabq/post_vep:v1' :
        'docker.io/yolandabq/post_vep:v1' }"

    input:
    tuple val(meta), path(vep_tsv), path(roh_automap)
    val maf
    val assembly

    output:
    tuple val(meta), path("*.SNV.INDEL.annotated.tsv"), emit: pvm_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def automap = roh_automap ? "--automap ${roh_automap}" : ''

    """
    postVEP_modification.R \\
    --input ${vep_tsv} \\
    --output ${prefix}.${assembly}.SNV.INDEL.annotated.tsv \\
    --maf ${maf} \\
    ${automap}

    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}