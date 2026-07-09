process POSTVEP {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/yolandabq/post_vep:v1' :
        'docker.io/yolandabq/post_vep:v1' }"

    input:
    tuple val(meta), path(vep_tsv), path(roh_automap)
    path(pvm_script)
    val maf
    val assembly
    path glowgenes_ranking
    path glowgenes_sgds
    path gene_list
    path extra_files_pvm 
    

    output:
    tuple val(meta), path("*.SNV.INDEL.annotated.tsv"), emit: pvm_tsv
    

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def automap = roh_automap ? "--automap '${roh_automap}'" : ''
    def glowgenes = glowgenes_ranking ? "--glowgenes ${glowgenes_ranking}" : ''
    def sgds = glowgenes_sgds ? "--SGDS ${glowgenes_sgds}" : ''
    def gene_filter = gene_list ? "--genefilter ${gene_list}" : ''

    """

    ${pvm_script} \\
    --input ${vep_tsv} \\
    --output ${prefix}.${assembly}.SNV.INDEL.annotated.tsv \\
    --maf ${maf} \\
    ${automap} \\
    ${glowgenes} \\
    ${sgds} \\
    ${gene_filter} \\
    ${args}

    """
}