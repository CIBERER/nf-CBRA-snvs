process POSTANNOTSV {
    tag "${meta}"
    label 'process_single'

    container "/mnt/genetica5/singularity_images/bioinfotools_2.0.0.sif"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'docker://docker.io/yolandabq/post_vep:v1' :
    //    'docker.io/yolandabq/post_vep:v1' }"

    input:
    tuple val(meta), path(annotated_cnv), path(colnames)
    path genefilter
    path glowgenes

    output:
    tuple val(meta), path("${prefix}.CNV.annotated.final.tsv"), emit: annotated_cnv
    //tuple val(meta), path("colnames.txt"), emit: colnames

    when:
    task.ext.when == null || task.ext.when

    script:

    prefix = task.ext.prefix ?: "${meta.id}"
		def genefilter_field = genefilter ? "--genefilter ${genefilter} " : ''
		def glowgenes_field  = glowgenes  ? "--glowgenes ${glowgenes} " : ''

    """
		postAnnotsv_modification.R \\
		--input ${annotated_cnv} \\
		--outputfile ${prefix}.CNV.annotated.final.tsv \\
		--extracolnames ${colnames} \\
		${genefilter_field} \\
		${glowgenes_field}

    """
}