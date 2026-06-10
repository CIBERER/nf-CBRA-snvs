process POSTANNOTSV {
    tag "${meta}"
    label 'process_single'
    errorStrategy 'ignore'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'docker://docker.io/tblabfjd/postannotsv:latest' :
       'docker.io/tblabfjd/postannotsv:latest' }"

       
    input:
    tuple val(meta), path(annotated_cnv), path(colnames)
    path genefilter
    path glowgenes

    output:
    tuple val(meta), path("${prefix}.CNV.annotated.final.tsv"), emit: annotated_cnv

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