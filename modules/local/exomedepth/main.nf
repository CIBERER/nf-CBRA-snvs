process EXOMEDEPTH {
    tag "${runname}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'docker://docker.io/yolandabq/exomedepth:latest' :
       'docker.io/yolandabq/exomedepth:latest' }"

    input:
    path bam
    path bai
    path bed
    val runname 

    output:
    tuple val(runname), path("exomedepth*"), emit: cnvs

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    exomeDepth.R -d . -o . -b ${bed} -n ${runname}
    
    """
}