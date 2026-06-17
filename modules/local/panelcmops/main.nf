process PANELCNMOPS {
    tag "${runname}"
    label 'process_single'
    errorStrategy 'ignore'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'docker://docker.io/yolandabq/panelcnmops:latest' :
       'docker.io/yolandabq/panelcnmops:latest' }"
    
    input:
    path bam
    path bai
    path bed
    val runname 

    output:
    tuple val(runname), path("panelcn.MOPS*"), emit: cnvs

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    panelcnMops.R -d . -o . -b ${bed} -n ${runname}
    """
}