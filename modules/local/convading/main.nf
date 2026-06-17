process CONVADING {
    tag "${runname}"
    label 'process_high'
    errorStrategy 'ignore'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'docker://docker.io/yolandabq/convading:latest' :
       'docker.io/yolandabq/convading:latest' }"

    input:
    path bam
    path bai
    path bed
    path fai
    val runname 


    output:
    tuple val(runname), path("CoNVaDING*"), emit: cnvs

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    cp ${moduleDir}/resources/usr/bin/CoNVaDING.py .
    cp ${moduleDir}/resources/usr/bin/CoNVaDING.pl .
    CoNVading_pipeline.py . ${bed} ./ ${runname} ${fai}
    
    """
}