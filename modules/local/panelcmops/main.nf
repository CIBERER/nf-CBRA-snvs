process PANELCNMOPS {
    tag "${runname}"
    label 'process_single'

    container "/mnt/genetica5/singularity_images/bioinfotools_2.0.0.sif"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'docker://docker.io/yolandabq/post_vep:v1' :
    //    'docker.io/yolandabq/post_vep:v1' }"

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