process CNV_RESULT_MIXER {
    tag "${meta}"
    label 'process_single'

    container "/mnt/genetica5/singularity_images/bioinfotools_2.0.0.sif"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'docker://docker.io/yolandabq/post_vep:v1' :
    //    'docker.io/yolandabq/post_vep:v1' }"

    input:
    tuple val(meta), path(cnvs)
    path samples2analyce

    output:
    tuple val(meta), path("${meta}.CNV.merged.bed"), emit: merged_bed

    when:
    task.ext.when == null || task.ext.when

    script:

    def samples2analyce_field   = samples2analyce ? "--samples ${samples2analyce} " : ""

    """
    CNV_result_mixer.R \
		--inputdir . \
		--outputfile ${meta}.CNV.merged.bed \
		${samples2analyce_field}
    """
}