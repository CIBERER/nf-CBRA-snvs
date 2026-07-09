process FILTER_MAF {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::htslib=1.21 conda-forge::gawk=5.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/htslib_gawk:1fe5bc0e7b4bb666' :
        'community.wave.seqera.io/library/htslib_gawk:1fe5bc0e7b4bb666' }"

    input:
    tuple val(meta), path(vep_tab)
    val(maf)

    output:
    tuple val(meta), path("*.filtered.tab.gz"), emit: maf_filtered

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    maf_filter.sh ${vep_tab} ${prefix}.filtered.tab.gz ${maf}
    """
}