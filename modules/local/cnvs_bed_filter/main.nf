process CNVS_BED_FILTER {
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0'
    : 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0'}"


    input:
    path cnvs_bed
    path fai
    val min_target
    val chromosomes

    output:
    path("*cnv.bed"), emit: cnvs_bed_filtered

    when:
    task.ext.when == null || task.ext.when

    script:
    def chrom_filter = ''

    if (chromosomes) {
    chrom_filter = chromosomes
        .split(',')
        .collect { "\$1!=\"" + it + "\"" }
        .join(' && ')
    }

    def awk_filter = chrom_filter ?
        "(\$3-\$2)>${min_target} && ${chrom_filter}" :
        "(\$3-\$2)>${min_target}"

    """
    panel="\$(basename ${cnvs_bed} .bed)"

    awk '{if(${awk_filter}){print \$0}}' ${cnvs_bed} \
    | bedtools sort -g ${fai} -i stdin > \${panel}.min${min_target}bp.cnv.bed
    
    """
}