process GLOWGENES {
    label 'process_single'
    errorStrategy 'ignore'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'docker://docker.io/tblabfjd/glowgenes:latest' :
       'docker.io/tblabfjd/glowgenes:latest' }"

    input:
    path gene_list


    output:
    path("GLOWgenes_prioritization_Random.txt"), emit: glow_ranking

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python /opt/GLOWgenes/GLOWgenes.py -i ${gene_list} -n /opt/GLOWgenesNets/GLOWgenesNets/networks_knowledgeCategories.cfg -o .
    """
}