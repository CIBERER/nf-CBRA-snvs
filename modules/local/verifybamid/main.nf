
process VERIFYBAMID {
    tag "$meta.id"
    label 'VERIFYBAMID'

    errorStrategy 'retry'
    maxRetries 1

    conda "${moduleDir}/environment.yml"
    publishDir "$params.outdir/$meta.id/verifyBamID/", mode: 'copy'

    input:    
    tuple val(meta), path(bam), path(bai)
    path(fasta)

    output:
    tuple val(meta), path("*.{selfSM,Ancestry}"), emit: files
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (task.attempt == 1) {
        """
        verifybamid2 --SVDPrefix ${params.SVDPrefix} --Reference $fasta --BamFile $bam --Output $prefix --NumThread $task.cpus

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            verifybamid2: \$(verifybamid2 --version 2>&1 | grep "Version" | sed 's/.*Version: //')
        END_VERSIONS
        """
    } 
    else{
        """  
        echo "VerifyBamID can't be executed." > ${prefix}.selfSM
        echo "VerifyBamID can't be executed." > ${prefix}.Ancestry

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            verifybamid2: \$(verifybamid2 --version 2>&1 | grep "Version" | sed 's/.*Version: //')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """  
    touch ${prefix}.selfSM
    touch ${prefix}.Ancestry

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verifybamid2: \$(verifybamid2 --version 2>&1 | grep "Version" | sed 's/.*Version: //')
    END_VERSIONS
    """
}
