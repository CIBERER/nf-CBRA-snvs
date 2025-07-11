process VCF_FILTER_MAF {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
      'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(vep_tab)
    val maf

    output:
    tuple val(meta), path("${prefix}.maf.filtered.tab"), emit: maf_filtered_tab

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    zcat ${vep_tab} | awk -v thresh=${maf} '
    BEGIN { FS = OFS = "\\t"; header_found = 0 }
    # Keep all header lines
    /^#/ {
        print
        # Detect column header (single # at beginning)
        if (!header_found && \$0 ~ /^#[^#]/) {
            header_found = 1
            for (i=1; i<=NF; i++) {
                if (\$i == "MAX_AF") {
                    max_af_col = i
                    break
                }
            }
            if (!max_af_col) {
                print "ERROR: MAX_AF column not found." > "/dev/stderr"
                exit 1
            }
        }
        next
    }
    # Filter data rows based on MAX_AF threshold
    {
        af = \$max_af_col == "." || \$max_af_col == "" ? 0 : \$max_af_col
        if (af < thresh) print
    }
    ' > ${prefix}.maf.filtered.tab

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.maf.filtered.tab

    """
}