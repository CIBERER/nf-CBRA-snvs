process VALIDATE_AUTOMAP {
    tag "$meta.id"
    input:
    tuple val(meta), path(vcf)
   
    output:
    tuple val(meta), path(vcf), env('VALIDATION_STATUS'), emit: validate_automap
   
    script:
    
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    line_count=\$(zgrep -v "#" ${vcf} | wc -l)
    if [ \$line_count -ge 1000 ]; then
        VALIDATION_STATUS="PASS"
        echo "File validation passed" > ${prefix}.automap.validation.log
    else
        VALIDATION_STATUS="FAIL"
        echo "File validation failed: only \$line_count lines" > ${prefix}.automap.validation.fail.log
    fi
    """
}
