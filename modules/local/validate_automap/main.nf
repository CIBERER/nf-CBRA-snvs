process VALIDATE_AUTOMAP {
    tag "$meta.id"
    input:
    tuple val(meta), path(vcf)
   
    output:
    tuple val(meta), path(vcf), env('VALIDATION_STATUS'), emit: validate_automap
   
    script:
    """
    line_count=\$(zgrep -v "#" ${vcf} | wc -l)
    if [ \$line_count -ge 1000 ]; then
        VALIDATION_STATUS="PASS"
        echo "File validation passed"
    else
        VALIDATION_STATUS="FAIL"
        echo "File validation failed: only \$line_count lines"
    fi
    """
}
