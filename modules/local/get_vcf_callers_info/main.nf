process GET_VCF_CALLERS_INFO {
  
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
      'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(gt_file)

    output:
    tuple val(meta), path("${prefix}_SF.txt"), emit: sf_file
    tuple val(meta), path("*.txt"), emit: column_files
    path  "versions.yml"          , emit: versions


    script:

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Extract the header line
    HEADER=\$(head -n 1 ${gt_file})

    # Count the number of columns
    COLS=\$(echo "\$HEADER" | awk '{print NF}')

    # Initialize an array to store file names
    FILES=()

    # Create files for each column based on the headers
    for ((i=1; i<=\$COLS; i++)); do
      # Get the header name for the current column
      COLUMN_HEADER=\$(echo "\$HEADER" | cut -f \$i)

      # Extract the column data (skip the header row)
      tail -n +2 ${gt_file} | cut -f \$i | awk -v header="\$COLUMN_HEADER" '{if (\$1 != "") print header; else print ""}' > "\${COLUMN_HEADER}.txt"

      # Add the file name to the array
      FILES+=("\${COLUMN_HEADER}.txt")
    done

    # Use paste to combine all the generated files and format the output with underscores
    paste -d "\t" "\${FILES[@]}" | perl -alne 'print join "_", @F' > ${prefix}_SF.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl --version | head -n 2 | tail -n 1 | cut -d' ' -f4)
    END_VERSIONS
    """
}