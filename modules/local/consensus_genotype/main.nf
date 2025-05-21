process CONSENSUS_GENOTYPE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9.5 conda-forge::pandas=1.3.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(gt_file)

    output:
    tuple val(meta), path("*_consensus_GT.txt"), emit: consensus_gt
    tuple val(meta), path("*_discordances.txt"), emit: discordances
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    import pandas as pd
    import sys

    df = pd.read_table("${gt_file}", sep="\t")

    # Count NAs in each column and sort
    na_counts = df.isna().sum()
    na_counts = na_counts.sort_values(ascending=False)
    # Convert Series to DataFrame
    df_na_counts = na_counts.reset_index()
    df_na_counts.columns = ['Column', 'NA_Count']

    # Each caller has a different weight so that in case of tie, the order of prioritization goes from the one that got the most calls
    # Add a new column with values starting from 2 and increasing by 1
    df_na_counts['weight'] = range(2, 2 + len(df_na_counts))

    df_new = pd.DataFrame()

    for program in list(df.columns):
        df_weighted = pd.concat([df[program]] * list(df_na_counts[df_na_counts["Column"] == program]["weight"])[0], axis=1, ignore_index=True)
        df_new = pd.concat([df_new, df_weighted], axis = 1)

    # Write consensus Genotype
    df_new.mode(axis=1).to_csv("${prefix}_consensus_GT.txt", sep='\t', header=False, index=False)

    # Calculate and write discordances in genotype
    discordance_list = list(map(lambda x: len(set(filter(lambda y: y == y, x))) - 1, df.values))
    with open("${prefix}_discordances.txt", 'w') as fp:
        for item in discordance_list:
            fp.write(f"{item}\\n")

    with open("versions.yml", "w") as f:
        f.write('''
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    ''')
    """
}