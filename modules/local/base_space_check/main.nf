process BS_CHECK {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    val meta
    val project
    val baseuser
    path samples
    val analysis

    output:
    tuple val(meta), path("projects.txt")        , emit: bsproyects
    tuple val(meta), path("controlsamples.txt")  , emit: controlsamples
    tuple val(meta), path("samples2analyce.txt") , emit: samples2analyce
    tuple val(meta), path("datasets.txt")        , emit: datasets
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def baseuser_config = baseuser ? "--config ${baseuser} " : ''

    if(samples && analysis.contains("CNV"))

        """

        // controlsamples: all samples need to be mapped for CNV calling 
	    // samples2analyce: results (SNVs and CNVs) are only reported for the specified sample(s)

	    // List all projects in the BaseSpace account
	    ${baseuser}bs list projects -f csv -F Name > projects.txt

	    // Check if the given project exist and if so, get all the sample names
        if grep -Fq ${project} projects.txt; then
            ${baseuser}bs list biosample --sort-by=BioSampleName -f csv -F BioSampleName --project-name=${project} | sort | uniq > controlsamples.txt    
	    else
	    	>&2 echo "ERROR: Project '${project}' does not exist in basespace\n"
	    	exit 1
	    fi

	    // Check that the specified samples exist inside the project 
	    for sample in \$(cat ${samples}); do
	    	if grep -q \${sample} controlsamples.txt; then
	    		grep \${sample} controlsamples.txt >> samples2analyce.txt
	    	else
	    		>&2 echo "ERROR: Sample '\${sample}' does not exist in the project '${project}'\n"
	    		exit 1
	    	fi
	    done


	    // Get the datasets id. 
	    // First the last appsession Id is retrieved. 
	    // Appsessions are the analysis done in a project. 
	    // We assume that these are basecalling and that the last one is the correct one.
	    appsession_id=\$(${baseuser}bs list appsession -f csv -F Id --project-name "${project}" | tail -n 1)

	    // List all the Output.Datasets (folders containing the reads per sample and per lane) and 
	    // filter to keep the ones containing the pattern *_L* to avoid duplicates.
	    // edit elby Graci on 21/05/2025: new approach for downloading: "grep _L ..." is not needed anymore, BUT for old cases, second line is the good one
	    ${baseuser}bs appsession property get -i "\${appsession_id}" --property-name="Output.Datasets" -f csv -F Id -F Name > datasets.txt
	    // ${baseuser}bs appsession property get -i "\${appsession_id}" --property-name="Output.Datasets" -f csv -F Id -F Name | grep "_L" | grep -v "Undetermined" > datasets.txt
        """

    else if(samples) 

        """
		// controlsamples: only the specifies sammple(s) is(are) mapped 
		// samples2analyce: results (SNVs) are only reported for the specified sample(s)

		// List all projects in the BaseSpace account
		${baseuser}bs list projects -f csv -F Name > projects.txt


		// Check if the given project exist and if so, get all the sample names
		if grep -Fq ${project} projects.txt; then
			${baseuser}bs list biosample --sort-by=BioSampleName -f csv -F BioSampleName --project-name=${project} | sort | uniq > controlsamples.txt    
		else
			>&2 echo "ERROR: Project '${project}' does not exist in basespace\n"
			exit 1
		fi

		// Check that the specified samples exist inside the project
		for sample in \$(cat ${samples}); do
			if grep -q \${sample} controlsamples.txt; then
				grep \${sample} controlsamples.txt >> samples2analyce.txt
			else
				>&2 echo "ERROR: Sample '\${sample}' does not exist in the project '${project}'\n"
				exit 1
			fi   
		done

		// controlsamples are the same ones as samples2analyce
		cat samples2analyce.txt > controlsamples.txt


		// Get the datasets id. 
		// First the last appsession Id is retrieved. 
		// Appsessions are the analysis done in a project. 
		// We assume that these are basecalling and that the last one is the correct one.
		appsession_id=\$(${baseuser}bs list appsession -f csv -F Id --project-name "${project}" | tail -n 1)

		// List all the Output.Datasets (folders containing the reads per sample and per lane) and 
		// filter to keep the ones containing the pattern *_L* to avoid duplicates.
		// edit el 21/05/2025: cambia la manera de descargar y ya no hay que hacer grep _L y tal, si se quieren descargar algunos antiguos igual si hace falta 
		${baseuser}bs appsession property get -i "\${appsession_id}" --property-name="Output.Datasets" -f csv -F Id -F Name > datasets.txt
		// ${baseuser}bs appsession property get -i "\${appsession_id}" --property-name="Output.Datasets" -f csv -F Id -F Name | grep "_L" | grep -v "Undetermined" > datasets.txt
		"""

    else 
		"""
		// controlsamples: all samples need to be mapped for SNV and CNV calling 
		// samples2analyce: results (SNVs and CNVs) are reported for all samples

		// List all projects in the BaseSpace account
		${baseuser}bs list projects -f csv -F Name > projects.txt

		// Check if the given project exist and if so, get all the sample names
		if grep -Fq ${project} projects.txt; then
			${baseuser}bs list biosample --sort-by=BioSampleName -f csv -F BioSampleName --project-name=${project} | sort | uniq > controlsamples.txt    
		else
			>&2 echo "ERROR: Project '${project}' does not exist in basespace\n"
			exit 1
		fi

		// samples2analyce are the same ones as controlsamples
		cat controlsamples.txt > samples2analyce.txt


		// Get the datasets id. 
		// First the last appsession Id is retrieved. 
		// Appsessions are the analysis done in a project. 
		// We assume that these are basecalling and that the last one is the correct one.
		appsession_id=\$(${baseuser}bs list appsession -f csv -F Id --project-name "${project}" | tail -n 1)

		// List all the Output.Datasets (folders containing the reads per sample and per lane) and 
		// filter to keep the ones containing the pattern *_L* to avoid duplicates.
		// edit el 21/05/2025: cambia la manera de descargar y ya no hay que hacer grep _L y tal, si se quieren descargar algunos antiguos igual si hace falta 
		${baseuser}bs appsession property get -i "\${appsession_id}" --property-name="Output.Datasets" -f csv -F Id -F Name > datasets.txt
		// ${baseuser}bs appsession property get -i "\${appsession_id}" --property-name="Output.Datasets" -f csv -F Id -F Name | grep "_L" | grep -v "Undetermined" > datasets.txt
		"""

    END_VERSIONS
    """

