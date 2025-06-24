process CREATE_AUTOMAP_LOG {
    label 'process_single'
    
    conda "conda-forge::python=3.9.5 conda-forge::pandas=1.3.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    val status_list

    output:
    path "automap_execution_log.txt", emit: log_file
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    
    from datetime import datetime
    import subprocess
    
    # Parse the status list string
    status_data_str = '''${status_list}'''
    
    # Create log file
    with open("automap_execution_log.txt", "w") as log:
        log.write("AUTOMAP Execution Log\\n")
        log.write("=" * 50 + "\\n")
        log.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")
        log.write("=" * 50 + "\\n\\n")
        
        try:
            # Remove outer brackets and split by '], ['
            entries = status_data_str.strip('[]').split('], [')
            
            for entry in entries:
                entry = entry.strip('[]')
                # Split by comma, but be careful with nested structures
                parts = entry.split(', ', 2)  # Split into max 3 parts
                
                if len(parts) >= 3:
                    meta = parts[0].strip()
                    status = parts[1].strip()
                    vcf = parts[2].strip()
                    
                    log.write(f"SAMPLE: {meta}\\n")
                    log.write(f"STATUS: {status}\\n")
                    log.write(f"VCF: {vcf}\\n\\n")
                    
        except Exception as e:
            log.write(f"Error parsing status list: {e}\\n")
            log.write(f"Raw data: {status_data_str}\\n")
        
        log.write(f"\\nLog generation completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")

    # Create versions file using Python
    try:
        python_version = subprocess.check_output(['python3', '--version'], 
                                                stderr=subprocess.STDOUT, 
                                                text=True).strip().replace('Python ', '')
    except:
        python_version = "unknown"
    
    with open("versions.yml", "w") as versions:
        versions.write(f'"${task.process}":\\n')
        versions.write(f'    python: {python_version}\\n')
    """
}