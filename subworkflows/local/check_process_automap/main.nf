include { AUTOMAP } from '../../../modules/local/automap/main'
include { CREATE_AUTOMAP_LOG } from '../../../modules/local/create_automap_log/main'

workflow CHECK_PROCESS_AUTOMAP {

    take:
    ch_vcfs        // channel (mandatory): tuple val(meta), path(vcfs), path(tbis)
    automap_assembly // channel (mandatory): val(assembly)        

    main:

    // Run AUTOMAP - failures will be ignored
    AUTOMAP (
        ch_vcfs,
        automap_assembly
    )
    
    // Create a channel with all samples and their status
    // Join input with output to determine success/failure
    status_channel = ch_vcfs
        .map { meta, vcf -> [meta, vcf] }  // Simplify to just meta and vcf
        .join(AUTOMAP.out.roh_automap_file, remainder: true)
        .map { meta, vcf, automap_file -> 
            def status = automap_file ? "SUCCESS" : "FAILED"
            [meta, status, vcf]
        }

    // Create log file
    CREATE_AUTOMAP_LOG (
        status_channel.collect()
    )

    // Emit only successful results
    successful_samples = AUTOMAP.out.roh_automap_file

    emit:
    automap = successful_samples
    log_file = CREATE_AUTOMAP_LOG.out.log_file
}
