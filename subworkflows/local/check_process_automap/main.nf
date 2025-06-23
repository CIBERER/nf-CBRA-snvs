include { AUTOMAP } from '../../../modules/local/automap/main'
include { VALIDATE_AUTOMAP } from '../../../modules/local/validate_automap/main'

workflow CHECK_PROCESS_AUTOMAP {

    take:
    ch_vcfs        // channel (mandatory): tuple val(meta), path(vcfs), path(tbis)
    automap_assembly // channel (mandatory): val(assembly)        

    main:

    VALIDATE_AUTOMAP (
        ch_vcfs
    )
    
    VALIDATE_AUTOMAP.out.validate_automap.view()

    // Filter only files that passed validation
    valid_files = VALIDATE_AUTOMAP.out.validate_automap
        | filter { meta, file, status -> status == "PASS" }
        | map { meta, file, status -> [meta, file] }

    AUTOMAP (
        valid_files,
        automap_assembly
    )

    automap = AUTOMAP.out.roh_automap_file

    emit:
    automap                                        // channel: [ [val(meta)], path(tsv)]

}
