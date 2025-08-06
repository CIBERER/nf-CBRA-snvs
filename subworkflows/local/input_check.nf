//
// Check input samplesheet and get read channels
//

include { samplesheetToList } from 'plugin/nf-schema'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    Channel.fromList(samplesheetToList(samplesheet, "assets/schema_input.json"))
        .map { create_fastq_channel(it) }
        .set { reads }

    Channel.fromList(samplesheetToList(samplesheet, "assets/schema_input.json"))
        .map { create_ped_channel(it) }
        .set { ped }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    ped                                       // channel: [ val(meta), path(ped) ]
    versions = Channel.empty() // SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(ArrayList row) {
    // gather meta
    def meta = row.get(0)
    // Remove family field if it's empty
    if (meta.family == [] || meta.family == null || meta.family == "") {
        meta.remove('family')
    }

    // meta.single_end depending on optional fastq_2 field
    meta.single_end = row.get(2)?.trim() ? false : true

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.get(1)).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.get(1)}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.get(1)) ] ]
    } else {
        if (!file(row.get(2)).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.get(2)}"
        }
        fastq_meta = [ meta, [ file(row.get(1)), file(row.get(2)) ] ]
    }
    return fastq_meta
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_ped_channel(ArrayList row) {
    // gather meta
    def meta = row.get(0)
    
    // Remove family field if it's empty
    // if (meta.family == [] || meta.family == null || meta.family == "") {
    //     meta.remove('family')
    // }

    // meta.single_end depending on optional fastq_2 field
    meta.single_end = row.get(2)?.trim() ? false : true

    // add path(s) of the fastq file(s) to the meta map
    def ped_meta = []

    if (file(row.get(3)).exists()) {
        //check family field if it's empty
        if (meta.family == [] || meta.family == null || meta.family == "") {
            exit 1, "ERROR: Please check input samplesheet -> Family field cannot be empty for trio analysis!\n${row.get(0)}"
        } else {
            def meta_family = [:]
            meta_family.id = meta.family
            ped_meta = [ meta_family,  file(row.get(3))  ]
        }
    }

    return ped_meta
}
