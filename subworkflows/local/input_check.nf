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

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = Channel.empty() // SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(ArrayList row) {
    // gather meta
    meta            = row.get(0)
    // meta.single_end depending on optional fastq_2 field
    meta.single_end = row.get(2)?.trim() ? true : false

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
