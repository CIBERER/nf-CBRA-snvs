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

    Channel.fromList(samplesheetToList(samplesheet, "assets/schema_input.json"))
        .map { create_bam_channel(it) }
        .set { bams }

    Channel.fromList(samplesheetToList(samplesheet, "assets/schema_input.json"))
        .map { create_vcf_channel(it) }
        .set { vcfs }

    emit:
    reads                       // channel: [ val(meta), [ reads ] ]
    bams                        // channel: [ val(meta), [ bam, bai ] ]
    vcfs                       // channel: [ val(meta), [ vcf, tbi ] ]
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

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []

    if (row.get(1)) {
        if (!file(row.get(1)).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.get(1)}"
        }
        
        // meta.single_end depending on optional fastq_2 field
        meta.single_end = row.get(2) ? false : true

        if (meta.single_end) {
            fastq_meta = [ meta, [ file(row.get(1)) ] ]
        } else {
            if (!file(row.get(2)).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.get(2)}"
            }
            fastq_meta = [ meta, [ file(row.get(1)), file(row.get(2)) ] ]
        }
    } else {
        fastq_meta = [ meta, [] ]
    }

    return fastq_meta
}

// Function to get list of [ meta, [ ped ] ]
// TODO: Try to avoid requiring the PED file to be included for all family members in the samplesheet.
def create_ped_channel(ArrayList row) {
    // gather meta
    def meta = row.get(0)
    
    // Remove family field if it's empty
    // The following block is commented out because removing the 'family' field here may interfere with downstream trio analysis.
    // It is preserved for future consideration as per the TODO above.
    // if (meta.family == [] || meta.family == null || meta.family == "") {
    //     meta.remove('family')
    // }

    // meta.single_end depending on optional fastq_2 field
    meta.single_end = row.get(2) ? false : true

    // add path(s) of the fastq file(s) to the meta map
    def ped_meta = []

if (row.get(7)) {
    if (file(row.get(7)).exists()) {
        //check family field if it's empty
        if (meta.family == [] || meta.family == null || meta.family == "") {
            exit 1, "ERROR: Please check input samplesheet -> Family field cannot be empty for trio analysis!\n${row.get(0)}"
        } else {
            def meta_family = [:]
            meta_family.id = meta.family
            ped_meta = [ meta_family,  file(row.get(7))  ]
        }
    }

    return ped_meta
    }
}

def create_bam_channel(ArrayList row) {
    // gather meta
    def meta = row.get(0)
    def bam_bai_meta = []

    if (meta.family == [] || meta.family == null || meta.family == "") {
        meta.remove('family')
    }

    if (row.get(3)) {
        if (file(row.get(3)).exists()) {
            if (file(row.get(4)).exists()) {
                bam_bai_meta = [ meta, file(row.get(3)), file(row.get(4)) ]
            } else {
                exit 1, "ERROR: Please check input samplesheet -> given bam file but not bai file"
            }
        } else {
            exit 1, "ERROR: Please check input samplesheet -> BAM file does not exist!\n${row.get(3)}"
        }
        return bam_bai_meta
    } else {
        bam_bai_meta = [ meta, [], [] ]
        return bam_bai_meta
    }
}

def create_vcf_channel(ArrayList row) {
    // gather meta
    def meta = row.get(0)
    def vcf_tbi_meta = []

    if (meta.family == [] || meta.family == null || meta.family == "") {
        meta.remove('family')
    }

    if (row.get(5)) {
        if (file(row.get(5)).exists()) {
            if (file(row.get(6)).exists()) {
                vcf_tbi_meta = [ meta, file(row.get(5)), file(row.get(6)) ]
            } else {
                exit 1, "ERROR: Please check input samplesheet -> given vcf file but not tbi file"
            }
        } else {
            exit 1, "ERROR: Please check input samplesheet -> VCF file does not exist!\n${row.get(5)}"
        }
        return vcf_tbi_meta
    } else {
        vcf_tbi_meta = [ meta, [], [] ]
        return vcf_tbi_meta
    }
}