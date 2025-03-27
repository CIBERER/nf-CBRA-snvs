include { DEEPVARIANT                                                    }      from '../../../subworkflows/nf-core/deepvariant/main'
include { BCFTOOLS_FILTER                                                }      from '../../../modules/nf-core/bcftools/filter/main'


workflow DEEP_VARIANT_VCF {

    take:
    ch_bam          // channel (mandatory): [ val(meta), path(bam), path(bai) ]
    ch_intervals    // channel (mandatory) : [ val(meta), path(bed) ]
    ch_fasta        // channel (mandatory) : [ val(meta2), path(fasta) ]
    ch_fai          // channel (mandatory) : [ val(meta2), path(fai) ]

    main:

    ch_versions = Channel.empty()


    ch_input = ch_bam
    .join(ch_intervals, by: 0)
    .map { meta, bam, bai, intervals -> 
        if (intervals instanceof List && intervals.isEmpty()) {
            return [ meta, bam, bai, [] ]
        } else {
            return [ meta, bam, bai, intervals ]
        }
    }


    DEEPVARIANT (
        ch_input,
        ch_fasta,
        ch_fai,
        Channel.empty(),
        Channel.empty()
    )
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions.first())

    BCFTOOLS_FILTER(
        DEEPVARIANT.out.vcf.join(DEEPVARIANT.out.vcf_index)
    )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())
 
    vcf = BCFTOOLS_FILTER.out.vcf.join(BCFTOOLS_FILTER.out.tbi)

    emit:
    vcf                                        // channel: [ val(meta), path(vcf), path(tbi)]
    versions = ch_versions                     // channel: [ versions.yml ]

}
