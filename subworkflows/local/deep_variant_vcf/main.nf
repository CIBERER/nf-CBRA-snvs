include { DEEPVARIANT                                                    }      from '../../../subworkflows/nf-core/deepvariant/main'
include { BCFTOOLS_FILTER                                                }      from '../../../modules/nf-core/bcftools/filter/main'
include { SPLITMULTIALLELIC                                                }      from '../../../modules/local/splitmultiallelic/main'

workflow DEEP_VARIANT_VCF {

    take:
    ch_bam          // channel (mandatory): [ val(meta), path(bam), path(bai) ]
    ch_intervals    // channel (mandatory) : [ val(meta), path(bed) ]
    ch_fasta        // channel (mandatory) : [ val(meta2), path(fasta) ]
    ch_fai          // channel (mandatory) : [ val(meta2), path(fai) ]
    ch_gzi          // channel: [ val(meta4), path(gzi) ]
    ch_par_bed      // channel: [ val(meta5), path(par_bed) ]

    main:

    ch_versions = Channel.empty()

    ch_input_for_deepvariant = ch_bam.join(ch_intervals)

    DEEPVARIANT (
        ch_input_for_deepvariant,
        ch_fasta,
        ch_fai,
        ch_gzi,
        ch_par_bed
    )
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions.first())

    BCFTOOLS_FILTER(
        DEEPVARIANT.out.vcf.join(DEEPVARIANT.out.vcf_index)
    )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())

    SPLITMULTIALLELIC (
        BCFTOOLS_FILTER.out.vcf.join(BCFTOOLS_FILTER.out.tbi),
        ch_fasta,
        "deepvariant"
    )
    ch_versions = ch_versions.mix(SPLITMULTIALLELIC.out.versions.first())

    vcf = SPLITMULTIALLELIC.out.biallelic_renamed_vcf

    emit:
    vcf                                        // channel: [ val(meta), path(vcf), path(tbi)]
    versions = ch_versions                     // channel: [ versions.yml ]

}
