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

   /* ch_input = ch_bam.join(ch_intervals, by: 0)
        .map { meta, bam, bai, intervals -> 
            [ meta, bam, bai, intervals.flatten() ]
        } */

    /*ch_input = ch_bam.join(ch_intervals, by: 0)
    .map { meta, bam, bai, intervals -> 
        def interval_path = intervals instanceof Path ? intervals : []
        [ meta, bam, bai, interval_path ]
    }*/
    
    /*ch_input = ch_bam.join(ch_intervals, by: 0)
    .map { meta, bam, bai, intervals -> 
        [ meta, bam, bai, intervals ?: [] ]
    }*/


    /*ch_input = ch_bam.join(ch_intervals, by: 0, remainder: true)
    .map { meta, bam, bai, intervals -> 
        [ meta, bam, bai, intervals ?: [] ]
    }*/

        // Combine ch_bam and ch_intervals
    /*ch_input = ch_bam.combine(ch_intervals, by: 0)
        .map { meta, bam, bai, intervals -> 
            [ meta, bam, bai, intervals ?: [] ]
    }*/

    // Add some debugging
    //ch_input.view { "DEBUG: ch_input: $it" }

    ch_bam.view()
    ch_intervals.view()


    ch_bam.join(ch_intervals).view()

    DEEPVARIANT (
        ch_bam.join(ch_intervals),
        ch_fasta,
        ch_fai,
        [[],[]],
        [[],[]]
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
