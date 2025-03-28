include { GATK4_CALIBRATEDRAGSTRMODEL                                       } from '../../../modules/nf-core/gatk4/calibratedragstrmodel/main'
include { GATK4_HAPLOTYPECALLER                                             }      from '../../../modules/nf-core/gatk4/haplotypecaller/main'

workflow DRAGEN_VCF {

    take:
    ch_bam         // channel (mandatory): [ val(meta), path(bam), path(bai) ]
    ch_fasta       // channel (mandatory): [ val(meta), path(fasta) ]
    ch_fai         // channel (mandatory): [ val(meta), path(fai) ]
    ch_refdict     // channel (mandatory): [ val(meta), path(dict) ]
    ch_ref_str     // channel (mandatory): [ val(meta), path(ref_str) ]
    ch_intervals    // channel (mandatory) : [ val(meta), path(bed) ]
    ch_dbsnp  // channel (mandatory) : [ val(meta3), path(vcf) ]
    ch_dbsnp_tbi  // channel (mandatory) : [ val(meta3), path(vcf) ]

    main:

    ch_versions = Channel.empty()
    
    GATK4_CALIBRATEDRAGSTRMODEL (
        ch_bam,
        ch_fasta.map { meta, fasta -> fasta },
        ch_fai.map { meta, fai -> fai },
        ch_refdict.map { meta, dict -> dict },
        ch_ref_str //.map { meta, ref_str -> ref_str }
    )
    ch_versions = ch_versions.mix(GATK4_CALIBRATEDRAGSTRMODEL.out.versions.first())

    ch_dragstr_model = GATK4_CALIBRATEDRAGSTRMODEL.out.dragstr_model
    ch_bam.view()
    ch_intervals.view()
    ch_dragstr_model.view()
    ch_bam.join(ch_intervals).join(ch_dragstr_model).view()

    GATK4_HAPLOTYPECALLER (
        ch_bam.join(ch_intervals).join(ch_dragstr_model),
        ch_fasta,
        ch_fai,
        ch_refdict,
        ch_dbsnp,
        ch_dbsnp_tbi
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())

    ch_dragstr_model = GATK4_HAPLOTYPECALLER.out.vcf

    emit:
    dragstr_model = ch_dragstr_model  // channel: [ val(meta), path(vcf) ]
    versions = ch_versions            // channel: [ versions.yml ]

}
