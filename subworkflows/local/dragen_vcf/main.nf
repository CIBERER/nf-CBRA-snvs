include { GATK4_CALIBRATEDRAGSTRMODEL } from '../../../modules/nf-core/gatk4/calibratedragstrmodel/main'

workflow DRAGEN_VCF {

    take:
    ch_bam         // channel (mandatory): [ val(meta), path(bam), path(bai) ]
    ch_fasta       // channel (mandatory): [ val(meta), path(fasta) ]
    ch_fai         // channel (mandatory): [ val(meta), path(fai) ]
    ch_refdict     // channel (mandatory): [ val(meta), path(dict) ]
    ch_ref_str     // channel (mandatory): [ val(meta), path(ref_str) ]

    main:

    ch_versions = Channel.empty()
    
    GATK4_CALIBRATEDRAGSTRMODEL (
        ch_bam,
        ch_fasta.map { meta, fasta -> fasta },
        ch_fai.map { meta, fai -> fai },
        ch_refdict.map { meta, dict -> dict },
        ch_ref_str.map { meta, ref_str -> ref_str }
    )
    ch_versions = ch_versions.mix(GATK4_CALIBRATEDRAGSTRMODEL.out.versions.first())

    ch_dragstr_model = GATK4_CALIBRATEDRAGSTRMODEL.out.dragstr_model

    emit:
    dragstr_model = ch_dragstr_model  // channel: [ val(meta), path(dragstr_model) ]
    versions = ch_versions            // channel: [ versions.yml ]

}
