include { SPLITMULTIALLELIC                 } from '../../../modules/local/splitmultiallelic/main'
include { BCFTOOLS_MERGE } from '../../../modules/nf-core/bcftools/merge/main'
include { TABIX_TABIX } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_MERGE_VARIANTCALLERS {

    take:
    ch_vcfs        // channel (mandatory): tuple val(meta), path(vcfs), path(tbis)
    ch_fasta        // channel (mandatory) : [ val(meta3), path(fasta) ]
    //ch_fai          // channel (mandatory) : [ val(meta3), path(fai) ]
    //ch_intervals    // channel (mandatory) : [ val(meta), path(bed) ]
    //ch_index        // channel (mandatory): [ val(meta2), path(index) ]

    main:

    ch_versions = Channel.empty()

    ch_vcfs_for_splitmultiallelic = ch_vcfs.transpose().map{ meta, vcf, tbi, program -> [meta, vcf, tbi, program]}.view()
    //he conseguido que el canal haga lo que quiero, pero los vcfs de prueba que estoy usando no funcionan. 
    
    SPLITMULTIALLELIC (
        ch_vcfs_for_splitmultiallelic,
        ch_fasta
    )
    
    ch_versions = ch_versions.mix(SPLITMULTIALLELIC.out.versions.first())

    vcf = SPLITMULTIALLELIC.out.biallelic_renamed_vcf

    emit:
    vcf  // channel: [ [val(meta)], path(bam), path(bai)]
    versions = ch_versions                     // channel: [ versions.yml ]

}
