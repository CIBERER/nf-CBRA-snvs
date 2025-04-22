include { SPLITMULTIALLELIC                 } from '../../../modules/local/splitmultiallelic/main'
include { BCFTOOLS_MERGE } from '../../../modules/nf-core/bcftools/merge/main'
include { TABIX_TABIX } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_MERGE_VARIANTCALLERS {

    take:
    ch_vcfs        // channel (mandatory): tuple val(meta), path(vcfs), path(tbis), val(program)
    ch_fasta        // channel (mandatory) : [ val(meta3), path(fasta) ]
    //ch_fai          // channel (mandatory) : [ val(meta3), path(fai) ]
    //ch_intervals    // channel (mandatory) : [ val(meta), path(bed) ]
    //ch_index        // channel (mandatory): [ val(meta2), path(index) ]

    main:

    ch_versions = Channel.empty()

    //ch_vcfs_for_splitmultiallelic = ch_vcfs.transpose().map{ meta, vcf, tbi, program -> [meta, vcf, tbi, program]}.view()
    //he conseguido que el canal haga lo que quiero, pero los vcfs de prueba que estoy usando no funcionan. 
    ch_vcfs.view()

    // Split the input channel by program
    ch_vcfs_for_splitmultiallelic = ch_vcfs.flatMap { meta, vcf, tbi, programs ->
        programs.collect { program ->
            def new_meta = meta + [program: program]
            [ new_meta, vcf, tbi, program ]
        }
    }


    ch_vcfs_for_splitmultiallelic.view()

    SPLITMULTIALLELIC (
        ch_vcfs_for_splitmultiallelic,
        ch_fasta
    )

    // ahora da error porque Process `VCF_MERGE_VARIANTCALLERS:SPLITMULTIALLELIC` input file name collision -- There are multiple input files for each of the following file names: test_sample.gatk.PASS.vcf.gz, test_sample.gatk.PASS.vcf.gz.tbi
    // probar cuando esté dragen con su vcf. 

    SPLITMULTIALLELIC.out.biallelic_renamed_vcf
        .view { meta, vcf, tbi -> 
            "Output from SPLITMULTIALLELIC: meta=${meta}, vcf=${vcf}, tbi=${tbi}"
        }
        .branch {
            gk: it[0].program == 'GK'
            dv: it[0].program == 'DV'
        }
        .set { split_output }

    split_output.gk
        .view { "GK output: $it" }
    split_output.dv
        .view { "DV output: $it" }

    ch_versions = ch_versions.mix(SPLITMULTIALLELIC.out.versions.first())

    vcf = SPLITMULTIALLELIC.out.biallelic_renamed_vcf

    emit:
    vcf  // channel: [ [val(meta)], path(vcf), path(tbi)]
    versions = ch_versions                     // channel: [ versions.yml ]

}
