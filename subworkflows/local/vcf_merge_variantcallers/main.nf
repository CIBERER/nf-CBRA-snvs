include { BCFTOOLS_MERGE } from '../../../modules/nf-core/bcftools/merge/main'
include { TABIX_TABIX } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TABIX_FINAL_VCF } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_QUERY_STATS } from '../../../modules/local/bcftools_query_stats/main'
include { CONSENSUS_GENOTYPE } from '../../../modules/local/consensus_genotype/main'
include { GET_VCF_CALLERS_INFO } from '../../../modules/local/get_vcf_callers_info/main'
include { CREATE_SAMPLE_INFO } from '../../../modules/local/create_sample_info/main'



workflow VCF_MERGE_VARIANTCALLERS {

    take:
    ch_vcfs        // channel (mandatory): tuple val(meta), path(vcfs), path(tbis)
    ch_fasta        // channel (mandatory) : [ val(meta3), path(fasta) ]
    ch_fai          // channel (mandatory) : [ val(meta3), path(fai) ]
    ch_intervals    // channel (mandatory) : [ val(meta), path(bed) ]
    ch_assembly

    main:

    ch_versions = Channel.empty()

   ch_vcfs.map { item ->
        def meta = item[0]
        def files = item[1..-1].collect { file(it) }
        def vcfs = files.findAll { it.name.endsWith('.vcf.gz') }
        def tbis = files.findAll { it.name.endsWith('.vcf.gz.tbi') }
        [meta, vcfs, tbis]
    }
    .set { ch_for_bcftoolsmerge }

    ch_for_bcftoolsmerge.view()

    BCFTOOLS_MERGE (
        ch_for_bcftoolsmerge,
        ch_fasta,
        ch_fai,
        ch_intervals
    )

    TABIX_TABIX(
        BCFTOOLS_MERGE.out.vcf
    )

    BCFTOOLS_QUERY_STATS (
        BCFTOOLS_MERGE.out.vcf.join(TABIX_TABIX.out.tbi)
    )

    CONSENSUS_GENOTYPE (
        BCFTOOLS_QUERY_STATS.out.gt
    )

    GET_VCF_CALLERS_INFO (
        BCFTOOLS_QUERY_STATS.out.gt
    )

    ch_stats = CONSENSUS_GENOTYPE.out.consensus_gt.join(CONSENSUS_GENOTYPE.out.discordances).join(BCFTOOLS_QUERY_STATS.out.ad_mean).join(BCFTOOLS_QUERY_STATS.out.dp_mean).join(BCFTOOLS_QUERY_STATS.out.rd_mean).join(BCFTOOLS_QUERY_STATS.out.vd_mean).join(BCFTOOLS_QUERY_STATS.out.vaf).join(GET_VCF_CALLERS_INFO.out.sf_file).join(BCFTOOLS_QUERY_STATS.out.programs)

    CREATE_SAMPLE_INFO (
        BCFTOOLS_MERGE.out.vcf.join(TABIX_TABIX.out.tbi).join(ch_stats),
        ch_assembly
    )

    TABIX_TABIX_FINAL_VCF (
        CREATE_SAMPLE_INFO.out.final_vcf
    )

    vcf = CREATE_SAMPLE_INFO.out.final_vcf.join(TABIX_TABIX_FINAL_VCF.out.tbi).view()

    emit:
    vcf                                        // channel: [ [val(meta)], path(vcf), path(tbi)]
    versions = ch_versions                     // channel: [ versions.yml ]

}
