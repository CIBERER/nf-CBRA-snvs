include { EXOMEDEPTH } from '../../../modules/local/exomedepth/main'
include { PANELCNMOPS } from '../../../modules/local/panelcmops/main'
include { CONVADING } from '../../../modules/local/convading/main'
include { CNVS_BED_FILTER } from '../../../modules/local/cnvs_bed_filter/main'
include { CNVS_RESULT_MIXER } from '../../../modules/local/cnvs_result_mixer/main'

include { ANNOTSV_ANNOTSV } from '../../../modules/nf-core/annotsv/annotsv/main'

include { POSTANNOTSV } from '../../../modules/local/postannotsv/main'

workflow CNVS_CALLING {

    take:
    bam_file_list        // channel (mandatory): [ path(bam)]
    bai_file_list        // channel (mandatory): [ path(bai)]
    ch_intervals_cnvs    // channel (mandatory) : [ path(bed) ]
    ch_fai          // channel (mandatory) : [ val(meta2), path(fai) ]
    runname // val(runname) from params.runname
    samples2analyce // path(samples2analyce) from params.samples2analyce
    annotations // channel (mandatory) : [ val(meta), path(annotations) ] from params.annotatio
    ch_small_variants // channel (optional) : [ path(small_variants) ] or [ [] ]
    ch_gene_transcripts // channel (optional) : [ val(meta), path(gene_transcripts) ] from params.gene_transcripts
    ch_candidate_genes // channel (optional) : [ val(meta), path(candidate_genes) ] from params.candidate_genes
    ch_false_positive_snv // channel (optional) : [ val(meta), path(false_positive_snv) ] from params.false_positive_snv
    ch_glowgenes_panel

    main:

    ch_versions = Channel.empty()
    
    CNVS_BED_FILTER (
        ch_intervals_cnvs,
        ch_fai.map{ meta, fai -> fai },
        params.min_target,
        params.chromosomes
    )

    // ─── Initialize empty channels for each CNV caller ───
    ch_convading_cnvs   = channel.empty()
    ch_panelcnmops_cnvs = channel.empty()
    ch_exomedepth_cnvs  = channel.empty()

    //CONVADING
    if (params.convading) {

        CONVADING(
            bam_file_list,
            bai_file_list,
            CNVS_BED_FILTER.out.cnvs_bed_filtered,
            ch_fai.map{ meta, fai -> fai },
            runname
        )

        ch_convading_cnvs = CONVADING.out.cnvs

    }

    //PANELCNMOPS

    if (params.panelcmops) {

        PANELCNMOPS(
            bam_file_list,
            bai_file_list,
            CNVS_BED_FILTER.out.cnvs_bed_filtered,
            runname
        )
        ch_panelcnmops_cnvs = PANELCNMOPS.out.cnvs
    }


    //EXOMEDEPTH

    if (params.exomedepth) {

        EXOMEDEPTH(
            bam_file_list,
            bai_file_list,
            CNVS_BED_FILTER.out.cnvs_bed_filtered,
            runname
        )
        EXOMEDEPTH.out.cnvs
        ch_exomedepth_cnvs = EXOMEDEPTH.out.cnvs

    }

   ch_cnvs = ch_exomedepth_cnvs
        .mix(ch_panelcnmops_cnvs, ch_convading_cnvs)
        .groupTuple().map { meta, file_lists -> [meta, file_lists.flatten()] }

    // Mix the CNV results from different callers and prepare the input for AnnotSV

    CNVS_RESULT_MIXER (
        ch_cnvs,
        samples2analyce
    )

    ch_for_annotsv = CNVS_RESULT_MIXER.out.merged_bed
        .combine(ch_small_variants)
        .map { meta, sv_vcf, candidate ->
            def sv_vcf_idx = []
            def cand = candidate.name == 'NO_FILE' ? [] : candidate
            [[id: meta], sv_vcf, sv_vcf_idx, cand]
        }

    
    ANNOTSV_ANNOTSV (
        ch_for_annotsv,
        annotations,
        ch_candidate_genes,
        ch_false_positive_snv, 
        ch_gene_transcripts
    )

    POSTANNOTSV (
        ANNOTSV_ANNOTSV.out.tsv.join(CNVS_RESULT_MIXER.out.colnames.map {meta, colnames -> [[id:meta], colnames] }),
        ch_candidate_genes.map{ meta, file -> file },
        ch_glowgenes_panel
    )

    emit:
    cnvs_annotated = POSTANNOTSV.out.annotated_cnv // channel: [ val(meta), path(tsv)]
    versions = ch_versions          // channel: [ versions.yml ]

}
