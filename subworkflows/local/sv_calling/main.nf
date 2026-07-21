include { MANTA_GERMLINE          } from '../../../modules/nf-core/manta/germline/main'
include { ANNOTSV_ANNOTSV         } from '../../../modules/nf-core/annotsv/annotsv/main'

workflow SV_CALLING {

    take:
    ch_bam                // channel (mandatory): [ val(meta), path(bam), path(bai) ]
    ch_fasta              // channel (mandatory): [ val(meta2), path(fasta) ]
    ch_fai                // channel (mandatory): [ val(meta3), path(fai) ]
    ch_manta_config       // channel (optional):  [ path(config) ]
    ch_annotsv_annotations // channel (mandatory): [ val(meta), path(annotations) ]
    ch_candidate_genes    // channel (optional):  [ val(meta), path(candidate_genes) ]
    ch_false_positive_snv // channel (optional):  [ val(meta), path(false_positive_snv) ]
    ch_gene_transcripts   // channel (optional):  [ val(meta), path(gene_transcripts) ]

    main:

    ch_versions = Channel.empty()

    //
    // Run Manta germline SV calling
    //
    // Manta expects: [ meta, bam, bai, target_bed, target_bed_tbi ]
    // For WGS there is no target BED, so pass empty files

    if (params.manta_joint) {
        log.info "Running Manta in joint calling mode"
            // Define the run name
        if (params.runname) { runname = params.runname }
        else { runname = new Date().format("yyyy-MM-dd_HH-mm") }
        println "Run name: $runname"

        bam_file_list = ch_bam
        .map{ meta, bam, bai -> bam }
        .collect()
        .map { files -> 
        def meta = runname
        [[id:meta], files.sort { it.name }] }
                
        
        bai_file_list = ch_bam
        .map{ meta, bam, bai -> bai }
        .collect().map { files -> 
        def meta = runname
        [[id:meta], files.sort { it.name }] }


        ch_manta_input = bam_file_list.join(bai_file_list)
        .map { meta, bam, bai ->
        [ meta, bam, bai, [], [] ]
        }

    } else {

    ch_manta_input = ch_bam.map { meta, bam, bai ->
        [ meta, bam, bai, [], [] ]
        }
    }

    MANTA_GERMLINE (
        ch_manta_input,
        ch_fasta,
        ch_fai,
        ch_manta_config
    )
    ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions.first())

    //
    // Prepare AnnotSV input from Manta diploid SV output
    //
    // MANTA_GERMLINE emits: diploid_sv_vcf [ meta, vcf.gz ] and diploid_sv_vcf_tbi [ meta, tbi ]
    ch_annotsv_input = MANTA_GERMLINE.out.diploid_sv_vcf
        .join(MANTA_GERMLINE.out.diploid_sv_vcf_tbi)
        .map { meta, vcf, tbi ->
            [ meta, vcf, tbi, [] ]  // no candidate_small_variants
        }

    //
    // Annotate structural variants with AnnotSV
    //
    ANNOTSV_ANNOTSV (
        ch_annotsv_input,
        ch_annotsv_annotations,
        ch_candidate_genes,
        ch_false_positive_snv,
        ch_gene_transcripts
    )
    ch_versions = ch_versions.mix(ANNOTSV_ANNOTSV.out.versions.first())

    emit:
    annotated_tsv     = ANNOTSV_ANNOTSV.out.tsv             // channel: [ val(meta), path(tsv) ]
    annotated_vcf     = ANNOTSV_ANNOTSV.out.vcf             // channel: [ val(meta), path(vcf) ] (optional, requires -vcf 1 in args)
    unannotated_tsv   = ANNOTSV_ANNOTSV.out.unannotated_tsv // channel: [ val(meta), path(tsv) ]
    diploid_sv_vcf    = MANTA_GERMLINE.out.diploid_sv_vcf   // channel: [ val(meta), path(vcf.gz) ]
    candidate_sv_vcf  = MANTA_GERMLINE.out.candidate_sv_vcf // channel: [ val(meta), path(vcf.gz) ]
    versions          = ch_versions                         // channel: [ versions.yml ]
}

