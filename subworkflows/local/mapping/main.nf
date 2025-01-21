include { BWA_MEM                 } from '../../../modules/nf-core/bwa/mem/main'
include { GATK4_FASTQTOSAM } from '../../../modules/nf-core/gatk4/fastqtosam/main'
include { GATK4_MERGEBAMALIGNMENT } from '../../../modules/nf-core/gatk4/mergebamalignment/main'
include { GATK4SPARK_MARKDUPLICATES } from '../../../modules/nf-core/gatk4spark/markduplicates/main'
include { PICARD_SORTSAM } from '../../../modules/nf-core/picard/sortsam/main'
include { GATK4_BASERECALIBRATOR } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { SAMTOOLS_SORT   } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX   } from '../../../modules/nf-core/samtools/index/main'


workflow MAPPING {

    take:
    ch_reads        // channel (mandatory): [ val(meta), [ path(reads) ] ]
    ch_index        // channel (mandatory): [ val(meta2), path(index) ]
    ch_fasta        // channel (mandatory) : [ val(meta3), path(fasta) ]
    ch_fai          // channel (mandatory) : [ val(meta3), path(fai) ]
    ch_refdict      // channel (mandatory) : [ val(meta3), path(dict) ]
    ch_intervals    // channel (mandatory) : [ val(meta), path(bed) ]
    ch_known_sites  // channel (mandatory) : [ val(meta4), path(vcf) ]
    ch_known_sites_tbi  // channel (mandatory) : [ val(meta4), path(vcf) ]

    main:

    ch_versions = Channel.empty()
    
    BWA_MEM (
        ch_reads,
        ch_index,
        ch_fasta,
        false
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    GATK4_FASTQTOSAM (
        ch_reads
    )
    ch_versions = ch_versions.mix(GATK4_FASTQTOSAM.out.versions.first())

    GATK4_MERGEBAMALIGNMENT (
        BWA_MEM.out.bam.join(GATK4_FASTQTOSAM.out.bam),
        ch_fasta,
        ch_refdict
    )
    ch_versions = ch_versions.mix(GATK4_MERGEBAMALIGNMENT.out.versions.first())

    GATK4SPARK_MARKDUPLICATES (
        GATK4_MERGEBAMALIGNMENT.out.bam,
        ch_fasta.map {meta , fasta -> [fasta]},
        ch_fai.map {meta, fai -> [fai]}, 
        ch_refdict.map {meta, dict -> [dict]}
    )
    ch_versions = ch_versions.mix(GATK4SPARK_MARKDUPLICATES.out.versions.first())

    GATK4_BASERECALIBRATOR (   
         GATK4SPARK_MARKDUPLICATES.out.output.join(GATK4SPARK_MARKDUPLICATES.out.bam_index.join(ch_intervals)),
         ch_fasta.map {meta , fasta -> [fasta]},
         ch_fai.map {meta, fai -> [fai]}, 
         ch_refdict.map {meta, dict -> [dict]},
         ch_known_sites,
         ch_known_sites_tbi
    )
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions.first())

    GATK4_APPLYBQSR (
        GATK4SPARK_MARKDUPLICATES.out.output.join(GATK4SPARK_MARKDUPLICATES.out.bam_index).join(GATK4_BASERECALIBRATOR.out.table).join(ch_intervals),
        ch_fasta.map {meta , fasta -> [fasta]},
        ch_fai.map {meta, fai -> [fai]}, 
        ch_refdict.map {meta, dict -> [dict]}
    )
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions.first())

    SAMTOOLS_INDEX (
        GATK4_APPLYBQSR.out.bam 
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    bam = GATK4_APPLYBQSR.out.bam.join(SAMTOOLS_INDEX.out.bai)

    emit:
    bam  // channel: [ [val(meta)], path(bam), path(bai)]
    versions = ch_versions                     // channel: [ versions.yml ]

}
