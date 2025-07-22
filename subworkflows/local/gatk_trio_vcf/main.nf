include { GATK4_HAPLOTYPECALLER                                          }      from '../../../modules/nf-core/gatk4/haplotypecaller/main'

include { BCFTOOLS_SORT                                                  }      from '../../../modules/nf-core/bcftools/sort/main'
include { GATK4_MERGEVCFS                                                }      from '../../../modules/nf-core/gatk4/mergevcfs/main'

include { TABIX_TABIX                                                    }      from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_FILTER                                                }      from '../../../modules/nf-core/bcftools/filter/main'
include { SPLITMULTIALLELIC                                                }      from '../../../modules/local/splitmultiallelic/main'

include { GATK4_SELECTVARIANTS  as  GATK4_SELECTVARIANTS_SNP             }      from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS  as  GATK4_SELECTVARIANTS_INDEL           }      from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS  as  GATK4_SELECTVARIANTS_MIX             }      from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_VARIANTFILTRATION  as  GATK4_VARIANTFILTRATION_SNV       }      from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_VARIANTFILTRATION  as  GATK4_VARIANTFILTRATION_INDEL     }      from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_VARIANTFILTRATION  as  GATK4_VARIANTFILTRATION_MIX       }      from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_VARIANTFILTRATION  as  GATK4_VARIANTFILTRATION_GENOTYPEPOSTERIOR       }      from '../../../modules/nf-core/gatk4/variantfiltration/main'

include { GATK4_GENOMICSDBIMPORT } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { GATK4_CALCULATEGENOTYPEPOSTERIORS } from '../../../modules/local/gatk4/calculategenotypeposteriors/main'

include { GATK4_VARIANTRECALIBRATOR as VARIANTRECALIBRATOR_INDEL} from '../../../modules/nf-core/gatk4/variantrecalibrator/main'
include { GATK4_VARIANTRECALIBRATOR as VARIANTRECALIBRATOR_SNP} from '../../../modules/nf-core/gatk4/variantrecalibrator/main'

// TODO: Create a module for CalculateGenotypePosteriors
// TODO: Create a module for ANNOTATEVARIANTS


workflow GATK_TRIO_VCF {

    take:
    ch_bam        // channel (mandatory): [ val(meta), path(bam), path(bai) ]
    ch_intervals    // channel (mandatory) : [ val(meta), path(bed) ]
    ch_fasta        // channel (mandatory) : [ val(meta2), path(fasta) ]
    ch_fai          // channel (mandatory) : [ val(meta2), path(fai) ]
    ch_refdict      // channel (mandatory) : [ val(meta2), path(dict) ]
    ch_dbsnp  // channel (mandatory) : [ val(meta3), path(vcf) ]
    ch_dbsnp_tbi  // channel (mandatory) : [ val(meta3), path(vcf) ]
    ch_intervals_genomicsdbimport
    no_intervals // channel (mandatory) : [ boolean ]
    ch_ped

    main:

    ch_versions = Channel.empty()
    
    ch_dragstr_model = ch_bam.map {meta, bam, bai -> tuple(meta, [])}


    GATK4_HAPLOTYPECALLER (
        ch_bam.join(ch_intervals).join(ch_dragstr_model),
        ch_fasta,
        ch_fai,
        ch_refdict,
        ch_dbsnp,
        ch_dbsnp_tbi
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())

    // Collect all VCFs and TBIs from all samples
    // First, prepare your joint VCF data
    joint_vcf_data = GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi, by: 0)
        .toList()
        .map { list_of_tuples ->
            def vcfs = []
            def tbis = []
            def sample_metas = []

            list_of_tuples.each { meta, vcf, tbi ->
                vcfs.add(vcf)
                tbis.add(tbi)
                sample_metas.add(meta)
            }
            
            def joint_meta = [
                id: 'joint_variant_calling',
                samples: sample_metas.collect { it.id },
                sample_count: sample_metas.size(),
            ]
            tuple(joint_meta, vcfs, tbis)
        }//.view()

    // Then combine with intervals
    joint_gvcf_ch = joint_vcf_data
        .combine(ch_intervals_genomicsdbimport)
        .map { joint_meta, vcfs, tbis, intervals ->
            tuple(joint_meta, vcfs, tbis, intervals, [], [])
        }
        //.view()

    // Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport
    GATK4_GENOMICSDBIMPORT(joint_gvcf_ch, false, false, false)
    //GATK4_GENOMICSDBIMPORT.out.genomicsdb.view()

    if (no_intervals) {
        // If no intervals are provided, we can use the whole genome
        genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{ meta, genomicsdb -> [ meta, genomicsdb, [], [], [] ] }
    } else {
        // Use the provided intervals
        genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.combine(ch_intervals.map { meta, bed -> bed }.first()).map{ meta, genomicsdb, intervals -> [ meta, genomicsdb, [], intervals, [] ] }
    }

    genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{ meta, genomicsdb -> [ meta, genomicsdb, [], [], [] ] }

    GATK4_GENOTYPEGVCFS (
        genotype_input,
        ch_fasta,
        ch_fai,
        ch_refdict,
        ch_dbsnp,
        ch_dbsnp_tbi
    )

    //GATK4_GENOTYPEGVCFS.out.vcf.view()
    //ch_intervals_genomicsdbimport.view()

    //GATK4_GENOTYPEGVCFS.out.vcf.join(GATK4_GENOTYPEGVCFS.out.tbi).combine(ch_intervals_genomicsdbimport)

    if (no_intervals) {
        // If no intervals are provided, we can use the whole genome
        ch_bed_for_selectvariants = GATK4_GENOTYPEGVCFS.out.vcf
        .join(GATK4_GENOTYPEGVCFS.out.tbi).map { meta, vcf, tbi -> [meta, vcf, tbi, []]
        }
    } else {
        // Use the provided intervals
        ch_bed_for_selectvariants = GATK4_GENOTYPEGVCFS.out.vcf
        .join(GATK4_GENOTYPEGVCFS.out.tbi)
        .combine(ch_intervals.map { meta, bed -> bed }.first())
    }

    ch_bed_for_selectvariants.view()


    // TODO: create subworkflows for HardFiltering
    
    GATK4_SELECTVARIANTS_SNP (
        ch_bed_for_selectvariants
    )
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS_SNP.out.versions.first())

    GATK4_SELECTVARIANTS_INDEL(
        ch_bed_for_selectvariants
    )
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS_INDEL.out.versions.first())

    GATK4_SELECTVARIANTS_MIX(
        ch_bed_for_selectvariants
    )
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS_MIX.out.versions.first())

    GATK4_VARIANTFILTRATION_SNV(
        GATK4_SELECTVARIANTS_SNP.out.vcf.join(GATK4_SELECTVARIANTS_SNP.out.tbi),
        ch_fasta,
        ch_fai,
        ch_refdict,
    ) 
    ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION_SNV.out.versions.first())

    GATK4_VARIANTFILTRATION_INDEL(
        GATK4_SELECTVARIANTS_INDEL.out.vcf.join(GATK4_SELECTVARIANTS_INDEL.out.tbi),
        ch_fasta,
        ch_fai,
        ch_refdict,
    )
    ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION_INDEL.out.versions.first())

    GATK4_VARIANTFILTRATION_MIX(
        GATK4_SELECTVARIANTS_MIX.out.vcf.join(GATK4_SELECTVARIANTS_MIX.out.tbi),
        ch_fasta,
        ch_fai,
        ch_refdict,
    )
    ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION_MIX.out.versions.first())

    BCFTOOLS_SORT(
        GATK4_VARIANTFILTRATION_SNV.out.vcf.concat( GATK4_VARIANTFILTRATION_INDEL.out.vcf, GATK4_VARIANTFILTRATION_MIX.out.vcf )
    )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())

    GATK4_MERGEVCFS(
        BCFTOOLS_SORT.out.vcf.groupTuple(),
        ch_refdict
    )
    ch_versions = ch_versions.mix(GATK4_MERGEVCFS.out.versions.first())

    BCFTOOLS_FILTER(
        GATK4_MERGEVCFS.out.vcf.join(GATK4_MERGEVCFS.out.tbi)
    )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())
    
    // SPLITMULTIALLELIC (
    //     BCFTOOLS_FILTER.out.vcf.join(BCFTOOLS_FILTER.out.tbi),
    //     ch_fasta,
    //     ""
    // )

    // TODO: add a module for CalculateGenotypePosteriors
    // TODO: add a module for ANNOTATEVARIANTS


    ch_ped_genotypeposterior = ch_ped
        .toList()
        .map { list_of_tuples ->
            def peds = []
            def sample_metas = []

            list_of_tuples.each { meta, ped ->
                peds.add(ped)
                sample_metas.add(meta)
            }
            
            def joint_meta = [
                id: 'joint_variant_calling',
                samples: sample_metas.collect { it.id },
                sample_count: sample_metas.size(),
            ]
            tuple(joint_meta, peds.first())
        }.view() 


    if (no_intervals) {
        // If no intervals are provided, we can use the whole genome
        ch_input_genotypeposterior = BCFTOOLS_FILTER.out.vcf.join(BCFTOOLS_FILTER.out.tbi).join(ch_ped_genotypeposterior).map { meta, vcf, ped -> [meta, vcf, [], ped]}
    } else {
        // Use the provided intervals
        ch_input_genotypeposterior = BCFTOOLS_FILTER.out.vcf.join(BCFTOOLS_FILTER.out.tbi).join(ch_ped_genotypeposterior).combine(ch_intervals.map { meta, bed -> bed }.first())
        map { meta, vcf, ped, interval -> [meta, vcf, interval, ped]}
    }

    //ch_input_genotypeposterior.view()

    GATK4_CALCULATEGENOTYPEPOSTERIORS (
        ch_input_genotypeposterior
    )

 // TODO: add variantfiltration. 
 
    vcf = BCFTOOLS_FILTER.out.vcf.join(BCFTOOLS_FILTER.out.tbi).view()

    emit:
    vcf // channel: [ val(meta), path(vcf), path(tbi)]
    versions = ch_versions                     // channel: [ versions.yml ]

}
