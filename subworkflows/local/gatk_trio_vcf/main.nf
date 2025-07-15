include { GATK4_HAPLOTYPECALLER                                          }      from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_SELECTVARIANTS  }      from '../../../modules/nf-core/gatk4/selectvariants/main'

include { GATK4_VARIANTFILTRATION  }      from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { BCFTOOLS_SORT                                                  }      from '../../../modules/nf-core/bcftools/sort/main'
include { GATK4_MERGEVCFS                                                }      from '../../../modules/nf-core/gatk4/mergevcfs/main'

include { TABIX_TABIX                                                    }      from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_FILTER                                                }      from '../../../modules/nf-core/bcftools/filter/main'
include { SPLITMULTIALLELIC                                                }      from '../../../modules/local/splitmultiallelic/main'


include { GATK4_GENOMICSDBIMPORT } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { GATK4_VARIANTRECALIBRATOR } from '../../../modules/nf-core/gatk4/variantrecalibrator/main'
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
    joint_gvcf_ch = GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi, by: 0)  // Join by meta (first element)
        .toList().view()
        .map { list_of_tuples ->
            def vcfs = []
            def tbis = []

            list_of_tuples.each { meta, vcf, tbi ->
                vcfs.add(vcf)
                tbis.add(tbi)
            }
            
            def joint_meta = [id: 'joint_variant_calling']
            
            tuple(joint_meta, vcfs, tbis, ch_intervals_genomicsdbimport, [], [])
        }
        .view()

    // Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport
    GATK4_GENOMICSDBIMPORT(joint_gvcf_ch, false, false, false)
    GATK4_GENOMICSDBIMPORT.out.genomicsdb.view()

    genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{ meta, genomicsdb -> [ meta, genomicsdb, [], [], [] ] }

    GATK4_GENOTYPEGVCFS (
        genotype_input,
        ch_fasta,
        ch_fai,
        ch_refdict,
        ch_dbsnp,
        ch_dbsnp_tbi
    )

    GATK4_GENOTYPEGVCFS.out.vcf.view()
 
    gvcf = GATK4_HAPLOTYPECALLER.out.vcf

    emit:
    gvcf // channel: [ val(meta), path(vcf), path(tbi)]
    versions = ch_versions                     // channel: [ versions.yml ]

}
