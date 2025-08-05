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
 
    vcf = GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi).view()

    emit:
    vcf // channel: [ val(meta), path(vcf), path(tbi)]
    versions = ch_versions                     // channel: [ versions.yml ]

}