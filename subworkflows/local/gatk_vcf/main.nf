include { GATK4_HAPLOTYPECALLER                                          }      from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_SELECTVARIANTS  as  GATK4_SELECTVARIANTS_SNP             }      from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS  as  GATK4_SELECTVARIANTS_INDEL           }      from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS  as  GATK4_SELECTVARIANTS_MIX             }      from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_VARIANTFILTRATION  as  GATK4_VARIANTFILTRATION_SNV       }      from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_VARIANTFILTRATION  as  GATK4_VARIANTFILTRATION_INDEL     }      from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_VARIANTFILTRATION  as  GATK4_VARIANTFILTRATION_MIX       }      from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { BCFTOOLS_SORT                                                  }      from '../../../modules/nf-core/bcftools/sort/main'
include { GATK4_MERGEVCFS                                                }      from '../../../modules/nf-core/gatk4/mergevcfs/main'

include { TABIX_TABIX                                                    }      from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_FILTER                                                }      from '../../../modules/nf-core/bcftools/filter/main'


workflow GATK_VCF {

    take:
    ch_bam        // channel (mandatory): [ val(meta), path(bam), path(bai) ]
    ch_intervals    // channel (mandatory) : [ val(meta), path(bed) ]
    ch_fasta        // channel (mandatory) : [ val(meta2), path(fasta) ]
    ch_fai          // channel (mandatory) : [ val(meta2), path(fai) ]
    ch_refdict      // channel (mandatory) : [ val(meta2), path(dict) ]
    ch_dbsnp  // channel (mandatory) : [ val(meta3), path(vcf) ]
    ch_dbsnp_tbi  // channel (mandatory) : [ val(meta3), path(vcf) ]

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



    GATK4_SELECTVARIANTS_SNP (
        GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi).join(ch_intervals)
    )
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS_SNP.out.versions.first())

    GATK4_SELECTVARIANTS_INDEL(
        GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi).join(ch_intervals)
    )
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS_INDEL.out.versions.first())

    GATK4_SELECTVARIANTS_MIX(
        GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi).join(ch_intervals)
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
 
    vcf = BCFTOOLS_FILTER.out.vcf.join(BCFTOOLS_FILTER.out.tbi)

    emit:
    vcf // channel: [ val(meta), path(bam)]
    versions = ch_versions                     // channel: [ versions.yml ]

}