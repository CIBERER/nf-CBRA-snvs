# GdTBioinfo-nf/snvs: Documentation about Variant Calling

## [GATK subworkflow](https://github.com/CIBERER/GdTBioinfo-nf-snvs/blob/gatk_subworkflow/subworkflows/local/gatk_vcf/main.nf)

This sub-workflow detects variants with [GATK Haplotypecaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) and selects PASS variants. 
By default, the pipeline considers the following criteria to mark a variant as PASS. The criteria are coded in the [modules.config file](https://github.com/CIBERER/GdTBioinfo-nf-snvs/blob/3ecd0886380c507d467ce2da78ce4a5868829c05/conf/modules.config#L69C1-L111C6) and can be modified with a user-defined nextflow config file.

**1. Filters for SNVs:** 
  * -filter "QD < 2.0" --filter-name "QD2" 
  * -filter "QUAL < 30.0" --filter-name "QUAL30" 
  * -filter "SOR > 3.0" --filter-name "SOR3" 
  * -filter "FS > 60.0" --filter-name "FS60" 
  * -filter "MQ < 40.0" --filter-name "MQ40" 
  * -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" 
  * -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

**2. Filters for INDELs:**
  * -filter "QD < 2.0" --filter-name "QD2" 
  * -filter "QUAL < 30.0" --filter-name "QUAL30" 
  * -filter "FS > 200.0" --filter-name "FS200" 
  * -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

**3. Filters for MIXED, MNP and SYMBOLIC:**
  * -filter "QD < 2.0" --filter-name "QD2" 
  * -filter "QUAL < 30.0" --filter-name "QUAL30"


More information about the [variant types detected by GATK tools](https://gatk.broadinstitute.org/hc/en-us/articles/360035530752-What-types-of-variants-can-GATK-tools-detect-or-handle). 

More information about [INFO fields in the vcf file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format).

