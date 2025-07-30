# nf-CBRA-snvs: Documentation about Variant Calling

## [GATK subworkflow](https://github.com/CIBERER/GdTBioinfo-nf-snvs/blob/gatk_subworkflow/subworkflows/local/gatk_vcf/main.nf)

This sub-workflow detects variants with [GATK Haplotypecaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller), apply [hard-filter](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering) based on GATK4 recomendations and selects PASS variants. 
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

After filtering PASS variants, [bcftools norm](https://samtools.github.io/bcftools/bcftools.html#norm) is used to split multiallelic variants. 

## [Dragen subworkflow](subworkflows/local/dragen_vcf/main.nf)

This sub-workflow detects variants with GATK Haplotypecaller with `--dragen-mode`. The parameters for the DRAGstr model for the input sample are estimated with [GATK4 CalibrateDragstrModel](https://gatk.broadinstitute.org/hc/en-us/articles/21905129631259-CalibrateDragstrModel). Hard Filtering is based on QUAL score in the [modules.config file](https://github.com/CIBERER/GdTBioinfo-nf-snvs/blob/3ecd0886380c507d467ce2da78ce4a5868829c05/conf/modules.config#L69C1-L111C6) 

´´´
gatk VariantFiltration \
      -V output_file.vcf \
      --filter-expression "QUAL < 10.4139" \
      --filter-name "DRAGENHardQUAL" \
      -O output_filtered.vcf
´´´

More information about running germline single sample short variant discovery in DRAGEN mode [here](https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939--How-to-Run-germline-single-sample-short-variant-discovery-in-DRAGEN-mode).

After filtering PASS variants, [bcftools norm](https://samtools.github.io/bcftools/bcftools.html#norm) is used to split multiallelic variants. 


## [DeepVariant subworkflow](subworkflows/local/deep_variant_vcf/main.nf)

This sub-workflow detects variants with [DeepVariant](https://github.com/google/deepvariant). It conteins three steps: 

* makeexamples: Converts the input alignment file to a tfrecord format suitable for the deep learning model
* callvariants: Call variants based on input tfrecords. The output is also in tfrecord format, and needs postprocessing to convert it to vcf.
* postprocessvariants: Convert variant calls from callvariants to VCF. 

See [DeepVariant subworkflow](https://github.com/nf-core/modules/tree/master/modules/nf-core/deepvariant) if you cant more information. 

After filtering PASS variants, [bcftools norm](https://samtools.github.io/bcftools/bcftools.html#norm) is used to split multiallelic variants. 



