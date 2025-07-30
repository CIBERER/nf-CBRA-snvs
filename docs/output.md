# nf-CBRA-snvs: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FastQC](#fastqc) - Raw read QC
- [Mapping](#Mapping) - Map reads to reference (BAW-MEM) and process bam file (`GATK MarkDuplicates`, `GATK BaseRecalibrator` and `GATK ApplyBQSR`)
- [Variant Calling](#Variant-Calling) - Detect variants with 3 tools:
  - [GATK4 Haplotypecaller](#GATK4-Haplotypecaller)
  - [Dragen](#Dragen)
  - [DeepVariant](#DeepVariant)
- [Merge and Integration](#Merge-and-Integration) - Merge and integrate the variants from the vcfs obtained with the different tools
- [Annotation](#Annotation) - Annotate the variants with [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) and add regions of homozygosity (ROHs) with [AUTOMAP](https://github.com/mquinodo/AutoMap) and other custom information. 
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

:::note
The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.
:::

### Mapping

<details markdown="1">
<summary>Output files</summary>

- `alignment/`
  - `*.bam`: bam file.
  - `*.bam.bai`: bam index file.

</details>

[BWA-MEM](https://github.com/lh3/bwa) map reads to reference genome. Resulting bam file is processed following [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery) with `GATK MarkDuplicates`, `GATK BaseRecalibrator` and `GATK ApplyBQSR`

### Variant Calling

<details markdown="1">
<summary>Output files</summary>

- `individual_callers_snvs/`
  - `*.gatk.PASS.vcf.gz`: vcf file obtained with gatk4 haplotypecaller.
  - `*.gatk.PASS.vcf.gz.tbi`: vcf index file obtained with gatk4 haplotypecaller.
  - `*.dragen.PASS.vcf.gz`: vcf file obtained with dragen.
  - `*.dragen.PASS.vcf.gz.tbi`: vcf index file obtained with dragen.
  - `*.deepvariant.PASS.vcf.gz`: vcf file obtained with deepvariant.
  - `*.deepvariant.PASS.vcf.gz.tbi`: vcf index file obtained with deepvariant.

</details>

See [variant calling](docs/variant_calling.md) section for more information. 

### Merge and Integration

<details markdown="1">
<summary>Output files</summary>

- `snvs/`
  - `*.final.vcf.gz`: final vcf file after merging and integrating all individual vcf files. 
  - `*.final.vcf.gz.tbi`: vcf index file. 

</details>

Vcf files are merged with [bcftools merge](https://samtools.github.io/bcftools/bcftools.html#merge). GT, DP and AD of detected variants with each variant caller are obtained with [BCFTOOLS_QUERY_STATS](modules/local/bcftools_query_stats/main.nf) module, calculating mean DP, mean AD and variant allele depth (VAF), and information about the variant callers that detect each variant is extracted in [GET_VCF_CALLERS_INFO](modules/local/get_vcf_callers_info/main.nf) module. A consensus genotype for each variant is established based on the most numerical genotype with the module [CONSENSUS_GENOTYPE](modules/local/consensus_genotype/main.nf). All this information is integrated in the step [CREATE_SAMPLE_INFO](modules/local/create_sample_info/main.nf), obtaing the final vcf file. 

### Annotation

<details markdown="1">
<summary>Output files</summary>

- `snvs/`
  - `*.SNV.INDEL.annotated.tsv`: final tsv file with annotated variants. 

</details>

The variants are annotated with [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) using the flag `--everything`, which includes the following options: `--sift b, --polyphen b, --ccds, --hgvs, --symbol, --numbers, --domains, --regulatory, --canonical, --protein, --biotype, --af, --af_1kg, --af_esp, --af_gnomade, --af_gnomadg, --max_af, --pubmed, --uniprot, --mane, --tsl, --appris, --variant_class, --gene_phenotype, --mirna`. See [this page](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html) for more information. `--custom` flag is used to include INFO field of the vcf file in the final annotated tsv file. 

[POSTVEP](modules/local/postvep/main.nf) step takes the VEP tab demilited output, filter variants by minor allele frequency (`--maf`) and add other custom annotations, as regions of homozigosity (ROHs) detected with [AUTOMAP](https://github.com/mquinodo/AutoMap) and [GLOWgenes](https://www.translationalbioinformaticslab.es/tblab-home-page/tools/glowgenes), a network-based algorithm developed to prioritize novel candidate genes associated with rare diseases.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
