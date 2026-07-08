# Introduction

**CBRA** (CIBERER Bioinformatics for Rare diseases Analysis - Small Nucleotide Variant) is a workflow optimized for the analysis of rare diseases, designed to detect SNVs and INDELs in targeted sequencing data (CES/WES) as well as whole genome sequencing (WGS).

This pipeline is developed using Nextflow, a workflow management system that enables an easy execution across various computing environments. It uses Docker or Singularity containers, simplifying setup and ensuring reproducibility of results. The pipeline assigns a container to each process, which simplifies the management and updating of software dependencies. When possible, processes are sourced from nf-core/modules, promoting reusability across all nf-core pipelines and contributing to the broader Nextflow community.

<p align="center">
  <img title="CBRA Workflow" src="docs/images/workflow.png" width=80%>
</p>

# Pipeline summary

The pipeline can perform the following steps:

- **Mapping** (`mapping = true`) of the reads to reference (BWA-MEM) 
- Process BAM file (`GATK MarkDuplicates`, `GATK BaseRecalibrator` and `GATK ApplyBQSR`)
- **Single Nucleotide Variant (SNVs) calling** (`variant_calling = true`) with the following tools:

  - GATK4 Haplotypecaller (`run_gatk = true`). This subworkflow includes:
    - **GATK4 Haplotypecaller**.
    - **Hard Filters** and **VarianFiltration** to mark PASS variants. More information [here](docs/variant_calling.md).
    - **Bcftools Filter** to keep PASS variants on chr1-22, X, Y.
    - **Split Multialletic**.
  - Dragen (`run_dragen = true`). This subworkflow includes:
    - **GATK4 Calibratedragstrmodel**
    - **GATK4 Haplotypecaller** with `--dragen-mode`.
    - **VarianFiltration** with `--filter-expression "QUAL < 10.4139" --filter-name "DRAGENHardQUAL"`to mark PASS variants. More information [here](https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939--How-to-Run-germline-single-sample-short-variant-discovery-in-DRAGEN-mode).
    - **Bcftools Filter** to keep PASS variants on chr1-22, X, Y.
    - **Split Multialletic**.
  - DeepVariant (`run_deepvariant = true`). This subworkflow includes:
    - **DeepVariant makeexamples**: Converts the input alignment file to a tfrecord format suitable for the deep learning model.
    - **DeepVariant callvariants**: Call variants based on input tfrecords. The output is also in tfrecord format, and needs postprocessing to convert it to vcf.
    - **DeepVariant postprocessvariants**: Convert variant calls from callvariants to VCF, and also create GVCF files based on genomic information from makeexamples. More information [here](https://github.com/nf-core/modules/tree/master/modules/nf-core/deepvariant).
    - **Bcftools Filter** to keep PASS variants on chr1-22, X, Y.
    - **Split Multialletic**.

  - In addition to these three variant callers, a trio analysis can be performed with option `trio_analysis = true`. Based on [GATK guides](https://gatk.broadinstitute.org/hc/en-us/articles/360035531432-Genotype-Refinement-workflow-for-germline-short-variants). This subworkflow includes: 
    - **GATK4 Haplotypecaller**.
    - **GATK4 Genomicsdbimport** to merge GVCFs from multiple samples
    - **GATK4 Genotypegvcfs** to perform joint genotyping
    - **Hard Filters** and **VarianFiltration** to mark PASS variants. More information [here](docs/variant_calling.md).
    - **Bcftools Filter** to keep PASS variants on chr1-22, X, Y.
    - **GATK4 Calculategenotypeposteriors** to calculate genotype posterior probabilities given the family 
    - **GATK4 Variantfiltration** based on genotypeposterior GQ<20
    - **GATK4 Variantannotator** to annotate possible de novo mutations in trios
    - **Filter proband ref**: filter variants that are REF in the proband
    - **Split Multialletic**.


- **Merge and integration** of the vcfs obtained with the different tools.
- **Annotation of SNVs** (`annotation = true`) of the variants:
  - Regions of homozygosity (ROHs) with [AUTOMAP](https://github.com/mquinodo/AutoMap)
  - Effect of the variants with [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) using the flag `--everything`, which includes the following options: `--sift b, --polyphen b, --ccds, --hgvs, --symbol, --numbers, --domains, --regulatory, --canonical, --protein, --biotype, --af, --af_1kg, --af_esp, --af_gnomade, --af_gnomadg, --max_af, --pubmed, --uniprot, --mane, --tsl, --appris, --variant_class, --gene_phenotype, --mirna`
  - Filtering by maf: If variants are only annotated with VEP's `--max_af` option, filtering is based on the MAX_AF field. Variants are retained if MAX_AF is missing or lower than the user-defined MAF threshold. Therefore, only variants with MAX_AF ≥ maf are removed. If variants are annotated using the gnomAD exome and genome annotation tables (using `--custom` option in vep), filtering uses the gnomADe_AF_grpmax, gnomADg_AF_grpmax, gnomADe_filt, and gnomADg_filt fields. A variant is removed only when its maximum population allele frequency (AF_grpmax) is greater than or equal to the MAF threshold and the corresponding gnomAD record has a PASS filter status. Variants with missing allele frequencies or non-PASS filter status are retained. Exome and genome annotations are evaluated independently, and a variant is discarded if it meets the removal criteria in either dataset.
  - POSTVEP is a module that runs a script that takes the table generated by VEP – which has been filtered by maf – and formats certain fields to facilitate analysis. Furthermore, if a list of genes has been added (`--gene_list`), it will be used to filter the variants  or to run GLOWgenes (if this has been enabled by setting `--glowgenes true`). In addition, this script can be used to add further annotations relating to the genes. 
  - You can enhance the annotation by incorporating gene rankings from [GLOWgenes](https://www.translationalbioinformaticslab.es/tblab-home-page/tools/glowgenes), a network-based algorithm developed to prioritize novel candidate genes associated with rare diseases. To add GLOWgenes, set the parameter `--glowgenes true`. GLOWgenes will run using the list of genes specified in the `--gen_list` parameter (an example [here](https://github.com/CIBERER/CBRA/blob/devel/assets/gene_list.csv)). If you already have a GLOWgenes ranking file, you can add it directly using `--glowgenes_ranking`. Precomputed rankings based on PanelApp gene panels are available [here](https://github.com/TBLabFJD/GLOWgenes/blob/master/precomputed_panelAPP/GLOWgenes_precomputed_panelAPP.tsv). To include a specific GLOWgenes ranking, use the option `--glowgenes_ranking (path to the panel.txt)`, for example: `--glowgenes_ranking https://raw.githubusercontent.com/TBLabFJD/GLOWgenes/refs/heads/master/precomputed_panelAPP/GLOWgenes_prioritization_Neurological_ciliopathies_GA.txt`. 
  
  *note*: If a list of genes is added using `--gene_list` and GLOWgenes is enabled, the genes in the list will be assigned position 0 (regardless of whether the GLOWgenes ranking was added using `--glowgenes_ranking` or was generated from the gene list). 
  - Additionally, you can include the Gene-Disease Specificity Score (SGDS) using: `--sgds`. This score ranges from 0 to 1, where 1 indicates a gene ranks highly for only a few specific diseases (high specificity), and 0 indicates the gene consistently ranks highly across many diseases (low specificity). 

- **Structural variants (SVs) analysis (`--svs = true and --ngs_type = wgs`):** For WGS, Manta germline is used for calling structural variants (SVs) from mapped paired-end sequencing reads. It also includes the following step:
  - **SVs Annotation**: [AnnotSV](https://lbgi.fr/AnnotSV/) is used to annotate the merged results. AnnotSV needs the annotations files. They can be downloaded using `annotsv_install_annotations = true`. The path to the notes folder can be specified using `--annotsv_annotations folder_path`. If `--annotsv_annotations` is not specified, the annotations files will be downloadad directly.

- **Copy number variants (CNVs) calling** (`--svs = true and --ngs_type = wes`), with the following steps for WES:
  - **Bed file filtering**: Module to filter the bed file used for targered sequencing, to keep only the regions with a length > `--min_target` (default 20) and to exclude the regions in `--chromosomes` (default 'chrX,X,chrY,Y,chrM,MT'). 
  - **Software for detecting CNVs**: These tools require a set of samples sequenced in the same batch in order to detect changes in coverage that indicate the presence of a CNV.
    - [ExomeDepth](https://github.com/vplagnol/ExomeDepth) (`exomedepth = true`). 
    - [panelcn.MOPS](https://github.com/bioinf-jku/panelcn.mops) (`panelcmops = true`).
    - [CoNVaDING](https://github.com/molgenis/CoNVaDING) (`convading = true`).
  - **Combining the results from the various CNVs caller**
  - **CNVs Annotation**: [AnnotSV](https://lbgi.fr/AnnotSV/) is used to annotate the merged results. AnnotSV needs the annotations files. They can be downloaded using `annotsv_install_annotations = true`. The path to the notes folder can be specified using `--annotsv_annotations folder_path`. If `--annotsv_annotations` is not specified, the annotations files will be downloadad directly. 

- **Short tandem repeats:** Expansion Hunter (`--run_expansionhunter true`) for targeted genotyping of short tandem repeats (STRs) and flanking variants. 

- **Additional analysis:** Mosdepth (`--run_mosdepth true`) for calculating genome-wide sequencing coverage. Mosdepth has different argument configurations that have been defined in 3 modes: 
    - **Fast**: When only the file .quantized.bed.gz is needed. It also decompress the file automatically. Arguments used: --quantize 10: -n -x (`--run_mosdepth true --mosdepth_mode fast`).
    - **Full**: When you want a complete analysis without including a bed file. Arguments used: --quantize 10: (`--run_mosdepth true --mosdepth_mode full`).
    - **Full_with_bed**: When you want a complete analysis and you include a bed file (for CNVs analysis). Arguments used: --quantize 10: --thresholds 1,5,10,30,50,100 (`--run_mosdepth true --mosdepth_bed /path/to/intervals.bed --mosdepth_mode full_with_bed`).


# Usage

First, prepare a samplesheet with your input data:

```
sample,fastq_1,fastq_2
SAMPLE_PAIRED_END,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a pair of paired end fastq files. 

You can run the pipeline using: 

```
nextflow run CBRA/main.nf \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

For more details and further functionality, please refer to the [usage](docs/usage.md) documentation.


# Pipeline output

For details about the output files and reports, please refer to the [output](docs/output.md) documentation.

# Credits

CBRA was developed within the framework of a call for intramural cooperative and complementary actions (ACCI) funded by CIBERER (Biomedical Research Network Centre for Rare Diseases).

**Main Developer**
- [Yolanda Benítez Quesada](https://github.com/yolandabq)

**Coordinator**
- [Carlos Ruiz Arenas](https://github.com/yocra3)

**Other contributors**
- [Graciela Uría Regojo](https://github.com/guriaregojo)
- [Cristina Arias Sardá](https://github.com/CrisAriasSarda)
- [Pedro Garrido Rodríguez](https://github.com/pedro-garridor)
- [Rafa Farias Varona](https://github.com/RafaFariasVarona)
- [Pablo Minguez](https://github.com/pminguez)
- [Daniel Lopez](https://github.com/dlopez-bioinfo)


