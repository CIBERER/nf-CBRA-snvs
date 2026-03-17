# CIBERER Groups Profile Parameters

This document describes the parameters for different CIBERER group profiles. Each profile requires specific annotation databases and resource files for variant annotation.

---


## tblab Profile

This profile runs the pipeline with the annotation databases and resources used by the Translational Bioinformatics Lab (tblab) at the Instituto de Investigación Sanitaria Fundación Jiménez Díaz (IIS-FJD) (date 06/11/2025). We have included a [template](https://github.com/CIBERER/nf-CBRA-snvs/blob/conf/CIBERER_templates/tblab_parameters.config) with the required annotation files used to run this profile.

### Required Parameters

#### VEP Configuration
- **`plugins_dir`** - VEP plugins directory

#### Splicing & Functional Predictions
- **`dbscSNV`** - Splicing consensus predictions (`.txt.gz`)
- **`dbNSFP`** - Non-synonymous functional predictions (`.vcf.gz`)
  - Includes: LRT, M-CAP, MetaLR, MetaSVM, MutationAssessor, MutationTaster, PROVEAN, FATHMM, MetaRNN, PrimateAI, DEOGEN2, BayesDel, ClinPred, LIST-S2, Aloft, fathmm-MKL, fathmm-XF, PolyPhen2, phyloP, phastCons, GERP++, Interpro, GTEx eQTL
- **`dbNSFP_gene`** - Gene-level annotations

#### Population Frequencies
- **`gnomadg`** - gnomAD genomes VCF (`.vcf.gz`)
- **`gnomadg_cov`** - gnomAD genomes coverage (`.vcf.gz`)
- **`gnomade`** - gnomAD exomes VCF (`.vcf.gz`)
- **`gnomade_cov`** - gnomAD exomes coverage (`.vcf.gz`)
- **`cSVS`** - Spanish population frequencies (optional)
- **`kaviar`** - Multi-source variant database (optional)
- **`maf_FJD_COHORT`** - Internal cohort MAF (optional)

#### Pathogenicity Predictions
- **`CADD_SNVS`** - CADD SNV scores (`.tsv.gz`)
- **`CADD_INDELS`** - CADD INDEL scores (`.tsv.gz`)
- **`mutscore`** - Mutation scores (optional)
- **`revel`** - REVEL pathogenicity scores (optional)
- **`loFtool`** - Loss-of-function tolerance (optional)
- **`exACpLI`** - ExAC pLI scores (optional)

#### Clinical & Disease Databases
- **`clinvar`** - ClinVar variants (`.vcf.gz`)
- **`omim`** - OMIM gene annotations (optional)
- **`cCRS_DB`** - Constrained coding regions (optional)
- **`dENOVO_DB`** - De novo variants (optional)

#### Splicing Predictions
- **`spliceAI_SNV`** - SpliceAI SNV scores (optional)
- **`spliceAI_INDEL`** - SpliceAI INDEL scores (optional)
- **`maxEntScan`** - MaxEntScan directory (optional)

#### Additional Annotations
- **`domino`** - Gene-disease associations (optional)
- **`tissue_expression`** - GTEx expression data (optional)
- **`region_dict`** - Custom region annotations (optional)
- **`gene_panels`** - Clinical gene panels (optional)

### TBI index Files

TBI index files for VCF based files are automatically derived from the corresponding VCF paths. If they tbi file is different from the default, you can set them manually. TBI files are defined for the following parameters: `dbscSNV_tbi`, `dbNSFP_tbi`, `gnomadg_tbi`, `gnomadg_cov_tbi`, `gnomade_tbi`, `gnomade_cov_tbi`, `cSVS_tbi`, `kaviar_tbi`, `maf_FJD_COHORT_tbi`, `CADD_SNVS_TBI`, `CADD_INDELS_TBI`, `mutscore_tbi`, `revel_tbi`, `clinvar_tbi`, `cCRS_DB_tbi`, `dENOVO_DB_tbi`, `spliceAI_SNV_tbi`, `spliceAI_INDEL_tbi`



---

## Other Profiles

### [Add other CIBERER group profiles here]

_To be documented as additional profiles are configured._

---
