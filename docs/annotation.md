# nf-CBRA-snvs: Documentation about the SNV/Indel annotation subworkflow

## [snv_annotation subworkflow](https://github.com/CIBERER/GdTBioinfo-nf-snvs/blob/annot_subworkflow_postvep/subworkflows/local/snv_annotation/main.nf)

This sub-workflow annotate variants from a vcf file with [Ensembl Variant Effect Predictor (Ensembl VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html). 
By default, VEP is run with the option "--everything", which is a shorcut flag to switch on all of the following:

```
--sift b, --polyphen b, --ccds, --hgvs, --symbol, --numbers, --domains, --regulatory, --canonical, --protein, --biotype, --af, --af_1kg, --af_esp, --af_gnomade, --af_gnomadg, --max_af, --pubmed, --uniprot, --mane, --tsl, --appris, --variant_class, --gene_phenotype, --mirna
```

For more information about the annotations, see the [VEP tutorial](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html)

**Cache**

Cache files path for VEP can be indicated with cache_path param. If not, cache can be downloaded automatically. It will be stored in the folder *vep_cache* inside the outdir. If interested in RefSeq transcripts you may download homo_sapiens_refseq cache file. For that, use the param *refseq = true*, which will download homo_sapiens_refseq cache and will add the flag *--refseq* to the *ENSEMBLVEP_VEP* process. 

**Extra annotations**

If you want to add more annotations, you can do it with --plugins and --custom options in VEP. The extra files can be added to work folder with the params --custom_extra_files and --extra_files. 

The annotation can be customized adding a configuration file. An example can be seen [here](https://github.com/CIBERER/GdTBioinfo-nf-snvs/blob/annot_subworkflow_postvep/subworkflows/local/snv_annotation/tests/nextflow.config). 
