#!/usr/bin/env Rscript

# Author: Gonzalo Nunez Moreno
# Date: 16/03/2021
## Reedited by: Yolanda Benítez Quesada
## Date: 23/05/2025

library(optparse)
library(data.table)

#############
# Arguments # 
#############
option_list = list(

  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="\t\t input TSV file (output from VEP)", metavar="character"),
  
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="\t\t Output file", metavar="character"),

  make_option(c("-a", "--automap"), type="character", default=NULL,
              help="\t\tAutomap output file (optional)", metavar="character"),
  
  make_option(c("-f", "--maf"), type="double", default=0.1,
               help="\t\tMinimum allele frequency to filter", metavar="character"),
  
  make_option(c("-s", "--SGDS"), type="character", default=NULL,
            help="\t\tGLOWgenes Score of Gene-Disease Specificity", metavar="character"),

  make_option(c("-d", "--dbNSFPgene"), type="character", default=NULL, 
              help="\t\tdbNSFP_gene file", metavar="character"),
  
  make_option(c("-r", "--regiondict"), type="character", default=NULL,
              help="\t\tRegion dictionary", metavar="character"),
  
  make_option(c("-m", "--omim"), type="character", default=NULL,
              help="\t\tOMIM information", metavar="character"),
  
  make_option(c("-D", "--domino"), type="character", default=NULL, 
              help="\t\tdomino file", metavar="character"),

  make_option(c("-e", "--expression"), type="character", default=NULL, 
            help="\t\ttissue expression file", metavar="character"),
  
  make_option(c("-g", "--genefilter"), type="character", default=NULL,
              help="\t\tGene list to filter the resutls", metavar="character"),
  
  make_option(c("-w", "--glowgenes"), type="character", default=NULL,
              help="\t\tGLOWgenes output file to annotate and srt the results", metavar="character"),

  make_option(c("-p", "--panel_annotation_file"), type="character", default=NULL,
              help="\t\tGene-Panel file to annotate", metavar="character")

)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

input = opt$input
output = opt$output
automap_path = opt$automap
maf = opt$maf
glowgenes_path = opt$glowgenes
SGDS_path = opt$SGDS

dbNSFPgenepath <- opt$dbNSFPgene
dominopath <- opt$domino
expression_path <- opt$expression
dict_region_path <- opt$regiondict
omim_path = opt$omim
genefilter_path = opt$genefilter
panels_path = opt$panel_annotation_file

################
# Data loading # 
################

##### VEP

print("Read VEP file")

#### Find the line number of the line starting with "#Uploaded_variation"
lines <- readLines(input)
start_line <- grep("^#Uploaded_variation", lines)

#### Read the file starting from the detected line
if (length(start_line) > 0) {
  vep <- read.delim(input, skip = start_line - 1, header = TRUE, stringsAsFactors = F, quote = "", check.names=F, colClasses = "character")
} else {
  stop("The line starting with '#Uploaded_variation' was not found.")
}


###############
## Filtering ##
###############


#### Filtering variants by MAF

print("Number of variants before filtering by MAF")
print(nrow(vep))

#vep <- vep[is.na(vep$MAX_AF) | (!is.na(vep$MAX_AF) & vep$MAX_AF < as.numeric(maf)), ]

vep$gnomADe_AF_grpmax = as.numeric(unlist(lapply(vep$gnomADe_AF_grpmax, function(x) strsplit(x, ",")[[1]][1])))
vep$gnomADg_AF_grpmax = as.numeric(unlist(lapply(vep$gnomADg_AF_grpmax, function(x) strsplit(x, ",")[[1]][1])))

vep = vep[is.na(vep$gnomADe_AF_grpmax) | as.numeric(vep$gnomADe_AF_grpmax) < as.numeric(maf) | vep$gnomADe_filt != "PASS",]
print(nrow(vep))
vep = vep[is.na(vep$gnomADg_AF_grpmax) | as.numeric(vep$gnomADg_AF_grpmax) < as.numeric(maf) | vep$gnomADg_filt != "PASS",]
print(nrow(vep))

print("Number of variants after filtering by MAF")
print(nrow(vep))

## Filtering variants by gene panel if included without GLOWgenes ranking

# genefilter
if (!is.null(genefilter_path)){
  genefilter = read.delim(genefilter_path, header = F, stringsAsFactors = F, quote = "", check.names=F)
}

# Gene Filter
if ((!is.null(genefilter_path)) & (is.null(glowgenes_path))){
  vep = vep[vep$SYMBOL %in% genefilter$V1,]
}

## Check if there are remaining variants
if (nrow(vep) == 0) {
  stop("There are no remaining variants")
}

#########################
## Include custom info ##
#########################

# dbNSFP gene
dbNSFP_gene = read.delim(dbNSFPgenepath, header = TRUE, stringsAsFactors = F, quote = "")
vep = merge(vep, dbNSFP_gene, by.x = "SYMBOL", by.y = "Gene_name", all.x = T)

#### include GLOWgenes and SGDS if included
if (!is.null(glowgenes_path)){
  print("Include GLOWgenes ranking")
  glowgenes = read.delim(glowgenes_path, header = F, stringsAsFactors = F, quote = "", check.names=F)
  colnames(glowgenes) = c("SYMBOL", "GLOWgenes")

  vep = merge(vep, glowgenes[c("SYMBOL", "GLOWgenes")], by= "SYMBOL", all.x = T)
}

if (!is.null(SGDS_path)) {
  print("Include GLOWgenes SGDS")
  SGDS <- read.delim(SGDS_path, sep = ",", header = TRUE, stringsAsFactors = FALSE, quote = "", check.names = FALSE)
  colnames(SGDS) = c("SYMBOL", "SGDS", "GLOWgenes_best_ranking", "GLOWgenes_median_ranking")
  vep = merge(vep, SGDS, by = "SYMBOL", all.x = TRUE)
  
}

#### OMIM
if (!is.null(omim_path)){
  print("Include OMIM")
  omim = read.delim(omim_path, header = F, stringsAsFactors = F, comment.char = "#", quote = "", check.names=F)
  colnames(omim) = c("Chromosome", "Genomic_Position_Start", "Genomic Position End", "Cyto_Location", "Computed_Cyto_Location", "MIM_Number",
    "Gene_Symbols", "Gene_Name",	"Approved_Gene_Symbol", "Entrez_Gene_ID", "Ensembl_Gene_ID", "Comments", "Phenotypes", "Mouse_Gene_Symbol-ID")
  vep = merge(vep, omim, by.x = "SYMBOL", by.y = "Approved_Gene_Symbol", all.x = T)
}

#### Region dictionary
if (!is.null(dict_region_path)){
  dict_region = read.csv(dict_region_path, header = F, sep = ",", stringsAsFactors = F)
  priority_list = c("SPLICING", "5UTR", "3UTR", "ncRNA", "regulatory", "UPSTREAM", "DOWNSTREAM", "EXONIC", "INTRONIC", "INTERGENIC", "-")
  dict_region$V2[is.na(dict_region$V2)] = "-"
  dict_region$V2 = factor(dict_region$V2, priority_list)
  rownames(dict_region) = dict_region$V1
}

# include domino
domino = read.delim(dominopath, header = TRUE, stringsAsFactors = F, quote = "")
vep = merge(vep, domino, by.x = "SYMBOL", by.y = "Gene_name", all.x = T)

# include tissue expression
expression = read.delim(expression_path, header = TRUE, stringsAsFactors = F, quote = "")
vep = merge(vep, expression, by.x = "SYMBOL", by.y = "Gene.name", all.x = T)


### create output dataframe
df_out  = data.frame(row.names = 1:nrow(vep), stringsAsFactors = F)

## Save the columns that will be added later to the output 

# Remove columns starting with "SAMPLE"
columns_to_remove <- grep("^SAMPLE", colnames(vep))
# Add "USED_REF" and "Allele" to the removal
columns_to_remove <- c(columns_to_remove, which(colnames(vep) %in% c("#Uploaded_variation","USED_REF", "Allele")))

# Subset the dataframe
vep_cleaned_columns <- colnames(vep[, -columns_to_remove])

# View the cleaned dataframe
print((vep_cleaned_columns))


#==================================#
# Basic information of the variant #
#==================================#
print("Basic information of the variant")

df_out$CHROM = unlist(lapply(vep$Location, function(x) strsplit(x, ":")[[1]][1]))
df_out$POS = as.numeric(unlist(lapply(vep$`#Uploaded_variation`, function(x) rev(strsplit(x, "_")[[1]])[2])))
df_out$REF = vep$USED_REF
df_out$ALT = vep$Allele

df_out$Location = vep$Location
df_out$SYMBOL = vep$SYMBOL
df_out$Gene_full_name = vep$Gene_full_name
if (!is.null(glowgenes_path)) df_out$GLOWgenes = vep$GLOWgenes
if ((!is.null(genefilter_path)) & (!is.null(glowgenes_path))) df_out$GLOWgenes[df_out$SYMBOL %in% genefilter$V1] = 0 # We assume the genes of the list are the genes from the panel
df_out$VARIANT_CLASS = vep$VARIANT_CLASS
df_out$Panels_name = vep$panels


#=====================#
# Add all the columns #
#=====================#

#df_out <- cbind(df_out,vep[vep_cleaned_columns])

#=====================#
# Feature information #
#=====================#
print("Feature information")

df_out$Existing_variation = vep$Existing_variation
if (!is.null(dict_region_path)) df_out$Genomic_region = unlist(lapply(vep$Consequence, function(x) as.character(dict_region[strsplit(x, ",")[[1]],2])[which.min(dict_region[strsplit(x, ",")[[1]],2])]))
df_out$CANONICAL = vep$CANONICAL
df_out$Feature = vep$Feature
df_out$Feature_type = vep$Feature_type
df_out$BIOTYPE = vep$BIOTYPE
df_out$Consequence = vep$Consequence
df_out$INTRON = vep$INTRON
df_out$EXON = vep$EXON
df_out$HGVSc = vep$HGVSc
df_out$HGVSp = vep$HGVSp
df_out$DISTANCE = as.numeric(vep$DISTANCE)
df_out$STRAND = vep$STRAND
df_out$Interpro_domain = vep$Interpro_domain
df_out$Interpro_domain = vep$Interpro_domain
df_out$Domino_Score = vep$Domino_Score


#===============#
# Pathogenicity #
#===============#
print("Pathogenicity")

df_out$CLNSIG = vep$ClinVar_CLNSIG
df_out$CLNREVSTAT = vep$ClinVar_CLNREVSTAT
df_out$CLNDN = vep$ClinVar_CLNDN
df_out$CLNSIGCONF = vep$ClinVar_CLNSIGCONF                                                                      
df_out$OMIM_phenotype = vep$Phenotypes
df_out$Orphanet_disorder = vep$Orphanet_disorder
df_out$Orphanet_association_type = vep$Orphanet_association_type
df_out$HPO_name = vep$HPO_name
df_out$PUBMED = vep$PUBMED





#=============#
# Frequencies #
#=============#
print("Frequencies")

df_out$gnomADg_AF =  as.numeric(unlist(lapply(vep$gnomADg_AF, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADg_AC = as.numeric(unlist(lapply(vep$gnomADg_AC, function(x) strsplit(x, ",")[[1]][1]))) 
df_out$gnomADg_AN = as.numeric(unlist(lapply(vep$gnomADg_AN, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADg_nhomalt = as.numeric(unlist(lapply(vep$gnomADg_nhomalt, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADg_cov_median = round(unlist(lapply(vep$gnomADg_cov_median, function(x) mean(as.numeric(strsplit(gsub("(?<![eE])-","0",x, perl = T), ",")[[1]])))))
df_out$gnomADg_cov_perc_20x = round(unlist(lapply(vep$gnomADg_cov_perc_20x, function(x) mean(as.numeric(strsplit(gsub("(?<![eE])-","0",x, perl = T), ",")[[1]])))),2)
df_out$gnomADg_filter = vep$gnomADg_filt
df_out$gnomADg_popmax = vep$gnomADg_grpmax
df_out$gnomADg_AF_popmax = vep$gnomADg_AF_grpmax
df_out$gnomADg_AC_popmax = as.numeric(unlist(lapply(vep$gnomADg_AC_grpmax, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADg_AF_nfe = as.numeric(unlist(lapply(vep$gnomADg_AF_nfe, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADg_AC_nfe = as.numeric(unlist(lapply(vep$gnomADg_AC_nfe, function(x) strsplit(x, ",")[[1]][1])))

df_out$gnomADe_AF = as.numeric(unlist(lapply(vep$gnomADe_AF, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADe_AC = as.numeric(unlist(lapply(vep$gnomADe_AC, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADe_AN = as.numeric(unlist(lapply(vep$gnomADe_AN, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADe_nhomalt = as.numeric(unlist(lapply(vep$gnomADe_nhomalt, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADe_cov_median = round(unlist(lapply(vep$gnomADe_cov_median, function(x) mean(as.numeric(strsplit(gsub("(?<![eE])-","0",x, perl = T), ",")[[1]])))))
df_out$gnomADe_cov_perc_20x = round(unlist(lapply(vep$gnomADe_cov_perc_20x, function(x) mean(as.numeric(strsplit(gsub("(?<![eE])-","0",x, perl = T), ",")[[1]])))),2)
df_out$gnomADe_filter = vep$gnomADe_filt
df_out$gnomADe_popmax = vep$gnomADe_grpmax
df_out$gnomADe_AF_popmax = vep$gnomADe_AF_grpmax
df_out$gnomADe_AC_popmax = as.numeric(unlist(lapply(vep$gnomADe_AC_grpmax, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADe_AF_nfe = as.numeric(unlist(lapply(vep$gnomADe_AF_nfe, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADe_AC_nfe = as.numeric(unlist(lapply(vep$gnomADe_AC_nfe, function(x) strsplit(x, ",")[[1]][1])))

df_out$kaviar_AF = vep$kaviar_AF
df_out$kaviar_AC = vep$kaviar_AC
df_out$CSVS_AF = as.numeric(unlist(lapply(vep$CSVS_AF, function(x) strsplit(x, ",")[[1]][1])))
df_out$CSVS_AC = as.numeric(unlist(lapply(vep$CSVS_AC, function(x) strsplit(x, ",")[[1]][1])))
df_out$FJD_MAF_AF = as.numeric(unlist(lapply(vep$FJD_MAF_AF, function(x) strsplit(x, ",")[[1]][1])))
df_out$FJD_MAF_AC = as.numeric(unlist(lapply(vep$FJD_MAF_AC, function(x) strsplit(x, ",")[[1]][1])))
##add new columns del MAF_FJD de DHR vs pseudocontroles (SON DE OJO LOS PSEUDOCONTROLES)
df_out$FJD_MAF_AF_DS_IRD = as.numeric(unlist(lapply(vep$FJD_MAF_AF_DS_irdt, function(x) strsplit(x, ",")[[1]][1])))
df_out$FJD_MAF_AC_DS_IRD = as.numeric(unlist(lapply(vep$FJD_MAF_AC_DS_irdt, function(x) strsplit(x, ",")[[1]][1])))
df_out$FJD_MAF_AF_P_IRD = as.numeric(unlist(lapply(vep$FJD_MAF_AF_P_eyeg, function(x) strsplit(x, ",")[[1]][1])))
df_out$FJD_MAF_AC_P_IRD = as.numeric(unlist(lapply(vep$FJD_MAF_AC_P_eyeg, function(x) strsplit(x, ",")[[1]][1])))
                                             
df_out$denovoVariants_SAMPLE_CT = vep$denovoVariants_SAMPLE_CT





#==========================#
# Pathogenicity prediction #
#==========================#
print("Pathogenicity prediction")

df_out$CADD_PHRED = as.numeric(vep$CADD_PHRED)
df_out$CADD_RAW = as.numeric(vep$CADD_RAW)
df_out$MutScore = as.numeric(vep$Mut_Score)
df_out$REVELScore = as.numeric(vep$REVEL_Score)

patho_norm_func = function(predictions){
  predictions = gsub(";", ",", predictions)
  predictions = tolower(predictions)
  unlist(lapply(predictions, function(x) {
    y = strsplit(x,",")[[1]]
    y = y[!y %in% c("-", ".", "U", "")]
    y[y %in% c("tolerated", "tolerated_low_confidence", "benign", "n", "l", "p", "t")] = "T"
    y[y %in% c("deleterious", "deleterious_low_confidence", "probably_damaging", 
               "possibly_damaging", "a", "m", "h", "d", "Dominant", "Recessive")] = "D"
    if ("D" %in% y) { return("D") }
    else if ("T" %in% y) { return("T") }
    else {return("")}
  }))
}  

df_pathogenic_predictors = data.frame(row.names = 1:nrow(vep))
df_pathogenic_predictors$SIFT = patho_norm_func(vep$SIFT)
df_pathogenic_predictors$PolyPhen = patho_norm_func(vep$PolyPhen)
df_pathogenic_predictors$Polyphen2_HDIV_pred = patho_norm_func(vep$Polyphen2_HDIV_pred)
df_pathogenic_predictors$Polyphen2_HVAR_pred = patho_norm_func(vep$Polyphen2_HVAR_pred)
df_pathogenic_predictors$LRT_pred = patho_norm_func(vep$LRT_pred)
df_pathogenic_predictors$`M-CAP_pred` = patho_norm_func(vep$`M-CAP_pred`)
df_pathogenic_predictors$MetaLR_pred = patho_norm_func(vep$MetaLR_pred)
df_pathogenic_predictors$MetaSVM_pred = patho_norm_func(vep$MetaSVM_pred)
df_pathogenic_predictors$MutationAssessor_pred = patho_norm_func(vep$MutationAssessor_pred)
df_pathogenic_predictors$MutationTaster_pred = patho_norm_func(vep$MutationTaster_pred)
df_pathogenic_predictors$PROVEAN_pred = patho_norm_func(vep$PROVEAN_pred)
df_pathogenic_predictors$FATHMM_pred = patho_norm_func(vep$FATHMM_pred)
df_pathogenic_predictors$MetaRNN_pred = patho_norm_func(vep$MetaRNN_pred)
df_pathogenic_predictors$PrimateAI_pred = patho_norm_func(vep$PrimateAI_pred)
df_pathogenic_predictors$DEOGEN2_pred = patho_norm_func(vep$DEOGEN2_pred)
df_pathogenic_predictors$BayesDel_addAF_pred = patho_norm_func(vep$BayesDel_addAF_pred)
df_pathogenic_predictors$BayesDel_noAF_pred = patho_norm_func(vep$BayesDel_noAF_pred)
df_pathogenic_predictors$ClinPred_pred = patho_norm_func(vep$ClinPred_pred)
df_pathogenic_predictors$`LIST-S2_pred` = patho_norm_func(vep$`LIST-S2_pred`)
df_pathogenic_predictors$Aloft_pred = patho_norm_func(vep$Aloft_pred)
df_pathogenic_predictors$`fathmm-MKL_coding_pred` = patho_norm_func(vep$`fathmm-MKL_coding_pred`)
df_pathogenic_predictors$`fathmm-XF_coding_pred` = patho_norm_func(vep$`fathmm-XF_coding_pred`)
  

df_out$N_Pathogenic_pred = apply(df_pathogenic_predictors, 1, function(x) table(x)["D"])
df_out$N_Pathogenic_pred[is.na(df_out$N_Pathogenic_pred)] = 0
df_out$N_Benign_pred = apply(df_pathogenic_predictors, 1, function(x) table(x)["T"])
df_out$N_Benign_pred[is.na(df_out$N_Benign_pred)] = 0
df_out$N_predictions = df_out$N_Pathogenic_pred + df_out$N_Benign_pred
df_out$Pathogenic_pred = apply(df_pathogenic_predictors, 1, function(x) paste(names(x)[which(x == "D")], collapse = ","))
df_out$Benign_pred = apply(df_pathogenic_predictors, 1, function(x) paste(names(x)[which(x == "T")], collapse = ","))





#=====================#
# Splicing predictors #
#=====================#
# print("Splicing predictors")

# # Select one splice prediction per row
# for (j in c("SpliceAI_SNV_SpliceAI", "SpliceAI_INDEL_SpliceAI")){
#   multi_gene_sites = grep(",",vep[,j])
#   for (i in multi_gene_sites){
#     splice_predictions = do.call("rbind",strsplit(strsplit(vep[i,j], ",")[[1]], "|",fixed = T))
#     if (vep$SYMBOL[i] %in% splice_predictions[,2]) {
#       vep[i,j] = paste(splice_predictions[splice_predictions[,2] == vep$SYMBOL[i],,drop = F][1,],collapse = "|")
#     } else {
#       max_value_row = which(splice_predictions[,3:6] == max(splice_predictions[,3:6]), arr.ind = TRUE)[1,1]
#       vep[i,j] = paste(splice_predictions[max_value_row,],collapse = "|")
#     }
#   }
# }
# # Merge SpliceAI predictions for INDELs and SNVs 
# vep$SpliceAI_INDEL_SpliceAI[vep$SpliceAI_INDEL_SpliceAI == "-"] = vep$SpliceAI_SNV_SpliceAI[vep$SpliceAI_INDEL_SpliceAI == "-"]
# # Create data.frame with separated SpliceAI values
# SpliceAI = data.frame(do.call("rbind", strsplit(vep$SpliceAI_INDEL_SpliceAI, "|", fixed = T)), stringsAsFactors = F)
# colnames(SpliceAI) = c("ALLELE", "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL")
# df_out$SpliceAI_SYMBOL = SpliceAI$SYMBOL
# df_out$SpliceAI_DS_AG = as.numeric(SpliceAI$DS_AG)
# df_out$SpliceAI_DS_AL = as.numeric(SpliceAI$DS_AL)
# df_out$SpliceAI_DS_DG = as.numeric(SpliceAI$DS_DG)
# df_out$SpliceAI_DS_DL = as.numeric(SpliceAI$DS_DL)
# df_out$SpliceAI_DS_Max = apply(df_out[c("SpliceAI_DS_AG", "SpliceAI_DS_AL", "SpliceAI_DS_DG", "SpliceAI_DS_DL")], 1, max)                           
# df_out$SpliceAI_DP_AG = as.numeric(SpliceAI$DP_AG)
# df_out$SpliceAI_DP_AL = as.numeric(SpliceAI$DP_AL)
# df_out$SpliceAI_DP_DG = as.numeric(SpliceAI$DP_DG)
# df_out$SpliceAI_DP_DL = as.numeric(SpliceAI$DP_DL)

# df_out$ada_score = as.numeric(vep$ada_score)
# df_out$rf_score = as.numeric(vep$rf_score)
# df_out$MaxEntScan_alt = as.numeric(vep$MaxEntScan_alt)
# df_out$MaxEntScan_diff = as.numeric(vep$MaxEntScan_diff)
# df_out$MaxEntScan_ref = as.numeric(vep$MaxEntScan_ref)

print("Splicing predictors")

# Select one splice prediction per row
for (j in c("SpliceAI_SNV_SpliceAI", "SpliceAI_INDEL_SpliceAI")){
  multi_gene_sites = grep(",", vep[,j])
  
  for (i in multi_gene_sites){
    splice_predictions = do.call("rbind",
                                 strsplit(strsplit(vep[i,j], ",")[[1]], "|", fixed = TRUE))
    
    if (vep$SYMBOL[i] %in% splice_predictions[,2]) {
      vep[i,j] = paste(
        splice_predictions[splice_predictions[,2] == vep$SYMBOL[i], , drop = FALSE][1,],
        collapse = "|"
      )
    } else {
      max_value_row = which(
        splice_predictions[,3:6] == max(splice_predictions[,3:6]),
        arr.ind = TRUE
      )[1,1]
      
      vep[i,j] = paste(splice_predictions[max_value_row,], collapse = "|")
    }
  }
}

# Merge SpliceAI predictions for INDELs and SNVs 
vep$SpliceAI_INDEL_SpliceAI[vep$SpliceAI_INDEL_SpliceAI == "-"] =
  vep$SpliceAI_SNV_SpliceAI[vep$SpliceAI_INDEL_SpliceAI == "-"]

# ----------- FIX CLAVE AQUÍ -----------

split_spliceai <- function(x) {
  if (is.na(x) || x == "-") {
    return(rep(NA, 10))
  }
  parts <- strsplit(x, "|", fixed = TRUE)[[1]]
  length(parts) <- 10
  return(parts)
}

SpliceAI = data.frame(
  do.call("rbind", lapply(vep$SpliceAI_INDEL_SpliceAI, split_spliceai)),
  stringsAsFactors = FALSE
)

colnames(SpliceAI) = c(
  "ALLELE", "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL",
  "DP_AG", "DP_AL", "DP_DG", "DP_DL"
)

# -------------------------------------

df_out$SpliceAI_SYMBOL = SpliceAI$SYMBOL
df_out$SpliceAI_DS_AG = as.numeric(SpliceAI$DS_AG)
df_out$SpliceAI_DS_AL = as.numeric(SpliceAI$DS_AL)
df_out$SpliceAI_DS_DG = as.numeric(SpliceAI$DS_DG)
df_out$SpliceAI_DS_DL = as.numeric(SpliceAI$DS_DL)

df_out$SpliceAI_DS_Max = apply(
  df_out[c("SpliceAI_DS_AG", "SpliceAI_DS_AL", "SpliceAI_DS_DG", "SpliceAI_DS_DL")],
  1,
  max,
  na.rm = TRUE
)

df_out$SpliceAI_DP_AG = as.numeric(SpliceAI$DP_AG)
df_out$SpliceAI_DP_AL = as.numeric(SpliceAI$DP_AL)
df_out$SpliceAI_DP_DG = as.numeric(SpliceAI$DP_DG)
df_out$SpliceAI_DP_DL = as.numeric(SpliceAI$DP_DL)

df_out$ada_score = as.numeric(vep$ada_score)
df_out$rf_score = as.numeric(vep$rf_score)
df_out$MaxEntScan_alt = as.numeric(vep$MaxEntScan_alt)
df_out$MaxEntScan_diff = as.numeric(vep$MaxEntScan_diff)
df_out$MaxEntScan_ref = as.numeric(vep$MaxEntScan_ref)



#============================#
# Conservation and phylogeny #
#============================#
print("Conservation and phylogeny")

df_out$LoFtool = as.numeric(vep$LoFtool)
#antiguo "ExACpLI ahora es pLI_gene_value
df_out$ExACpLI = as.numeric(vep$pLI_gene_value)
df_out$gnomAD_exomes_CCR = vep$gnomAD_exomes_CCR
df_out$phastCons30way_mammalian = as.numeric(vep$phastCons470way_mammalian)
df_out$phyloP30way_mammalian = as.numeric(vep$phyloP470way_mammalian)
df_out$MGI_mouse_phenotype = vep$MGI_mouse_phenotype_filt





#=====================================================#
#Expression, process, route, function and interaction #
#=====================================================#
print("Expression, process, route, function and interaction")

df_out$GTEx_V8_gene = vep$GTEx_V8_eQTL_gene
df_out$GTEx_V8_tissue = vep$GTEx_V8_eQTL_tissue
df_out$`Expression_GNF-Atlas` = vep$Expression.GNF.Atlas.
df_out$Pathway_KEGG = vep$Pathway.KEGG._full
df_out$GO_biological_process = vep$GO_biological_process
df_out$GO_cellular_component = vep$GO_cellular_component
df_out$GO_molecular_function = vep$GO_molecular_function
df_out$Interactions_IntAct = vep$Interactions.IntAct.

df_out$retina_RNA_tissue_consensus = round(vep$retina,2)
df_out$testis_RNA_tissue_consensus = round(vep$testis,2)
df_out$kidney_RNA_tissue_consensus = round(vep$kidney,2)
df_out$brain_max_RNA_tissue_consensus = round(vep$brain_max,2)
df_out$glands_max_RNA_tissue_consensus = round(vep$glands_max,2)
df_out$digestive_max_RNA_tissue_consensus = round(vep$digestive_max,2)
df_out$heart_RNA_tissue_consensus = round(vep$heart.muscle,2)
df_out$liver_RNA_tissue_consensus = round(vep$liver,2)
df_out$lung_RNA_tissue_consensus = round(vep$lung,2)
df_out$pancreas_RNA_tissue_consensus = round(vep$pancreas,2)
df_out$skel_muscle_RNA_tissue_consensus = round(vep$skeletal.muscle,2)
df_out$skin_RNA_tissue_consensus = round(vep$skin,2)
df_out$mean_expression_RNA_tissue_consensus = round(vep$mean_exp,2)
df_out$retina_ratio_exp_RNA_tissue_consensus = round(vep$retina_ratio, 2)
                           
df_out$Original_pos = vep$SAMPLE_Original_pos
df_out$variant_id = vep$SAMPLE_variant_id



#====================#
# Sample information #
#====================#
print("Sample information")

samples = unique(gsub("_.*$", "", gsub("^SAMPLE_", "", colnames(vep)[grepl(".*_GT$", colnames(vep), perl = T)])))
for (sample in samples){
  for (field in c("GT", "VAF", "AD", "DP", "SF", "GD", "GQ", "FT")){
    tryCatch(
      {
        print(paste0(sample, "_", field))
        df_out[,paste0(sample, "_", field)] = vep[,paste0("SAMPLE_", sample, "_", field)]
      },
      error=function(e) print(paste0("There is no ", field, " information of the sample ", sample)),
      warning=function(e) print(paste0("There is no ", field, " information of the sample ", sample))
      )
  }
  df_out[,paste0(sample, "_ROH")] <- "NaN"
  # Sacar del output de autopmap
  tryCatch(
    {
      automap = read.delim(automap_path, header = F, comment.char = "#", stringsAsFactors = F)
      df_out[,paste0(sample,"_ROH")] = "False"
      for (i in 1:nrow(automap)) {
        df_out[df_out$POS >= automap$V2[i] & df_out$POS <= automap$V3[i] & gsub("chr","",df_out$CHROM) == gsub("chr","",automap$V1[i]), paste0(sample,"_ROH")] = "True"
      }
    },
    error = function(e) {
    # Handle errors: File not found or other issues
    print(paste0("There is no AutoMap information for the sample ", sample))
    },
    warning = function(w) {
    # Handle warnings
    print(paste0("Warning encountered for the sample ", sample, ": ", conditionMessage(w)))
    }
  )
}

df_out$hiConfDeNovo = vep$SAMPLE_hiConfDeNovo
df_out$loConfDeNovo = vep$SAMPLE_loConfDeNovo


#==============================================#
#Extra sample information (individual callers) #
#==============================================#

# Extract column names that match SAMPLE_* pattern
sample_columns <- grep("^SAMPLE_*", colnames(vep), value = TRUE)
program_suffixes <- gsub("SAMPLE_", "", sample_columns)

for (sample in samples) {
  program_suffixes <-  gsub(paste0(sample, "_"), "", program_suffixes)
}
# Extract the {program}_{suffix} part only if it exists after {samplename}
program_suffixes <- program_suffixes[program_suffixes != "variant_id" & program_suffixes != "Original_pos"]
program_suffixes_field <- grep("_", program_suffixes, value = TRUE)

for (sample in samples){
  for (program_field in program_suffixes_field){
    tryCatch(
      {
        print(paste0(sample, "_", program_field))
        df_out[,paste0(sample, "_", program_field)] = vep[,paste0("SAMPLE_", sample, "_", program_field)]
      },
      error=function(e) print(paste0("There is no ", program_field, " information of the sample ", sample)),
      warning=function(e) print(paste0("There is no ", program_field, " information of the sample ", sample))
    )
  }
} 

df_out$Original_pos = vep$SAMPLE_Original_pos
df_out$variant_id = vep$SAMPLE_variant_id


# #======== #
# # Sorting #
# #======== #

## Sort the output
if (!is.null(glowgenes_path)){
  df_out = df_out[order(df_out$GLOWgenes, df_out$POS),]
} else {
  df_out = df_out[order(df_out$CHROM, df_out$POS),]
}


#==============#
# Write output #
#==============#
df_out[df_out=="-"] = NA
write.table(df_out, output, sep = "\t", col.names = T, row.names = F, quote = F, na = "")