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

#### Find the line number where the header starts, using shell grep
start_line <- as.integer(
  system(paste0("zgrep -n -m1 '^#Uploaded_variation' ", shQuote(input), " | cut -d: -f1"),
         intern = TRUE)
)

if (length(start_line) == 0 || is.na(start_line)) {
  stop("The line starting with '#Uploaded_variation' was not found.")
}

# Fichero temporal en el mismo directorio de trabajo (con espacio garantizado),
# en vez de depender de /tmp
tmp_file <- file.path(dirname(input), "vep_body_tmp.tsv")
system(paste0("zcat ", shQuote(input), " | tail -n +", start_line, " > ", shQuote(tmp_file)))

vep <- fread(tmp_file, header = TRUE, sep = "\t",
             colClasses = "character", quote = "",
             data.table = TRUE, na.strings = c("", "NA", "-"))

file.remove(tmp_file)

if ("Uploaded_variation" %in% colnames(vep) && !("#Uploaded_variation" %in% colnames(vep))) {
  setnames(vep, "Uploaded_variation", "#Uploaded_variation")
}

# genefilter
if (!is.null(genefilter_path)){
  genefilter = read.delim(genefilter_path, header = F, stringsAsFactors = F, quote = "", check.names=F)
}

#### include GLOWgenes and SGDS 
if (!is.null(glowgenes_path)){
  
  glowgenes = read.delim(glowgenes_path, header = F, stringsAsFactors = F, quote = "", check.names=F)
  colnames(glowgenes) = c("SYMBOL", "score", "GLOWgenes")
  
  # Add 0 to the genes used to run GLOWgenes but that are not in the output file of GLOWgenes
  if (!is.null(genefilter_path)){
    genefilter$score = NA
    genefilter$GLOWgenes = 0
    colnames(genefilter) = c("SYMBOL", "score", "GLOWgenes")
  
    glowgenes = rbind(genefilter, glowgenes)
  }

  vep = merge(vep, glowgenes[c("SYMBOL", "GLOWgenes")], by= "SYMBOL", all.x = T)
}

# Gene Filter
if ((!is.null(genefilter_path)) & (is.null(glowgenes_path))){
  vep = vep[vep$SYMBOL %in% genefilter$V1,]
}

if (!is.null(SGDS_path)) {

  SGDS <- read.delim(SGDS_path, sep = ",", header = TRUE, stringsAsFactors = FALSE, quote = "", check.names = FALSE)
  colnames(SGDS) = c("SYMBOL", "SGDS", "GLOWgenes_best_ranking", "GLOWgenes_median_ranking")
  vep = merge(vep, SGDS, by = "SYMBOL", all.x = TRUE)
  
}

#### OMIM
if (!is.null(omim_path)){
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


df_out  = data.frame(row.names = 1:nrow(vep), stringsAsFactors = F)

## Save the columns that will be added later to the output 

# Remove columns starting with "SAMPLE"
columns_to_remove <- grep("^SAMPLE", colnames(vep))
# Add "USED_REF" and "Allele" to the removal
columns_to_remove <- c(columns_to_remove, which(colnames(vep) %in% c("#Uploaded_variation","USED_REF", "Allele", "SYMBOL", "Location", "VARIANT_CLASS")))

# Subset the dataframe
vep_cleaned_columns <- vep[, !columns_to_remove, with = FALSE]

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
df_out$VARIANT_CLASS = vep$VARIANT_CLASS
df_out$Panels_name = vep$panels


#=====================#
# Add all the columns #
#=====================#

df_out <- cbind(df_out, vep_cleaned_columns)

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


# #===============#
# # Pathogenicity #
# #===============#

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