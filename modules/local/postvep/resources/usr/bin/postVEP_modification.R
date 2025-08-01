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

  make_option(c("-w", "--glowgenes"), type="character", default=NULL,
            help="\t\tGLOWgenes output file to annotate and sort the results", metavar="character"),
  
  make_option(c("-s", "--SGDS"), type="character", default=NULL,
            help="\t\tGLOWgenes Score of Gene-Disease Specificity", metavar="character")

)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

input = opt$input
output = opt$output
automap_path = opt$automap
maf = opt$maf
glowgenes_path = opt$glowgenes
SGDS_path = opt$SGDS

################
# Data loading # 
################

# VEP

print("Read VEP file")

# Find the line number of the line starting with "#Uploaded_variation"
lines <- readLines(input)
start_line <- grep("^#Uploaded_variation", lines)

# Read the file starting from the detected line
if (length(start_line) > 0) {
  vep <- read.delim(input, skip = start_line - 1, header = TRUE, stringsAsFactors = F, quote = "", check.names=F, colClasses = "character")
} else {
  stop("The line starting with '#Uploaded_variation' was not found.")
}

#head(vep)
# Filtering variants by MAF

print("Number of variants before filtering by MAF")
print(nrow(vep))

vep <- vep[is.na(vep$MAX_AF) | (!is.na(vep$MAX_AF) & vep$MAX_AF < as.numeric(maf)), ]

print("Number of variants after filtering by MAF")
print(nrow(vep))


# GLOWgenes
if (!is.null(glowgenes_path)){
  
  glowgenes = read.delim(glowgenes_path, header = F, stringsAsFactors = F, quote = "", check.names=F)
  colnames(glowgenes) = c("SYMBOL", "GLOWgenes")

  vep = merge(vep, glowgenes[c("SYMBOL", "GLOWgenes")], by= "SYMBOL", all.x = T)
}

if (!is.null(SGDS_path)) {

  SGDS <- read.delim(SGDS_path, sep = ",", header = TRUE, stringsAsFactors = FALSE, quote = "", check.names = FALSE)
  colnames(SGDS) = c("SYMBOL", "SGDS", "GLOWgenes_best_ranking", "GLOWgenes_median_ranking")
  vep = merge(vep, SGDS, by = "SYMBOL", all.x = TRUE)
  
}


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

#=====================#
# Add all the columns #
#=====================#

df_out <- cbind(df_out,vep[vep_cleaned_columns])


#====================#
# Sample information #
#====================#
print("Sample information")

samples = unique(gsub("_.*$", "", gsub("^SAMPLE_", "", colnames(vep)[grepl(".*_GT$", colnames(vep), perl = T)])))
for (sample in samples){
  for (field in c("GT", "VAF", "AD", "DP", "SF", "GD")){
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