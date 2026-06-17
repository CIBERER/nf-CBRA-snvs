#!/usr/bin/env Rscript

### CNV analysis 
### Author: Gonzalo Núñez Moreno



rm(list=ls())


#**********************#
# package requirements #
#**********************#

library(optparse)
library(bedr)
library(plyr)

#==============#
#   Arguments  #
#==============#
option_list=list(
  make_option(c('-i','--input'),type="character", help="Input directory."),
  make_option(c('-o','--outputfile'),type="character", help="Output directory."),
  make_option(c('-e','--extracolnames'), type="character", default=NULL, help="Extra_colnames"),
  make_option(c("-g", "--genefilter"), type="character", default=NULL, help="\t\tGene list to filter the resutls", metavar="character"),
  make_option(c("-w", "--glowgenes"), type="character", default=NULL,help="\t\tGLOWgenes output file to annotate and srt the results", metavar="character")
)

opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser) #list of the args

input <- opt$input
outputfile <- opt$outputfile
extracolnames <- opt$extracolnames
genefilter_path <- opt$genefilter
glowgenes_path <- opt$glowgenes



#==================#
#   Input loading  #
#==================#
annotated_tsv = read.delim(input, header = T, stringsAsFactors = F, quote = "")
extracolnames_tsv = read.delim(extracolnames, header = F, stringsAsFactors = F, quote = "")




#=================#
#   Change order  #
#=================#
extracolnames_positions = which(colnames(annotated_tsv) %in% extracolnames_tsv$V1)
first_cols  = match(c("ACMG_class", "Samples_ID", "SV_chrom", "SV_start", "SV_end", "SV_length", "SV_type", "N_PROGRAMS"),colnames(annotated_tsv))
annotated_tsv$ACMG_class = gsub("full=", "", annotated_tsv$ACMG_class)

new_order = c(first_cols, (max(extracolnames_positions)+1):(ncol(annotated_tsv)-1), extracolnames_positions)

annotated_tsv_ordered = annotated_tsv[,new_order]




#======================#
#   Gene Priorization  #
#======================#
if ((!is.null(genefilter_path)) & (!is.null(glowgenes_path))){
  genefilter = read.delim(genefilter_path, header = F, stringsAsFactors = F, quote = "", check.names=F)
  colnames(genefilter) = "SYMBOL"
  genefilter$GLOWgenes = 0
  glowgenes = read.delim(glowgenes_path, header = F, stringsAsFactors = F, quote = "", check.names=F)[,c(1,3)]
  colnames(glowgenes) = c("SYMBOL", "GLOWgenes")
  
  glowgenes = rbind(genefilter, glowgenes)
  rownames(glowgenes) = glowgenes$SYMBOL
  
  GLOWgenes = data.frame(row.names = 1:nrow(annotated_tsv_ordered))
  
  for(i in nrow(glowgenes):1){
    aa = glowgenes[i,]
    GLOWgenes[grep(aa$SYMBOL, annotated_tsv_ordered$Gene_name),1] = aa$GLOWgenes
  }
  
  annotated_tsv_ordered = cbind(GLOWgenes, annotated_tsv_ordered)
  colnames(annotated_tsv_ordered)[1]  = "GLOWgenes"
} else if ((!is.null(genefilter_path))) {
  genefilter = read.delim(genefilter_path, header = F, stringsAsFactors = F, quote = "", check.names=F)
  annotated_tsv_ordered = annotated_tsv_ordered[grep(paste(genefilter$V1, collapse = "|"), annotated_tsv_ordered$Gene_name),]
  
}




#=================#
#   Write output  #
#=================#
write.table(annotated_tsv_ordered, outputfile, col.names = T, row.names = F, quote = F, sep = "\t", na = "")

library(openxlsx)
write.xlsx(annotated_tsv_ordered, paste0(outputfile, ".xlsx"), colNames = TRUE)


print("Finish")