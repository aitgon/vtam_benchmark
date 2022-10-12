#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#install_github("tobiasgf/lulu")
#install.packages("dplyr")
#install.packages("lulu")


if (length(args)!=3) {
  stop("3 arguments must be supplied", call.=FALSE)
}

suppressMessages(suppressWarnings(library("dplyr")))
suppressMessages(suppressWarnings(library("lulu")))

#setwd("~/vtam_benchmark/dalu_bat/dada")

#wdir_path = "out/dalu_bat/dada"
#asv_table_names_csv_path = file.path(wdir_path, "asv_table_names.csv")
asv_table_names_csv_path = args[1]
#asv_auto_blast_tsv_path = file.path(wdir_path, "asv_auto_blast.tsv")
asv_auto_blast_tsv_path = args[2]
#asv_table_lulu_csv_path = file.path(wdir_path, "asv_table_lulu.csv")
asv_table_lulu_csv_path = args[3]

otutab <- read.csv(asv_table_names_csv_path, sep=',', header=TRUE, as.is=TRUE, row.names = 1, check.names=FALSE)
#colnames(otutab)
matchlist <- read.table(asv_auto_blast_tsv_path, sep='\t', header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)
lulu_results <- lulu(otutab, matchlist)

write.csv(lulu_results$curated_table, file=asv_table_lulu_csv_path)

q()


