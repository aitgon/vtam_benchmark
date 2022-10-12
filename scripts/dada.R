#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#BiocManager::install("dada2")
#install.packages("devtools")
#library(devtools)
#install.packages("dplyr")

if (length(args)!=8) {
  stop("8 arguments must be supplied", call.=FALSE)
}

suppressMessages(suppressWarnings(library("dada2")))
suppressMessages(suppressWarnings(library("dplyr")))
suppressMessages(suppressWarnings(library("ShortRead")))


# path <- "~/vtam_benchmark/dalu_bat/fastq"
#setwd("~/vtam_benchmark/dalu_bat/dada")
# indir_path <- "out/dalu_bat/fastq"
indir_path = args[1]
#outdir_path = "out/dalu_bat/dada"
#track_csv_path = file.path(outdir_path, "track.csv")
track_csv_path = args[2]
#asv_fasta_path = file.path(outdir_path, "asv.fasta")
asv_fasta_path = args[3]
#asv_table_names_csv_path = file.path(outdir_path, "asv_table_names.csv")
asv_table_names_csv_path = args[4]
truncLen0 = as.numeric(args[5]) # 120 bat, 170 mfzr, 150 zfzr
# seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 175:190] for MFZR
# seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 154:169] for ZFZR
# seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 124:163] for bat
seq_length_range0 = as.numeric(args[6])
seq_length_range1 = as.numeric(args[7])
dataset = args[8] # bat or fish


dir.create(dirname(track_csv_path), showWarnings = FALSE)
dir.create(dirname(asv_fasta_path), showWarnings = FALSE)
dir.create(dirname(asv_table_names_csv_path), showWarnings = FALSE)
#list.files(path)

# read filenames
fnFs <- sort(list.files(indir_path, pattern="_fw.fastq", full.names = TRUE))
fnRs <- sort(list.files(indir_path, pattern="_rv.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_fw"), `[`, 1)
#sample.names

# make quality plots
#png(filename="plotQualityProfile_F.png")
#plotQualityProfile (fnFs[1:6])
#dev.off()
#png(filename="plotQualityProfile_R.png")
#plotQualityProfile (fnRs[1:6])
#dev.off()

# Filter and trim
filtFs <- file.path(indir_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(indir_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(truncLen0, truncLen0), maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE,matchID=TRUE)
#dim(out)

if (dataset=="bat") {
    print(dataset)
    # delete samples from the list with 0 reads
    filtFs <- filtFs[!(names(filtFs) %in% c("COI-NC_index2-3"))]
    filtRs <- filtRs[!(names(filtRs) %in% c("COI-NC_index2-3"))]
    length(filtFs)
    length(filtRs)
}

# make quality plots after filtering
#png(filename="plotQualityProfile_after_filter_F.png")
#plotQualityProfile (filtFs[1:6])
#dev.off()
#png(filename="plotQualityProfile_after_filter_R.png")
#plotQualityProfile (filtRs[1:6])
#dev.off()

# learn errors
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# plot error
#png(filename="plotError_F.png")
#plotErrors(errF, nominalQ=TRUE)
#dev.off()
#png(filename="plotError_R.png")
#plotErrors(errR, nominalQ=TRUE)
#dev.off()


# Infer ASV
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#dadaFs[[1]]

# merge
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
#head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
#dim(seqtab)
table(nchar(getSequences(seqtab)))

#select sequence length range
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq_length_range0:seq_length_range1]

# Delete chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
#dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))
sum(seqtab.nochim)/sum(seqtab)

# get read counts
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
#head(track)
# write track info 
#write.csv(track, file="track.csv")
write.csv(track, file=track_csv_path)

# write ASV table 
# write.csv(t(seqtab.nochim), file="asv_table.csv")


###
### Make files for LULU
###
# Write fasta
sq <- getSequences(seqtab.nochim)
id <- paste0(1:length(sq))
names(sq) <- id
#writeFasta(sq, file="asv.fasta")
writeFasta(sq, file=asv_fasta_path)

# add seqID instead of sequences to ASV table
transposed_seqtab.nochim <- t(seqtab.nochim)
rownames(transposed_seqtab.nochim) <- id
#write.csv(transposed_seqtab.nochim, file="asv_table_names.csv")
write.csv(transposed_seqtab.nochim, file=asv_table_names_csv_path)


q()


