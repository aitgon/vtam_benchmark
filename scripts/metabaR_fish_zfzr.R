args = commandArgs(trailingOnly=TRUE)

#install.packages("devtools")
#devtools::install_github("metabaRfactory/metabaR")

#install.packages("BiocManager")
#BiocManager::install("biomformat")
#install.packages("remotes")
#remotes::install_github("metabaRfactory/metabaR")
#install.packages("prettydoc")

if (length(args)!=5) {
  stop("5 arguments must be supplied", call.=FALSE)
}

file_reads = args[1]
file_motus = args[2]
file_pcrs = args[3]
file_samples = args[4]
obibar_base = args[5]

#library(help="metabaR")
library(metabaR) 
library(ggplot2)
library(reshape2)
library(cowplot)
library(igraph)

#setwd("~/vtam_benchmark/obibar_fish/zfzr_obiout/")
###  import
#mblist <- tabfiles_to_metabarlist(file_reads = "reads_obi.tsv",
#                                    file_motus = "motus_obi_taxa2.tsv",
#                                    file_pcrs = "~/vtam_benchmark/metafiles/pcrs_adapted_to_metabaR_fish_mfzr.tsv",
#                                    file_samples = "~/vtam_benchmark/metafiles/samples_wo_control_fish.tsv")

#mblist <- tabfiles_to_metabarlist(file_reads = "out/obibar_fish/zfzr_obiout/reads_obi.tsv",
#                                    file_motus = "out/obibar_fish/zfzr_obiout/motus_obi_taxa2.tsv",
#                                    file_pcrs = "metafiles/pcrs_adapted_to_metabaR_fish_zfzr.tsv",
#                                    file_samples = "metafiles/samples_wo_control_fish.tsv")


mblist <- tabfiles_to_metabarlist(file_reads = file_reads,
                                    file_motus = file_motus,
                                    file_pcrs = file_pcrs,
                                    file_samples = file_samples)


#setwd("~/vtam_benchmark/obibar_fish/zfzr_metabarout/")

###  Compute the number of reads per pcr
mblist$pcrs$nb_reads <- rowSums(mblist$reads)
### Compute the number of motus per pcr
mblist$pcrs$nb_motus <- rowSums(mblist$reads>0)

############
# Flagging spurious signal
############
###
## Identifying extraction contaminants
###

# A contaminant should be preferentially amplified in negative controls since 
#there is no competing DNA. The function contaslayer relies on this assumption 
# and detects MOTUs whose relative abundance across the whole dataset is highest
#in negative controls.

mblist <- contaslayer(mblist, 
                        control_types = "extraction",
                        output_col = "not_an_extraction_conta")

# Compute relative abundance of all pcr contaminants together 
a <- data.frame(conta.relab = rowSums(mblist$reads[,!mblist$motus$not_an_extraction_conta]) / 
                  rowSums(mblist$reads))
# Add information on control types
a$control_type <- mblist$pcrs$control_type[match(rownames(a), rownames(mblist$pcrs))]

#png(file="out/obibar_fish/zfzr_metabarout/relative_abundance_external_contaminants.png")
#ggplot(a, aes(x=control_type, y=conta.relab, color=control_type)) + 
#  geom_boxplot() + geom_jitter(alpha=0.5) +
#  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
#  labs(x=NULL, y="Prop. Reads (log10)") + 
#  theme_bw() + 
#  scale_y_log10()
#dev.off()

# eliminate pcrs with > 10% of conta
mblist$pcrs$low_contamination_level <- 
  ifelse(a$conta.relab[match(rownames(mblist$pcrs), rownames(a))]>1e-1,  F, T)

###
## Flagging spurious or non-target MOTUs
###

#Flag MOTUs corresponding to target (TRUE) vs. non-target (FALSE) taxa 
mblist$motus$target_taxon <- grepl("Eukaryota", mblist$motus$domain)

###
# Identify MOTUs whose sequence is too dissimilar from references
# Flag not degraded (TRUE) vs. potentially degraded sequences (FALSE)
mblist$motus$not_degraded <-
  ifelse(mblist$motus$identity < 80, F, T)
# if no %id for taxassign, the value is NA => change this to FALSE
mblist$motus$not_degraded <- ifelse(is.na(mblist$motus$not_degraded), FALSE, mblist$motus$not_degraded)

###
## Detecting PCR outliers
###

mblist$pcrs$nb_reads <- rowSums(mblist$reads)
mblist$pcrs$nb_motus <- rowSums(mblist$reads>0)
#png(file="out/obibar_fish/zfzr_metabarout/sequencing_depth.png")
#ggplot(mblist$pcrs, aes(nb_reads)) +
#  geom_histogram(bins=40, color="grey", fill="white") + 
#  geom_vline(xintercept = 5e2, lty=2, color="orange") + # threshold
#  scale_x_log10() + 
#  labs(x="# Reads (with all OTUs and PCRs)", 
#       y="# PCRs") +
#  theme_bw() + 
#  theme(panel.grid = element_blank())
#dev.off()

# Flag pcrs with an acceptable sequencing depth (TRUE) or inacceptable one (FALSE)
mblist$pcrs$seqdepth_ok <- ifelse(mblist$pcrs$nb_reads < 5e2, F, T)

###
## A second way to evaluate the PCR quality is to assess their reproducibility.
###
# Subsetting the metabarlist to eliminate pcrs with 0 reads and negative controls
mblist_sub <- subset_metabarlist(mblist, 
                                   table = "pcrs", 
                                   indices = mblist$pcrs$nb_reads>0 & (
                                     is.na(mblist$pcrs$control_type) |
                                       mblist$pcrs$control_type=="positive"))
comp1 = pcr_within_between(mblist_sub)
#png(file="out/obibar_fish/zfzr_metabarout/pcr_slayer.png")
#check_pcr_thresh(comp1)
#dev.off()

mblist_sub <- pcrslayer(mblist_sub, output_col = "replicating_pcr", plot = F)

# Distinguish between pcrs obtained from samples from positive controls
#png(file="out/obibar_fish/zfzr_metabarout/pcr_distances.png")
#mds = check_pcr_repl(mblist_sub, 
#                     groups = mblist_sub$pcrs$type, 
#                     funcpcr = mblist_sub$pcrs$replicating_pcr)
#mds + labs(color = "pcr type") + scale_color_manual(values = c("cyan4", "gray"))
#dev.off()

# Now report the flagging in the initial metabarlist
mblist$pcrs$replicating_pcr <- TRUE
mblist$pcrs[rownames(mblist_sub$pcrs),"replicating_pcr"] <- mblist_sub$pcrs$replicating_pcr

###
# Lowering tag-jumps
###

# Define a vector of thresholds to test
thresholds <- c(0,1e-4,1e-3, 5e-3, 1e-2, 3e-2, 5e-2) 
# Run the tests and stores the results in a list
tests <- lapply(thresholds, function(x) tagjumpslayer(mblist,x))
names(tests) <- paste("t_", thresholds, sep="")

# Format the data for ggplot with amount of reads at each threshold
tmp <- melt(as.matrix(do.call("rbind", lapply(tests, function(x) rowSums(x$reads)))))
colnames(tmp) <- c("threshold", "sample", "abundance")

# Add richness in MOTUs at each threshold
tmp$richness <-
  melt(as.matrix(do.call("rbind", lapply(tests, function(x) {
    rowSums(x$reads > 0)
  }))))$value

# Add control type information on pcrs and make data curation threshold numeric
tmp$controls <- mblist$pcrs$control_type[match(tmp$sample, rownames(mblist$pcrs))]
tmp$threshold <- as.numeric(gsub("t_", "", tmp$threshold))

# New table formatting for ggplot
tmp2 <- melt(tmp, id.vars=colnames(tmp)[-grep("abundance|richness", colnames(tmp))])
#png(file="out/obibar_fish/zfzr_metabarout/tagjumpslayer.png")
#ggplot(tmp2, aes(x=as.factor(threshold), y=value)) + 
#  geom_boxplot(color="grey40") + 
#  geom_vline(xintercept = which(levels(as.factor(tmp2$threshold)) == "0.01"), col="orange", lty=2) + 
#  geom_jitter(aes(color=controls), width = 0.2, alpha=0.5) + 
#  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
#  facet_wrap(~variable+controls, scale="free_y", ncol=5) + 
#  theme_bw() + 
#  scale_y_log10() +
#  labs(x="MOTU pcr : total abundance filtering threshold", y="# Reads/MOTUs") + 
#  theme(panel.grid = element_blank(), 
#        strip.background = element_blank(), 
#        axis.text.x = element_text(angle=40, h=1), 
#        legend.position = "none")
#dev.off()


# Use tag-jump corrected metabarlist with the threshold identified above
# 0.005 for bat and fish
clean <- tests[["t_0.005"]]

###
# Eliminate motus and pcrs did not pass filtering
###
# eliminate motus that did not pass filtering
clean_womotus <- subset_metabarlist(clean, "motus", 
                          indices = rowSums(clean$motus[,c("not_an_extraction_conta", "target_taxon",
                                                         "not_degraded")]) == 3)
# eliminate pcrs that did not pass filtering
clean_womotus_wopcrs <- subset_metabarlist(clean_womotus, "pcrs", 
                           indices = rowSums(clean_womotus$pcrs[,c("low_contamination_level", 
                                                            "seqdepth_ok", "replicating_pcr")]) == 3)
summary_metabarlist(clean_womotus_wopcrs)

###
# print out asv table
###
reads_clean <- t(extract_table(clean_womotus_wopcrs, "reads"))
#write.csv2(reads_clean, file="out/obibar_fish/zfzr_metabarout/obibar_base.csv")
write.csv2(reads_clean, file=obibar_base)

