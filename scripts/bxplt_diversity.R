#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#install.packages("vegan")
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("reshape2")
#install.packages("harrietr")
#install.packages("ggpubr")
#install.packages("readxl")
#install.packages("cowplot")

print(length(args))
if (length(args)!=8) {
  stop("8 arguments must be supplied", call.=FALSE)
}


library(vegan)
library(ggplot2)
library(dplyr)
library(reshape2)
#library(harrietr)
library(ggpubr)
library(readxl)
library(cowplot)

#wd="out"
#dcb_path = file.path(wd, 'dalu_bat/dalu_bat_final_taxa.tsv')
#dcf_path = file.path(wd, 'dalu_fish/dalu_fish_final_taxa.tsv')
#ocb_path = file.path(wd, 'obibar_bat/obibar_bat_final_taxa.tsv')
#ocf_path = file.path(wd, 'obibar_fish/obibar_fish_final_taxa.tsv')
#vtamb_path = file.path(wd, 'vtam_bat/vtam_bat_final_taxa.tsv')
#vtamf_path = file.path(wd, 'vtam_fish/vtam_fish_final_taxa.tsv')

dcb_path = args[1]
dcf_path = args[2]
ocb_path = args[3]
ocf_path = args[4]
vtamb_path = args[5]
vtamf_path = args[6]

#bxplt_asvrichness_path = file.path(wd, "bxplt_asvrichness.png")
#bxplt_betadiversity_path = file.path(wd, "bxplt_betadiversity.png")

bxplt_asvrichness_path = args[7]
bxplt_betadiversity_path = args[8]

#####2022/07/14
##Calcul of beta diversity of different pipelines and comparison to of DALu and OBIbaR to VTAM


dcb<-read.table(file = dcb_path, sep = '\t', header = TRUE)
dcb<-dcb[, c(1,14:370)]
dcf<-read.table(file = dcf_path, sep = '\t', header = TRUE)
dcf<-dcf[, c(1,2:36)]


ocb<-read.table(file = ocb_path, sep = '\t', header = TRUE)
ocb<-ocb[, c(1,14:370)]
ocf<-read.table(file = ocf_path, sep = '\t', header = TRUE)
ocf<-ocf[, c(1,2:36)]


vtamb<-read.table(file = vtamb_path, sep = '\t', header = TRUE)
vvb<-vtamb[,c(3,17:373)]
vtamf<-read.table(file = vtamf_path, sep = '\t', header = TRUE)
vvf<-vtamf[,c(1,13:47)]

#setwd("~/vtam_benchmark/comparaison")

#remise en t forme des df + calcul de l'alpha diversit? des ASVs par somme des lignes
forme<-function(x,v){
  y<-x[-1]
  rownames(y)<-x[,1]
  y<-t(y)
  y[y>0]<-1
  y<-data.frame(y)
  y$Alpha<-apply(y,1,sum)
  y$Condition<-as.character(v)
  y$Pipeline<-substring(as.character(v), first=1, last=1)
  y$MaxMin<-substring(as.character(v), first=2, last=2)
  y$Dataset<- substring(as.character(v), first=3, last=3)
  y<-select(y,"Alpha","Condition","Dataset","Pipeline","MaxMin")
  y
  }


#####################################################################################
################alphadiversity########################################################
#####################################################################################
AlphaDF<-rbind(
  forme(dcb,"dcb"),
  forme(dcf,"dcf"),
  
  forme(ocb,"ocb"),
  forme(ocf,"ocf"),
  
  forme(vvb,"vvb"),
  forme(vvf,"vvf")
)

AlphaDF$Pipeline<-paste(AlphaDF$Pipeline, AlphaDF$MaxMin,"")
AlphaDF$Pipeline[AlphaDF$Pipeline=="d c "]<-"DALU"
AlphaDF$Pipeline[AlphaDF$Pipeline=="o c "]<-"OBIbaR"
AlphaDF$Pipeline[AlphaDF$Pipeline=="v v "]<-"VTAM"

AlphaDF$Dataset[AlphaDF$Dataset=="b"]<-"Bat"
AlphaDF$Dataset[AlphaDF$Dataset=="f"]<-"Fish"

#-------------------------------- boxplot
AlphaDF$Pipeline=factor(AlphaDF$Pipeline, levels=c("VTAM", "DALU", "OBIbaR"))
title = "Comparison of ASV richnesses"
ylab = "ASV richness"
label_font_size = 12
ylim_max = 120
stat_label_y = 110
p = ggplot(AlphaDF, aes(x=Pipeline, y=Alpha, fill=Pipeline))
p = p + geom_boxplot()
p = p + facet_grid(. ~ Dataset)
p = p + ggtitle(title) # for the main title
p = p + theme_light()
# p = p + theme_linedraw()

p = p + theme(axis.text.x = element_text(size = label_font_size, angle = 45, hjust = 1))
p = p + theme(axis.title.x = element_blank())
p = p + theme(axis.title.y = element_text(size = label_font_size))
p = p + theme(legend.position = "None")
p = p + theme(panel.grid.major.x = element_blank())
p = p + theme(plot.title = element_text(size=label_font_size))
p = p + theme(plot.title = element_text(size = label_font_size, hjust = 0.5))
p = p + theme(strip.text.x = element_text(size = label_font_size))
p = p + ylab(ylab)
p = p + ylim(0, ylim_max)
p = p + theme(strip.background = element_rect(fill = "black"))
p = p + theme(strip.text = element_text(colour = 'white'))

my_comparisons <- list( c("VTAM", "DALU"), c("VTAM", "OBIbaR"))
p = p + stat_compare_means(comparisons = my_comparisons, label.y = c(60, 75, 105, 85), label = "p.signif", size=6)

#p
#dir.create("out", showWarnings = F)
ggsave(bxplt_asvrichness_path, width = 12, height = 12, units = "cm")

#####################################################################################
################betadiversity########################################################
#####################################################################################
#repartir des fichiers gener?s plus haut

#fonction creation
Beta<-function(x,v){
  y<-x[-1]
  rownames(y)<-x[,1]
  y<-t(y)
  y[y>0]<-1
  y<-data.frame(y)
  B<-data.frame(Beta=c(vegdist(y, method="jaccard")))
  B$Condition<-as.character(v)
  B$Filtering<-substring(as.character(v), first=1, last=1)
  B$MaxMin<-substring(as.character(v), first=2, last=2)
  
  B$Dataset<- substring(as.character(v), first=3, last=3)
  B<-select(B,"Beta","Condition","Dataset","Filtering","MaxMin")
  B
}


Beta(dcb,"dcb")

###Calcul Beta diversit?
BetaDF<-rbind(
  Beta(dcb,"dcb"),
  Beta(dcf,"dcf"),
  
  Beta(ocb,"ocb"),
  Beta(ocf,"ocf"),
  
  Beta(vvb,"vvb"),
  Beta(vvf,"vvf")
)

BetaDF$Pipeline<-paste(BetaDF$Filtering, BetaDF$MaxMin,"")
BetaDF$Pipeline[BetaDF$Pipeline=="d c "]<-"DALU"
BetaDF$Pipeline[BetaDF$Pipeline=="o c "]<-"OBIbaR"
BetaDF$Pipeline[BetaDF$Pipeline=="v v "]<-"VTAM"

BetaDF$Dataset[BetaDF$Dataset=="b"]<-"Bat"
BetaDF$Dataset[BetaDF$Dataset=="f"]<-"Fish"

#-------------------------------- boxplot
BetaDF$Pipeline=factor(BetaDF$Pipeline, levels=c("VTAM", "DALU", "OBIbaR"))
title = "Comparison of Beta diversities"
ylab = "Beta diversity"
xlab = "Pipelines"
label_font_size = 12
ylim_max = 1.25
p = ggplot(BetaDF, aes(x=Pipeline, y=Beta, fill=Pipeline))

p = p + geom_boxplot()
p = p + facet_grid(. ~ Dataset)
p = p + ggtitle(title) # for the main title
p = p + theme_light()
# p = p + theme_linedraw()

p = p + theme(axis.text.x = element_text(size = label_font_size, angle = 45, hjust = 1))
p = p + theme(axis.title.x = element_blank())
p = p + theme(axis.title.y = element_text(size = label_font_size))
p = p + theme(legend.position = "None")
p = p + theme(panel.grid.major.x = element_blank())
p = p + theme(plot.title = element_text(size=label_font_size))
p = p + theme(plot.title = element_text(size = label_font_size, hjust = 0.5))
p = p + theme(strip.text.x = element_text(size = label_font_size))
p = p + ylab(ylab)
p = p + ylim(0, ylim_max)
p = p + theme(strip.background = element_rect(fill = "black"))
p = p + theme(strip.text = element_text(colour = 'white'))

my_comparisons <- list( c("VTAM", "DALU"), c("VTAM", "OBIbaR"))
p = p + stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(1.0, 1.05, 1.1, 1.15), label = "p.signif", size=6)
#p

#dir.create("out", showWarnings = F)
ggsave(bxplt_betadiversity_path, width = 12, height = 12, units = "cm")

