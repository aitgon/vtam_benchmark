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
if (length(args)!=9) {
  stop("9 arguments must be supplied", call.=FALSE)
}

library(vegan)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(readxl)
library(cowplot)


if(0){
#setwd("/home/meglecz/vtam_benchmark_local/092022_bxplt")
wd = "/home/meglecz/vtam_benchmark_local/092022_bxplt"
species_bat_path = file.path(wd, 'bat_sp.tsv')
dcb_path = file.path(wd, 'dalu_bat_final_taxa.tsv')
dcf_path = file.path(wd, 'dalu_fish_final_taxa.tsv')
ocb_path = file.path(wd, 'obibar_bat_final_taxa.tsv')
ocf_path = file.path(wd, 'obibar_fish_final_taxa.tsv')
vtamb_path = file.path(wd, 'vtam_bat_final_taxa.tsv')
vtamf_path = file.path(wd, 'vtam_fish_final_taxa.tsv')
bxplt_asvrichness_path = "bxplt_asvrichness.png"
bxplt_betadiversity_path = "bxplt_betadiversity.png"
} else {
species_bat_path = args[1]
dcb_path = args[2]
dcf_path = args[3]
ocb_path = args[4]
ocf_path = args[5]
vtamb_path = args[6]
vtamf_path = args[7]
bxplt_asvrichness_path = args[8]
bxplt_betadiversity_path = args[9]
}


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


#dcb<-read.table(file = "dalu_bat_final_taxa.tsv", sep = '\t', header = TRUE)
#dcb<-dcb[, c(1,14:370)]
#dcf<-read.table(file = "dalu_fish_final_taxa.tsv", sep = '\t', header = TRUE)
#dcf<-dcf[, c(1,2:36)]

#ocb<-read.table(file = "obibar_bat_final_taxa.tsv", sep = '\t', header = TRUE)
#ocb<-ocb[, c(1,14:370)]
#ocf<-read.table(file = "obibar_fish_final_taxa.tsv", sep = '\t', header = TRUE)
#ocf<-ocf[, c(1,2:36)]

#vtamb<-read.table(file = "vtam_bat_final_taxa.tsv", sep = '\t', header = TRUE)
#vvb<-vtamb[,c(3,17:373)]
#vtamf<-read.table(file = "vtam_fish_final_taxa.tsv", sep = '\t', header = TRUE)
#vvf<-vtamf[,c(1,13:47)]

species<-read.csv(file = species_bat_path, sep= '\t', header = TRUE)

# format df + calcul diversity alpha (sum of lines)
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
label_font_size = 11
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

my_comparisons <- list( c("VTAM", "DALU"), c("VTAM", "OBIbaR"), c("DALU", "OBIbaR"))
p = p + stat_compare_means(comparisons = my_comparisons, label.y = c(60, 75, 105, 85), label = "p.signif", size=6)

#p
#dir.create("out", showWarnings = F)
ggsave(bxplt_asvrichness_path, width = 12, height = 12, units = "cm")
#####################################################################################
################betadiversity########################################################
#####################################################################################
# format df
#fonction creation
Beta<-function(x,v){
  y<-x[-1]
  rownames(y)<-x[,1]
  y<-t(y)
  y[y>0]<-1
  y<-data.frame(y)
  m<-as.matrix(vegdist(y, method="jaccard"))
  m2 <- melt(m)[melt(upper.tri(m))$value,]
  names(m2) <- c("c1", "c2", "distance")
  B<-m2
  B$Condition<-as.character(v)
  B$Filtering<-substring(as.character(v), first=1, last=1)
  B$MaxMin<-substring(as.character(v), first=2, last=2)
  
  B$Dataset<- substring(as.character(v), first=3, last=3)
  B<-select(B,"c1","c2","distance","Condition","Dataset","Filtering","MaxMin")
  B
}


###Calcul Beta diversity
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


##Dimension intra vs inter secies 
BetaDF<-merge(BetaDF,species,by.x="c1",by.y="Sample.ID",all.x=TRUE)
BetaDF<-merge(BetaDF,species,by.x="c2",by.y="Sample.ID",all.x=TRUE)


BetaDF<-BetaDF %>% 
  mutate(comparison = if_else(
    Host_species.x == Host_species.y, "Intra-Species", "Inter-Species"))
  
BetaDF$comparison[BetaDF$Dataset=="Fish"]<-"Intra-Species"

BetaDF<-subset(BetaDF,BetaDF$Host_species.x!="NU" & BetaDF$Host_species.y!="NU"  |  (is.na(BetaDF$Host_species.x) & is.na(BetaDF$Host_species.y)))  #####removing pairwise comparison with at least one NU specimens

BetaDF$FiltComp<-paste(BetaDF$Pipeline,BetaDF$comparison,sep="_")
BetaDF$Panel<-paste(BetaDF$Dataset,BetaDF$comparison,sep=" ")

#-------------------------------- boxplot
BetaDF$Pipeline=factor(BetaDF$Pipeline, levels=c("VTAM", "DALU", "OBIbaR"))
BetaDF$Panel=factor(BetaDF$Panel, levels=c("Bat Inter-Species", "Bat Intra-Species", "Fish Intra-Species"))
title = "Comparison of Beta diversities"
ylab = "Beta diversity"
xlab = "Pipelines"
label_font_size = 11
ylim_max = 1.25
p = ggplot(BetaDF, aes(x=Pipeline, y=distance, fill=Pipeline))

p = p + geom_boxplot()
p = p + facet_grid(. ~ Panel)

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

my_comparisons <- list( c("VTAM", "DALU"), c("VTAM", "OBIbaR"), c("DALU", "OBIbaR"))

p = p + stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(1.0, 1.05, 1.1, 1.15), label = "p.signif", size=6)
#p

#dir.create("out", showWarnings = F)
ggsave(bxplt_betadiversity_path, width = 12, height = 12, units = "cm")

if(0){
##### Pairwise tests 28/09/2022
#Intra-Species Bat
DistIntrBats<-subset(BetaDF,BetaDF$Dataset=="Bat" & BetaDF$comparison=="Intra-Species")
kruskal.test(DistIntrBats$distance~DistIntrBats$Pipeline)
SubDistIntrBats<-subset(DistIntrBats,DistIntrBats$Pipeline=="VTAM" | DistIntrBats$Pipeline=="OBIbaR")   ##select the pipelines to be compared
kruskal.test(SubDistIntrBats$distance~SubDistIntrBats$Pipeline)
aggregate(distance~Pipeline,data=DistIntrBats,FUN=mean)   # mean distances for each group

#Inter-Species Bat
DistInterBats<-subset(BetaDF,BetaDF$Dataset=="Bat" & BetaDF$comparison=="Inter-Species")
kruskal.test(DistInterBats$distance~DistInterBats$Pipeline)
SubDistInterBats<-subset(DistInterBats,DistInterBats$Pipeline=="VTAM" | DistInterBats$Pipeline=="OBIbaR") ##select the pipelines to be compared
kruskal.test(SubDistInterBats$distance~SubDistInterBats$Pipeline)

#Intra-Species Fish
DistIntrFish<-subset(BetaDF,BetaDF$Dataset=="Fish" & BetaDF$comparison=="Intra-Species")
kruskal.test(DistIntrFish$distance~DistIntrFish$Pipeline)
SubDistIntrFish<-subset(DistIntrFish,DistIntrFish$Pipeline=="VTAM" | DistIntrFish$Pipeline=="OBIbaR") ##select the pipelines to be compared
kruskal.test(SubDistIntrFish$distance~SubDistIntrFish$Pipeline)

###Calcul de mean and SD for each group
cbind(
  data.frame(MeanIntraBats=aggregate(distance~Pipeline,data=DistIntrBats,FUN=mean)$distance),   # mean distances for each group
  data.frame(cvIntraBats=aggregate(distance~Pipeline,data=DistIntrBats,FUN=sd)$distance),   # sd and distances for each group
  data.frame(MeanInterBats=aggregate(distance~Pipeline,data=DistInterBats,FUN=mean)$distance),   
  data.frame(cvInterBats=aggregate(distance~Pipeline,data=DistInterBats,FUN=sd)$distance),
  data.frame(MeanIntraFish=aggregate(distance~Pipeline,data=DistIntrFish,FUN=mean) $distance), 
  data.frame(cvIntraFish=aggregate(distance~Pipeline,data=DistIntrFish,FUN=sd) $distance) ) 
}
