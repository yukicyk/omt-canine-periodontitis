# Canine MiSeq Analysis 12 June 2018
# 20181011 DeSeq2

setwd("~/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/4/")

#import data
load("~/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/Canine_ps.RData")
library(phyloseq)
library(dplyr)
# ps    original dada2 output dataset
# ps1   removed sample MT-E-V-2
# ps0   removed unassigned phyla
# ps2   filtered low prevalence taxa (sum of prevalence <10)
# ps3   filtered prevalence at threshold at 5% of total sample number
sample_data(ps1)$Treatment<- gsub("Recipeint", "Recipient", sample_data(ps1)$Treatment) #saved to R.data
rank_names(ps1)
# [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
# Create table, number of features for each phyla
ps1_spname <-factor(paste0(tax_table(ps1)[, "Genus"]," ", tax_table(ps1)[,"Species"]))
ps1_spname <-gsub(" NA", " sp.", ps1_spname)
ps1_spname <-gsub("NA sp.", "Unknown sp.", ps1_spname)
tax_table(ps1)<-cbind(tax_table(ps1), "Species2"=ps1_spname)
write.table(table(tax_table(ps1)[,"Species2"], exclude = NULL),"Species_ASV.txt", sep = "\t")
# Create table, number of ASV for each phyla in each Treatment_Timepoint
TT<-factor(paste0(sample_data(ps1)$Treatment,"_",sample_data(ps1)$Timepoint))
sample_data(ps1)<-cbind(sample_data(ps1),"TT"=TT)
levels(TT)
#[1] "Control_T1"   "Control_T2"   "Control_T3"   "Control_T4"   "Donor_T1"     "Donor_T2"     "Recipient_T1" "Recipient_T2"
#[9] "Recipient_T3" "Recipient_T4"

subset_samples(ps1, TT=="Control_T1") %>% table(tax_table(.)[,"Phylum"], exclude = NULL)


ps1g <- tax_glom(ps1, "Genus")
ntaxa(ps1g) #220
ps1c <- tax_glom(ps1, "Class")
ntaxa(ps1c) #32
tax_table(ps1c)
ps1p <- tax_glom(ps1, "Phylum")

#phyloseq heatmap-ASV
library(scales)
library(vegan)
order_of_samples <- c("MT-B-1","MT-D-1","MT-G-1","MT-I-1","MT-K-1","MT-M-1","MT-O-1","MT-R-1","MT-S-1","MT-B-2","MT-D-2","MT-G-2","MT-I-2","MT-K-2","MT-M-2","MT-O-2","MT-R-2","MT-S-2","MT-B-3","MT-D-3","MT-G-3","MT-I-3","MT-K-3","MT-M-3","MT-O-3","MT-R-3","MT-S-3","MT-B-4","MT-D-4","MT-G-4","MT-I-4","MT-K-4","MT-M-4","MT-O-4","MT-R-4","MT-S-4","Transplant-A","Transplant-C","Transplant-F","Transplant-H","Transplant-J","Transplant-L","Transplant-N","Transplant-P","Transplant-Q","MT-E-I-1","MT-E-II-1","MT-E-III-1","MT-E-IV-1","MT-E-V-1","MT-A-1","MT-C-1","MT-F-1","MT-H-1","MT-J-1","MT-L-1","MT-N-1","MT-P-1","MT-Q-1","MT-A-2","MT-C-2","MT-F-2","MT-H-2","MT-J-2","MT-L-2","MT-N-2","MT-P-2","MT-Q-2","MT-A-3","MT-C-3","MT-F-3","MT-H-3","MT-J-3","MT-L-3","MT-N-3","MT-P-3","MT-Q-3","MT-A-4","MT-C-4","MT-F-4","MT-H-4","MT-J-4","MT-L-4","MT-N-4","MT-P-4","MT-Q-4")
genus_heatmap <- plot_heatmap(ps1g, method ="NMDS", distance ="bray",na.value = "black",
                              trans = log_trans(10), sample.order = order_of_samples, 
                              title="Heatmap showing log10-transformed abundances of microbiome at the Genus rank", 
                              taxa.label="Genus")
class_heatmap <- plot_heatmap(ps1c, method ="NMDS", distance ="bray",na.value = "black",
                              trans = log_trans(10), sample.order = order_of_samples,
                              title="Heatmap showing log10-transformed abundances of microbiome at the Class rank", 
                              taxa.label="Class")
phylum_heatmap <- plot_heatmap(ps1p, method ="NMDS", distance ="bray",na.value = "black",
                              trans = log_trans(10), sample.order = order_of_samples,
                              title="Heatmap showing log10-transformed abundances of microbiome at the Phylum rank", 
                              taxa.label="Phylum")

#DESEQ2
library(DESeq2)
#function to calculate geometric means before size factor estimation
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# remove Donor T1 samples with ambiguous dog
ps01 <-subset_samples(ps1, ((!Dog =="DogA,C") & (!Dog =="DogF,H") & (!Dog=="DogJ,L")))
#remove donor samples from the dataset
ps02 <- subset_samples(ps1, !Treatment =="Donor") #saved to R.data

#simply put everything together
levels(sample_data(ps02)$Treatment)
#[1] "Control"   "Recipient"
levels(sample_data(ps02)$Timepoint)
#[1] "T1" "T2" "T3" "T4"
dds=phyloseq_to_deseq2(ps02, ~Treatment+Timepoint+Treatment:Timepoint)
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)
dds = DESeq(dds, fitType="local")
resultsNames(dds)
#[1] "Intercept"                      "Treatment_Recipeint_vs_Control" "Timepoint_T2_vs_T1"            
#[4] "Timepoint_T3_vs_T1"             "Timepoint_T4_vs_T1"             "TreatmentRecipeint.TimepointT2"
#[7] "TreatmentRecipeint.TimepointT3" "TreatmentRecipeint.TimepointT4"
plotDispEsts(dds)
library("IHW")
resT2T1<- results(dds, name="Timepoint_T2_vs_T1", filterFun=ihw, alpha=0.05)
summary(resT2T1) #6,8
resT3T2<- results(dds, contrast=c("Timepoint","T3","T2"), filterFun=ihw, alpha=0.05)
summary(resT3T2) #2,3
resT4T2<- results(dds, contrast=c("Timepoint","T4","T2"), filterFun=ihw, alpha=0.05)
summary(resT4T2) #6,6
resT4T3<- results(dds, contrast=c("Timepoint","T4","T3"), filterFun=ihw, alpha=0.05)
summary(resT4T3) #10,6

#DESEQ2
#20181011
library(phyloseq)
library(plyr)
library(dplyr)
library("IHW")
library(DESeq2)
#function to calculate geometric means before size factor estimation
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
#merge the treatment and timepoint and create a new group variable
sample_data(ps02)$group<- factor(paste0(sample_data(ps02)$Treatment,sample_data(ps02)$Timepoint))
levels(sample_data(ps02)$group)
#[1] "ControlT1"   "ControlT2"   "ControlT3"   "ControlT4"   "RecipeintT1" "RecipeintT2" "RecipeintT3" "RecipeintT4"
dds2=phyloseq_to_deseq2(ps02, ~0+group)
set.seed(20181011)
geoMeans2 = apply(counts(dds2), 1, gm_mean)
dds2 = estimateSizeFactors(dds2, geoMeans = geoMeans2)
dds2 = DESeq(dds2, fitType="local")
plotDispEsts(dds2)
resultsNames(dds2)
#[1] "groupControlT1"   "groupControlT2"   "groupControlT3"   "groupControlT4"   "groupRecipeintT1" "groupRecipeintT2"
#[7] "groupRecipeintT3" "groupRecipeintT4"
resRT3T2<- results(dds2, contrast=list("groupRecipientT2","groupRecipientT3"), filterFun=ihw, alpha=0.05)
summary(resRT3T2) #8,2
resRT4T3<- results(dds2, contrast=list("groupRecipientT4","groupRecipientT3"), filterFun=ihw, alpha=0.05)
summary(resRT4T3) #22,4
resRT4T2<- results(dds2, contrast=list("groupRecipientT4","groupRecipientT2"), filterFun=ihw, alpha=0.05)
summary(resRT4T2) #10,1
# same applied to control as well
resCT3T2<- results(dds2, contrast=list("groupControlT2","groupControlT3"), filterFun=ihw, alpha=0.05)
summary(resCT3T2) #4,2
resCT4T3<- results(dds2, contrast=list("groupControlT4","groupControlT3"), filterFun=ihw, alpha=0.05)
summary(resCT4T3) #11,6
resCT4T2<- results(dds2, contrast=list("groupControlT4","groupControlT2"), filterFun=ihw, alpha=0.05)
summary(resCT4T2) #7,5

#formatting the result table
resRT3T2 = resRT3T2[order(resRT3T2$log2FoldChange, na.last=NA),]
sigtabRT3T2 = resRT3T2[(resRT3T2$padj < 0.05),]
sigtabRT3T2 = cbind(as(sigtabRT3T2, "data.frame"), as(tax_table(ps02)[rownames(sigtabRT3T2), ], "matrix"))
head(sigtabRT3T2)
#replace NA species to sp.
sigtabRT3T2$Species = factor(sigtabRT3T2$Species, levels=c(levels(sigtabRT3T2$Species), "sp."))
#convert all NA's to sp.
sigtabRT3T2$Species[is.na(sigtabRT3T2$Species)] = "sp."
#create new taxa name
sigtabRT3T2$taxa <- factor(paste0(sigtabRT3T2$Genus," ", sigtabRT3T2$Species))
sigtabRT3T2$taxa <- gsub("_", " ",sigtabRT3T2$taxa)
write.table(sigtabRT3T2, "sigtabRT3T2.txt", sep="\t")
sigtabRT3T2<-rownames_to_column(sigtabRT3T2)


resRT4T3 = resRT4T3[order(resRT4T3$log2FoldChange, na.last=NA),]
sigtabRT4T3 = resRT4T3[(resRT4T3$padj < 0.05),]
sigtabRT4T3 = cbind(as(sigtabRT4T3, "data.frame"), as(tax_table(ps02)[rownames(sigtabRT4T3), ], "matrix"))
head(sigtabRT4T3)
sigtabRT4T3$Species = factor(sigtabRT4T3$Species, levels=c(levels(sigtabRT4T3$Species), "sp."))
sigtabRT4T3$Species[is.na(sigtabRT4T3$Species)] = "sp."
sigtabRT4T3$taxa <- factor(paste0(sigtabRT4T3$Genus," ", sigtabRT4T3$Species))
sigtabRT4T3$taxa<- gsub("_", " ",sigtabRT4T3$taxa)
write.table(sigtabRT4T3, "sigtabRT4T3.txt", sep="\t")
sigtabRT4T3<-rownames_to_column(sigtabRT4T3)

resRT4T2 = resRT4T2[order(resRT4T2$log2FoldChange, na.last=NA),]
sigtabRT4T2 = resRT4T2[(resRT4T2$padj < 0.05),]
sigtabRT4T2 = cbind(as(sigtabRT4T2, "data.frame"), as(tax_table(ps02)[rownames(sigtabRT4T2), ], "matrix"))
head(sigtabRT4T2)
sigtabRT4T2$Species = factor(sigtabRT4T2$Species, levels=c(levels(sigtabRT4T2$Species), "sp."))
sigtabRT4T2$Species[is.na(sigtabRT4T2$Species)] = "sp."
sigtabRT4T2$taxa <- factor(paste0(sigtabRT4T2$Genus," ", sigtabRT4T2$Species))
revalue(sigtabRT4T2$taxa, c("NA sp." = "Candidatus Endomicrobium sp.")) -> sigtabRT4T2$taxa
sigtabRT4T2$taxa <- gsub("_", " ", sigtabRT4T2$taxa)
write.table(sigtabRT4T2, "sigtabRT4T2.txt", sep="\t")
sigtabRT4T2<-rownames_to_column(sigtabRT4T2)

#Formatting table for Control sets
resCT3T2 = resCT3T2[order(resCT3T2$log2FoldChange, na.last=NA),]
sigtabCT3T2 = resCT3T2[(resCT3T2$padj < 0.05),]
sigtabCT3T2 = cbind(as(sigtabCT3T2, "data.frame"), as(tax_table(ps02)[rownames(sigtabCT3T2), ], "matrix"))
head(sigtabCT3T2)
sigtabCT3T2$Species = factor(sigtabCT3T2$Species, levels=c(levels(sigtabCT3T2$Species), "sp."))
sigtabCT3T2$Species[is.na(sigtabCT3T2$Species)] = "sp."
sigtabCT3T2$taxa <- factor(paste0(sigtabCT3T2$Genus," ", sigtabCT3T2$Species))
sigtabCT3T2$taxa <- gsub("_", " ", sigtabCT3T2$taxa)
write.table(sigtabCT3T2, "sigtabCT3T2.txt", sep="\t")
sigtabCT3T2<-rownames_to_column(sigtabCT3T2)

resCT4T3 = resCT4T3[order(resCT4T3$log2FoldChange, na.last=NA),]
sigtabCT4T3 = resCT4T3[(resCT4T3$padj < 0.05),]
sigtabCT4T3 = cbind(as(sigtabCT4T3, "data.frame"), as(tax_table(ps02)[rownames(sigtabCT4T3), ], "matrix"))
head(sigtabCT4T3)
sigtabCT4T3$Species = factor(sigtabCT4T3$Species, levels=c(levels(sigtabCT4T3$Species), "sp."))
sigtabCT4T3$Species[is.na(sigtabCT4T3$Species)] = "sp."
sigtabCT4T3$taxa <- factor(paste0(sigtabCT4T3$Genus," ", sigtabCT4T3$Species))
revalue(sigtabCT4T3$taxa, c("NA sp." = "Candidatus Endomicrobium sp.")) -> sigtabCT4T3$taxa
write.table(sigtabCT4T3, "sigtabCT4T3.txt", sep="\t")
sigtabCT4T3$taxa <- gsub("_", " ", sigtabCT4T3$taxa)
sigtabCT4T3<-rownames_to_column(sigtabCT4T3)

resCT4T2 = resCT4T2[order(resCT4T2$log2FoldChange, na.last=NA),]
sigtabCT4T2 = resCT4T2[(resCT4T2$padj < 0.05),]
sigtabCT4T2 = cbind(as(sigtabCT4T2, "data.frame"), as(tax_table(ps02)[rownames(sigtabCT4T2), ], "matrix"))
head(sigtabCT4T2)
sigtabCT4T2$Species = factor(sigtabCT4T2$Species, levels=c(levels(sigtabCT4T2$Species), "sp."))
sigtabCT4T2$Species[is.na(sigtabCT4T2$Species)] = "sp."
sigtabCT4T2$taxa <- factor(paste0(sigtabCT4T2$Genus," ", sigtabCT4T2$Species))
revalue(sigtabCT4T2$taxa, c("NA sp." = "Candidatus Endomicrobium sp.")) -> sigtabCT4T2$taxa
sigtabCT4T2$taxa <- gsub("_", " ", sigtabCT4T2$taxa)
sigtabCT4T2<-rownames_to_column(sigtabCT4T2)
write.table(sigtabCT4T2, "sigtabCT4T2.txt", sep="\t")


library(tidyverse)
library("ggplot2")
library(ggrepel)

theme_set(theme_bw())
RT3T2<- ggplot(sigtabRT3T2 %>%
                 rownames_to_column(), aes(y=taxa, x=log2FoldChange, color=Phylum, label = rowname)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.3) +
  geom_point(size=6) + 
  geom_text_repel(color="black") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  ggtitle("Recipient Timepoint 2 vs 3")

RT4T3<- ggplot(sigtabRT4T3 %>%
                 rownames_to_column(), aes(y=taxa, x=log2FoldChange, color=Phylum, label = rowname)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.3) +
  geom_point(size=6) + 
  geom_text_repel(color="black") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  ggtitle("Recipient Timepoint 4 vs 3")
RT4T2<- ggplot(sigtabRT4T2 %>%
                 rownames_to_column(), aes(y=taxa, x=log2FoldChange, color=Phylum, label = rowname)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.3) +
  geom_point(size=6) + 
  geom_text_repel(color="black") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  ggtitle("Recipient Timepoint 4 vs 2")

#horizontal bar plots
#creat custom color scale for phyla
library(RColorBrewer)
myColors <- brewer.pal(11,"Set3")
names(myColors) <- c("Firmicutes","Proteobacteria","Fusobacteria","Bacteroidetes","Actinobacteria","Synergistetes","Saccharibacteria(TM7)","Spirochaetes","Chlorobi", "Elusimicrobia","Gracilibacteria(GN02)")
colScale <- scale_colour_manual(name = "phylum",values = myColors)

theme_set(theme_bw(base_size=14))
  ggplot(sigtabCT3T2, aes(x=reorder(rowname, -log2FoldChange), y=log2FoldChange, fill=Phylum, label=taxa)) + 
  geom_hline(yintercept = 0.0, color = "black", size = 0.3) +
  geom_bar(stat="identity") + coord_flip() +
  geom_text(color="black", aes(y = ifelse(log2FoldChange > 0, -0.5, 0.5), hjust = ifelse(log2FoldChange > 0, 1, 0)))+
  scale_x_discrete(name="Taxa", breaks=NULL) +
  colScale + 
  scale_fill_manual(name = "Phylum", values = myColors)+
  labs(title = "Control Timepoint 2 vs 3") +
    geom_text(aes(label = rowname),position = position_stack(0.8),hjust = ifelse(sigtabCT3T2$log2FoldChange > 0, 0, 1))

theme_set(theme_bw(base_size=14))
  ggplot(sigtabCT4T2, aes(x=reorder(rowname, -log2FoldChange), y=log2FoldChange, fill=Phylum, label=taxa)) + 
  geom_hline(yintercept = 0.0, color = "black", size = 0.3) +
  geom_bar(stat="identity") + coord_flip() +
  geom_text(color="black", aes(y = ifelse(log2FoldChange > 0, -0.5, 0.5), hjust = ifelse(log2FoldChange > 0, 1, 0)))+
  scale_x_discrete(name="Taxa", breaks=NULL) +
  colScale + 
  scale_fill_manual(name = "Phylum", values = myColors)+
  labs(title = "Control Timepoint 2 vs 4") +
    geom_text(aes(label = rowname),position = position_stack(0.8),hjust = ifelse(sigtabCT4T2$log2FoldChange > 0, 0, 1))
  
theme_set(theme_bw(base_size=14))
  ggplot(sigtabCT4T3 , aes(x=reorder(rowname, -log2FoldChange), y=log2FoldChange, fill=Phylum, label=taxa)) + 
  geom_hline(yintercept = 0.0, color = "black", size = 0.3) +
  geom_bar(stat="identity") + coord_flip() +
  geom_text(color="black", aes(y = ifelse(log2FoldChange > 0, -0.5, 0.5), hjust = ifelse(log2FoldChange > 0, 1, 0)))+
  scale_x_discrete(name="Taxa", breaks=NULL) +
  colScale + 
  scale_fill_manual(name = "Phylum", values = myColors)+
  labs(title = "Control Timepoint 3 vs 4") +
    geom_text(aes(label = rowname),position = position_stack(0.8),hjust = ifelse(sigtabCT4T3$log2FoldChange > 0, 0, 1))

theme_set(theme_bw(base_size=14))
  ggplot(sigtabRT3T2, aes(x=reorder(rowname, -log2FoldChange), y=log2FoldChange, fill=Phylum, label=taxa)) + 
  geom_hline(yintercept = 0.0, color = "black", size = 0.3) +
  geom_bar(stat="identity") + coord_flip() +
  geom_text(color="black", aes(y = ifelse(log2FoldChange > 0, -0.5, 0.5), hjust = ifelse(log2FoldChange > 0, 1, 0)))+
  scale_x_discrete(name="Taxa", breaks=NULL) +
  colScale + 
  scale_fill_manual(name = "Phylum", values = myColors)+
  labs(title = "Recipient Timepoint 2 vs 3") +
    geom_text(aes(label = rowname),position = position_stack(0.8),hjust = ifelse(sigtabRT3T2$log2FoldChange > 0, 0, 1))

theme_set(theme_bw(base_size=14))
  ggplot(sigtabRT4T2, aes(x=reorder(rowname, -log2FoldChange), y=log2FoldChange, fill=Phylum, label=taxa)) + 
  geom_hline(yintercept = 0.0, color = "black", size = 0.3) +
  geom_bar(stat="identity") + coord_flip() +
  geom_text(color="black", aes(y = ifelse(log2FoldChange > 0, -0.5, 0.5), hjust = ifelse(log2FoldChange > 0, 1, 0)))+
  scale_x_discrete(name="Taxa", breaks=NULL) +
  colScale + 
  scale_fill_manual(name = "Phylum", values = myColors)+
    labs(title = "Recipient Timepoint 2 vs 4") + expand_limits(y = -20) +
      geom_text(aes(label = rowname),position = position_stack(0.8),hjust = ifelse(sigtabRT4T2$log2FoldChange > 0, 0, 1))

#need manual edit RT4T3 for multiple NA taxa
sigtabRT4T3<- read.delim("sigtabRT4T3.txt", sep="\t", header=TRUE)
sigtabRT4T3<-arrange(sigtabRT4T3, log2FoldChange)
sigtabRT4T3$taxa <- gsub("_", " ", sigtabRT4T3$taxa)
theme_set(theme_bw(base_size=16))

  ggplot(sigtabRT4T3, aes(x=reorder(rowname, -log2FoldChange), y=log2FoldChange, fill=Phylum, label=taxa)) + 
  geom_hline(yintercept = 0.0, color = "black", size = 0.3) +
  geom_bar(stat="identity") + coord_flip() +
  geom_text(color="black", aes(y = ifelse(log2FoldChange > 0, -0.5, 0.5), hjust = ifelse(log2FoldChange > 0, 1, 0)))+
  scale_x_discrete(name="Taxa", breaks=NULL) +
  colScale + 
  scale_fill_manual(name = "Phylum", values = myColors)+
  labs(title = "Recipient Timepoint 3 vs 4")  + expand_limits(y = -20) +
    geom_text(aes(label = rowname), position = position_stack(0.8),hjust = ifelse(sigtabRT4T3$log2FoldChange > 0, 0, 1))

  
#DESEQ2 time series analysis > see differences between recipient and control over time from T2-T3-T4
#https://support.bioconductor.org/p/62684/
  #design(dds) <- ~ time + treat + time:treat
  #dds <- DESeq(dds, test="LRT", reduced = ~ time + treat)
levels(sample_data(ps02)$Treatment)
#[1] "Control"   "Recipient"
levels(sample_data(ps02)$Timepoint)
#[1] "T1" "T2" "T3" "T4"
ps03 <- subset_samples(ps02, !Timepoint =="T1")
dds3 =phyloseq_to_deseq2(ps03, ~Timepoint+Treatment+Timepoint:Treatment)
dds3 <- DESeq(dds3, test="LRT", reduced = ~ Treatment+Timepoint )
resultsNames(dds3)
#[1] "Intercept"                      "Timepoint_T3_vs_T2"             "Timepoint_T4_vs_T2"            
#[4] "Treatment_Recipient_vs_Control" "TimepointT3.TreatmentRecipient" "TimepointT4.TreatmentRecipient"
plotDispEsts(dds3)
betas <- coef(dds3)
colnames(dds3)
res3 <- results(dds3) 
sum(!is.na(res3$pvalue)) #1345
sum(is.na(res3$pvalue)) #405
sum(res3$pvalue<0.05, na.rm=TRUE) #5
sum(res3$padj<0.1, na.rm=TRUE) #0
sum(res3$padj <1, na.rm=TRUE) #0, adjusted p-value all =1, nothing passed the multiple-corrections
res3 = res3[order(res3$pvalue, na.last=NA),]
head(res3[order(res3$pvalue, na.last=NA),],20)

res3
sigtab = res3[(res3$pvalue < 0.05),]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps03)[rownames(sigtab), ], "matrix"))
head(sigtab)
betas<- na.omit(betas)
topGenes <- head(rownames(res3),20)
mat <- betas[topGenes, -c(1,1)]
thr <-3
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
library("pheatmap")
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)

#try the microbiome R package
biocLite("microbiome")
library(microbiome)
summarize_phyloseq(ps1)
  
#Alternative package for differntial abundance analysis (seems more statistically sound)
# A differential abundance analysis for the comparison of two or more conditions. 
# Useful for analyzing data from standard RNA-seq or meta-RNA-seq assays 
# as well as selected and unselected values from in-vitro sequence selections. 
# Uses a Dirichlet-multinomial model to infer abundance from counts, 
# optimized for three or more experimental replicates. 
# The method infers biological and sampling variation to calculate the expected false discovery rate (FDR),
# given the variation, based on a Wilcox rank test or Welch t-test (via aldex.ttest), 
# or a glm and Kruskal-Wallis test (via aldex.glm). Reports p-values and Benjamini-Hochberg corrected p-values.
source("https://bioconductor.org/biocLite.R")
biocLite("ALDEx2")


