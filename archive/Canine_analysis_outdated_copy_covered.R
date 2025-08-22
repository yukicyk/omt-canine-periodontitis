# Canine MiSeq Analysis Part 2

setwd("~/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/2/")

#import all the files for phyloseq
samdat <- read.delim("mappingV2.txt", header=TRUE)
sampledata <-samdat[,2:ncol(samdat)]
rownames(sampledata)<-samdat[,1]
seqtab<- readRDS("seqtab_rename.rds")
taxa_RDP_species_singles<- readRDS("taxa_RDP_species_single.rds")
ASVcode <- paste0("ASV", 1:ncol(seqtab))
rownames(taxa_RDP_species_singles) <- ASVcode

#import MLtree
library("ape")
phylo.tree <- read.tree("ASV.nwk")

#import into phyloseq
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(sampledata), 
               tax_table(taxa_RDP_species_singles), phy_tree(phylo.tree))
### Remove sample MT-E-V-2, and update the mapping file
ps1<- subset_samples(ps, sample_names(ps) !="MT-E-V-2")
# Supervised Prevalence Filtering
#Remove the features (ASVs) with ambigurous phylum annotation
ps0 <- subset_taxa(ps1, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))
#Compute the total and average prevalences of the features in each phylum.
library(dplyr)
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# Define phyla to filter (sum of prevalence <10)
filterPhyla = c("Acidobacteria", "Chloroflexi", "Cyanobacteria/Chloroplast")
ps2 = subset_taxa(ps0, !Phylum %in% filterPhyla)
# Unsupervised prevalence filtering at prevalence threshold at 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps0)
keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
ps3 = prune_taxa(keepTaxa, ps2)

#phyloseq heatmap-ASV
library(scales)
order_of_samples <- c("MT-B-1","MT-D-1","MT-G-1","MT-I-1","MT-K-1","MT-M-1","MT-O-1","MT-R-1","MT-S-1","MT-B-2","MT-D-2","MT-G-2","MT-I-2","MT-K-2","MT-M-2","MT-O-2","MT-R-2","MT-S-2","MT-B-3","MT-D-3","MT-G-3","MT-I-3","MT-K-3","MT-M-3","MT-O-3","MT-R-3","MT-S-3","MT-B-4","MT-D-4","MT-G-4","MT-I-4","MT-K-4","MT-M-4","MT-O-4","MT-R-4","MT-S-4","Transplant-A","Transplant-C","Transplant-F","Transplant-H","Transplant-J","Transplant-L","Transplant-N","Transplant-P","Transplant-Q","MT-E-I-1","MT-E-II-1","MT-E-III-1","MT-E-IV-1","MT-E-V-1","MT-A-1","MT-C-1","MT-F-1","MT-H-1","MT-J-1","MT-L-1","MT-N-1","MT-P-1","MT-Q-1","MT-A-2","MT-C-2","MT-F-2","MT-H-2","MT-J-2","MT-L-2","MT-N-2","MT-P-2","MT-Q-2","MT-A-3","MT-C-3","MT-F-3","MT-H-3","MT-J-3","MT-L-3","MT-N-3","MT-P-3","MT-Q-3","MT-A-4","MT-C-4","MT-F-4","MT-H-4","MT-J-4","MT-L-4","MT-N-4","MT-P-4","MT-Q-4")
plot_heatmap(ps3, sample.order = order_of_samples, method ="NMDS", distance ="bray",na.value = "black",trans = log_trans(2))

#R package heatmap function
heatmap(otu_table(ps3))

#Try the plot tree function in phyloseq
#look at the top 50 ASVs
physeq = prune_taxa(taxa_names(ps3)[1:50], ps3)
plot_tree(physeq, nodelabf=nodeplotblank, label.tips="taxa_names", ladderize="left", color="Treatment", shape="Timepoint", base.spacing=0.03)
