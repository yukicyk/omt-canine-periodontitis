#Rstudio dada2 pipeline
#Using R studio Version 1.1.383; R version 3.4.2 (2017-09-28); Biocondcutor version 3.6 (2017-10-31); DADA2 version 1.6

library(dada2); packageVersion("dada2")

#####Filter and process HKMiSeq1
path1 <- "/Users/ODPC/Desktop/Canine_transplant_MiSeq/Clean_Data/cut_adapt/" 
list.files(path1)
#Filter and Trim
# Sort ensures forward/reverse reads are in same order
fnF1s <- sort(list.files(path1, pattern="_1.fq.gz.out", full.names=TRUE))
fnR1s <- sort(list.files(path1, pattern="_2.fq.gz.out", full.names=TRUE))
sample1.names <-gsub("(.*)_1.fq.gz.out", "\\1", basename(fnF1s))
# Assign the filenames for the filtered fastq.gz files.
filt_path1 <- file.path(path1, "filtered") # Place filtered files in filtered/ subdirectory
filtF1s <- file.path(filt_path1, paste0(sample1.names, "_F_filt.fastq.gz"))
filtR1s <- file.path(filt_path1, paste0(sample1.names, "_R_filt.fastq.gz"))

#Examine quality profiles of forward and reverse reads and export as pdf file
pdf("Canine_transplant_MiSeq_R1_QS.pdf")
plotQualityProfile(fnF1s[1:12])
plotQualityProfile(fnF1s[13:24])
plotQualityProfile(fnF1s[25:36])
plotQualityProfile(fnF1s[37:48])
plotQualityProfile(fnF1s[49:60])
plotQualityProfile(fnF1s[61:72])
plotQualityProfile(fnF1s[73:84])
plotQualityProfile(fnF1s[85:87])
dev.off()
pdf("Canine_transplant_MiSeq_R2_QS.pdf")
plotQualityProfile(fnR1s[1:12])
plotQualityProfile(fnR1s[13:24])
plotQualityProfile(fnR1s[25:36])
plotQualityProfile(fnR1s[37:48])
plotQualityProfile(fnR1s[49:60])
plotQualityProfile(fnR1s[61:72])
plotQualityProfile(fnR1s[73:84])
plotQualityProfile(fnR1s[85:87])
dev.off()

#Perform filtering and trimming
out1 <- filterAndTrim(fnF1s, filtF1s, fnR1s, filtR1s, truncLen=c(250,230),
                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE, matchIDs=TRUE) 
#Learn the Error Rates
errF1 <- learnErrors(filtF1s, multithread=TRUE)
errR1 <- learnErrors(filtR1s, multithread=TRUE)
#visualise the estimated error rates 
pdf("Error_plots.pdf")
plotErrors(errF1, nominalQ=TRUE)
plotErrors(errR1, nominalQ=TRUE)
dev.off()
#Dereplication 
derepF1s <- derepFastq(filtF1s, verbose=FALSE)
derepR1s <- derepFastq(filtR1s, verbose=FALSE)
# Name the derep-class objects by the sample names <-run this and belows
names(derepF1s) <- sample1.names
names(derepR1s) <- sample1.names
#Sample Inference: Infer the sequence variants in each sample
dadaF1s <- dada(derepF1s, err=errF1, multithread=TRUE)
dadaR1s <- dada(derepR1s, err=errR1, multithread=TRUE)
# Merge the denoised forward and reverse reads
mergers1 <- mergePairs(dadaF1s, derepF1s, dadaR1s, derepR1s, verbose=FALSE)
#Construct sequence table
seqtab <- makeSequenceTable(mergers1)
dim(seqtab)
table(nchar(getSequences(seqtab)))
#Remove short sequence <400bp (only keep those between 401-450)
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(401,450)]
#Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
#Evaluate the fraction of chimeras 
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)
#Track filtering of reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out1, sapply(dadaF1s, getN), sapply(mergers1, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
track
seqtab <- seqtab.nochim


### Update classification using RDP + eHOMD +COT (not using silva database due to confusion of the taxonomy systems used;
### in oarticular the superphylum, Patescibacteria that includes candidate phyla Gracilibacteria, Microgenomates, Parcubacteria, and Saccharibacteria


seqtab <- readRDS(file.choose())

# Assign taxonomy
taxa_RDP <- assignTaxonomy(seqtab, "/Users/OB/Desktop/temp/Taxonomy_files/2018/COT_HOMD_RDP_train_set_v2.txt", multithread=TRUE)
unname(head(taxa_RDP))

# Species level assignments
taxa_RDP_species<- addSpecies(taxa_RDP, "/Users/OB/Desktop/temp/Taxonomy_files/2018/COT_HOMD_RDP_spec_assign_v2.txt", allowMultiple=TRUE)
taxa_RDP_species_singles<- addSpecies(taxa_RDP, "/Users/OB/Desktop/temp/Taxonomy_files/2018/COT_HOMD_RDP_spec_assign_v2.txt", allowMultiple=FALSE)

#re-formatting the ASV table and taxonomy data
#Renaming the samples
#export the sequence in fasta format, under the ASV# code names (ASV1-ASVn)
library(seqinr)
ASVcode <- paste0("ASV", 1:ncol(seqtab))
write.fasta(as.list(colnames(seqtab)), ASVcode, "ASV.fasta")
#rename the ASVs using code ASV# (ASV1-ASVn)
colnames(seqtab) <- ASVcode
rownames(taxa_RDP) <- ASVcode
rownames(taxa_RDP_species) <- ASVcode
rownames(taxa_RDP_species_singles) <- ASVcode
rownames(seqtab) = gsub(pattern = "_F_filt.fastq.gz", replacement = "", x = rownames(seqtab))

# Write to disk
saveRDS(seqtab, "seqtab_rename.rds")
saveRDS(taxa_RDP, "taxa_RDP.rds")
saveRDS(taxa_RDP_species_singles, "taxa_RDP_species_single.rds")
saveRDS(taxa_RDP_species, "taxa_RDP_species.rds")

#export as csv text file
write.csv(seqtab, "seqtab_renamed.csv", col.names=NA)
write.csv(taxa_RDP, "taxa_RDP_renamed.csv", col.names=NA)
write.csv(taxa_RDP_species_singles, "species_RDP_renamed.csv", col.names=NA)
taxa_RDP_species<- readRDS("~/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/taxa_RDP_species.rds")
write.csv(taxa_RDP_species, "species_RDP_multiple_renamed.csv", col.names=NA)

#import the metadata mapping file
samdat <- read.delim("mappingV2.txt", header=TRUE)
sampledata <-samdat[,2:ncol(samdat)]
rownames(sampledata)<-samdat[,1]

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
ps1

# Plot Shannon index
plot_richness(ps1, x="Dog", measures=c("Shannon"), color="Treatment", shape="Timepoint") + 
  geom_point(size=7, alpha=0.75)+ ggtitle("Alpha Diversity measure Shannon Index")
plot_richness(ps1, x="Dog", measures=c("Observed"), color="Treatment", shape="Timepoint") + 
  geom_point(size=7, alpha=0.75)+ ggtitle("Alpha Diversity measure Observed ASVs")
theme_set(theme_bw())

#UniFrac weighted and unweighted distance calculation and ordination
wuf = ordinate(ps1, "PCoA", "unifrac", weighted=TRUE)
uwuf = ordinate(ps1, "PCoA", "unifrac", weighted=FALSE)
bc =ordinate(ps1, "PCoA", "bray")
plot_ordination(ps1, bc, color="Treatment", shape="Timepoint", label = "Dog") +
    geom_point(size=3) +
    ggtitle("MDS/PCoA on Bray-Curtis dissimilarity")


library(tibble)
#F1000 workflow Filtering of data based on prevalence at the phylum level

#export phyla data manual check of data
ps1p <- tax_glom(ps1, "Phylum")
write.table(otu_table(ps1p),"phylum.txt", sep="\t", col.names = tax_table(ps1p)[, "Phylum"])
write.table(taxa_sums(ps1p),"phylum_sums.txt", sep="\t", row.names = tax_table(ps1p)[, "Phylum"])
# Create table, number of features (ASVs) for each phylum
table(tax_table(ps1)[, "Phylum"], exclude = NULL) 

#Remove the features (ASVs) with ambigurous phylum annotation
ps0 <- subset_taxa(ps1, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
#Can we further calssify the ambiguorous sequences?
psna <- subset_taxa(ps1, is.na(Phylum))
all_ASV_seq <- read.fasta(file.choose())
na_seq_id <- taxa_names(psna)
names(all_ASV_seq) %in% na_seq_id
na_seq <- all_ASV_seq[c(which(names(all_ASV_seq) %in% na_seq_id))]
write.fasta(sequences=na_seq, names=names(na_seq), file.out="na_seq.fasta")

# Define prevalence of each taxa
# (in how many samples did each taxa appear at least once)
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))
#Compute the total and average prevalences of the features in each phylum.
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Supervised Prevalence Filtering
# Define phyla to filter (sum of prevalence <10)
filterPhyla = c("Acidobacteria", "Chloroflexi", "Cyanobacteria/Chloroplast")
# Filter entries with unidentified Phylum.
ps2 = subset_taxa(ps0, !Phylum %in% filterPhyla)
ps2 #1739 taxa

# Unsupervised prevalence filtering
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps2, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps0),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps0)
prevalenceThreshold
## [1] 4.35

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps3 = prune_taxa(keepTaxa, ps0)
ps3 #811 taxa

# Export table
ps3p <- tax_glom(ps3, "Phylum")
write.table(otu_table(ps3p),"ps3.txt", sep="\t", col.names = tax_table(ps3p)[, "Phylum"])
write.table(taxa_sums(ps3p),"phylum_sums_ps3.txt", sep="\t", row.names = tax_table(ps3p)[, "Phylum"])
# Plot Shannon index and Observed ASVs dot plots
plot_richness(ps3, x="Dog", measures=c("Shannon"), color="Treatment", shape="Timepoint") + 
  geom_point(size=7, alpha=0.75)+ ggtitle("Alpha Diversity measure Shannon Index")
plot_richness(ps3, x="Dog", measures=c("Observed"), color="Treatment", shape="Timepoint") + 
  geom_point(size=7, alpha=0.75)+ ggtitle("Alpha Diversity measure Observed ASVs")



# ordination plot samples
ord = ordinate(ps3, "NMDS", "unifrac", weighted=TRUE)
ord2 = ordinate(ps3, "NMDS", "unifrac", weighted=FALSE)
ord3 =ordinate(ps3, "NMDS", "bray")
ord4 = ordinate(ps3, "MDS","dpcoa")
plot_ordination(ps3, ord4, color="Treatment", shape="Timepoint", label = "Dog") +
  geom_point(size=3, alpha=0.5) +
  ggtitle("MDS on DPCoA, prevalence filtered data")
# Test various distance method for MDS ordination on ps3 data
dist_methods <- unlist(distanceMethodList)
# Remove the user-defined distance
dist_methods = dist_methods[-which(dist_methods=="ANY")]
dist_methods = dist_methods[-which(dist_methods=="manhattan")]
#Loop through each distance method, save each plot to a list, called plist
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  distance(ps3s, method=i) %>%
  ordinate(ps3s, "MDS", distance=.) %>%
  # Create plot, store as temp variable, p
  plot_ordination(ps3s, ., color="Treatment", shape="Timepoint") -> p
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}
# Combine results
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=Treatment, shape=Timepoint))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for prevalence filtered data (Species)")
p
#obtain the axis information (variance%)
df2 = laply(plist, function(x) x$labels)

# ordination plot taxa
theme_set(theme_bw())
phyla.ord <- ordinate(ps3, "NMDS", "bray")
plot_ordination(ps3, phyla.ord, type="taxa", color="Phylum", title="Taxa NMDS on Bray-Curtis Dissimilartiy, prevalence filtered data") + 
  facet_wrap(~Phylum, 6) 

ps2p <- tax_glom(ps2, "Phylum")
ps2c <- tax_glom(ps2, "Class")
ps2o <- tax_glom(ps2, "Order")
ps2f <- tax_glom(ps2, "Family")
ps2g <- tax_glom(ps2, "Genus")
ps2s <- tax_glom(ps2, "Species")

ps3p <- tax_glom(ps3, "Phylum")
ps3c <- tax_glom(ps3, "Class")
ps3o <- tax_glom(ps3, "Order")
ps3f <- tax_glom(ps3, "Family")
ps3g <- tax_glom(ps3, "Genus")
ps3s <- tax_glom(ps3, "Species")

#Heatmap
#phyloseq heatmap
order_of_samples <- c("MT-B-1","MT-D-1","MT-G-1","MT-I-1","MT-K-1","MT-M-1","MT-O-1","MT-R-1","MT-S-1","MT-B-2","MT-D-2","MT-G-2","MT-I-2","MT-K-2","MT-M-2","MT-O-2","MT-R-2","MT-S-2","MT-B-3","MT-D-3","MT-G-3","MT-I-3","MT-K-3","MT-M-3","MT-O-3","MT-R-3","MT-S-3","MT-B-4","MT-D-4","MT-G-4","MT-I-4","MT-K-4","MT-M-4","MT-O-4","MT-R-4","MT-S-4","Transplant-A","Transplant-C","Transplant-F","Transplant-H","Transplant-J","Transplant-L","Transplant-N","Transplant-P","Transplant-Q","MT-E-I-1","MT-E-II-1","MT-E-III-1","MT-E-IV-1","MT-E-V-1","MT-A-1","MT-C-1","MT-F-1","MT-H-1","MT-J-1","MT-L-1","MT-N-1","MT-P-1","MT-Q-1","MT-A-2","MT-C-2","MT-F-2","MT-H-2","MT-J-2","MT-L-2","MT-N-2","MT-P-2","MT-Q-2","MT-A-3","MT-C-3","MT-F-3","MT-H-3","MT-J-3","MT-L-3","MT-N-3","MT-P-3","MT-Q-3","MT-A-4","MT-C-4","MT-F-4","MT-H-4","MT-J-4","MT-L-4","MT-N-4","MT-P-4","MT-Q-4")
plot_heatmap(ps2p, taxa.label = "Phylum", low="#4575b4", high="#d73027", sample.label = "Treatment", sample.order = order_of_samples)
p = plot_heatmap(ps3, sample.label = "Treatment", sample.order = order_of_samples)
# p + geom_vline (aes(xintercept =9 ,color="white"))
heatmap(otu_table(ps3))

#superheat customized heatmap
library(superheat)
ps2_phyla <- read.csv(file.choose(), header=TRUE)
phyla_heatmap <-data.matrix(ps2_phyla[,2:ncol(ps2_phyla)])
rownames(phyla_heatmap) <- ps2_phyla[,1]
#transform to relative abundance (proportional transformation)
phyla_heatmapR = prop.table(phyla_heatmap, margin = 2)
#prepare mapping labels for the samples
Treatment <- rep("Control",times=36)
Treatment <-append(Treatment, rep("Donor", times=14))
Treatment <-append(Treatment, rep("Recipient", times=36))
Timepoint <- rep(c("T1","T2","T3","T4","T1"),each=9)
Timepoint <- append(Timepoint, rep(c("T2"),each=5))
Timepoint <- append(Timepoint, rep(c("T1","T2","T3","T4"),each=9))
sample_reads_sum <- c(21326,21381,21324,21666,23246,20485,22120,19874,25030,24055,20221,19090,20327,22043,20505,19698,20636,20441,23569,20919,19538,21091,23436,21587,20178,20535,20127,22875,20368,22783,22214,20961,19562,28986,20178,19579,19070,19030,29158,18614,21314,21476,20902,20442,19784,20725,19876,19000,19416,19545,20445,20350,23681,21837,22715,19864,20586,20873,21553,20798,25740,19841,21461,21650,20978,19635,20557,24177,24512,23427,20756,20464,25361,27041,19366,19983,20853,19980,21563,20552,22596,24029,21734,24779,20067,19333)
phyla_reads_sum <- c(350,870,9593,15039,25251,37809,43700,45025,48369,93099,104684,317726,494525,607403)

superheat(phyla_heatmapR, scale=FALSE, 
          heat.pal = c("#4575b4", "#91bfdb", "#91cf60","#fee090", "#fc8d59", "#d73027" ),  heat.lim = c(0, 1), #define color scale
          heat.pal.values = c(0, 0.1,0.2, 0.4, 0.6, 1),
          left.label.text.alignment = "right",  
          left.label.text.size = 4, grid.hline = FALSE, 
          left.label.col = "white", bottom.label = "none",
          row.title = "Phylum", 
          grid.vline.col = "white", 
          yt= colSums(phyla_heatmap), yt.plot.type = "bar",  
          yt.axis.name = "No. of reads\nper sample",
          yr = phyla_reads_sum  / 1000, yr.axis.name = "Total no. of reads \n(in thousands)", 
          yr.plot.type = "bar")

#Transform 
transform_sample_counts(GP.chl, function(OTU) OTU/sum(OTU) )
