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

seqtab <- readRDS(file.choose())

# Assign taxonomy
taxa_RDP <- assignTaxonomy(seqtab, "/Users/OB/Desktop/temp/Taxonomy_files/2018/COT_HOMD_RDP_train_set.fa", multithread=TRUE)
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
               tax_table(taxa_RDP_species), phy_tree(phylo.tree))


### Remove sample MT-E-V-2, and update the mapping file
### Update classification using RDP + eHOMD +COT (not using silva database due to confusion of the taxonomy systems used;
### in oarticular the superphylum, Patescibacteria that includes candidate phyla Gracilibacteria, Microgenomates, Parcubacteria, and Saccharibacteria

# Plot Shannon index
plot_richness(ps, x="Dog", measures=c("Shannon"), color="Treatment", shape="Timepoint") + 
  geom_point(size=7, alpha=0.75)+ ggtitle("Alpha Diversity measure Shannon Index")
plot_richness(ps, x="Dog", measures=c("Observed"), color="Treatment", shape="Timepoint") + 
  geom_point(size=7, alpha=0.75)+ ggtitle("Alpha Diversity measure Observed ASVs")
theme_set(theme_bw())

#plot bar plots
#All samples, phylum
p= plot_bar(ps, "Phylum", fill="Phylum", facet_grid=Timepoint~Treatment)
p  + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
#top 10 phyla
pps <- tax_glom(ps,"Phylum")
top10p <- names(sort(taxa_sums(pps), decreasing=TRUE))[1:10]
pps.top10 <- transform_sample_counts(pps, function(ASV) ASV/sum(ASV))
pps.top10 <- prune_taxa(top10p, pps.top10)
p = plot_bar(pps.top10, x="Dog", fill="Phylum") + facet_grid(Timepoint~Treatment)
#export phyla data
phylum <- tax_glom(ps, "Phylum")
write.table(otu_table(phylum),"phylum.txt", sep="\t")

#Plot abundance
plot_abundance = function(physeq, ylabn = "",
                          Facet = "Order",
                          Color = "Phylum"){
  mphyseq = psmelt(physeq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq,
         mapping = aes_string(x = "Dog", y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + ylab(ylabn) +
    scale_y_log10()
}
plot_abundance(ps)



#UniFrac weighted and unweighted distance calculation and ordination
ordu = ordinate(ps, "PCoA", "unifrac", weighted=TRUE)
ordu2 = ordinate(ps, "PCoA", "unifrac", weighted=FALSE)
ordu3 =ordinate(ps, "PCoA", "bray")
plot_ordination(ps, ordu, color="Treatment", shape="Timepoint", label = "Dog") +
  geom_point(size=2) +
  ggtitle("MDS/PCoA on weighted-UniFrac distance, GlobalPatterns")



ps %>%
  subset_samples(select= Timepoint == "T1") %>%
  ordinate(method = "PCoA", distance = "bray") %>%
  plot_ordination(ps, ., color = "Treatment", label = "Dog", title="PCoA-BrayCurtis") +
  geom_point(size=2) +
  geom_line()


ps %>%
  prune_samples(Timepoint =="T0") %>%
  ordinate(method = "PCoA", distance = "unifrac") %>%
  plot_ordination(ps, ., color = "Treatment", label = "Dog", title="PCoA-BrayCurtis-T1-4")

ps2 %>%
  ordinate(method = "PCoA", distance = "wunifrac") ->ord

ps2 %>%
  sample_data()

ps2 %>%
  sample_data() %>%
  as("matrix") %>%
  as.data.frame() -> samTab


plotDat <- ord$vectors[,c("Axis.1", "Axis.2")] %>%
  .[match(rownames(sampledata), rownames(.)),] %>%
  .[c(19,5,1,2,3,4,6,7,8,9,19,14,10,11,12,13,15,16,17,18,20,29,25,26,27,28,30,31,32,33,20,38,34,35,36,37,39,40,41,42,21,47,43,44,45,46,48,49,50,51,21,56,52,53,54,55,57,58,59,60,22,65,61,62,63,64,66,67,68,69,23,74,70,71,72,73,24,79,75,76,77,78,80,81,82,83,84,85,86,87),]

plotDat2 <- plotDat
plotDat2 <- cbind(plotDat2, rownames(plotDat))

plotDat3 <- plotDat2 %>%
  `rownames<-`(NULL) %>%
  as.data.frame() %>%
  mutate(dog = c(rep("a", 6), rep("b", 4), rep("c", 6), rep("d", 4), rep("f", 6), rep("g", 4), rep("h", 6), rep("i", 4), rep("j", 6), rep("k", 4), rep("l", 6), rep("m", 4), rep("n", 6), rep("o", 4), rep('p', 6), rep("q", 6), rep("r", 4), rep("s", 4)))

samTab %>%
  rownames_to_column("sample") -> samTab

left_join(plotDat3, samTab, by=c("V3" = "sample")) -> plotTab4

ggplot(plotTab4) +
  geom_line(aes(x = as.numeric(as.character(Axis.1)), y = as.numeric(as.character(Axis.2)), group = dog, color = Treatment)) +
  geom_point(aes(x = as.numeric(as.character(Axis.1)), y = as.numeric(as.character(Axis.2)), shape = Timepoint), size=3 )+
  labs(title="MDS/PCoA on weighted-UniFrac distance, GlobalPatterns", x="PCoA Axis 1 (35.1% total variation)" , y="PCoA Axis 2 (15.3% total variation)")



plotDat2
ggplot(plotDat2) +
  geom_line(aes(x = Axis.1, y = Axis.2, group = dog, color = dog))


ps %>%
  tax_glom("Genus") -> gps
ps %>%
  tax_glom("Class") -> cps
ps %>%
  tax_glom("Family") -> fps

library(tibble)

#F1000 workflow Filtering of data based on prevalence
# Show available ranks in the dataset
rank_names(ps)
# Create table, number of features (ASVs) for each phylum
table(tax_table(ps)[, "Phylum"], exclude = NULL)
#Remove the features (ASVs) with ambigurous phylum annotation
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
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
filterPhyla = c("Acidobacteria", "Chloroflexi", "Cyanobacteria", "Planctomycetes")
# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps0, !Phylum %in% filterPhyla)
ps1 #1736 taxa

# Unsupervised prevalence filtering
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
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
ps2 = prune_taxa(keepTaxa, ps0)
ps2 #809 taxa


# Plot Shannon index
plot_richness(ps2, x="Dog", measures=c("Shannon"), color="Treatment", shape="Timepoint") + 
  geom_point(size=7, alpha=0.75)+ ggtitle("Alpha Diversity measure Shannon Index")
plot_richness(ps2, x="Dog", measures=c("Observed"), color="Treatment", shape="Timepoint") + 
  geom_point(size=7, alpha=0.75)+ ggtitle("Alpha Diversity measure Observed ASVs")

#ordination
ordu = ordinate(ps2, "PCoA", "unifrac", weighted=TRUE)
ordu2 = ordinate(ps2, "PCoA", "unifrac", weighted=FALSE)
ordu3 =ordinate(ps2, "PCoA", "bray")
plot_ordination(ps2, ordu3, color="Treatment", shape="Timepoint", label = "Dog") +
  geom_point(size=2) +
  ggtitle("PCoA Bray-Curtis, prevalence filtered data")

