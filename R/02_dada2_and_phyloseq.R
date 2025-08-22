# Step 2: ASV Calling, Taxonomy Assignment, and Phyloseq Object Creation
## Script: `R/02_dada2_and_phyloseq.R`
## Source Files: `Canine_transplant_MiSeq_v2.R`
## Description: This is the core script for processing the trimmed reads. It performs quality filtering, error rate learning, dereplication, ASV inference, merging of paired-end reads, and chimera removal. It then assigns taxonomy using the custom reference files and combines the ASV table, taxonomy table, metadata, and phylogenetic tree into a single `phyloseq` object.
## Inputs:
###   Trimmed FASTQ files from Step 1.
###   `mappingV2.txt` metadata file.
###   Taxonomy reference files from Step 0.
###   A phylogenetic tree (`ASV.nwk`), which was built externally using `mafft` and `MEGA7`.
## Outputs:
###   `seqtab_final.rds`: The final ASV table.
###   `taxa_final.rds`: The final taxonomy table.
###   `phyloseq_object.rds`: The final, comprehensive phyloseq object for downstream analysis.

# Rstudio dada2 pipeline
# Using R studio Version 1.1.383; R version 3.4.2 (2017-09-28); Biocondcutor version 3.6 (2017-10-31); DADA2 version 1.6

library(dada2); packageVersion("dada2")

#### Filter and process MiSeq fastq files
path1 <- "/Users/ODPC/Desktop/Canine_transplant_MiSeq/Clean_Data/cut_adapt/" 
list.files(path1)

# Filter and Trim
# Sort ensures forward/reverse reads are in same order
fnF1s <- sort(list.files(path1, pattern="_1.fq.gz.out", full.names=TRUE))
fnR1s <- sort(list.files(path1, pattern="_2.fq.gz.out", full.names=TRUE))
sample1.names <-gsub("(.*)_1.fq.gz.out", "\\1", basename(fnF1s))

# Assign the filenames for the filtered fastq.gz files.
filt_path1 <- file.path(path1, "filtered") # Place filtered files in filtered/ subdirectory
filtF1s <- file.path(filt_path1, paste0(sample1.names, "_F_filt.fastq.gz"))
filtR1s <- file.path(filt_path1, paste0(sample1.names, "_R_filt.fastq.gz"))

# Examine quality profiles of forward and reverse reads and export as pdf file
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

# Perform filtering and trimming
out1 <- filterAndTrim(fnF1s, filtF1s, fnR1s, filtR1s, truncLen=c(250,230),
                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE, matchIDs=TRUE) 

# Learn the Error Rates
errF1 <- learnErrors(filtF1s, multithread=TRUE)
errR1 <- learnErrors(filtR1s, multithread=TRUE)

# Visualise the estimated error rates 
pdf("Error_plots.pdf")
plotErrors(errF1, nominalQ=TRUE)
plotErrors(errR1, nominalQ=TRUE)
dev.off()

# Dereplication 
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

# Construct sequence table
seqtab <- makeSequenceTable(mergers1)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# Remove short sequence <400bp (only keep those between 401-450)
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(401,450)]

# Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)

# Evaluate the fraction of chimeras 
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)

# Track filtering of reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out1, sapply(dadaF1s, getN), sapply(mergers1, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
track
seqtab <- seqtab.nochim


### Update classification using RDP + eHOMD +COT (not using silva database due to confusion of the taxonomy systems used;
### in particular the superphylum, Patescibacteria includes candidate phyla Gracilibacteria, Microgenomates, Parcubacteria, and Saccharibacteria
## Requires Taxonomy reference files from Step 0.

seqtab <- readRDS(file.choose())

# Assign taxonomy
taxa_RDP <- assignTaxonomy(seqtab, "/Users/OB/Desktop/temp/Taxonomy_files/2018/COT_HOMD_RDP_train_set_v2.txt", multithread=TRUE)
unname(head(taxa_RDP))

# Species level assignments
taxa_RDP_species<- addSpecies(taxa_RDP, "/Users/OB/Desktop/temp/Taxonomy_files/2018/COT_HOMD_RDP_spec_assign_v2.txt", allowMultiple=TRUE)
taxa_RDP_species_singles<- addSpecies(taxa_RDP, "/Users/OB/Desktop/temp/Taxonomy_files/2018/COT_HOMD_RDP_spec_assign_v2.txt", allowMultiple=FALSE)

# Re-formatting the ASV table and taxonomy data
# Renaming the samples
# Export the sequence in fasta format, under the ASV# code names (ASV1-ASVn)
library(seqinr)
ASVcode <- paste0("ASV", 1:ncol(seqtab))
write.fasta(as.list(colnames(seqtab)), ASVcode, "ASV.fasta")

# Rename the ASVs using code ASV# (ASV1-ASVn)
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

# Export as csv text file
write.csv(seqtab, "seqtab_renamed.csv", col.names=NA)
write.csv(taxa_RDP, "taxa_RDP_renamed.csv", col.names=NA)
write.csv(taxa_RDP_species_singles, "species_RDP_renamed.csv", col.names=NA)
taxa_RDP_species<- readRDS("~/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/taxa_RDP_species.rds")
write.csv(taxa_RDP_species, "species_RDP_multiple_renamed.csv", col.names=NA)

# Import the metadata mapping file
samdat <- read.delim("mappingV2.txt", header=TRUE)
sampledata <-samdat[,2:ncol(samdat)]
rownames(sampledata)<-samdat[,1]

# Import MLtree
library("ape")
phylo.tree <- read.tree("ASV.nwk")

# Import into phyloseq
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(sampledata), 
               tax_table(taxa_RDP_species_singles), phy_tree(phylo.tree))

### Remove obsolete sample MT-E-V-2, and update the mapping file
ps1<- subset_samples(ps, sample_names(ps) !="MT-E-V-2")
ps1