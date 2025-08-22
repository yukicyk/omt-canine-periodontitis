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
seqtab <- seqtabe.nochim

# Assign taxonomy
taxa <- assignTaxonomy(seqtab, "~/Desktop/Xiaolin/Shanghai/Dada2/Taxonomy_Reference/COT_HOMD_sliva_train_set.fa", multithread=TRUE)
unname(head(taxa))

# Species level assignments
taxa_species<- addSpecies(taxa, "~/Desktop/Xiaolin/Shanghai/Dada2/Taxonomy_Reference/COT_HOMD_sliva_spec_assign.fa", allowMultiple=FALSE)

# Write to disk
saveRDS(seqtab, "seqtab_final.rds")
saveRDS(taxa, "taxa.rds")
saveRDS(taxa_species, "species.rds")

#export as csv text file
write.csv(seqtab, "seqtab_final.csv", col.names=NA)

#import into phyloseq #shiny phyloseq
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
shiny::runGitHub("shiny-phyloseq","joey711")
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

#re-formatting the ASV table and taxonomy data
#Renaming the samples
library(dplyr)
rownames(seqtab) %>%
  as.data.frame() %>%
  apply(., 1, function(x){substr(x, 1, nchar(x) - 16)}) 
#export the sequence in fasta
seqs <- getSequences(seqtab)
#rename the ASVs using code ASV# instead of the sequences

library("DECIPHER")
seqs <- getSequences(seqtab)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
