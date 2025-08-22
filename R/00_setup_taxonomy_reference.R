# **Step 0: Setup Taxonomy Reference Database**
## Script: `R/00_setup_taxonomy_reference.R`
## Source Files: `tidying_up_taxonomy_ref.R`, `tidying_up_taxonomy_species_ref.R`
### Description: This script prepares the customized taxonomy reference files for oral bacteria with the RDP trainset, and the Canine Oral Microbiome (COT) and Human Oral Microbiome Database (HOMD). This is a one-time setup step required before running the main DADA2 pipeline.
### Purpose: The customized taxonomy file enhance bacterial classification to more relevant oral species/strains/phylotypes.
### Inputs: Raw FASTA and taxonomy files from RDP, HOMD, and COT.
### Outputs: Formatted `COT_HOMD_RDP_train_set.fa` and `COT_HOMD_RDP_spec_assign_v2.txt` files for DADA2.
### Reference: Revising_taxonomy_ref.txt; 
### Reference: COT - 416 full-length sequences, GenBank: JN713151–JN713566, (Dewhirst et al. 2012) 
### Reference: HOMD - eHOMD 16S rRNA RefSeq Version 15.1 (Chen et al. 2010, last updated on 2017-11-16)
### Reference: RDP - RDP trainset 16, derived from RDP Release 11.5 that consists of 3,356,809 aligned and annotated 16S rRNA gene sequences, last updated on 2016-09-30
### Reference: Greengene - RDP is no longer updated, for future usage, please look for GG reference trainset

# oral reference sequence sets were downloaded and store in tsv and ran through this tidying_up_taxonomy_ref.R by Raymond Tong
library(dplyr)
library(readr)
library(tibble)

cFasta <- read_tsv(file.choose())
cFasta <- as.data.frame(cFasta)
startLines <- cFasta %>%
  as.data.frame() %>%
  .[, 1] %>%
  grep(">", .)

seqStart <- startLines + 1
seqEnd <- c(startLines[-1] - 1, nrow(cFasta) )
seqPos <- data.frame(startLines, seqStart, seqEnd)

seqPos <- cFasta[seqPos$startLines, 1] %>%
  as.data.frame() %>%
  cbind(seqPos)

sequ <- apply(seqPos, 1, function(x){
  paste(cFasta[x[3] : x[4], 1], collapse = "")
})

seqPos <- cbind(seqPos, sequ)
seqTab <- seqPos[, c(1, 5)] %>%
  `colnames<-`(c("id", "sequence"))

cTax <- read_tsv(file.choose(), col_names = FALSE)
cTax <- as.data.frame(cTax)
cTax$X1 <- paste0(">",cTax$X1)

seqTab2 <- left_join(seqTab, cTax, by = c("id" = "X1"))
seqTab2$X2 <- paste0(">", seqTab2$X2)

seqTab2$level <- seqTab2$X2 %>%
  tokenizers::tokenize_regex(";") %>%
  lapply(., function(x){
    paste(x[1:6], collapse = ";")
  }) %>%
  unlist()

tSet <- seqTab2 %>%
  select(level, sequence)

tSet$level <- paste0(tSet$level, ";")

tSet <- as.vector(t(tSet)) %>%
  as.data.frame()

write_tsv(tSet, "COT_train_set.fa", col_names = FALSE)

seqTab3 <- seqTab2

seqTab3$level2 <- seqTab3$X2 %>%
  tokenizers::tokenize_regex(";") %>%
  lapply(., function(x){
    paste(x[6:7], collapse = " ")
  }) %>%
  unlist()

seqTab3$level2 <- paste(seqTab3$id, seqTab3$level2)

specAssign <- seqTab3 %>%
  select(level2, sequence)

specAssign <- as.vector(t(specAssign)) %>%
  as.data.frame()

write_tsv(specAssign, "COT_species_assignment.fa", col_names = FALSE)


# tidying_up_taxonnomy_species_ref.R by Raymond Tong
hFasta <- read_tsv(file.choose(), col_names = FALSE)

hTax <- read_tsv(file.choose())

hFastaTab <- hFasta %>%
  as.data.frame() %>%
  .[, 1] %>%
  matrix(ncol = 2, byrow = TRUE) %>%
  as.data.frame()

hFastaTab$id <- hFastaTab$V1 %>%
  as.character() %>%
  tokenizers::tokenize_regex(" \\| ") %>%
  lapply(., function(x){x[1]}) %>%
  unlist()

hFastaTab$Species <- hFastaTab$V1 %>%
  as.character() %>%
  tokenizers::tokenize_regex(" \\| ") %>%
  lapply(., function(x){x[2]}) %>%
  unlist()

hFastaTab$hmt <- hFastaTab$V1 %>%
  as.character() %>%
  tokenizers::tokenize_regex(" \\| ") %>%
  lapply(., function(x){x[3]}) %>%
  unlist() %>%
  substring(5)

hTax2 <- hTax %>%
  .[, 1:8] %>%
  mutate(level = paste(Domain, Phylum, Class, Order, Family, Genus, sep = ";"))

hTax2$level <- paste0(hTax2$level, ";")

hFastaTab2 <- left_join(hFastaTab, hTax2, by = c("hmt" = "HMT_ID"))

hFastaTab2 <- hFastaTab2 %>%
  mutate(level2 = paste(id, Species.x))

hFastaTab2$level <- paste0(">", hFastaTab2$level)

hTSet <- hFastaTab2 %>%
  select(level, V2) %>%
  as.data.frame() %>%
  t() %>%
  as.vector() %>%
  as.data.frame()

write_tsv(hTSet, "HOMD_train_set.fa", col_names = FALSE)

hSpecAssign <- hFastaTab2 %>%
  select(level2, V2) %>%
  as.data.frame() %>%
  t() %>%
  as.vector() %>%
  as.data.frame()

write_tsv(hSpecAssign, "HOMD_species_assignment.fa", col_names = FALSE)

