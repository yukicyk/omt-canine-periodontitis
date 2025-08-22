#### **Step 5: Functional Prediction Analysis**
## Script: `R/05_functional_prediction_analysis.R`
## Source Files: `dada2_picrust.R`, `koAnalysis.R`
## Description: This script details the workflow for predicting the functional potential of the microbial communities. It first prepares the DADA2 output for use with PICRUSt (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States). It then analyzes the resulting KEGG Orthology (KO) abundance tables to find functional pathways that differ between groups and timepoints.
## Inputs: `seqtab_final.rds`, Greengenes reference database.
## Outputs: Predicted KEGG pathway abundance tables and plots comparing functional profiles.

# DADA2 to PICRUSt 
# Canine MiSeq for canine transplant dataset
# 20180621 YC on ODPC iMac
# source of code https://github.com/vmaffei/dada2_to_picrust

# Set working directory where files are saved
setwd("~/Desktop/Canine_transplant_MiSeq/Analysis/dada2_to_picrust")
# Dependencies: ShortRead & biom
install.packages("devtools")
source("https://bioconductor.org/biocLite.R")
biocLite("rhdf5")
biocLite("ShortRead")
install_github("joey711/biom")
library(devtools)
library(ShortRead)
library(biom) # note: use Joey's biom latest dev version;
# 1) Make study db
# grab study seqs
seqtab.nochim <- readRDS(file=file.choose())
seqs_study <- colnames(seqtab.nochim)
ids_study <- paste("ASV", 1:ncol(seqtab.nochim), sep = "_")
# merge db and study seqs
db_out <- data.frame(ids=ids_study,seqs=seqs_study,count=colSums(seqtab.nochim))
fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))
# write study fasta for filtering
writeFasta(fasta, file = "~/Desktop/Canine_transplant_MiSeq/Analysis/dada2_to_picrust/genome_prediction/gg_13_5_study_db.fasta.pre")
# filter sequences that diverge from gg_13_5 by 97%
# depending on how well greengenes covers your study sequences, consider reducing 97% to 70 or 50%
system('vsearch --usearch_global genome_prediction/gg_13_5_study_db.fasta.pre --db gg_13_5.fasta --matched genome_prediction/gg_13_5_study_db.fasta --id 0.97
') ##Run this command in terminal instead 

setwd("~/Desktop/Canine_transplant_MiSeq/Analysis/dada2_to_picrust/genome_prediction")

id_filtered <- as.character(id(readFasta("gg_13_5_study_db.fasta")))
db_out_filt <- db_out[db_out$ids%in%id_filtered,]
seqtab_biom <- t(seqtab.nochim)
# 2) output seq variant count data as biom;
# subset seqtab and output sample count biom
seqtab_biom <- seqtab_biom[rownames(seqtab_biom)%in%db_out_filt$seqs,]
rownames(seqtab_biom) <- db_out_filt[db_out_filt$seqs%in%rownames(seqtab_biom),"ids"]
biom_object <- biom::make_biom(data = seqtab_biom)
biom::write_biom(biom_object, biom_file = "sample_counts.biom")
# create final study db
system('cat gg_13_5.fasta >> gg_13_5_study_db.fasta') #using full dataset
system('cat gg_97_otus.fasta >> gg_13_5_study_db.fasta') #using smaller pre-cluster dataset

# biom files were analyzed with picrust after nomalization normalize_by_copy_number.py and predict_metagenomes.py to predict 
# Output is a table of function counts (e.g. KEGG KOs) by sample ids, export in tab-delimted text

# ko prediction results are sorted
library(dplyr)
library(readr)
library(tidyr)
library(biomformat)
library(tibble)
library(ggplot2)

l2 <- read_tsv("predicted_metagenomes.L2.txt", skip = 1)
l3 <- read_tsv("predicted_metagenomes.L3.txt", skip = 1)
bm <- read_biom(file.choose())

l2 <- l2 %>%
  gather("sample", "value", -`#OTU ID`, -KEGG_Pathways) %>%
  `colnames<-`(c("ko", "pathways", "sample", "value"))

l3 <- l3 %>%
  gather("sample", "value", -`#OTU ID`, -KEGG_Pathways) %>%
  `colnames<-`(c("ko", "pathways", "sample", "value"))

mapDat <- read_tsv(file.choose())

l2 <- l2 %>%
  mutate(sample = substring(sample, 1, nchar(sample) - 16)) %>%
  left_join(mapDat, by = c("sample" = "Sample"))

l3 <- l3 %>%
  mutate(sample = substring(sample, 1, nchar(sample) - 16)) %>%
  left_join(mapDat, by = c("sample" = "Sample"))

bm <- bm %>%
  biom_data() %>%
  as("matrix") %>%
  as.data.frame() %>%
  rownames_to_column("ko") %>%
  gather("sample", "value", -ko)

bm <- bm %>%
  mutate(sample = substring(sample, 1, nchar(sample) - 16)) %>%
  left_join(mapDat, by = c("sample" = "Sample"))

# l2

rec_l2_t3_t2 <- l2 %>%
  filter(Timepoint %in% c("T2", "T3")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Recipeint") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T2_mean = mean(unlist(T2), na.rm = TRUE), T3_mean = mean(unlist(T3), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T2), unlist(T3), paired = TRUE)$p.value) 

rec_l2_t3_t2_c <- l2 %>%
  filter(Timepoint %in% c("T2", "T3")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Control") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T2_mean = mean(unlist(T2), na.rm = TRUE), T3_mean = mean(unlist(T3), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T2), unlist(T3), paired = TRUE)$p.value) 

rec_l2_t4_t2 <- l2 %>%
  filter(Timepoint %in% c("T2", "T4")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Recipeint") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T2_mean = mean(unlist(T2), na.rm = TRUE), T4_mean = mean(unlist(T4), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T2), unlist(T4), paired = TRUE)$p.value) 

rec_l2_t4_t2_c <- l2 %>%
  filter(Timepoint %in% c("T2", "T4")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Control") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T2_mean = mean(unlist(T2), na.rm = TRUE), T4_mean = mean(unlist(T4), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T2), unlist(T4), paired = TRUE)$p.value) 

rec_l2_t4_t3 <- l2 %>%
  filter(Timepoint %in% c("T4", "T3")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Recipeint") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T4_mean = mean(unlist(T4), na.rm = TRUE), T3_mean = mean(unlist(T3), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T4), unlist(T3), paired = TRUE)$p.value) 

rec_l2_t4_t3_c <- l2 %>%
  filter(Timepoint %in% c("T4", "T3")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Control") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T4_mean = mean(unlist(T4), na.rm = TRUE), T3_mean = mean(unlist(T3), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T4), unlist(T3), paired = TRUE)$p.value)

# l3

rec_l3_t3_t2 <- l3 %>%
  filter(Timepoint %in% c("T2", "T3")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Recipeint") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T2_mean = mean(unlist(T2), na.rm = TRUE), T3_mean = mean(unlist(T3), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T2), unlist(T3), paired = TRUE)$p.value) 

rec_l3_t3_t2_c <- l3 %>%
  filter(Timepoint %in% c("T2", "T3")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Control") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T2_mean = mean(unlist(T2), na.rm = TRUE), T3_mean = mean(unlist(T3), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T2), unlist(T3), paired = TRUE)$p.value) 

rec_l3_t4_t2 <- l3 %>%
  filter(Timepoint %in% c("T2", "T4")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Recipeint") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T2_mean = mean(unlist(T2), na.rm = TRUE), T4_mean = mean(unlist(T4), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T2), unlist(T4), paired = TRUE)$p.value) 

rec_l3_t4_t2_c <- l3 %>%
  filter(Timepoint %in% c("T2", "T4")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Control") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T2_mean = mean(unlist(T2), na.rm = TRUE), T4_mean = mean(unlist(T4), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T2), unlist(T4), paired = TRUE)$p.value) 

rec_l3_t4_t3 <- l3 %>%
  filter(Timepoint %in% c("T4", "T3")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Recipeint") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T4_mean = mean(unlist(T4), na.rm = TRUE), T3_mean = mean(unlist(T3), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T4), unlist(T3), paired = TRUE)$p.value) 

rec_l3_t4_t3_c <- l3 %>%
  filter(Timepoint %in% c("T4", "T3")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Control") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T4_mean = mean(unlist(T4), na.rm = TRUE), T3_mean = mean(unlist(T3), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T4), unlist(T3), paired = TRUE)$p.value)

# bm

rec_bm_t3_t2 <- bm %>%
  filter(Timepoint %in% c("T2", "T3")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Recipeint") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T2_mean = mean(unlist(T2), na.rm = TRUE), T3_mean = mean(unlist(T3), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T2), unlist(T3), paired = TRUE)$p.value) 

rec_bm_t3_t2_c <- bm %>%
  filter(Timepoint %in% c("T2", "T3")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Control") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T2_mean = mean(unlist(T2), na.rm = TRUE), T3_mean = mean(unlist(T3), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T2), unlist(T3), paired = TRUE)$p.value) 

rec_bm_t4_t2 <- bm %>%
  filter(Timepoint %in% c("T2", "T4")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Recipeint") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T2_mean = mean(unlist(T2), na.rm = TRUE), T4_mean = mean(unlist(T4), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T2), unlist(T4), paired = TRUE)$p.value) 

rec_bm_t4_t2_c <- bm %>%
  filter(Timepoint %in% c("T2", "T4")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Control") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T2_mean = mean(unlist(T2), na.rm = TRUE), T4_mean = mean(unlist(T4), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T2), unlist(T4), paired = TRUE)$p.value)

rec_bm_t4_t3 <- bm %>%
  filter(Timepoint %in% c("T4", "T3")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Recipeint") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T4_mean = mean(unlist(T4), na.rm = TRUE), T3_mean = mean(unlist(T3), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T4), unlist(T3), paired = TRUE)$p.value) 

rec_bm_t4_t3_c <- bm %>%
  filter(Timepoint %in% c("T4", "T3")) %>%
  filter(!grepl("Transplant", sample)) %>%
  filter(Treatment == "Control") %>%
  select(ko, Timepoint, value) %>%
  group_by(Timepoint, ko) %>%
  summarise(value = list(value)) %>%
  spread(Timepoint, value) %>%
  group_by(ko) %>%
  mutate(T4_mean = mean(unlist(T4), na.rm = TRUE), T3_mean = mean(unlist(T3), na.rm = TRUE)) %>%
  mutate(p_value = t.test(unlist(T4), unlist(T3), paired = TRUE)$p.value) 

# summary

summaryTab <- rec_l2_t3_t2 %>%
  filter(p_value <= 0.05) %>%
  mutate(T2 = 1, T3 = 1, Change = ifelse(T3_mean > T2_mean, "Up", "Down"), Treatment = "Recipient", Level = "l2") %>%
  select(ko, Level, T2, T3, Treatment, Change, p_value) %>%
  merge(rec_l2_t3_t2_c %>%
          filter(p_value <= 0.05) %>%
          mutate(T2 = 1, T3 = 1, Change = ifelse(T3_mean > T2_mean, "Up", "Down"), Treatment = "Control", Level = "l2") %>%
          select(ko, Level, T2, T3, Treatment, Change, p_value), all = TRUE) %>%
  merge(rec_l2_t4_t2 %>%
          filter(p_value <= 0.05) %>%
          mutate(T2 = 1, T4 = 1, Change = ifelse(T4_mean > T2_mean, "Up", "Down"), Treatment = "Recipient", Level = "l2") %>%
          select(ko, Level, T2, T4, Treatment, Change, p_value), all = TRUE) %>%
  merge(rec_l2_t4_t2_c %>%
          filter(p_value <= 0.05) %>%
          mutate(T2 = 1, T4 = 1, Change = ifelse(T4_mean > T2_mean, "Up", "Down"), Treatment = "Control", Level = "l2") %>%
          select(ko, Level, T2, T4, Treatment, Change, p_value), all = TRUE) %>%
  merge(rec_l2_t4_t3 %>%
          filter(p_value <= 0.05) %>%
          mutate(T3 = 1, T4 = 1, Change = ifelse(T4_mean > T3_mean, "Up", "Down"), Treatment = "Recipient", Level = "l2") %>%
          select(ko, Level, T3, T4, Treatment, Change, p_value), all = TRUE) %>%
  merge(rec_l2_t4_t3_c %>%
          filter(p_value <= 0.05) %>%
          mutate(T3 = 1, T4 = 1, Change = ifelse(T4_mean > T3_mean, "Up", "Down"), Treatment = "Control", Level = "l2") %>%
          select(ko, Level, T3, T4, Treatment, Change, p_value), all = TRUE) %>%
  merge(rec_l3_t3_t2 %>%
          filter(p_value <= 0.05) %>%
          mutate(T3 = 1, T2 = 1, Change = ifelse(T3_mean > T2_mean, "Up", "Down"), Treatment = "Recipient", Level = "l3") %>%
          select(ko, Level, T3, T2, Treatment, Change, p_value), all = TRUE) %>%
  merge(rec_l3_t3_t2_c %>%
          filter(p_value <= 0.05) %>%
          mutate(T3 = 1, T2 = 1, Change = ifelse(T3_mean > T2_mean, "Up", "Down"), Treatment = "Control", Level = "l3") %>%
          select(ko, Level, T3, T2, Treatment, Change, p_value), all = TRUE) %>%
  merge(rec_l3_t4_t2 %>%
          filter(p_value <= 0.05) %>%
          mutate(T4 = 1, T2 = 1, Change = ifelse(T4_mean > T2_mean, "Up", "Down"), Treatment = "Recipient", Level = "l3") %>%
          select(ko, Level, T4, T2, Treatment, Change, p_value), all = TRUE) %>%
  merge(rec_l3_t4_t2_c %>%
          filter(p_value <= 0.05) %>%
          mutate(T4 = 1, T2 = 1, Change = ifelse(T4_mean > T2_mean, "Up", "Down"), Treatment = "Control", Level = "l3") %>%
          select(ko, Level, T4, T2, Treatment, Change, p_value), all = TRUE) %>%
  merge(rec_l3_t4_t3 %>%
          filter(p_value <= 0.05) %>%
          mutate(T4 = 1, T3 = 1, Change = ifelse(T4_mean > T3_mean, "Up", "Down"), Treatment = "Recipient", Level = "l3") %>%
          select(ko, Level, T4, T3, Treatment, Change, p_value), all = TRUE) %>%
  merge(rec_l3_t4_t3_c %>%
          filter(p_value <= 0.05) %>%
          mutate(T4 = 1, T3 = 1, Change = ifelse(T4_mean > T3_mean, "Up", "Down"), Treatment = "Control", Level = "l3") %>%
          select(ko, Level, T4, T3, Treatment, Change, p_value), all = TRUE) %>%
  merge(rec_bm_t3_t2 %>%
          filter(p_value <= 0.05) %>%
          mutate(T2 = 1, T3 = 1, Change = ifelse(T3_mean > T2_mean, "Up", "Down"), Treatment = "Recipient", Level = "bm") %>%
          select(ko, Level, T2, T3, Treatment, Change, p_value), all = TRUE) %>%
  merge(rec_bm_t3_t2_c %>%
          filter(p_value <= 0.05) %>%
          mutate(T2 = 1, T3 = 1, Change = ifelse(T3_mean > T2_mean, "Up", "Down"), Treatment = "Control", Level = "bm") %>%
          select(ko, Level, T2, T3, Treatment, Change, p_value), all = TRUE) %>%
  merge(rec_bm_t4_t2 %>%
          filter(p_value <= 0.05) %>%
          mutate(T2 = 1, T4 = 1, Change = ifelse(T4_mean > T2_mean, "Up", "Down"), Treatment = "Recipient", Level = "bm") %>%
          select(ko, Level, T2, T4, Treatment, Change, p_value), all = TRUE) %>%
  merge(rec_bm_t4_t2_c %>%
          filter(p_value <= 0.05) %>%
          mutate(T2 = 1, T4 = 1, Change = ifelse(T4_mean > T2_mean, "Up", "Down"), Treatment = "Control", Level = "bm") %>%
          select(ko, Level, T2, T4, Treatment, Change, p_value), all = TRUE) %>%
  merge(rec_bm_t4_t3 %>%
          filter(p_value <= 0.05) %>%
          mutate(T3 = 1, T4 = 1, Change = ifelse(T4_mean > T3_mean, "Up", "Down"), Treatment = "Recipient", Level = "bm") %>%
          select(ko, Level, T3, T4, Treatment, Change, p_value), all = TRUE) %>%
  merge(rec_bm_t4_t3_c %>%
          filter(p_value <= 0.05) %>%
          mutate(T3 = 1, T4 = 1, Change = ifelse(T4_mean > T3_mean, "Up", "Down"), Treatment = "Control", Level = "bm") %>%
          select(ko, Level, T3, T4, Treatment, Change, p_value), all = TRUE)

# plots

plot_L2_T1_All <- l2 %>%
  filter(Treatment %in% c("Recipeint", "Control", "Donor")) %>%
  filter(Timepoint == "T1") %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, color = Treatment)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L2_T1_All")

ggsave("L2_T1_all.pdf", plot_L2_T1_All, "pdf")

plot_L2_T2_All <- l2 %>%
  filter(Treatment %in% c("Recipeint", "Control", "Donor")) %>%
  filter(Timepoint == "T2") %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, color = Treatment)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L2_T2_All")

ggsave("L2_T2_all.pdf", plot_L2_T2_All, "pdf")

plot_L2_T3_All <- l2 %>%
  filter(Treatment %in% c("Recipeint", "Control")) %>%
  filter(Timepoint == "T3") %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, color = Treatment)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L2_T3_All")

ggsave("L2_T3_all.pdf", plot_L2_T3_All, "pdf")

plot_L2_T4_All <- l2 %>%
  filter(Treatment %in% c("Recipeint", "Control")) %>%
  filter(Timepoint == "T4") %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, color = Treatment)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L2_T4_All")

ggsave("L2_T4_all.pdf", plot_L2_T4_All, "pdf")

plot_L2_Sig_all <- l2 %>%
  filter(Treatment %in% c("Recipeint", "Control", "Donor")) %>%
  filter(ko %in% (summaryTab %>%
                    filter(Level == "l2") %>%
                    .$ko)) %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, fill = Treatment, color = Timepoint)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L2_Significant") +
  scale_color_manual(values = rep("black", 4)) +
  guides(color = FALSE)

ggsave("L2_Sig_all.pdf", plot_L2_Sig_all, "pdf")

plot_l2_Sig_high <- l2 %>%
  filter(Treatment %in% c("Recipeint", "Control", "Donor")) %>%
  filter(ko %in% (summaryTab %>%
                    filter(Level == "l2") %>%
                    .$ko %>%
                    unique() %>%
                    .[c(1,4,10,21,24,25,26,28)])) %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, fill = Treatment, color = Timepoint)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L2_Significant_High") +
  scale_color_manual(values = rep("black", 4)) +
  guides(color = FALSE)

ggsave("L2_Sig_high.pdf", plot_l2_Sig_high, "pdf")

plot_l2_Sig_Med <- l2 %>%
  filter(Treatment %in% c("Recipeint", "Control", "Donor")) %>%
  filter(ko %in% (summaryTab %>%
                    filter(Level == "l2") %>%
                    .$ko %>%
                    unique() %>%
                    .[c(11,12,13,14,18,20,22,23,30)])) %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, fill = Treatment, color = Timepoint)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L2_Significant_Med") +
  scale_color_manual(values = rep("black", 4)) +
  guides(color = FALSE)

ggsave("L2_Sig_Med.pdf", plot_l2_Sig_Med, "pdf")

plot_l2_Sig_low <- l2 %>%
  filter(Treatment %in% c("Recipeint", "Control", "Donor")) %>%
  filter(ko %in% (summaryTab %>%
                    filter(Level == "l2") %>%
                    .$ko %>%
                    unique() %>%
                    .[c(2,3,5,6,7,8,9,15,16,17,19,27,29)])) %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, fill = Treatment, color = Timepoint)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L2_Significant_Low") +
  scale_color_manual(values = rep("black", 4)) +
  guides(color = FALSE)

ggsave("L2_Sig_Low.pdf", plot_l2_Sig_low, "pdf")

# l3

plot_L3_T2_All <- l3 %>%
  filter(Treatment %in% c("Recipeint", "Control")) %>%
  filter(Timepoint == "T2") %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, color = Treatment)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L3_T2_All")

ggsave("L3_T2_all.pdf", plot_L3_T2_All, "pdf")

plot_L3_T3_All <- l3 %>%
  filter(Treatment %in% c("Recipeint", "Control")) %>%
  filter(Timepoint == "T3") %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, color = Treatment)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L3_T3_All")

ggsave("L3_T3_all.pdf", plot_L3_T3_All, "pdf")

plot_L3_T4_All <- l3 %>%
  filter(Treatment %in% c("Recipeint", "Control")) %>%
  filter(Timepoint == "T4") %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, color = Treatment)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L3_T4_All")

ggsave("L3_T4_all.pdf", plot_L3_T4_All, "pdf")


l3 %>%
  filter(ko %in% (summaryTab %>%
                    filter(Level == "l3") %>%
                    .$ko %>%
                    unique())) %>%
  group_by(ko) %>%
  summarise(min = min(value)) %>%
  arrange(desc(min)) -> l3MinTab

cutoffPos <- seq(1, nrow(l3MinTab), length.out = 7)

plot_L3_Sig_1 <- l3 %>%
  filter(Treatment %in% c("Recipeint", "Control", "Donor")) %>%
  filter(ko %in% l3MinTab$ko[cutoffPos[1] : cutoffPos[2]]) %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, fill = Treatment, color = Timepoint)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L3_Significant_1") +
  scale_color_manual(values = rep("black", 4)) +
  guides(color = FALSE)

ggsave("L3_Sig_1.pdf", plot_L3_Sig_1, "pdf")

plot_L3_Sig_2 <- l3 %>%
  filter(Treatment %in% c("Recipeint", "Control", "Donor")) %>%
  filter(ko %in% l3MinTab$ko[(cutoffPos[2] + 1) : cutoffPos[3]]) %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, fill = Treatment, color = Timepoint)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L3_Significant_2") +
  scale_color_manual(values = rep("black", 4)) +
  guides(color = FALSE)

ggsave("L3_Sig_2.pdf", plot_L3_Sig_2, "pdf")

plot_L3_Sig_3 <- l3 %>%
  filter(Treatment %in% c("Recipeint", "Control", "Donor")) %>%
  filter(ko %in% l3MinTab$ko[(cutoffPos[3] + 1) : cutoffPos[4]]) %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, fill = Treatment, color = Timepoint)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L3_Significant_3") +
  scale_color_manual(values = rep("black", 4)) +
  guides(color = FALSE)

ggsave("L3_Sig_3.pdf", plot_L3_Sig_3, "pdf")

plot_L3_Sig_4 <- l3 %>%
  filter(Treatment %in% c("Recipeint", "Control", "Donor")) %>%
  filter(ko %in% l3MinTab$ko[(cutoffPos[4] + 1) : cutoffPos[5]]) %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, fill = Treatment, color = Timepoint)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L3_Significant_4") +
  scale_color_manual(values = rep("black", 4)) +
  guides(color = FALSE)

ggsave("L3_Sig_4.pdf", plot_L3_Sig_4, "pdf")

plot_L3_Sig_5 <- l3 %>%
  filter(Treatment %in% c("Recipeint", "Control", "Donor")) %>%
  filter(ko %in% l3MinTab$ko[(cutoffPos[5] + 1) : cutoffPos[6]]) %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, fill = Treatment, color = Timepoint)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L3_Significant_5") +
  scale_color_manual(values = rep("black", 4)) +
  guides(color = FALSE)

ggsave("L3_Sig_5.pdf", plot_L3_Sig_5, "pdf")

plot_L3_Sig_6 <- l3 %>%
  filter(Treatment %in% c("Recipeint", "Control", "Donor")) %>%
  filter(ko %in% l3MinTab$ko[(cutoffPos[6] + 1) : cutoffPos[7]]) %>%
  ggplot() +
  geom_boxplot(aes(x = ko, y = value, fill = Treatment, color = Timepoint)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("L3_Significant_6") +
  scale_color_manual(values = rep("black", 4)) +
  guides(color = FALSE)

ggsave("L3_Sig_6.pdf", plot_L3_Sig_6, "pdf")

summaryTab







