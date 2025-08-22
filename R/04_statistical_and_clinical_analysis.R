# **Step 4: Advanced Statistical and Clinical Analysis**
## Script: `R/04_statistical_and_clinical_analysis.R`
## Source Files: `PRC.docx`, `power.docx`, `correlationAnalysis.R`, `cforestAnalysisWithMeta.R`
## Description: This script brings together several advanced analyses:
### - Principal Response Curve (PRC): To visualize the effect of the OMT on the community over time relative to the control (Figure 4 in the paper).
### - PERMANOVA Power Analysis: A retrospective analysis to confirm the statistical power of the study design.
### - Correlation Analysis: Uses Pearson correlation to find significant associations between microbial taxa (ASVs/Species) and clinical parameters (BOP, PPD).
### - Predictive Modeling: Employs Random Forest to identify which taxa are most predictive of clinical outcomes.
## Inputs:`phyloseq_object.rds`, clinical data files (`Ergebnisse_BOP_PPD_summarized_BOP_PPD_3.txt`).
## Outputs:PRC plots, correlation heatmaps, variable importance plots from Random Forest, and associated statistical tables.

library("phyloseq")
library("ggplot2")
library("dplyr")

#Principal Response Curve (PRC) was conducted using the genepiper GUI-R pack

# Constrained Ordinations
# how environmental variables (BOP and PPD) are associated with the changes in community composition
# add the BOP and PPD to the phyloseq data
sample_data(ps1) %>%
  mutate(Sample = rownames(.))%>%
  full_join(BOP_dat, by =  "Sample")%>%
  left_join(PPD_dat, by =  "Sample") %>%
  left_join(plaque, by =  "Sample") %>%
  column_to_rownames(var = "Sample") ->sample_data(ps1) 
# Remove data points with missing metadata; remove donor samples; split control and recipient
ps.c <- ps1 %>%  subset_samples(Treatment=="Control")%>% prune_taxa(taxa_sums(.) > 0, .)
ps.r <- ps1 %>%  subset_samples(Treatment=="Recipient")%>% prune_taxa(taxa_sums(.) > 0, .)
bray.c <- phyloseq::distance(physeq = ps.c, method = "bray")
bray.r <- phyloseq::distance(physeq = ps.r, method = "bray")
# CAP ordinate
cap_ord.c <- ordinate(physeq = ps.c, method = "CAP", distance = bray.c, formula = ~ BOP+PPD+plaque)
cap_ord.r <- ordinate(physeq = ps.r, method = "CAP", distance = bray.r, formula = ~ BOP+PPD+plaque)

# CAP plot
cap_plot.c <- plot_ordination(physeq = ps.c, ordination = cap_ord.c, color = "Timepoint", axes = c(1,2)) + 
  geom_point(aes(colour = Timepoint), alpha = 0.7, size = 4) +
  geom_text(aes(label=sample_names(ps.c)),nudge_y = 0.1)

cap_plot.r <- plot_ordination(physeq = ps.r, ordination = cap_ord.r, color = "Timepoint", axes = c(1,2)) + 
  geom_point(aes(colour = Timepoint), alpha = 0.7, size = 4) +
  geom_text(aes(label=sample_names(ps.r)),nudge_y = 0.1)

arrowmat.r <- vegan::scores(cap_ord.r, display = "bp")
arrowdf.r <- data.frame(labels = rownames(arrowmat.r), arrowmat.r)

arrowmat.c <- vegan::scores(cap_ord.c, display = "bp")
arrowdf.c <- data.frame(labels = rownames(arrowmat.c), arrowmat.c)



arrow_map <- aes(xend = CAP1, yend = CAP2, x = 0, y = 0, shape = NULL, color = NULL, label = labels)
label_map <- aes(x = 1.3 * CAP1, y = 1.3 * CAP2, shape = NULL, color = NULL, label = labels)
arrowhead = arrow(length = unit(0.02, "npc"))

cap_plot.c + geom_segment(mapping = arrow_map, size = .5, data = arrowdf.c, color = "black", arrow = arrowhead) + 
  geom_text(mapping = label_map, size = 4,  data = arrowdf.c, show.legend = FALSE )+
  ggtitle("Canonical analysis of principal coordinates for Control samples")

cap_plot.r + geom_segment(mapping = arrow_map, size = .5, data = arrowdf.r, color = "black", arrow = arrowhead) + 
  geom_text(mapping = label_map, size = 4,  data = arrowdf.r, show.legend = FALSE )+
  ggtitle("Canonical analysis of principal coordinates for Recipient samples")

#PERMANOVA
ps.wuf<-phyloseq::distance(ps.ra, method="wunifrac")
metadata <- as(sample_data(ps.ra), "data.frame")
metadata$TT= paste0(metadata$Treatment, "_", metadata$Timepoint)
colnames(metadata)
#[1] "Dog"       "Timepoint" "Treatment" "Date"      "Gender"    "Week"      "TT"       
library("vegan")
long.factor=c("TT","Treatment","Timepoint","Gender")
results_wuf <- lapply(long.factor, function(x){
  form <- as.formula(paste("ps.wuf", x, sep="~")) 
  z <- adonis(form, data = metadata, permutations=9999)
  return(as.data.frame(z$aov.tab)) })
names(results_wuf) <- long.factor
results_wuf <- do.call(rbind, results_wuf)
write.csv(results_wuf, "~/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/adonis_wuf.csv")
pairwise.adonis.dm(ps.wuf,metadata$TT)

adonis2(ps.wuf ~ Treatment+Timepoint+Gender, data = metadata, by = "margin")

adonis(ps.r.wuf ~ Timepoint,data = (as(sample_data(ps.r.ra), "data.frame")), permutations=9999)
results_wuf.r <-pairwise.adonis.dm(ps.r.wuf,(as(sample_data(ps.r.ra), "data.frame")$Timepoint))
write.csv(results_wuf.r, "~/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/adonis_wuf_r.csv")

adonis(ps.c.wuf ~ Timepoint,data = (as(sample_data(ps.c.ra), "data.frame")), permutations=9999)
results_wuf.c <-pairwise.adonis.dm(ps.c.wuf,(as(sample_data(ps.c.ra), "data.frame")$Timepoint))
write.csv(results_wuf.c, "~/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/adonis_wuf_c.csv")

# remove self-comparisons
wu.m =  melt(as.matrix(ps.r.wuf)) %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)
# get sample data (S4 error OK and expected)
sd <-metadata %>%
  tibble::rownames_to_column(var="Sample") %>%
  filter(Treatment=="Recipient") %>%
  select(Sample, Timepoint) %>%
  mutate_if(is.factor,as.character)

# combined distances with sample data
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")

colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

# plot
ggplot(wu.sd, aes(x = Type2, y = value)) +
  theme_bw() +
  geom_point() +
  geom_boxplot(aes(color = ifelse(Type1 == Type2, "red", "black"))) +
  scale_color_identity() +
  facet_wrap(~ Type1, scales = "free_x") +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  ggtitle(paste0("Weighted Unifrac distance metric"))

# Correlation analysis with the clinical metadata
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(party)
library(partykit)
library(caret)
library(randomForest)
library(RColorBrewer)
library(Hmisc)

# Import BOP data
bopDat <- read_tsv("Ergebnisse_BOP_PPD_summarized_BOP_PPD_3.txt") %>%
  filter(Tooth == "BOP") %>%
  dplyr::rename(sampleID = Dog) %>%
  select(sampleID, BOP)

# Import PPD data
ppdDat <- read_tsv("Ergebnisse_BOP_PPD_summarized_BOP_PPD_3.txt") %>%
  filter(Tooth == "PPD") %>%
  dplyr::rename(sampleID = Dog) %>%
  select(sampleID, PPD)

# Import the microbiota data (OTU|ASV abundance table)
otuFilePath <- "ps1_ASVtable.txt"
otuTab <- read_tsv(otuFilePath)

# Reformatting the microbiota abundance table
otuTab <- otuTab %>%
  as.data.frame() %>%
  .[-c(1,2),] %>%
  `colnames<-`(.[1,]) %>%
  .[-1,] %>%
  remove_rownames() %>%
  column_to_rownames("ASV") %>%
  t()

# Combining with the clinical data (BOP & PPD)
allDf <- otuTab %>%
  apply(., c(1, 2), as.numeric) %>%
  data.frame(stringsAsFactors = FALSE) %>%
  t() %>%
  data.frame(stringsAsFactors = FALSE) %>%
  rownames_to_column("AVS") %>%
  mutate_if(is.numeric, funs(. / sum(.))) %>%
  gather(key = "sampleID", value = "value", -AVS) %>%
  spread(key = "AVS", value = "value") %>%
  .[1:77,] %>%
  left_join(bopDat %>%
              mutate(sampleID = make.names(sampleID))) %>%
  left_join(ppdDat %>%
              mutate(sampleID = make.names(sampleID))) %>%
  select(-sampleID)

### Pearson Correlation Analysis at ASVs 
allDf %>%
  mutate_all(funs(as.numeric)) %>%
  mutate_all(funs(ifelse(. == 0, NA, .))) %>%
  as.matrix() %>%
  Hmisc::rcorr(type = "pearson") -> corr

# Extracting correlation analysis results
combn(colnames(corr$r), 2) %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y")) -> combTab
# Collecting the correlations r values
corr$r %>%
  as.data.frame() %>%
  rownames_to_column("x") %>%
  gather(key = "y", value = "r", -x) -> corrR
# Collecting the number of observations (n) used in analyzing each pair of variables
corr$n %>%
  as.data.frame() %>%
  rownames_to_column("x") %>%
  gather(key = "y", value = "n", -x) -> corrN
# Collecting the asymptotic P-values
corr$P %>%
  as.data.frame() %>%
  rownames_to_column("x") %>%
  gather(key = "y", value = "p", -x) -> corrP
# Combining the correlation results (r-n-P) into one table
combTab %>%
  left_join(corrR, by = c("x" = "x", "y" = "y")) %>%
  left_join(corrP, by = c("x" = "x", "y" = "y")) %>%
  left_join(corrN, by = c("x" = "x", "y" = "y")) -> corrTab

# Filtering for significant correlation results 
# Bonferroni correction for multivariate correlation??
corrTab %>%
  filter(n > 10) %>%
  filter(p <= 0.05) %>%
  filter(abs(r) >= 0.5) %>%
  filter(x %in% c("BOP", "PPD") | y %in% c("BOP", "PPD")) -> selAsvTab

# Plotting the correlation results
labY <- unique(selAsvTab[, "y"])
labX <- unique(selAsvTab[, "x"])[-4]

corrTab %>%
  filter(x %in% labX) %>%
  filter(y %in% labY) %>%
  ggplot() +
  geom_tile(aes(x = x, y = y), fill = "white", color = "black") +
  coord_equal() +
  geom_point(aes(x = x, y = y, color = r, size = r))


### Pearson Correlation Analysis at the Species level 
# agglomerate species

taxTab3 <- read_csv("species_RDP_multiple_renamed.csv")

taxTab4 <- apply(taxTab3, 1, function(x){
  if(any(is.na(x))){
    naPos <- which(is.na(x))
    minNaPos <- min(naPos)
    lastTaxa <- x[minNaPos - 1]
    
    fillTaxa <- paste0("unknown_", lastTaxa)
    x[as.numeric(naPos)] <- fillTaxa
  }
  
  if(!grepl("unknown_", x[8])){
    x[8] <- paste(x[7], x[8], sep = "_")
  }
  
  x[8] <- strsplit(x[8], "/")[[1]][1]
  
  x
}) %>%
  t() %>%
  data.frame(stringsAsFactors = FALSE)

aggTab2 <- otuTab %>%
  t() %>%
  data.frame(stringsAsFactors = FALSE) %>%
  rownames_to_column("ASV") %>%
  left_join(taxTab4 %>%
              select(X1, Species), by = c("ASV" = "X1")) %>%
  group_by(Species) %>%
  select(-ASV) %>%
  summarise_all(funs(sum(as.numeric(.)))) %>%
  ungroup() %>%
  mutate_if(is.numeric, funs(. / sum(.))) %>%
  as.data.frame() %>%
  `rownames<-`(make.names(.$Species)) %>%
  select(-Species) %>%
  t()

speciesDf <- aggTab2 %>%
  as.data.frame() %>%
  rownames_to_column("sampleID") %>%
  .[1:77,] %>%
  left_join(bopDat %>%
              mutate(sampleID = make.names(sampleID))) %>%
  left_join(ppdDat %>%
              mutate(sampleID = make.names(sampleID))) %>%
  select(-sampleID)

speciesDf %>%
  mutate_all(funs(as.numeric)) %>%
  mutate_all(funs(ifelse(. == 0, NA, .))) %>%
  as.matrix() %>%
  Hmisc::rcorr(type = "pearson") -> corrS

combn(colnames(corrS$r), 2) %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y")) -> combTabS

corrS$r %>%
  as.data.frame() %>%
  rownames_to_column("x") %>%
  gather(key = "y", value = "r", -x) -> corrRS

corrS$n %>%
  as.data.frame() %>%
  rownames_to_column("x") %>%
  gather(key = "y", value = "n", -x) -> corrNS

corrS$P %>%
  as.data.frame() %>%
  rownames_to_column("x") %>%
  gather(key = "y", value = "p", -x) -> corrPS

combTabS %>%
  left_join(corrRS, by = c("x" = "x", "y" = "y")) %>%
  left_join(corrPS, by = c("x" = "x", "y" = "y")) %>%
  left_join(corrNS, by = c("x" = "x", "y" = "y")) -> corrTabS

corrTabS %>%
  filter(n > 10) %>%
  filter(p <= 0.05) %>%
  filter(abs(r) >= 0.5) %>%
  filter(x %in% c("BOP", "PPD") | y %in% c("BOP", "PPD")) -> selAsvTabS

labYS <- unique(selAsvTabS[, "y"])
labXS <- unique(selAsvTabS[, "x"])[-2]

corrTabS %>%
  filter(x %in% labXS) %>%
  filter(y %in% labYS) %>%
  ggplot() +
  geom_tile(aes(x = x, y = y), fill = "white", color = "black") +
  coord_equal() +
  geom_point(aes(x = x, y = y, color = r, size = r))

# Machine learning model building, Random forest analysis at ASVs
# Building a RF model to predict periodontal treatment outcomes (health/ disease in future timepoints) based on baseline (current) microbiota profile.
# Import OTU table
otuFilePath <- "ps1_ASVtable.txt"
otuTab <- read_tsv(otuFilePath)

otuTab <- otuTab %>%
  as.data.frame() %>%
  .[-c(1,2),] %>%
  `colnames<-`(.[1,]) %>%
  .[-1,] %>%
  remove_rownames() %>%
  column_to_rownames("ASV") %>%
  t()

# Import Metadata table
metaTab <- rownames(otuTab) %>%
  strsplit("-") %>%
  lapply(`[[`, 2) %>%
  unlist() %>%
  cbind(rownames(otuTab), .) %>%
  as.data.frame() %>%
  mutate(time = c(rep(c(0, 2, 4, 14), 4), rep(0, 5), rep(c(0, 2, 4, 14), 14), rep(2, 9)))

write.table(metaTab, "meta.tab", sep = ",", col.names = FALSE, row.names = FALSE)

# Create a Sample table
subTab <- levels(metaTab[,2]) %>%
  cbind(c(3,1,3,1,2,3,1,3,1,3,1,3,1,3,1,3,3,1,1))

subTab <- as.data.frame(subTab)

subTab[, 2] <- factor(subTab[, 2])

levels(subTab[, 2]) <- c("Control", "Donor", "Recepient")

colnames(subTab) <- c("subjectID", "role")

selSam <- metaTab[metaTab$time %in% c(0, 2), "V1"]

# Prepare train and test dataset
otuTab[rownames(otuTab) %in% selSam, ] %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(funs(as.numeric(as.character(.)))) %>%
  mutate_all(funs(. / sum(.))) %>%
  t() %>%
  `colnames<-`(colnames(otuTab)) %>%
  as.data.frame() %>%
  mutate(state = c(rep("D", 8), rep("H", 5), rep("D", 28), rep("H", 9))) -> dat

dat$state <- factor(dat$state)

set.seed(123)

smp_size <- floor(0.7 * nrow(dat))
smp_idx <- sample(1 : nrow(dat), smp_size)

trainSet <- dat[smp_idx,]
testSet <- dat[-smp_idx,]

# Training random forest model with trainset
forest1 <- cforest(state ~ ., data = trainSet, ntree = 500, mtry = 5)
pred1 <- Predict(forest1, newdata = testSet)

confusMat1 <- confusionMatrix(pred1, testSet$state)

varimp(forest1)

forestTest <- function(formula, train, test, ntrees, mtrys){
  res <- c()
  
  set.seed(123)
  
  for(i in ntrees){
    for(j in mtrys){
      forest <- cforest(as.formula(formula), data = train, ntree = i, mtry = j)
      pred <- Predict(forest, newdata = test)
      confusMat <- confusionMatrix(pred, test$state)
      res <- c(res, i, j, confusMat$overall[[1]])
    }
  }
  matrix(res, ncol = 3, byrow = TRUE) %>%
    `colnames<-`(c("ntree", "mtry", "accurancy"))
}


forTest1 <- forestTest("state ~ .", trainSet, testSet, c(100, 200, 500, 1000, 2000, 5000), c(3, 5, 10, 15, 30, 100))

matrix(forTest1, ncol = 3, byrow = TRUE)

set.seed(123)
forest2 <- cforest(state ~ ., data = trainSet, ntree = 500, mtry = 15) # set ntree 500 and mtry 30 with acurrancy 80%
impTab2 <- varimp(forest2)  # get importances

### Trial 3, pre-filter the ASV abundance table by eliminating the ASVs with importance < 1
selAsv3 <- names(impTab2)[abs(impTab2) >= 1] # select ASV with importance over 1

formula3 <- "state ~" %>%
  paste(paste(selAsv3, collapse = "+"))

set.seed(123)
forest3 <- cforest(as.formula(formula3), data = trainSet, ntree = 500, mtry = 15)
pred3 <- Predict(forest3, newdata = testSet)

confusMat3 <- confusionMatrix(pred3, testSet$state) # accurancy 100%

impTab3 <- varimp(forest3)

### Trial 4, progressive removal of ASVs with importance < 1
selAsv4 <- names(impTab3)[abs(impTab3) >= 1]

formula4 <- "state ~" %>%
  paste(paste(selAsv4, collapse = "+"))

set.seed(123)
forest4 <- cforest(as.formula(formula4), data = trainSet, ntree = 500, mtry = 15)
pred4 <- Predict(forest4, newdata = testSet)

confusMat4 <- confusionMatrix(pred4, testSet$state) # Prediction accuracy: 93%
varimp(forest4)

### Utilizing Trial 3 RF model for prediction
taxTab <- read_tsv("ps1_taxtable.txt")

impTaxTab <- names(impTab3) %>%
  cbind(impTab3) %>%
  as.tibble() %>%
  `colnames<-`(c("ID", "Importance")) %>%
  left_join(taxTab, by = c("ID" = "Kingdom"))

# Predict 14 week data
selSam2 <- metaTab[metaTab$time %in% c(14), "V1"]

dat_14w <- otuTab[rownames(otuTab) %in% selSam2,] %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(funs(as.numeric(as.character(.)))) %>%
  mutate_all(funs(. / sum(.))) %>%
  t() %>%
  `colnames<-`(colnames(otuTab)) %>%
  as.data.frame()

Predict(forest3, newdata = dat_14w)

# Predict 4 week data
selSam3 <- metaTab[metaTab$time %in% 4, "V1"]

dat_4w <- otuTab[rownames(otuTab) %in% selSam3,] %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(funs(as.numeric(as.character(.)))) %>%
  mutate_all(funs(. / sum(.))) %>%
  t() %>%
  `colnames<-`(colnames(otuTab)) %>%
  as.data.frame()

Predict(forest3, newdata = dat_4w)

# Predict 2 week data
selSam4 <- metaTab[metaTab$time %in% 2, "V1"]

dat_2w <- otuTab[rownames(otuTab) %in% selSam4,] %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(funs(as.numeric(as.character(.)))) %>%
  mutate_all(funs(. / sum(.))) %>%
  t() %>%
  `colnames<-`(colnames(otuTab)) %>%
  as.data.frame()

Predict(forest3, newdata = dat_2w)

### Plot individual ASV
samMetaTab <- metaTab %>%
  `colnames<-`(c("ID", "Dog", "time")) %>%
  filter(Dog != "E") %>%
  filter(time != 0) %>%
  .[1:54,]

plotTab <- otuTab[rownames(otuTab) %in% samMetaTab$ID,] %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(funs(as.numeric(as.character(.)))) %>%
  mutate_all(funs(. / sum(.))) %>%
  t() %>%
  `colnames<-`(colnames(otuTab)) %>%
  as.data.frame() %>%
  select(one_of(impTaxTab$ID)) %>%
  rownames_to_column("ID") %>%
  left_join(samMetaTab) %>%
  select(-ID) %>%
  left_join(subTab, by = c("Dog" = "subjectID")) %>%
  gather(key = "ASV", value = "value", -Dog, -time, -role)

gg <- plotTab %>%
  ggplot() +
  geom_line(aes(x = time, y = value, group = Dog, color = role)) +
  facet_wrap(~ ASV, scales = "free_y") +
  ylab("Relative Abundance") +
  xlab("Time (weeks)")

ggsave("impAsv.pdf", gg, device = "pdf")

# Machine learning model building, Random forest analysis at the Genus level
##### agglomerate to genus levels

aggTab <- otuTab %>%
  t() %>%
  data.frame(stringsAsFactors = FALSE) %>%
  rownames_to_column("ASV") %>%
  left_join(taxTab %>%
              select(Kingdom, Species), by = c("ASV" = "Kingdom")) %>%
  group_by(Species) %>%
  select(-ASV) %>%
  summarise_all(funs(sum(as.numeric(.)))) %>%
  ungroup() %>%
  mutate_if(is.numeric, funs(. / sum(.))) %>%
  as.data.frame() %>%
  `rownames<-`(make.names(.$Species)) %>%
  select(-Species) %>%
  t()

selSam <- metaTab[metaTab$time %in% c(0, 2), "V1"]

aggTab[rownames(aggTab) %in% make.names(selSam), ] %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(funs(as.numeric(as.character(.)))) %>%
  mutate_all(funs(. / sum(.))) %>%
  t() %>%
  `colnames<-`(colnames(aggTab)) %>%
  as.data.frame() %>%
  mutate(state = c(rep("D", 8), rep("H", 5), rep("D", 28), rep("H", 9))) -> dat

dat$state <- factor(dat$state)

set.seed(123)

smp_size <- floor(0.7 * nrow(dat))
smp_idx <- sample(1 : nrow(dat), smp_size)

trainSet <- dat[smp_idx,]
testSet <- dat[-smp_idx,]

forTest2 <- forestTest("state ~ . - NA.", trainSet, testSet, c(100, 200, 500, 1000, 2000, 5000), c(3, 5, 10, 15, 30, 100))

forest5 <- cforest(state ~ . - NA., trainSet, ntree = 500, mtry = 15)
impTab5 <- varimp(forest5)

selAsv6 <- names(impTab5)[abs(impTab5) >= 1]

formula6 <- "state ~" %>%
  paste(paste(selAsv6, collapse = "+"))

set.seed(123)
forest6 <- cforest(as.formula(formula6), data = trainSet, ntree = 500, mtry = 15)
pred6 <- Predict(forest6, newdata = testSet)

confusMat6 <- confusionMatrix(pred6, testSet$state) # accurancy 93%
impTab6 <- varimp(forest6)

taxTab2 <- taxTab
taxTab2$Species <- make.names(taxTab2$Species)


impTaxTab6 <- cbind(names(impTab6), impTab6) %>%
  `colnames<-`(c("Species", "Importance")) %>%
  as.data.frame() %>%
  left_join(taxTab2) %>%
  select(Kingdom, Phylum, Class, Order, Family, Genus, Species, Importance)

samMetaTab6 <- metaTab %>%
  `colnames<-`(c("ID", "Dog", "time")) %>%
  filter(Dog != "E") %>%
  filter(time != 0) %>%
  .[1:54,]

samMetaTab6$ID <- samMetaTab6$ID %>%
  as.character() %>%
  make.names()

plotTab6 <- aggTab[rownames(aggTab) %in% samMetaTab6$ID,] %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(funs(as.numeric(as.character(.)))) %>%
  mutate_all(funs(. / sum(.))) %>%
  t() %>%
  `colnames<-`(colnames(aggTab)) %>%
  as.data.frame() %>%
  select(one_of(names(impTab6))) %>%
  rownames_to_column("ID") %>%
  left_join(samMetaTab6) %>%
  select(-ID) %>%
  left_join(subTab, by = c("Dog" = "subjectID")) %>%
  gather(key = "Species", value = "value", -Dog, -time, -role)

gg6 <- plotTab6 %>%
  ggplot() +
  geom_line(aes(x = time, y = value, group = Dog, color = role)) +
  facet_wrap(~ Species, scales = "free_y") +
  ylab("Relative Abundance") +
  xlab("Time (weeks)")

ggsave("impSpecies.pdf", gg6, device = "pdf")

# Machine learning model building, Random forest analysis at the Species level
### agglomerate to species levels
taxTab3 <- read_csv("species_RDP_multiple_renamed.csv")

taxTab4 <- apply(taxTab3, 1, function(x){
  if(any(is.na(x))){
    naPos <- which(is.na(x))
    minNaPos <- min(naPos)
    lastTaxa <- x[minNaPos - 1]
    
    fillTaxa <- paste0("unknown_", lastTaxa)
    x[as.numeric(naPos)] <- fillTaxa
  }
  
  if(!grepl("unknown_", x[8])){
    x[8] <- paste(x[7], x[8], sep = "_")
  }
  
  x[8] <- strsplit(x[8], "/")[[1]][1]
  
  x
}) %>%
  t() %>%
  data.frame(stringsAsFactors = FALSE)

aggTab2 <- otuTab %>%
  t() %>%
  data.frame(stringsAsFactors = FALSE) %>%
  rownames_to_column("ASV") %>%
  left_join(taxTab4 %>%
              select(X1, Species), by = c("ASV" = "X1")) %>%
  group_by(Species) %>%
  select(-ASV) %>%
  summarise_all(funs(sum(as.numeric(.)))) %>%
  ungroup() %>%
  mutate_if(is.numeric, funs(. / sum(.))) %>%
  as.data.frame() %>%
  `rownames<-`(make.names(.$Species)) %>%
  select(-Species) %>%
  t()

selSam <- metaTab[metaTab$time %in% c(0, 2), "V1"]

aggTab2[rownames(aggTab2) %in% make.names(selSam), ] %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(funs(as.numeric(as.character(.)))) %>%
  mutate_all(funs(. / sum(.))) %>%
  t() %>%
  `colnames<-`(colnames(aggTab2)) %>%
  as.data.frame() %>%
  mutate(state = c(rep("D", 8), rep("H", 5), rep("D", 28), rep("H", 9))) -> dat

dat$state <- factor(dat$state)

set.seed(123)

smp_size <- floor(0.7 * nrow(dat))
smp_idx <- sample(1 : nrow(dat), smp_size)

trainSet <- dat[smp_idx,]
testSet <- dat[-smp_idx,]

forTest3 <- forestTest("state ~ . - NA.", trainSet, testSet, c(100, 200, 500, 1000, 2000, 5000), c(3, 5, 10, 15, 30, 100))

forest7 <- cforest(state ~ ., trainSet, ntree = 500, mtry = 10)
impTab7 <- varimp(forest7)

selAsv7 <- names(impTab7)[abs(impTab7) >= 1]

formula7 <- "state ~" %>%
  paste(paste(selAsv7, collapse = "+"))

set.seed(123)
forest8 <- cforest(as.formula(formula7), data = trainSet, ntree = 500, mtry = 10)
pred8 <- Predict(forest8, newdata = testSet)

confusMat8 <- confusionMatrix(pred8, testSet$state) # accurancy 100%
impTab8 <- varimp(forest8)

selAsv9 <- names(impTab8)[abs(impTab8) >= 1]

formula9 <- "state ~" %>%
  paste(paste(selAsv9, collapse = "+"))

set.seed(123)
forest9 <- cforest(as.formula(formula9), data = trainSet, ntree = 500, mtry = 10)
pred9 <- Predict(forest9, newdata = testSet)

confusMat9 <- confusionMatrix(pred9, testSet$state) # accurancy 100%
impTab9 <- varimp(forest9)

samMetaTab6 <- metaTab %>%
  `colnames<-`(c("ID", "Dog", "time")) %>%
  filter(Dog != "E") %>%
  filter(time != 0) %>%
  .[1:54,]

samMetaTab6$ID <- samMetaTab6$ID %>%
  as.character() %>%
  make.names()

plotTab9 <- aggTab2[rownames(aggTab2) %in% samMetaTab6$ID,] %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(funs(as.numeric(as.character(.)))) %>%
  mutate_all(funs(. / sum(.))) %>%
  t() %>%
  `colnames<-`(colnames(aggTab2)) %>%
  as.data.frame() %>%
  select(one_of(names(impTab9))) %>%
  rownames_to_column("ID") %>%
  left_join(samMetaTab6) %>%
  select(-ID) %>%
  left_join(subTab, by = c("Dog" = "subjectID")) %>%
  gather(key = "Species", value = "value", -Dog, -time, -role)

gg9 <- plotTab9 %>%
  ggplot() +
  geom_line(aes(x = time, y = value, group = Dog, color = role)) +
  facet_wrap(~ Species, scales = "free_y") +
  ylab("Relative Abundance") +
  xlab("Time (weeks)")

ggsave("impSpecies.pdf", gg9, device = "pdf")

##### random forest package
set.seed(123)
rf1 <- randomForest(state ~ ., data = dat, subset = smp_idx)
plot(rf1)

selSpecies <- rownames(importance(rf1))[importance(rf1) != 0]
rfFormula2 <- "state ~" %>%
  paste(paste(selSpecies, collapse = "+"))

rf2 <- randomForest(as.formula(rfFormula2), data = dat, subset = smp_idx)
plot(rf2)

selSpecies2 <- rownames(importance(rf2))[importance(rf2) != 0]
rfFormula3 <- "state ~" %>%
  paste(paste(selSpecies2, collapse = "+"))

rf3 <- randomForest(as.formula(rfFormula3), data = dat, subset = smp_idx)
plot(rf3)

selSpecies3 <- rownames(importance(rf3))[importance(rf3) != 0]
rfFormula4 <- "state ~" %>%
  paste(paste(selSpecies3, collapse = "+"))

rf4 <- randomForest(as.formula(rfFormula4), data = dat, subset = smp_idx)
plot(rf4)

selSpecies4 <- rownames(importance(rf4))[importance(rf4) != 0]
rfFormula5 <- "state ~" %>%
  paste(paste(selSpecies4, collapse = "+"))

rf5 <- randomForest(as.formula(rfFormula5), data = dat, subset = smp_idx)
plot(rf5)

selSpecies5 <- rownames(importance(rf5))[importance(rf5) != 0]
rfFormula6 <- "state ~" %>%
  paste(paste(selSpecies5, collapse = "+"))

rf6 <- randomForest(as.formula(rfFormula6), data = dat, subset = smp_idx)
plot(rf6)

selSpecies6 <- rownames(importance(rf6))[importance(rf6) != 0]
rfFormula7 <- "state ~" %>%
  paste(paste(selSpecies6, collapse = "+"))

rf7 <- randomForest(as.formula(rfFormula7), data = dat, subset = smp_idx)
plot(rf7)

selSpecies7 <- rownames(importance(rf7))[importance(rf7) != 0]
rfFormula8 <- "state ~" %>%
  paste(paste(selSpecies7, collapse = "+"))

set.seed(123)
rf8 <- randomForest(as.formula(rfFormula8), data = dat, subset = smp_idx)
plot(rf8)

selSpecies8 <- rownames(importance(rf8))[importance(rf8) != 0]
rfFormula9 <- "state ~" %>%
  paste(paste(selSpecies8, collapse = "+"))

set.seed(123)
rf9 <- randomForest(as.formula(rfFormula9), data = dat, subset = smp_idx)
plot(rf9)

selSpecies9 <- rownames(importance(rf9))[importance(rf9) != 0]
rfFormula10 <- "state ~" %>%
  paste(paste(selSpecies9, collapse = "+"))

set.seed(123)
rf10 <- randomForest(as.formula(rfFormula10), data = dat, subset = smp_idx)
plot(rf10)

selSpecies10 <- rownames(importance(rf10))[importance(rf10) != 0]
rfFormula11 <- "state ~" %>%
  paste(paste(selSpecies10, collapse = "+"))

set.seed(123)
rf11 <- randomForest(as.formula(rfFormula11), data = dat, subset = smp_idx)
plot(rf11)


rfTest <- function(dat, subset, ntree = 500, seed = 123){
  rfFormula <- "state ~ ."
  
  rfs <- list()
  
  set.seed(seed)
  rf <- randomForest(as.formula(rfFormula), data = dat, subset = subset, ntree = ntree)
  rfs[[1]] <- rf
  
  err <- round(sum(rf$confusion[5:6]),4)
  lastErr <- err
  
  i <- 2
  while(TRUE){
    
    selSpecies <- rownames(importance(rf))[importance(rf) != 0]
    rfFormula <- "state ~" %>%
      paste(paste(selSpecies, collapse = "+"))
    
    set.seed(seed)
    rf <- randomForest(as.formula(rfFormula), data = dat, subset = subset, ntree = ntree)
    err <- round(sum(rf$confusion[5:6]),4)
    cat(err)
    
    if(round(err, 4) <= lastErr && sum(importance(rf) == 0) != 0){
      lastErr <- err
      rfs[[i]] <- rf
      i <- i + 1
    }else{
      break()
    }
  }
  rfs
}

rft <- rfTest(dat, smp_idx)

impTabs <- lapply(rft, function(x){
  importance(x) %>%
    as.data.frame() %>%
    rownames_to_column("Species")
})

for(i in 1 : length(impTabs)){
  impTabs[[i]] <- impTabs[[i]] %>%
    mutate(n = i)
}

impSum <- do.call(rbind, impTabs)

impSum %>%
  group_by(Species) %>%
  summarise(sumMDG = sum(MeanDecreaseGini)) %>%
  top_n(20, sumMDG) %>%
  arrange(desc(sumMDG)) %>%
  .$Species -> topSpecies

impVarPlot1 <- impSum %>%
  filter(Species %in% topSpecies) %>%
  as.data.frame() %>%
  left_join(dat %>%
              summarise_if(is.numeric, funs(mean)) %>%
              t() %>%
              as.data.frame() %>%
              rownames_to_column("Species") %>%
              rename(mean_abund = V1)) %>%
  ggplot() +
  geom_point(aes(x = MeanDecreaseGini, y = factor(Species, levels = rev(topSpecies)), color = factor(n), size = mean_abund)) +
  scale_color_manual(values = rev(brewer.pal(4, "Spectral")), name = "Rounds of test") +
  ylab("Top 20 Species") +
  scale_size_continuous(name = "Mean of Abundance") +
  ggtitle("Importance of Variables (ntree = 500)")


ggsave("impVarPlot_500.pdf", device = "pdf")

##### ntree = 1000

rft <- rfTest(dat, smp_idx, 1000)

impTabs <- lapply(rft, function(x){
  importance(x) %>%
    as.data.frame() %>%
    rownames_to_column("Species")
})

for(i in 1 : length(impTabs)){
  impTabs[[i]] <- impTabs[[i]] %>%
    mutate(n = i)
}

impSum <- do.call(rbind, impTabs)

impSum %>%
  group_by(Species) %>%
  summarise(sumMDG = sum(MeanDecreaseGini)) %>%
  top_n(20, sumMDG) %>%
  arrange(desc(sumMDG)) %>%
  .$Species -> topSpecies

impVarPlot1 <- impSum %>%
  filter(Species %in% topSpecies) %>%
  as.data.frame() %>%
  left_join(dat %>%
              summarise_if(is.numeric, funs(mean)) %>%
              t() %>%
              as.data.frame() %>%
              rownames_to_column("Species") %>%
              rename(mean_abund = V1)) %>%
  ggplot() +
  geom_point(aes(x = MeanDecreaseGini, y = factor(Species, levels = rev(topSpecies)), color = n, size = mean_abund)) +
  scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")), name = "Rounds of Test") +
  ylab("Top 20 Species") +
  scale_size_continuous(name = "Mean of Abundance") +
  ggtitle("Importance of Variables (ntree = 1000)")


ggsave("impVarPlot_1000.pdf", device = "pdf")


##### seed = 1434

rft <- rfTest(dat, smp_idx, 1000, 1434)

impTabs <- lapply(rft, function(x){
  importance(x) %>%
    as.data.frame() %>%
    rownames_to_column("Species")
})

for(i in 1 : length(impTabs)){
  impTabs[[i]] <- impTabs[[i]] %>%
    mutate(n = i)
}

impSum <- do.call(rbind, impTabs)

impSum %>%
  group_by(Species) %>%
  summarise(sumMDG = sum(MeanDecreaseGini)) %>%
  top_n(20, sumMDG) %>%
  arrange(desc(sumMDG)) %>%
  .$Species -> topSpecies

impVarPlot1 <- impSum %>%
  filter(Species %in% topSpecies) %>%
  as.data.frame() %>%
  left_join(dat %>%
              summarise_if(is.numeric, funs(mean)) %>%
              t() %>%
              as.data.frame() %>%
              rownames_to_column("Species") %>%
              rename(mean_abund = V1)) %>%
  ggplot() +
  geom_point(aes(x = MeanDecreaseGini, y = factor(Species, levels = rev(topSpecies)), color = n, size = mean_abund)) +
  scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")), name = "Rounds of Test") +
  ylab("Top 20 Species") +
  scale_size_continuous(name = "Mean of Abundance") +
  ggtitle("Importance of Variables (ntree = 1000, seed = 1434)")


ggsave(filename = "impVarPlot_1000_1434.pdf", impVarPlot1, device = "pdf")


##### ntree = 2000

rft <- rfTest(dat, smp_idx, 2000)

impTabs <- lapply(rft, function(x){
  importance(x) %>%
    as.data.frame() %>%
    rownames_to_column("Species")
})

for(i in 1 : length(impTabs)){
  impTabs[[i]] <- impTabs[[i]] %>%
    mutate(n = i)
}

impSum <- do.call(rbind, impTabs)

impSum %>%
  group_by(Species) %>%
  summarise(sumMDG = sum(MeanDecreaseGini)) %>%
  top_n(20, sumMDG) %>%
  arrange(desc(sumMDG)) %>%
  .$Species -> topSpecies

impVarPlot1 <- impSum %>%
  filter(Species %in% topSpecies) %>%
  as.data.frame() %>%
  left_join(dat %>%
              summarise_if(is.numeric, funs(mean)) %>%
              t() %>%
              as.data.frame() %>%
              rownames_to_column("Species") %>%
              rename(mean_abund = V1)) %>%
  ggplot() +
  geom_point(aes(x = MeanDecreaseGini, y = factor(Species, levels = rev(topSpecies)), color = n, size = mean_abund)) +
  scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")), name = "Rounds of Test") +
  ylab("Top 20 Species") +
  scale_size_continuous(name = "Mean of Abundance") +
  ggtitle("Importance of Variables (ntree = 2000)")


ggsave(filename = "impVarPlot_2000.pdf", impVarPlot1, device = "pdf")


##### 4 = disease , 14 = health

selSam <- metaTab %>%
  filter(time %in% c(4, 14)) %>%
  .$V1

otuTab[rownames(otuTab) %in% selSam,] %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  t() %>%
  data.frame(stringsAsFactors = FALSE) %>%
  rownames_to_column("ASV") %>%
  left_join(taxTab4 %>%
              select(X1, Species), by = c("ASV" = "X1")) %>%
  group_by(Species) %>%
  select(-ASV) %>%
  summarise_all(funs(sum(as.numeric(.)))) %>%
  ungroup() %>%
  mutate_if(is.numeric, funs(. / sum(.))) %>%
  as.data.frame() %>%
  `rownames<-`(make.names(.$Species)) %>%
  select(-Species) %>%
  t()