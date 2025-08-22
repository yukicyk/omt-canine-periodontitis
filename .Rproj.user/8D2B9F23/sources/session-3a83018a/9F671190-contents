library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(party)
library(partykit)
library(caret)
library(randomForest)
library(RColorBrewer)

# OTU table
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

# Meta table
metaTab <- rownames(otuTab) %>%
  strsplit("-") %>%
  lapply(`[[`, 2) %>%
  unlist() %>%
  cbind(rownames(otuTab), .) %>%
  as.data.frame() %>%
  mutate(time = c(rep(c(0, 2, 4, 14), 4), rep(0, 5), rep(c(0, 2, 4, 14), 14), rep(2, 9)))

write.table(metaTab, "meta.tab", sep = ",", col.names = FALSE, row.names = FALSE)

# Sample table
subTab <- levels(metaTab[,2]) %>%
  cbind(c(3,1,3,1,2,3,1,3,1,3,1,3,1,3,1,3,3,1,1))

subTab <- as.data.frame(subTab)

subTab[, 2] <- factor(subTab[, 2])

levels(subTab[, 2]) <- c("Control", "Donor", "Recepient")

colnames(subTab) <- c("subjectID", "role")

selSam <- metaTab[metaTab$time %in% c(0, 2), "V1"]

# prepare train and test data
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

### try 3 with eliminate ASV importance lower than 1
selAsv3 <- names(impTab2)[abs(impTab2) >= 1] # select ASV with importance over 1

formula3 <- "state ~" %>%
  paste(paste(selAsv3, collapse = "+"))

set.seed(123)
forest3 <- cforest(as.formula(formula3), data = trainSet, ntree = 500, mtry = 15)
pred3 <- Predict(forest3, newdata = testSet)

confusMat3 <- confusionMatrix(pred3, testSet$state) # accurancy 100%

impTab3 <- varimp(forest3)

### try 4 with further eliminate ASV importance lower than 1
selAsv4 <- names(impTab3)[abs(impTab3) >= 1]

formula4 <- "state ~" %>%
  paste(paste(selAsv4, collapse = "+"))

set.seed(123)
forest4 <- cforest(as.formula(formula4), data = trainSet, ntree = 500, mtry = 15)
pred4 <- Predict(forest4, newdata = testSet)

confusMat4 <- confusionMatrix(pred4, testSet$state) # accurancy 93%
varimp(forest4)

### according to try 3
taxTab <- read_tsv("ps1_taxtable.txt")

impTaxTab <- names(impTab3) %>%
  cbind(impTab3) %>%
  as.tibble() %>%
  `colnames<-`(c("ID", "Importance")) %>%
  left_join(taxTab, by = c("ID" = "Kingdom"))

# predict 14 week data
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

# predict 4 week data
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

# predict 2 week data
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

### plot individual ASV
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

















