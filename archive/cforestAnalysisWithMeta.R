library(tidyverse)
library(ROCR)
library(randomForest)

bopDat <- read_tsv("Ergebnisse_BOP_PPD_summarized_BOP_PPD_3.txt") %>%
  filter(Tooth == "BOP") %>%
  dplyr::rename(sampleID = Dog) %>%
  select(sampleID, BOP)

ppdDat <- read_tsv("Ergebnisse_BOP_PPD_summarized_BOP_PPD_3.txt") %>%
  filter(Tooth == "PPD") %>%
  dplyr::rename(sampleID = Dog) %>%
  select(sampleID, PPD)

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

### BOP
bopDf <- otuTab %>%
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
  select(-sampleID)

smpSize <- floor(0.7 * nrow(bopDf))
smpIdx <- sample(1 : nrow(bopDf), smpSize)

trainSet <- bopDf[smpIdx,]
testSet <- bopDf[-smpIdx,]

rfTune <- tuneRF(trainSet[, 1:1759], trainSet[, 1759])

set.seed(123)

forestBop1 <- randomForest(BOP ~ ., data = trainSet, ntree = 500, mtry = 1759)
predBop1 <- predict(forestBop1, testSet)

testSet$BOP %>%
  cbind(predBop1) %>%
  `colnames<-`(c("True_value", "Predict_value")) %>%
  as.tibble() %>%
  ggplot() +
  geom_point(aes(x = Predict_value, y = True_value)) +
  geom_abline(intercept = 0, slope = 1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

varImpPlot(forestBop1)

RMSEBop1 <- sqrt(mean((testSet$BOP - predBop1)^2))
MAEBop1 <- mean(abs(testSet$BOP - predBop1))

# Round 2
trainSet2 <- trainSet[, c(importance(forestBop1) > 0, TRUE)]
testSet2 <- testSet[, c(importance(forestBop1) > 0, TRUE)]

set.seed(123)

rfTune2 <- tuneRF(trainSet2[, 1:838], trainSet2[, 839], ntreeTry = 500)

set.seed(123)

forestBop2 <- randomForest(BOP ~ ., data = trainSet2, ntree = 500, mtry = 558)
predBop2 <- predict(forestBop2, testSet2)

testSet2$BOP %>%
  cbind(predBop2) %>%
  `colnames<-`(c("True_value", "Predict_value")) %>%
  as.tibble() %>%
  ggplot() +
  geom_point(aes(x = Predict_value, y = True_value)) +
  geom_abline(intercept = 0, slope = 1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

varImpPlot(forestBop2)

RMSEBop2 <- sqrt(mean((testSet2$BOP - predBop2)^2))
MAEBop2 <- mean(abs(testSet2$BOP - predBop2))

importance(forestBop2) %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  top_n(20, IncNodePurity) %>%
  arrange(desc(IncNodePurity)) %>%
  left_join(trainSet2 %>%
              summarise_all(funs(mean)) %>%
              as.data.frame() %>%
              t() %>%
              as.data.frame() %>%
              rownames_to_column("species")) %>%
  mutate(species = factor(species, levels = rev(species))) %>%
  ggplot() +
  geom_point(aes(x = IncNodePurity, y = species, size = V1)) +
  labs(y = "Species", x = "Importance") +
  guides(size = FALSE) +
  ggtitle("Importance of random forest on BOP (tree = 500, round = 2)")

ggsave("raymond_result/impBopAsvR2.pdf", device = "pdf")

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

bopSpeciesDf <- aggTab2 %>%
  as.data.frame() %>%
  rownames_to_column("sampleID") %>%
  .[1:77,] %>%
  left_join(bopDat %>%
              mutate(sampleID = make.names(sampleID))) %>%
  select(-sampleID)

trainSet3 <- bopSpeciesDf[smpIdx,]
testSet3 <- bopSpeciesDf[-smpIdx,]

rfTune <- tuneRF(trainSet3[, 1:419], trainSet3[, 420], ntreeTry = 500)

set.seed(123)

forestBop3 <- randomForest(BOP ~ ., data = trainSet3, ntree = 500, mtry = 70)
predBop3 <- predict(forestBop3, testSet3)

testSet3$BOP %>%
  cbind(predBop3) %>%
  `colnames<-`(c("True_value", "Predict_value")) %>%
  as.tibble() %>%
  ggplot() +
  geom_point(aes(x = Predict_value, y = True_value)) +
  geom_abline(intercept = 0, slope = 1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

varImpPlot(forestBop3)

RMSEBop3 <- sqrt(mean((testSet3$BOP - predBop3)^2))
MAEBop3 <- mean(abs(testSet3$BOP - predBop3))

# Round 2
trainSet4 <- trainSet3[, c(importance(forestBop3) > 0, TRUE)]
testSet4 <- testSet3[, c(importance(forestBop3) > 0, TRUE)]

set.seed(123)

rfTune2 <- tuneRF(trainSet4[, 1:349], trainSet2[, 350], ntreeTry = 500)

set.seed(123)

forestBop4 <- randomForest(BOP ~ ., data = trainSet4, ntree = 500, mtry = 349)
predBop4 <- predict(forestBop4, testSet4)

testSet4$BOP %>%
  cbind(predBop4) %>%
  `colnames<-`(c("True_value", "Predict_value")) %>%
  as.tibble() %>%
  ggplot() +
  geom_point(aes(x = Predict_value, y = True_value)) +
  geom_abline(intercept = 0, slope = 1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

varImpPlot(forestBop4)

RMSEBop4 <- sqrt(mean((testSet4$BOP - predBop4)^2))
MAEBop4 <- mean(abs(testSet4$BOP - predBop4))

# Round 3
trainSet5 <- trainSet4[, c(importance(forestBop4) > 0, TRUE)]
testSet5 <- testSet4[, c(importance(forestBop4) > 0, TRUE)]

set.seed(123)

rfTune2 <- tuneRF(trainSet5[, 1:316], trainSet5[, 317], ntreeTry = 500)

set.seed(123)

forestBop5 <- randomForest(BOP ~ ., data = trainSet5, ntree = 500, mtry = 105)
predBop5 <- predict(forestBop5, testSet5)

testSet5$BOP %>%
  cbind(predBop5) %>%
  `colnames<-`(c("True_value", "Predict_value")) %>%
  as.tibble() %>%
  ggplot() +
  geom_point(aes(x = Predict_value, y = True_value)) +
  geom_abline(intercept = 0, slope = 1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

varImpPlot(forestBop5)

RMSEBop5 <- sqrt(mean((testSet5$BOP - predBop5)^2))
MAEBop5 <- mean(abs(testSet5$BOP - predBop5))

importance(forestBop5) %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  top_n(20, IncNodePurity) %>%
  arrange(desc(IncNodePurity)) %>%
  left_join(trainSet5 %>%
              summarise_all(funs(mean)) %>%
              as.data.frame() %>%
              t() %>%
              as.data.frame() %>%
              rownames_to_column("species")) %>%
  mutate(species = factor(species, levels = rev(species))) %>%
  ggplot() +
  geom_point(aes(x = IncNodePurity, y = species, size = V1)) +
  labs(y = "Species", x = "Importance") +
  guides(size = FALSE) +
  ggtitle("Importance of random forest on BOP (tree = 500, round = 3)")

ggsave("raymond_result/impBopR3.pdf", device = "pdf")

### PPD
ppdDf <- otuTab %>%
  apply(., c(1, 2), as.numeric) %>%
  data.frame(stringsAsFactors = FALSE) %>%
  t() %>%
  data.frame(stringsAsFactors = FALSE) %>%
  rownames_to_column("AVS") %>%
  mutate_if(is.numeric, funs(. / sum(.))) %>%
  gather(key = "sampleID", value = "value", -AVS) %>%
  spread(key = "AVS", value = "value") %>%
  .[1:77,] %>%
  left_join(ppdDat %>%
              mutate(sampleID = make.names(sampleID))) %>%
  select(-sampleID)

pTrainSet <- ppdDf[smpIdx,]
pTestSet <- ppdDf[-smpIdx,]

tuneRF(pTrainSet[, 1:1759], pTrainSet[, 1759])

set.seed(123)

forestPpd1 <- randomForest(PPD ~ ., data = pTrainSet, ntree = 500, mtry = 1759)
predPpd1 <- predict(forestPpd1, pTestSet)

pTestSet$PPD %>%
  cbind(predPpd1) %>%
  `colnames<-`(c("True_value", "Predict_value")) %>%
  as.tibble() %>%
  ggplot() +
  geom_point(aes(x = Predict_value, y = True_value)) +
  geom_abline(intercept = 0, slope = 1) +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 5))

varImpPlot(forestPpd1)

RMSEPpd1 <- sqrt(mean((pTestSet$PPD - predPpd1)^2))
MAEPpd1 <- mean(abs(pTestSet$PPD - predPpd1))

# Round 2
pTrainSet2 <- pTrainSet[, c(importance(forestPpd1) > 0, TRUE)]
pTestSet2 <- pTestSet[, c(importance(forestPpd1) > 0, TRUE)]

set.seed(123)

tuneRF(pTrainSet2[, 1:840], pTrainSet2[, 841], ntreeTry = 500)

set.seed(123)

forestPpd2 <- randomForest(PPD ~ ., data = pTrainSet2, ntree = 500, mtry = 280)
predPpd2 <- predict(forestPpd2, pTestSet2)

pTestSet2$PPD %>%
  cbind(predPpd2) %>%
  `colnames<-`(c("True_value", "Predict_value")) %>%
  as.tibble() %>%
  ggplot() +
  geom_point(aes(x = Predict_value, y = True_value)) +
  geom_abline(intercept = 0, slope = 1) +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 5))

varImpPlot(forestPpd2)

RMSEPpd2 <- sqrt(mean((pTestSet2$PPD - predPpd2)^2))
MAEPpd2 <- mean(abs(pTestSet2$PPD - predPpd2))

importance(forestPpd2) %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  top_n(20, IncNodePurity) %>%
  arrange(desc(IncNodePurity)) %>%
  left_join(pTrainSet2 %>%
              summarise_all(funs(mean)) %>%
              as.data.frame() %>%
              t() %>%
              as.data.frame() %>%
              rownames_to_column("species")) %>%
  mutate(species = factor(species, levels = rev(species))) %>%
  ggplot() +
  geom_point(aes(x = IncNodePurity, y = species, size = V1)) +
  labs(y = "Species", x = "Importance") +
  guides(size = FALSE) +
  ggtitle("Importance of random forest on PPD (tree = 500, round = 2)")

ggsave("raymond_result/impPpdAsvR2.pdf", device = "pdf")

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

ppdSpeciesDf <- aggTab2 %>%
  as.data.frame() %>%
  rownames_to_column("sampleID") %>%
  .[1:77,] %>%
  left_join(ppdDat %>%
              mutate(sampleID = make.names(sampleID))) %>%
  select(-sampleID)

pTrainSet3 <- ppdSpeciesDf[smpIdx,]
pTestSet3 <- ppdSpeciesDf[-smpIdx,]

set.seed(123)

tuneRF(trainSet3[, 1:419], trainSet3[, 420], ntreeTry = 500)

set.seed(123)

forestPpd3 <- randomForest(PPD ~ ., data = pTrainSet3, ntree = 500, mtry = 70)
predPpd3 <- predict(forestPpd3, pTestSet3)

pTestSet3$PPD %>%
  cbind(predPpd3) %>%
  `colnames<-`(c("True_value", "Predict_value")) %>%
  as.tibble() %>%
  ggplot() +
  geom_point(aes(x = Predict_value, y = True_value)) +
  geom_abline(intercept = 0, slope = 1) +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 5))

varImpPlot(forestPpd3)

RMSEPpd3 <- sqrt(mean((pTestSet3$PPD - predPpd3)^2))
MAEPpd3 <- mean(abs(pTestSet3$PPD - predPpd3))

# Round 2
pTrainSet4 <- pTrainSet3[, c(importance(forestPpd3) > 0, TRUE)]
pTestSet4 <- pTestSet3[, c(importance(forestPpd3) > 0, TRUE)]

set.seed(123)

tuneRF(trainSet4[, 1:348], trainSet2[, 349], ntreeTry = 500)

set.seed(123)

forestPpd4 <- randomForest(PPD ~ ., data = pTrainSet4, ntree = 500, mtry = 232)
predPpd4 <- predict(forestPpd4, pTestSet4)

pTestSet4$PPD %>%
  cbind(predPpd4) %>%
  `colnames<-`(c("True_value", "Predict_value")) %>%
  as.tibble() %>%
  ggplot() +
  geom_point(aes(x = Predict_value, y = True_value)) +
  geom_abline(intercept = 0, slope = 1) +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 5))

varImpPlot(forestPpd4)

RMSEPpd4 <- sqrt(mean((pTestSet4$PPD - predPpd4)^2))
MAEPpd4 <- mean(abs(pTestSet4$PPD - predPpd4))

# Round 3
pTrainSet5 <- pTrainSet4[, c(importance(forestPpd4) > 0, TRUE)]
pTestSet5 <- pTestSet4[, c(importance(forestPpd4) > 0, TRUE)]

set.seed(123)

tuneRF(pTrainSet5[, 1:325], pTrainSet5[, 326], ntreeTry = 500)

set.seed(123)

forestPpd5 <- randomForest(PPD ~ ., data = pTrainSet5, ntree = 500, mtry = 54)
predPpd5 <- predict(forestPpd5, pTestSet5)

pTestSet5$PPD %>%
  cbind(predPpd5) %>%
  `colnames<-`(c("True_value", "Predict_value")) %>%
  as.tibble() %>%
  ggplot() +
  geom_point(aes(x = Predict_value, y = True_value)) +
  geom_abline(intercept = 0, slope = 1) +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 5))

varImpPlot(forestPpd5)

RMSEPpd5 <- sqrt(mean((pTestSet5$PPD - predPpd5)^2))
MAEPpd5 <- mean(abs(pTestSet5$PPD - predPpd5))

importance(forestPpd5) %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  top_n(20, IncNodePurity) %>%
  arrange(desc(IncNodePurity)) %>%
  left_join(pTrainSet5 %>%
              summarise_all(funs(mean)) %>%
              as.data.frame() %>%
              t() %>%
              as.data.frame() %>%
              rownames_to_column("species")) %>%
  mutate(species = factor(species, levels = rev(species))) %>%
  ggplot() +
  geom_point(aes(x = IncNodePurity, y = species, size = V1)) +
  labs(y = "Species", x = "Importance") +
  guides(size = FALSE) +
  ggtitle("Importance of random forest on PPD (tree = 500, round = 3)")

ggsave("raymond_result/impPpdSpeciesR3.pdf", device = "pdf")











