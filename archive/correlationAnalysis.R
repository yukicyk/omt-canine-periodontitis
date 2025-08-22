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

### Pearson
allDf %>%
  mutate_all(funs(as.numeric)) %>%
  mutate_all(funs(ifelse(. == 0, NA, .))) %>%
  as.matrix() %>%
  Hmisc::rcorr(type = "pearson") -> corr

combn(colnames(corr$r), 2) %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y")) -> combTab

corr$r %>%
  as.data.frame() %>%
  rownames_to_column("x") %>%
  gather(key = "y", value = "r", -x) -> corrR

corr$n %>%
  as.data.frame() %>%
  rownames_to_column("x") %>%
  gather(key = "y", value = "n", -x) -> corrN

corr$P %>%
  as.data.frame() %>%
  rownames_to_column("x") %>%
  gather(key = "y", value = "p", -x) -> corrP

combTab %>%
  left_join(corrR, by = c("x" = "x", "y" = "y")) %>%
  left_join(corrP, by = c("x" = "x", "y" = "y")) %>%
  left_join(corrN, by = c("x" = "x", "y" = "y")) -> corrTab

corrTab %>%
  filter(n > 10) %>%
  filter(p <= 0.05) %>%
  filter(abs(r) >= 0.5) %>%
  filter(x %in% c("BOP", "PPD") | y %in% c("BOP", "PPD")) -> selAsvTab

labY <- unique(selAsvTab[, "y"])
labX <- unique(selAsvTab[, "x"])[-4]

corrTab %>%
  filter(x %in% labX) %>%
  filter(y %in% labY) %>%
  ggplot() +
  geom_tile(aes(x = x, y = y), fill = "white", color = "black") +
  coord_equal() +
  geom_point(aes(x = x, y = y, color = r, size = r))


### agglomerate species

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
