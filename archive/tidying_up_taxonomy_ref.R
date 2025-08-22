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









