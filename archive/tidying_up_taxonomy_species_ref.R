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
