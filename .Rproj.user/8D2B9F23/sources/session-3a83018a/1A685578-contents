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







