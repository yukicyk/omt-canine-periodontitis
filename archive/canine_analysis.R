#Canine microbiome transplant phyloseq, microbiome and other downstream analysis
#20181220

setwd(dir = "~/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/")
library("phyloseq","ggplot2","dplyr")
#creating the finalize dataset
# load("~/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/history/Canine_ps.RData")
#import the updated metadata mapping file (with dates)
samdat <- read.delim("~/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/mappingV4.txt", header=TRUE)
sampledata <-samdat[,2:ncol(samdat)]
rownames(sampledata)<-samdat[,1]
sample_data(ps1)<-sampledata
# ps1 = prune_taxa(taxa_sums(ps1) > 0, ps1)
# sample_data(ps1)$Treatment = gsub("Recipeint", "Recipient", sample_data(ps1)$Treatment)
saveRDS(ps1,"~/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/ps1.rds")

#load ps1 data
ps1<- readRDS("~/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/ps1.rds")

#prepare dataframe table of data
otutab<-as.data.frame(otu_table(ps1)) %>% tibble::rownames_to_column(.)
colnames(otutab)[colnames(otutab)=="rowname"] <- "Sample"
metadata <-read.delim("~/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/mappingV4.txt", header=TRUE)

#From Raymond's script
tax_dat <- read_tsv("//Users/OB/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/model_picrust/ps1_taxtable.txt") %>%
  rowwise() %>%
  mutate(
    Phylum = ifelse(is.na(Phylum), Kingdom, Phylum),
    Class = ifelse(is.na(Class), Phylum, Class),
    Order = ifelse(is.na(Order), Class, Order),
    Family = ifelse(is.na(Family), Order, Family),
    Genus = ifelse(is.na(Genus), Family, Genus),
    Species = ifelse(is.na(Species), Genus, Species)
  ) %>%
  mutate(Species = ifelse(Species == Genus, paste(Species, "sp.", sep = "_"), paste(Genus, Species, sep = "_")))

#data stats counting number of sequences
summary(sample_sums(ps1))
subset_samples(ps1, sample_data(ps1)$Treatment == "Control") %>% 
  sample_sums(.) %>% summary(.)

#plot number of sequences and number of ASVs
#From Raymond's script
otu_tab <- read_tsv("//Users/OB/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/model_picrust/otu.tab")
getColor <- function(n){
  
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1 : n]
  
}

col <- getColor(3) %>%
{c(.[2], .[1], .[3])}

num_ASV_sample <- otu_tab %>%
  gather(Sample, value, -`#OTU ID`) %>%
  mutate(value = value != 0) %>%
  group_by(Sample) %>%
  summarise(`Number of ASV` = sum(value)) %>%
  left_join(metadata) %>%
  mutate(Treatment = factor(Treatment, levels = c("Donor", "Control", "Recipient")))
#plot number of ASVs for each sample
num_ASV_sample %>%
  rename(Num = `Number of ASV`) %>%
  ggplot(aes(Dog, Num, fill = Treatment, group = Sample)) +
  geom_col(position = "dodge") +
  labs(y = "Number of ASV") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = col)
#counting number of ASVs
subset_samples(ps1, sample_data(ps1)$Treatment == "Recipient") %>% 
  prune_taxa(taxa_sums(.) > 0, .) %>% ntaxa(.)
num_ASV_sample %>%
  filter(., Treatment =="Control") %>%
  summary(.)
#plot number of sequences for each sample
num_seq_sample <- otu_tab %>%
  gather(Sample, value, -`#OTU ID`) %>%
  group_by(Sample) %>%
  summarise(`Number of seq` = sum(value)) %>%
  left_join(map_dat) %>%
  mutate(Treatment = factor(Treatment, levels = c("Donor", "Control", "Recipient")))
num_seq_sample %>%
  rename(Num = `Number of seq`) %>%
  ggplot(aes(Dog, Num, fill = Treatment, group = Sample)) +
  geom_col(position = "dodge") +
  labs(y = "Number of reads") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = col)

#plot number of ASV for each treatment group
num_ASV_tr <- otu_tab %>%
  gather(Sample, value, -`#OTU ID`) %>%
  left_join(map_dat) %>%
  group_by(Treatment, `#OTU ID`) %>%
  summarise(value = sum(value)) %>%
  mutate(value = value != 0) %>%
  group_by(Treatment) %>%
  summarise(`Number of ASV` = sum(value))

num_ASV_tr %>%
  as.data.frame()

num_ASV_tr %>%
  mutate(Treatment = factor(Treatment, levels = c("Donor", "Control", "Recipient"))) %>%
  ggplot(aes(Treatment, `Number of ASV`, fill = Treatment)) +
  geom_col() +
  scale_fill_manual(values = col)

#Plot rarefraction curves for all samples, and for the 6 clinical categories
#load the phyloseq extension richness.R functions
library("vegan", "ggplot2")
p.rare <- ggrare(ps1, step = 1000, color = "Treatment",label="Sample", se = FALSE)
p.rare + xlim(0, 25000)
p.all <- p.rare + facet_wrap(~Treatment)+ xlim(0, 25000)
p.all


#venn diagram showing number of ASVs shared among 3 treatment group
otu_tab %>%
  gather(Sample, value, -`#OTU ID`) %>%
  left_join(map_dat) %>%
  group_by(Treatment, `#OTU ID`) %>%
  summarise(value = sum(value)) %>%
  filter(value != 0) %>%
  group_by(Treatment) %>%
  summarise(ASV = list(`#OTU ID`)) %>%
  {`names<-`(.$ASV, .$Treatment)} %>%
  venn.diagram(filename = NULL,main.fontfamily = "serif", sub.fontfamily = "serif",fill = c(col[2], col[1], col[3])) %>%
  grid.draw()

# rarefy_even_depth: Resample an OTU table such that all samples have the same sampling depth
ps.rf = rarefy_even_depth(ps1, sample.size = 18614,
                             rngseed = 20181221, replace = TRUE, trimOTUs = FALSE, verbose = TRUE)



# Plot Shannon index and ASVs across dogs
plot_richness(ps.rf, x="Dog", measures=c("Shannon"), color="Treatment", shape="Timepoint") + 
  geom_point(size=7, alpha=0.75)+ ggtitle("Alpha Diversity: Shannon Entropy across samples") + ylim(0,6)
  
plot_richness(ps.rf, x="Dog", measures=c("Observed"), color="Treatment", shape="Timepoint") + 
  geom_point(size=7, alpha=0.75)+ ggtitle("Alpha Diversity: Observed ASVs across samples") 
theme_set(theme_bw())

# Plot Shannon index across time
donor.rf<-subset_samples(ps.rf, Treatment=="Donor") %>% prune_taxa(taxa_sums(.) > 0, .)
recipient.rf<-subset_samples(ps.rf, Treatment=="Recipient")%>% prune_taxa(taxa_sums(.) > 0, .)
control.rf<-subset_samples(ps.rf, Treatment=="Control")%>% prune_taxa(taxa_sums(.) > 0, .)
 
p<-plot_richness(recipient.rf, x="Week", measures=c("Shannon"), color="Dog") + 
  geom_point(size=7, alpha=0.75) +ggtitle("Shannon Entropies across timepoints in the Recipient dogs")+ylim(0,6)
p+geom_line(data=p$data, aes(x=Week, y=value, group=Dog))+ theme(axis.text.x= element_text(size=12,angle=0))+scale_x_continuous(limits=c(0,14), breaks=seq(0,14,2))

p<-plot_richness(control.rf, x="Week", measures=c("Shannon"), color="Dog") + 
  geom_point(size=7, alpha=0.75) +ggtitle("Shannon Entropies across timepoints in the Control dogs")+ylim(0,6)
p+geom_line(data=p$data, aes(x=Week, y=value, group=Dog))+ theme(axis.text.x= element_text(size=12,angle=0))+scale_x_continuous(limits=c(0,14), breaks=seq(0,14,2))

p<-plot_richness(ps.rf%>%subset_samples(Treatment=="Control"|Treatment=="Recipient"), x="Week", measures=c("Shannon"), color="Dog") + 
  geom_point(size=7, alpha=0.75) +ggtitle("Diversity changes across timepoints in the Control and Recipient dogs")
p+geom_line(data=p$data, aes(x=Week, y=value, group=Dog))+ theme(axis.text.x= element_text(size=12,angle=0))+scale_x_continuous(limits=c(0,14), breaks=seq(0,14,2))+
facet_wrap(~Treatment, ncol = 2)+theme(legend.position="bottom")+ylab("Shannon Entropies")+expand_limits(y=0)


library("ggplot2")
theme_set(theme_bw(base_size=14))
p<- plot_richness(donor.rf, x="Date", measures=c("Shannon"), color="Timepoint") + 
  geom_point(size=7, alpha=0.75) +ggtitle("Shannon Entropies changes along time in Donor dog E")+ scale_color_manual(labels = c("T1: Baseline", "T2: Transplantation"), values = c("blue", "red")) + theme(axis.text.x= element_text(size=10),legend.text= element_text(size=10))


p$data %<>%
  mutate(Date=as.Date(Date, format="%d/%m/%Y")) %>%
  mutate(Date = gsub("0016", "2016", Date)) %>%
  mutate(Date = gsub("0017", "2017", Date))
p$data %>%
  mutate(samples = gsub("Transplant","MT", samples)) ->p$data
p +geom_line(data=p$data, aes(x=Date, y=value, group=Dog), colour="gray") + ylim(0,5) +
  geom_text(data=p$data, aes(label=samples), nudge_y = -0.4)
#adding segment of lines instead
p$data$value
#[1] 4.836509 4.083073 4.369847 4.650906 4.098328 4.342519 4.202698 4.368164 3.891975 4.156283 4.211308 3.831847 3.974388
#[14] 4.142916
p$data$Date
#[1] "2016-10-07" "2017-02-01" "2017-03-01" "2017-06-06" "2017-09-06" "2016-10-20" "2016-10-21" "2017-02-15" "2017-02-16"
#[10] "2017-03-21" "2017-03-24" "2017-06-20" "2017-09-20" "2017-09-21"
p+geom_segment(aes(x = "2016-10-07", y = 4.836509, xend = "2016-10-20", yend = 4.342519), colour ="gray")+
  geom_segment(aes(x = "2016-10-20", y = 4.342519, xend = "2016-10-21", yend = 4.202698), colour ="gray")+
  geom_segment(aes(x = "2017-02-01", y = 4.083073, xend = "2017-02-15", yend = 4.368164), colour ="gray")+
  geom_segment(aes(x = "2017-02-15", y = 4.368164, xend = "2017-02-16", yend = 3.891975), colour ="gray")+
  geom_segment(aes(x = "2017-03-01", y = 4.369847, xend = "2017-03-21", yend = 4.156283), colour ="gray")+
  geom_segment(aes(x = "2017-03-21", y = 4.156283, xend = "2017-03-24", yend = 4.211308), colour ="gray")+
  geom_segment(aes(x = "2017-06-06", y = 4.650906, xend = "2017-06-20", yend = 3.831847), colour ="gray")+
  geom_segment(aes(x = "2017-09-06", y = 4.098328, xend = "2017-09-20", yend = 3.974388), colour ="gray")+
  geom_segment(aes(x = "2017-09-20", y = 3.974388, xend = "2017-09-21", yend = 4.142916), colour ="gray")

#compare alpha-diversity at baseline (T1), and across different timepoint
library("ggsignif")
alpha_meas = c("Observed", "Shannon")
theme_set(theme_bw())
plot_richness(ps.rf%>% subset_samples(Timepoint=="T1"), x="Treatment", measures=alpha_meas, color="Treatment")+
  geom_boxplot(aes(x=Treatment, y=value, color=Treatment), alpha=0.1) +
    theme(strip.background = element_rect(colour = "black", fill = "white"), text=element_text(size=12))+
  geom_signif(comparisons = list(c("Control","Donor"),c("Recipient","Control"),c("Recipient","Donor")), test=wilcox.test, map_signif_level=TRUE, color='black', step_increase=0.1)+
  ggtitle("Comparing alpha-diversities of microbiome at baseline (T1)")+expand_limits(y=0)
plot_richness(ps.rf%>% subset_samples(Timepoint=="T2"), x="Treatment", measures=alpha_meas, color="Treatment")+
  geom_boxplot(aes(x=Treatment, y=value, color=Treatment), alpha=0.1) +
  theme(strip.background = element_rect(colour = "black", fill = "white"), text=element_text(size=12))+
  geom_signif(comparisons = list(c("Control","Donor"),c("Recipient","Control"),c("Recipient","Donor")), test=wilcox.test, map_signif_level=TRUE, color='black', step_increase=0.1)+
  ggtitle("Comparing alpha-diversities of microbiome at before transplant (T2)")+expand_limits(y=0)
plot_richness(ps.rf%>% subset_samples(Timepoint=="T3"), x="Treatment", measures=alpha_meas, color="Treatment")+
  geom_boxplot(aes(x=Treatment, y=value, color=Treatment), alpha=0.1) +
  theme(strip.background = element_rect(colour = "black", fill = "white"), text=element_text(size=12))+
  geom_signif(comparisons = list(c("Recipient","Control")), test=wilcox.test, map_signif_level=TRUE, color='black')+
  ggtitle("Comparing alpha-diversities of microbiome 2 weeks after transplant (T3)")+expand_limits(y=0)
plot_richness(ps.rf%>% subset_samples(Timepoint=="T4"), x="Treatment", measures=alpha_meas, color="Treatment")+
  geom_boxplot(aes(x=Treatment, y=value, color=Treatment), alpha=0.1) +
  theme(strip.background = element_rect(colour = "black", fill = "white"), text=element_text(size=12))+
  geom_signif(comparisons = list(c("Recipient","Control")), test=wilcox.test, map_signif_level=TRUE, color='black')+
  ggtitle("Comparing alpha-diversities of microbiome 12 weeks after transplant (T4)")+expand_limits(y=0)

plot_richness(ps.rf%>% subset_samples(Treatment=="Recipient"), x="Timepoint", measures=alpha_meas, color="Timepoint")+
  geom_boxplot(aes(x=Timepoint, y=value, color=Timepoint), alpha=0.1) +
  theme(strip.background = element_rect(colour = "black", fill = "white"), text=element_text(size=12))+
  geom_signif(comparisons = list(c("T1","T2"),c("T2","T3"),c("T2","T4"),c("T1","T4")), test=wilcox.test, map_signif_level=TRUE, color='black', step_increase=0.1)+
  ggtitle("Comparing alpha-diversities of microbiome in recipient")+expand_limits(y=0)

plot_richness(ps.rf%>% subset_samples(Treatment=="Control"), x="Timepoint", measures=alpha_meas, color="Timepoint")+
  geom_boxplot(aes(x=Timepoint, y=value, color=Timepoint), alpha=0.1) +
  theme(strip.background = element_rect(colour = "black", fill = "white"), text=element_text(size=12))+
  geom_signif(comparisons = list(c("T1","T2"),c("T2","T3"),c("T2","T4"),c("T1","T4")), test=wilcox.test, map_signif_level=TRUE, color='black', step_increase=0.1)+
  ggtitle("Comparing alpha-diversities of microbiome in control")+expand_limits(y=0)

plot_richness(ps.rf%>% subset_samples(Treatment=="Control"|Treatment=="Recipient"), x="Timepoint", measures="Shannon", color="Treatment")+
  geom_boxplot(aes(x=Timepoint, y=value, color=Treatment), alpha=0.1) +
  theme(strip.background = element_rect(colour = "black", fill = "white"), text=element_text(size=12))+
  geom_signif(comparisons = list(c("T1","T2"),c("T2","T3"),c("T2","T4"),c("T1","T4")), test=wilcox.test, map_signif_level=TRUE, color='black', step_increase=0.1)+
  ggtitle("Comparing alpha-diversities of microbiome in Control and Recipient across timepoints")+facet_wrap(~Treatment, ncol = 2)+ylab("Shannon entropies")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


#transform to relative abundance
ps.ra = transform_sample_counts(ps1, function(x){x / sum(x)})
donor.ra<-subset_samples(ps.ra, Treatment=="Donor") %>% prune_taxa(taxa_sums(.) > 0, .)
recipient.ra<-subset_samples(ps.ra, Treatment=="Recipient")%>% prune_taxa(taxa_sums(.) > 0, .)
control.ra<-subset_samples(ps.ra, Treatment=="Control")%>% prune_taxa(taxa_sums(.) > 0, .)
#nMDS weighted Unifrac
ps.wuf<-phyloseq::distance(ps.ra, method="wunifrac")
nmds.wuf <-ordinate(ps.ra, "NMDS", ps.wuf)
plot_ordination(ps.ra, nmds.wuf, "samples", color="Treatment", shape="Timepoint") +
  geom_point(size = 7, alpha=0.7) +
  theme_light(base_size = 14) +
  labs(x="NMDS1", y="NMDS2", title="nMDS ordination of weighted UniFrac distance")

plot_ordination(donor.ra, ordinate(donor.ra,method = "NMDS", distance="wunifrac"), "samples", color="Timepoint") +
  geom_point(size = 7, alpha=0.7) +
  theme_light(base_size = 14) +
  labs(x="NMDS1", y="NMDS2", title="nMDS ordination of weighted UniFrac distance")

plot_ordination(control.ra, ordinate(control.ra,method = "NMDS", distance="wunifrac"), "samples", color="Timepoint" ) +
  geom_point(size = 7, alpha=0.7) +
  theme_light(base_size = 14) +
  labs(x="NMDS1", y="NMDS2", title="nMDS ordination of weighted UniFrac distance")

plot_ordination(recipient.ra, ordinate(recipient.ra,method = "NMDS", distance="wunifrac"), "samples", color="Timepoint" ) +
  geom_point(size = 7, alpha=0.7) +
  theme_light(base_size = 14) +
  labs(x="NMDS1", y="NMDS2", title="nMDS ordination of weighted UniFrac distance")

library("vegan")
adonis(ps.wuf ~ Treatment, data = metadata)


#phyloseq heatmap
plot_heatmap(physeq=ps.ra)

# 100%stack barplot for composition
# prune out phyla below 2% in each sample for total dataset
# split donor, Control and recipient into 3 seperate plots
ps.phylum <- ps1 %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum


tax_table(ps1) = gsub("\\[", "",tax_table(ps1))
tax_table(ps1) = gsub("\\]", "_",tax_table(ps1))
tax_table(ps1)= gsub("_$", "",tax_table(ps1))


ps.genus <- ps1 %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.001) %>%                         # Filter out low abundance taxa
  arrange(Genus)                                      # Sort data frame alphabetically by genus

unique(ps.genus$Genus) #208 levels -> some repetitive genus with []
#204 levels 
#genus plot not informative
library("ggplot2")
theme_set(theme_bw(base_size=14))
ggplot(ps.genus%>%filter(Treatment=="Recipient"), aes(x = Timepoint, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity") + ylim(0,1) +
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genus > 0.1%) \n") +
  ggtitle("Genus Composition of Bacterial Communities in Recipient Dogs") +theme(legend.position="bottom")

# Set colors for plotting
unique(ps.phylum$Phylum)
colourCount =length(unique(ps.phylum$Phylum))
library(RColorBrewer)
myColors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928")
names(myColors) <- c(as.character(unique(ps.phylum$Phylum)))
colScale <- scale_fill_manual(name = "Phylum",values = myColors)

# Plot
library("ggplot2")
theme_set(theme_bw(base_size=14))
ggplot(ps.phylum%>%filter(Treatment=="Recipient"), aes(x = Timepoint, y = Abundance, fill = Phylum)) + 
  facet_grid(rows=vars(Dog)) +
  geom_bar(stat = "identity") + ylim(0,1) +
  theme(axis.title.x = element_blank()) + 
  colScale +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition of Bacterial Communities in Recipient Dogs") 
p<-ggplot(ps.phylum%>%filter(Treatment=="Donor"), aes(x = Date, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") + ylim(0,1) + 
  theme(axis.title.x = element_blank()) + 
  colScale +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition of Bacterial Communities in the Donor Dog") 
p$data %<>%
  mutate(Date=as.Date(Date, format="%d/%m/%Y")) %>%
  mutate(Date = gsub("0016", "2016", Date)) %>%
  mutate(Date = gsub("0017", "2017", Date))
p$data %>%  mutate(Sample = gsub("Transplant","MT", Sample)) ->p$data
p+ theme(axis.text.x = element_text(angle=90))+
  geom_text(data=p$data%>% 
              filter(Phylum=="Bacteroidetes"), 
                     aes(label=Sample), angle=90,position = position_stack(vjust = 0.5))

ggplot(ps.phylum%>%filter(Treatment=="Recipient"), aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") + ylim(0,1) +
  theme(axis.title.x = element_blank()) + 
  colScale +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition of Bacterial Communities in Recipient Dogs") + theme(axis.text.x = element_text(angle=90))

#nMDS BrayCurtis
#ordination
ps.bc<-phyloseq::distance(ps.ra, method="bray")
nmds.bc <-ordinate(ps.ra, "NMDS", ps.bc)
plot_ordination(ps.ra, nmds.bc, type="taxa", color="Phylum", title="nMDS ordination showing the bacterial taxa")+geom_abline(intercept = 0.5, slope = -2.5, colour ="grey")+ facet_wrap(~Phylum,6)
plot_ordination(ps.ra, nmds.bc, "samples", color="Treatment", shape="Timepoint") +
  geom_point(size = 7, alpha=0.7) +
  theme_light(base_size = 14) +
  scale_x_continuous(breaks= waiver()) +
  #geom_segment(aes(x = -0.1, y = 0.5, xend = 0.3, yend = -0.25), colour ="gray")+
  geom_abline(intercept = 0.5, slope = -2.5, colour ="grey") +
  labs(x="NMDS1", y="NMDS2", title="nMDS ordination of Bray-Curtis dissimilarity")
plot_ordination(ps.ra, nmds.bc, type="biplot", color="Phylum", shape="Treatment", label="Timepoint", title="nMDS ordination showing both Sample and Taxa") 

#remove taxa irrelevant to transplant (absent in all recipient samples or absent in donor samples)
# Define phyla to filter
filterPhyla = c("Acidobacteria", "Chloroflexi", "Cyanobacteria/Chloroplast","Elusimicrobia","Synergistetes", "WPS-2","NA")
#Remove the features (ASVs) with ambigurous phylum annotation and unidentified Phylum
ps2 <- subset_taxa(ps1, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) %>%
       subset_taxa(!Phylum %in% filterPhyla) #1687 taxa

#ordination
ps2.ra = transform_sample_counts(ps2, function(x){x / sum(x)})
ps2.bc<-phyloseq::distance(ps2.ra, method="bray")
nmds2.bc <-ordinate(ps2.ra, "NMDS", ps2.bc)
plot_ordination(ps2.ra, nmds2.bc, type="taxa", color="Phylum", title="nMDS ordination showing the bacterial taxa")+geom_abline(intercept = 0.5, slope = -2.5, colour ="grey")+ facet_wrap(~Phylum,4)
plot_ordination(ps2.ra, nmds2.bc, "samples", color="Treatment", shape="Timepoint") +
  geom_point(size = 7, alpha=0.7) +
  theme_light(base_size = 14) +
  scale_x_continuous(breaks= waiver()) +
  #geom_segment(aes(x = -0.1, y = 0.5, xend = 0.3, yend = -0.25), colour ="gray")+
  geom_abline(intercept = 0.5, slope = -2.5, colour ="grey") +
  labs(x="NMDS1", y="NMDS2", title="nMDS ordination of Bray-Curtis dissimilarity")

#from Raymond's script
BOP_dat <-read_tsv("//Users/OB/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/model_picrust/bop_ppd.txt") %>%
  filter(Tooth == "BOP") %>%
  rename(SampleID = Sample) %>%
  select(SampleID, BOP)

PPD_dat <- read_tsv("//Users/OB/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/model_picrust/bop_ppd.txt") %>%
  filter(Tooth == "PPD") %>%
  rename(SampleID = Sample) %>%
  select(SampleID, PPD)

plaque <- read_tsv("//Users/OB/Desktop/Yuki_OBMac/Canine_MiSeq/Analysis/model_picrust/plaque.txt") %>%
  select(SampleID, plaque)

#BOPplot
BOP_dat %>%
  left_join(map_dat, by = c("SampleID" = "Sample")) %>%
  filter(Treatment != "Donor") %>%
  mutate(Timepoint = factor(Timepoint)) %>%
  mutate(Timepoint = `levels<-`(Timepoint, c(0, 2, 4, 14))) %>%
  mutate(Timepoint = as.numeric(as.character(Timepoint))) %>%
  {
    ggplot(., aes(Timepoint, BOP, color = Treatment, group = Dog)) +
      geom_line() +
      geom_text_repel(data = filter(., Timepoint == 14), aes(Timepoint, BOP, label = Dog), segment.color = "black") +
      facet_grid(cols = vars(Treatment)) +
      scale_color_manual(values = col[-1]) +
      scale_x_continuous(limits=c(0,14), breaks=seq(0,14,2)) + 
      labs(x = "Week", y = "% BOP") +
      guides(colour = "none")
  }
#PPD plot
PPD_dat %>%
  left_join(map_dat, by = c("SampleID" = "Sample")) %>%
  filter(Treatment != "Donor") %>%
  mutate(Timepoint = factor(Timepoint)) %>%
  mutate(Timepoint = `levels<-`(Timepoint, c(0, 2, 4, 14))) %>%
  mutate(Timepoint = as.numeric(as.character(Timepoint))) %>%
  {
    ggplot(., aes(Timepoint, PPD, color = Treatment, group = Dog)) +
      geom_line() +
      facet_grid(cols = vars(Treatment)) +
      geom_text_repel(data = filter(., Timepoint == 14), aes(Timepoint, PPD, label = Dog), segment.color = "black") +
      scale_color_manual(values = col[-1]) + 
      scale_x_continuous(limits=c(0,14), breaks=seq(0,14,2)) +
      labs(x = "Week", y = "mean PPD") + 
      guides(colour = "none")
  }
#Plaque% plot
plaque %>%
  left_join(map_dat, by = c("SampleID" = "Sample")) %>%
  filter(Treatment != "Donor") %>%
  mutate(Timepoint = factor(Timepoint)) %>%
  mutate(Timepoint = `levels<-`(Timepoint, c(0, 2, 4, 14))) %>%
  mutate(Timepoint = as.numeric(as.character(Timepoint))) %>%
  {
    ggplot(., aes(Timepoint, plaque, color = Treatment, group = Dog)) +
      geom_line() +
      facet_grid(cols = vars(Treatment)) +
      geom_text_repel(data = filter(., Timepoint == 14), aes(Timepoint, plaque, label = Dog), segment.color = "black") +
      scale_color_manual(values = col[-1]) + 
      scale_x_continuous(limits=c(0,14), breaks=seq(0,14,2)) +
      labs(x = "Week", y = "% plaque") + 
      guides(colour = "none")
  }
# Constrained Ordinations
# how environmental variables (BOP and PPD) are associated with the changes in community composition
# add the BOP and PPD to the phyloseq data
sample_data(ps2) %>%
  mutate(SampleID = rownames(.))%>%
  full_join(BOP_dat, by =  "SampleID")%>%
  left_join(PPD_dat, by =  "SampleID") %>%
  left_join(plaque, by =  "SampleID") %>%
  column_to_rownames(var = "SampleID") ->sample_data(ps2) 
# Remove data points with missing metadata; remove T1 samples
ps2_not_na <- ps2 %>%
  subset_samples(!is.na(BOP) & !is.na(PPD) & Timepoint!="T1")
bray_not_na <- phyloseq::distance(physeq = ps2_not_na, method = "bray")

# CAP ordinate
cap_ord <- ordinate(physeq = ps2_not_na, method = "CAP", distance = bray_not_na, formula = ~ BOP+PPD+plaque)

# CAP plot
cap_plot <- plot_ordination(physeq = ps2_not_na, ordination = cap_ord, color = "Treatment", axes = c(1,2)) + 
  aes(shape = Timepoint) + scale_color_manual(values = col[-1]) +
  geom_point(aes(colour = Treatment), alpha = 0.7, size = 4) +geom_text(aes(label=sample_names(ps2_not_na)),nudge_y = 0.1)
# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, yend = CAP2, x = 0, y = 0, shape = NULL, color = NULL, label = labels)
label_map <- aes(x = 1.3 * CAP1, y = 1.3 * CAP2, shape = NULL, color = NULL, label = labels)
arrowhead = arrow(length = unit(0.02, "npc"))
# Make a new graphic
cap_plot + geom_segment(mapping = arrow_map, size = .5, data = arrowdf, color = "black", arrow = arrowhead) + 
  geom_text(mapping = label_map, size = 4,  data = arrowdf, show.legend = FALSE )

