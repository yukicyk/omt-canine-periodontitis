#### **Step 3: Diversity and Ordination Analysis**
## Script: `R/03_diversity_and_ordination_analysis.R`
## Source Files: `canine_analysisv2.R`
## Description: This script handles the core ecological analyses presented in the paper. It calculates alpha diversity metrics (Shannon, Observed ASVs) and performs beta diversity analysis using Bray-Curtis and weighted/unweighted UniFrac distances. It generates ordination plots (PCoA, NMDS) to visualize community differences (Figures 2 & 3 in the paper) and performs PERMANOVA tests to assess statistical significance.
## Inputs: `phyloseq_object.rds`.
## Outputs: Alpha and beta diversity plots, ordination plots, and statistical results tables.
# OTU: Operational Taxonomic Unit, ASV: Amplicon Sequence Variant

# Import phyloseq file from step 02.
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")

#Plot rarefraction curves for all samples
#load the phyloseq extension richness.R functions
library("vegan", "ggplot2")
p.rare <- ggrare(ps1, step = 1000, color = "Treatment",label="Sample", se = FALSE)
p.rare + xlim(0, 25000)
p.all <- p.rare + facet_wrap(~Treatment)+ xlim(0, 25000)
p.all

#venn diagram showing number of ASVs shared among 3 groups
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


# Plot Shannon index
plot_richness(ps1, x="Dog", measures=c("Shannon"), color="Treatment", shape="Timepoint") + 
  geom_point(size=7, alpha=0.75)+ ggtitle("Alpha Diversity measure Shannon Index")
plot_richness(ps1, x="Dog", measures=c("Observed"), color="Treatment", shape="Timepoint") + 
  geom_point(size=7, alpha=0.75)+ ggtitle("Alpha Diversity measure Observed ASVs")
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
  ggtitle("Comparing alpha-diversities of microbiome at T2 before transplant")+expand_limits(y=0)
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


# UniFrac weighted and unweighted distance calculation and ordination
wuf = ordinate(ps1, "PCoA", "unifrac", weighted=TRUE)
uwuf = ordinate(ps1, "PCoA", "unifrac", weighted=FALSE)
bc =ordinate(ps1, "PCoA", "bray")
plot_ordination(ps1, bc, color="Treatment", shape="Timepoint", label = "Dog") +
  geom_point(size=3) +
  ggtitle("MDS/PCoA on Bray-Curtis dissimilarity")


library(tibble)
# F1000 workflow Filtering of data based on prevalence at the phylum level

# Export phyla data for manual check of data
ps1p <- tax_glom(ps1, "Phylum")
write.table(otu_table(ps1p),"phylum.txt", sep="\t", col.names = tax_table(ps1p)[, "Phylum"])
write.table(taxa_sums(ps1p),"phylum_sums.txt", sep="\t", row.names = tax_table(ps1p)[, "Phylum"])

# Create table, number of features (ASVs) for each phylum
table(tax_table(ps1)[, "Phylum"], exclude = NULL) 

# Remove the features (ASVs) with ambigurous phylum annotation
ps0 <- subset_taxa(ps1, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#Can we further classify the ambiguorous sequences?
psna <- subset_taxa(ps1, is.na(Phylum))
all_ASV_seq <- read.fasta(file.choose())
na_seq_id <- taxa_names(psna)
names(all_ASV_seq) %in% na_seq_id
na_seq <- all_ASV_seq[c(which(names(all_ASV_seq) %in% na_seq_id))]
write.fasta(sequences=na_seq, names=names(na_seq), file.out="na_seq.fasta")

# Define prevalence of each taxa
# (in how many samples does each taxa appear at least once)
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))

# Compute the total and average prevalences of the features in each phylum.
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Supervised Prevalence Filtering
# Define phyla to filter (sum of prevalence <10)
filterPhyla = c("Acidobacteria", "Chloroflexi", "Cyanobacteria/Chloroplast")

# Filter entries with unidentified Phylum.
ps2 = subset_taxa(ps0, !Phylum %in% filterPhyla)
ps2 #1739 taxa

# Unsupervised prevalence filtering
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps2, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps0),color=Phylum)) +
  
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps0)
prevalenceThreshold
## [1] 4.35

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps3 = prune_taxa(keepTaxa, ps0)
ps3 #811 taxa

# Export table
ps3p <- tax_glom(ps3, "Phylum")
write.table(otu_table(ps3p),"ps3.txt", sep="\t", col.names = tax_table(ps3p)[, "Phylum"])
write.table(taxa_sums(ps3p),"phylum_sums_ps3.txt", sep="\t", row.names = tax_table(ps3p)[, "Phylum"])

# Plot Shannon index and Observed ASVs dot plots
plot_richness(ps3, x="Dog", measures=c("Shannon"), color="Treatment", shape="Timepoint") + 
  geom_point(size=7, alpha=0.75)+ ggtitle("Alpha Diversity measure Shannon Index")
plot_richness(ps3, x="Dog", measures=c("Observed"), color="Treatment", shape="Timepoint") + 
  geom_point(size=7, alpha=0.75)+ ggtitle("Alpha Diversity measure Observed ASVs")

# ordination plot samples
ord = ordinate(ps3, "NMDS", "unifrac", weighted=TRUE)
ord2 = ordinate(ps3, "NMDS", "unifrac", weighted=FALSE)
ord3 =ordinate(ps3, "NMDS", "bray")
ord4 = ordinate(ps3, "MDS","dpcoa")
plot_ordination(ps3, ord4, color="Treatment", shape="Timepoint", label = "Dog") +
  geom_point(size=3, alpha=0.5) +
  ggtitle("MDS on DPCoA, prevalence filtered data")

# Test various distance method for MDS ordination on ps3 data
dist_methods <- unlist(distanceMethodList)

# Remove the user-defined distance
dist_methods = dist_methods[-which(dist_methods=="ANY")]
dist_methods = dist_methods[-which(dist_methods=="manhattan")]

# Loop through each distance method, save each plot to a list, called plist
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  distance(ps3s, method=i) %>%
    ordinate(ps3s, "MDS", distance=.) %>%
    # Create plot, store as temp variable, p
    plot_ordination(ps3s, ., color="Treatment", shape="Timepoint") -> p
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

# Combine results
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=Treatment, shape=Timepoint))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for prevalence filtered data (Species)")
p

# Obtain the axis information (variance%)
df2 = laply(plist, function(x) x$labels)

# Ordination plot taxa
theme_set(theme_bw())
phyla.ord <- ordinate(ps3, "NMDS", "bray")
plot_ordination(ps3, phyla.ord, type="taxa", color="Phylum", title="Taxa NMDS on Bray-Curtis Dissimilartiy, prevalence filtered data") + 
  facet_wrap(~Phylum, 6) 

ps2p <- tax_glom(ps2, "Phylum")
ps2c <- tax_glom(ps2, "Class")
ps2o <- tax_glom(ps2, "Order")
ps2f <- tax_glom(ps2, "Family")
ps2g <- tax_glom(ps2, "Genus")
ps2s <- tax_glom(ps2, "Species")

ps3p <- tax_glom(ps3, "Phylum")
ps3c <- tax_glom(ps3, "Class")
ps3o <- tax_glom(ps3, "Order")
ps3f <- tax_glom(ps3, "Family")
ps3g <- tax_glom(ps3, "Genus")
ps3s <- tax_glom(ps3, "Species")

# Heatmap
# Phyloseq heatmap
order_of_samples <- c("MT-B-1","MT-D-1","MT-G-1","MT-I-1","MT-K-1","MT-M-1","MT-O-1","MT-R-1","MT-S-1","MT-B-2","MT-D-2","MT-G-2","MT-I-2","MT-K-2","MT-M-2","MT-O-2","MT-R-2","MT-S-2","MT-B-3","MT-D-3","MT-G-3","MT-I-3","MT-K-3","MT-M-3","MT-O-3","MT-R-3","MT-S-3","MT-B-4","MT-D-4","MT-G-4","MT-I-4","MT-K-4","MT-M-4","MT-O-4","MT-R-4","MT-S-4","Transplant-A","Transplant-C","Transplant-F","Transplant-H","Transplant-J","Transplant-L","Transplant-N","Transplant-P","Transplant-Q","MT-E-I-1","MT-E-II-1","MT-E-III-1","MT-E-IV-1","MT-E-V-1","MT-A-1","MT-C-1","MT-F-1","MT-H-1","MT-J-1","MT-L-1","MT-N-1","MT-P-1","MT-Q-1","MT-A-2","MT-C-2","MT-F-2","MT-H-2","MT-J-2","MT-L-2","MT-N-2","MT-P-2","MT-Q-2","MT-A-3","MT-C-3","MT-F-3","MT-H-3","MT-J-3","MT-L-3","MT-N-3","MT-P-3","MT-Q-3","MT-A-4","MT-C-4","MT-F-4","MT-H-4","MT-J-4","MT-L-4","MT-N-4","MT-P-4","MT-Q-4")
plot_heatmap(ps2p, taxa.label = "Phylum", low="#4575b4", high="#d73027", sample.label = "Treatment", sample.order = order_of_samples)
p = plot_heatmap(ps3, sample.label = "Treatment", sample.order = order_of_samples)
# p + geom_vline (aes(xintercept =9 ,color="white"))
heatmap(otu_table(ps3))

# Superheat customized heatmap
library(superheat)
ps2_phyla <- read.csv(file.choose(), header=TRUE)
phyla_heatmap <-data.matrix(ps2_phyla[,2:ncol(ps2_phyla)])
rownames(phyla_heatmap) <- ps2_phyla[,1]
# Transform to relative abundance (proportional transformation)
phyla_heatmapR = prop.table(phyla_heatmap, margin = 2)
# Prepare mapping labels for the samples
Treatment <- rep("Control",times=36)
Treatment <-append(Treatment, rep("Donor", times=14))
Treatment <-append(Treatment, rep("Recipient", times=36))
Timepoint <- rep(c("T1","T2","T3","T4","T1"),each=9)
Timepoint <- append(Timepoint, rep(c("T2"),each=5))
Timepoint <- append(Timepoint, rep(c("T1","T2","T3","T4"),each=9))
sample_reads_sum <- c(21326,21381,21324,21666,23246,20485,22120,19874,25030,24055,20221,19090,20327,22043,20505,19698,20636,20441,23569,20919,19538,21091,23436,21587,20178,20535,20127,22875,20368,22783,22214,20961,19562,28986,20178,19579,19070,19030,29158,18614,21314,21476,20902,20442,19784,20725,19876,19000,19416,19545,20445,20350,23681,21837,22715,19864,20586,20873,21553,20798,25740,19841,21461,21650,20978,19635,20557,24177,24512,23427,20756,20464,25361,27041,19366,19983,20853,19980,21563,20552,22596,24029,21734,24779,20067,19333)
phyla_reads_sum <- c(350,870,9593,15039,25251,37809,43700,45025,48369,93099,104684,317726,494525,607403)

superheat(phyla_heatmapR, scale=FALSE, 
          heat.pal = c("#4575b4", "#91bfdb", "#91cf60","#fee090", "#fc8d59", "#d73027" ),  heat.lim = c(0, 1), #define color scale
          heat.pal.values = c(0, 0.1,0.2, 0.4, 0.6, 1),
          left.label.text.alignment = "right",  
          left.label.text.size = 4, grid.hline = FALSE, 
          left.label.col = "white", bottom.label = "none",
          row.title = "Phylum", 
          grid.vline.col = "white", 
          yt= colSums(phyla_heatmap), yt.plot.type = "bar",  
          yt.axis.name = "No. of reads\nper sample",
          yr = phyla_reads_sum  / 1000, yr.axis.name = "Total no. of reads \n(in thousands)", 
          yr.plot.type = "bar")

#transform to relative abundance
ps.ra = transform_sample_counts(ps1, function(x){x / sum(x)})

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

#MDS/pcoA on Bray-Curtis
plot_ordination(control.ra, ordinate(control.ra,method = "MDS", distance="bray"), "samples", color="Timepoint" ) +
  geom_point(size = 7, alpha=0.3) +
  theme_light(base_size = 14) +
  labs(x="MDS1", y="MDS2", title="MDS/PCoA ordination of Bray-Curtis Dissimilarities")+
  stat_ellipse(type = "t", linetype = 2,size=1, show.legend = FALSE)+
  ggrepel::geom_text_repel(aes(label=sample_names(control.ra)), show.legend = FALSE)

plot_ordination(recipient.ra, ordinate(recipient.ra,method = "MDS", distance="bray"), "samples", color="Timepoint" ) +
  geom_point(size = 7, alpha=0.3) +
  theme_light(base_size = 14) +
  labs(x="MDS1", y="MDS2", title="MDS/PCoA ordination of Bray-Curtis Dissimilarities")+
  stat_ellipse(type = "t", linetype = 2, size=1,show.legend = FALSE)+
  ggrepel::geom_text_repel(aes(label=sample_names(recipient.ra)), show.legend = FALSE)

#plot only donor and recipients
dr<-subset_samples(ps.ra, Treatment!="Control") %>% prune_taxa(taxa_sums(.) > 0, .)
plot_ordination(dr, ordinate(dr,method = "MDS", distance="bray"), "samples", color="Timepoint", shape="Treatment") +
  geom_point(size = 7, alpha=0.7) +
  theme_light(base_size = 16) +
  labs(x="MDS1", y="MDS2", title="MDS/PCoA ordination of Bray-Curtis Dissimilarities")+
  stat_ellipse(type = "t", linetype = 2, show.legend = FALSE)+
  ggrepel::geom_text_repel(aes(label=sample_names(dr)), show.legend = FALSE)
#plot only donor and controls
dc<-subset_samples(ps.ra, Treatment!="Recipient") %>% prune_taxa(taxa_sums(.) > 0, .)
dc.bc<-phyloseq::distance(dc, method="bray")
dc.bc.pcoa <-ordinate(dc, "MDS", dc.bc)
plot_ordination(dc, dc.bc.pcoa, "samples", color="Timepoint", shape="Treatment") +
  geom_point(size = 7, alpha=0.7) +
  theme_light(base_size = 16) +
  labs(x="MDS1", y="MDS2", title="MDS/PCoA ordination of Bray-Curtis Dissimilarities")+
  stat_ellipse(type = "t", linetype = 2, show.legend = FALSE)+
  ggrepel::geom_text_repel(aes(label=sample_names(dc)), show.legend = FALSE)
#plot recipients and controls at T1 & T2
rc12<-subset_samples(ps.ra, Treatment!="Donor" & (Timepoint=="T1"| Timepoint=="T2")) %>% prune_taxa(taxa_sums(.) > 0, .)
rc12.bc<-phyloseq::distance(rc12, method="bray")
rc12.bc.pcoa <-ordinate(rc12, "MDS", rc12.bc)
plot_ordination(rc12, rc12.bc.pcoa, "samples", color="Treatment", shape="Timepoint") +
  geom_point(size = 7, alpha=0.7) +
  theme_light(base_size = 14) +
  labs(x="MDS1", y="MDS2", title="MDS/PCoA ordination of Bray-Curtis Dissimilarities")+
  stat_ellipse(type = "t", linetype = 2, show.legend = FALSE)+
  ggrepel::geom_text_repel(aes(label=sample_names(rc12)), show.legend = FALSE)

#plot recipients and controls at T3&T4
rc34<-subset_samples(ps.ra, Treatment!="Donor" & (Timepoint=="T3"| Timepoint=="T4")) %>% prune_taxa(taxa_sums(.) > 0, .)
rc34.bc<-phyloseq::distance(rc34, method="bray")
rc34.bc.pcoa <-ordinate(rc34, "MDS", rc34.bc)
plot_ordination(rc34, rc34.bc.pcoa, "samples", color="Treatment", shape="Timepoint") +
  geom_point(size = 7, alpha=0.7) +
  theme_light(base_size = 14) +
  labs(x="MDS1", y="MDS2", title="MDS/PCoA ordination of Bray-Curtis Dissimilarities")+
  stat_ellipse(type = "t", linetype = 2, show.legend = FALSE)+
  ggrepel::geom_text_repel(aes(label=sample_names(rc34)), show.legend = FALSE)

