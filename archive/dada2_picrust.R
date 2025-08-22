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