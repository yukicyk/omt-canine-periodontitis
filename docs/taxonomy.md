# Taxonomy Assignment

(*Excerpt from internal report dated 2018-08-25, with modern updates noted*)

The DADA2 package provides a native implementation of the naive Bayesian classifier method (Wang et al., 2007) for taxonomy assignment. It takes as input a set of sequences to be classified and a training set of reference sequences with known taxonomy, and it outputs taxonomic assignments with at least `minBoot` (default 50) bootstrap confidence. The `dada2` package also implements a method to make species-level assignments based on exact matching between Amplicon Sequence Variants (ASVs) and sequenced reference strains.

We used a curated reference sequence database that comprised of the `RDP trainset 16` (Cole et al., 2014), `eHOMD 16S rRNA RefSeq Version 15.1` (Chen et al., 2010), and the `full set of canine oral microbiome 16S rRNA gene sequences` (Dewhirst et al., 2012). `Mothur` (Schloss et al., 2009) was used to process the reference sequences into a DADA2-compatible format.

Following the genus-level classification, species-level assignments were made based on the principle of exact matching (100% identity) between ASVs and the species-level reference sequences from the combined database. This stringent approach is the recommended standard for assigning species to 16S gene fragments (Edgar, 2018). For this high level of species-strain identification and matching, the source of the reference sequences in the database is critical. Since the project collects oral plaque samples from beagle dogs, the most relevant references are from the Canine Oral Taxa (COT) set. The Human Oral Microbiome Database (eHOMD, or "HOT") set, being the gold standard for dentistry, was also included for its high relevance to the oral microbiota context. During classification, there was a preference in assignation: when an ASV had tied exact matches across databases, the assignment was given in the order of COT > eHOMD > RDP.

Notably, the SILVA database was deliberately not used in this workflow in 2019 due to (1) its incompatible taxonomic hierarchy (e.g., the superphylum Patescibacteria) with the other reference databases, which ensured consistency in classification across the combined dataset, and (2) the fact that most of its sequences originated from irrelevant environments, while the oral niche databases were more pertinent to this project.

> **_Update Note (added 2025, not part of the original 2018 document):_**
> It is now acknowledged that the RDP database is obsolete. For any future or reproduced studies, the huge and well-maintained SILVA database provides the most robust reference for better taxonomic resolution. It is recommended to use an updated SILVA dataset, as mentioned in the [`data/data_README.md`](./data/data_README.md), for any further work. For maximum accuracy, specialized databases like COT and eHOMD (and any other new, well-documented canine oral 16S rRNA sequences) should also be formatted and combined with the SILVA dataset.

### References

-   Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. *Applied and environmental microbiology*, 73(16), 5261-5267.
-   Schloss, P.D., et al. (2009). Introducing mothur: open-source, platform-independent, community-supported software for describing and comparing microbial communities. *Applied and environmental microbiology*, 75(23), pp.7537-7541.
-   Edgar, R. C. (2018). Updating the 97% identity threshold for 16S ribosomal RNA OTUs. *Bioinformatics*, 34(14), 2371-2375.
-   Cole, J. R., et al. (2014). Ribosomal Database Project: data and tools for high throughput rRNA analysis. *Nucleic Acids Research*, 42(Database issue):D633-D642.
-   Chen, T., et al. (2010). The Human Oral Microbiome Database: a web accessible resource for investigating oral microbe taxonomic and genomic information. *Database*, Vol. 2010, Article ID baq013.
-   Dewhirst, F. E., et al. (2012). The canine oral microbiome. *PloS one*, 7(4), e36067.