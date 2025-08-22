#==============================
# 01_preprocessing.sh
#==============================
# Step 1: Raw Sequence Pre-processing
# Source Files: `Canine_MiSeq_LOG copy.txt`, `Cutadapt.log.txt`
# Description: This shell script uses `cutadapt` to trim the 341F and 806R 16S rRNA gene primers from the raw, paired-end `.fq.gz` files.
# Inputs: Raw `.fq.gz` files for forward and reverse reads.
# Outputs: Trimmed `.fq.gz.out` files for each sample, ready for DADA2.
# Requirement: cutadapt 1.14; Python 3.6.2rc1

# Trim forward reads (R1)
find . -name "*_1.fq.gz" -exec cutadapt -g ACTCCTACGGGAGGCAGCAG -o {}.out {} \;

# Trim reverse reads (R2)
find . -name "*_2.fq.gz" -exec cutadapt -g GGACTACHVGGGTWTCTAAT -o {}.out {} \;