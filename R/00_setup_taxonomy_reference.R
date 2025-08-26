#
# SCRIPT: 00_setup_taxonomy_reference.R
#
# PURPOSE: Prepare a customized taxonomy reference database for dada2. This script
#          parses, cleans, and combines three sources:
#          1. Canine Oral Taxa (COT)
#          2. Human Oral Microbiome Database (eHOMD)
#          3. RDP Trainset 16
#          This is a one-time setup step.
#
# INPUTS:
#   - Raw FASTA and taxonomy files from RDP, HOMD, and COT, located in:
#     - `data/raw_taxonomy/COT.fasta`
#     - `data/raw_taxonomy/COT_taxonomy.tsv`
#     - `data/raw_taxonomy/HOMD.fasta`
#     - `data/raw_taxonomy/HOMD_taxonomy.tsv`
#     - `data/raw_taxonomy/rdp_train_set_16.fa`
#     - `data/raw_taxonomy/rdp_species_assignment_16.fa`
#
# OUTPUTS:
#   - `data/taxonomy/custom_train_set.fa`: Combined training set for assignTaxonomy.
#   - `data/taxonomy/custom_species_assignment.fa`: Combined set for addSpecies.
#

# --- 1. Load Libraries and Set Paths ---

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, readr, tibble, tidyr)

# Define relative paths for clarity and reproducibility
raw_tax_path <- "data/raw_taxonomy"
output_path <- "data/taxonomy"

# Create output directory if it doesn't exist
if (!dir.exists(output_path)) dir.create(output_path)

# --- 2. Process Canine Oral Taxa (COT) ---

print("Processing Canine Oral Taxa (COT)...")

# Load the raw multi-line FASTA and taxonomy mapping file
cot_fasta_raw <- read_lines(file.path(raw_tax_path, "COT.fasta"))
cot_tax_map <- read_tsv(file.path(raw_tax_path, "COT_taxonomy.tsv"), col_names = FALSE)

# -- Parse the multi-line FASTA format --
# Find the line numbers for headers (lines starting with '>')
header_lines <- grep(">", cot_fasta_raw)
# Create a dataframe to define the start and end of each sequence
seq_pos <- data.frame(
  header = cot_fasta_raw[header_lines],
  start = header_lines + 1,
  end = c(header_lines[-1] - 1, length(cot_fasta_raw))
)
# Use the positions to paste the sequence lines together
sequences <- apply(seq_pos, 1, function(row) {
  paste(cot_fasta_raw[row[["start"]]:row[["end"]]], collapse = "")
})
# Create a clean FASTA dataframe
cot_fasta_df <- data.frame(id = seq_pos$header, sequence = sequences)

# -- Combine FASTA with taxonomy and format for dada2 --
cot_combined <- cot_fasta_df %>%
  # Join with the taxonomy mapping file
  left_join(cot_tax_map, by = c("id" = "X1")) %>%
  rename(taxonomy_string = X2) %>%
  # Ensure the taxonomy string starts with the root
  mutate(taxonomy_string = paste0("Root;", taxonomy_string))

# Create the training set (6 levels: Kingdom to Genus)
cot_train_set <- cot_combined %>%
  mutate(dada2_header = sub("^(.*?;.*?;.*?;.*?;.*?;.*?);.*", "\\1;", taxonomy_string)) %>%
  select(dada2_header, sequence)

# Create the species assignment set (Genus species)
cot_species_assign <- cot_combined %>%
  mutate(dada2_header = sub(".*?;(.*?);(.*?)$", "\\1 \\2", taxonomy_string)) %>%
  select(dada2_header, sequence)


# --- 3. Process Human Oral Microbiome Database (HOMD) ---

print("Processing Human Oral Microbiome Database (HOMD)...")

# Load the raw two-line FASTA and taxonomy mapping file
homd_fasta_raw <- read_lines(file.path(raw_tax_path, "HOMD.fasta"))
homd_tax_map <- read_tsv(file.path(raw_tax_path, "HOMD_taxonomy.tsv"))

# -- Parse the two-line FASTA format --
# The matrix trick converts the vector of lines into a 2-column df
homd_fasta_df <- matrix(homd_fasta_raw, ncol = 2, byrow = TRUE) %>%
  as.data.frame() %>%
  rename(header = V1, sequence = V2) %>%
  # Extract the HMT-ID from the complex header for joining
  mutate(HMT_ID = sub(".*\\|HMT-([0-9]+).*", "\\1", header))

# -- Combine FASTA with taxonomy and format for dada2 --
homd_combined <- homd_fasta_df %>%
  left_join(homd_tax_map, by = "HMT_ID") %>%
  # Construct the full taxonomy string required by dada2
  mutate(
    taxonomy_string = paste("Root", Domain, Phylum, Class, Order, Family, Genus, sep = ";"),
    species = Species.x # Use the species name from the FASTA header
  )

# Create the training set (Genus level)
homd_train_set <- homd_combined %>%
  mutate(dada2_header = paste0(taxonomy_string, ";")) %>%
  select(dada2_header, sequence)

# Create the species assignment set
homd_species_assign <- homd_combined %>%
  mutate(dada2_header = paste(Genus, species)) %>%
  select(dada2_header, sequence)


# --- 4. Process RDP Reference Files ---

print("Processing RDP reference files...")

# RDP files are often already in the correct dada2 format. We just need to load them.
# Function to parse a standard dada2-formatted reference FASTA
parse_dada2_fasta <- function(filepath) {
  raw_lines <- read_lines(filepath)
  fasta_df <- matrix(raw_lines, ncol = 2, byrow = TRUE) %>%
    as.data.frame() %>%
    rename(dada2_header = V1, sequence = V2)
  return(fasta_df)
}

rdp_train_set <- parse_dada2_fasta(file.path(raw_tax_path, "rdp_train_set_16.fa"))
rdp_species_assign <- parse_dada2_fasta(file.path(raw_tax_path, "rdp_species_assignment_16.fa"))

# Remove the ">" from the headers for consistency before merging
rdp_train_set$dada2_header <- sub(">", "", rdp_train_set$dada2_header)
rdp_species_assign$dada2_header <- sub(">", "", rdp_species_assign$dada2_header)


# --- 5. Combine All Sources and Write Final Files ---

print("Combining all sources and writing final reference files...")

# Combine the training sets from all three sources
full_train_set <- bind_rows(cot_train_set, homd_train_set, rdp_train_set) %>%
  # Remove any duplicate sequences, keeping the first instance
  distinct(sequence, .keep_all = TRUE) %>%
  # Add the ">" back to the header for the final FASTA format
  mutate(dada2_header = paste0(">", dada2_header))

# Combine the species assignment sets
full_species_assign <- bind_rows(cot_species_assign, homd_species_assign, rdp_species_assign) %>%
  distinct(sequence, .keep_all = TRUE) %>%
  mutate(dada2_header = paste0(">", dada2_header))

# -- Write the final files in the correct two-line FASTA format --
# We interleave the header and sequence rows to create the final output

# Write the training set
train_set_output <- full_train_set %>%
  pivot_longer(everything(), names_to = "type", values_to = "line") %>%
  select(line)
write.table(train_set_output,
            file = file.path(output_path, "custom_train_set.fa"),
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# Write the species assignment set
species_assign_output <- full_species_assign %>%
  pivot_longer(everything(), names_to = "type", values_to = "line") %>%
  select(line)
write.table(species_assign_output,
            file = file.path(output_path, "custom_species_assignment.fa"),
            col.names = FALSE, row.names = FALSE, quote = FALSE)

print("---------------------------------------------------------")
print("Taxonomy reference setup complete.")
print(paste("Final files saved in:", output_path))
print("---------------------------------------------------------")