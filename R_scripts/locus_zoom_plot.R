# Add dplyr
library(dplyr)

# Clone repository with LocusZoom functions (https://github.com/Geeketics/LocusZooms)
# The directory called LocusZooms-master will be created
# Then add necessary function via source:
source("/LocusZooms-master/functions/locus_zoom.R")

# Set path to files for analysis:
current_path = "your_path"
setwd(current_path)

# Upload FINEMAP .cred1 output file - this file contains 95% credible set of SNPs
# We skip first 5 rows because this is a specific header
finemap_result <- read.table("your_file.cred1", sep = " ", header = T, skip = 5) 

# Upload GWAS sumstats
# GWAS sumstats file should contain only columns:
# "CHR", "SNP", "BP", "P" - cromosome, SNP id, base pair, p-value
GWAS <- read.table("your_file.tsv", sep = "\t", header = T)

# Upload LD reference matrix 
LD <- read.table("your_file.ld", stringsAsFactors = FALSE, header = TRUE)

# Upload a table with genes coordinates in GRCh37
# This file is provided with LocusZoom function
UCSC_GRCh37_Genes_UniqueList.txt <- read.delim("UCSC_GRCh37_Genes_UniqueList.txt", stringsAsFactors = FALSE, header = TRUE)

# Select from GWAS data only SNPs that are present in LD matrix
# This step is optional, in this particular case we calculated LD matrix from different cohort (1000 Genomes), so not all the SNPs are common between GWAS data and LD matrix
GWAS <- GWAS %>% 
  dplyr::filter(SNP %in% LD$SNP_B | SNP %in% LD$SNP_A)

# Create locus zoom plot, SNPs in 95% credible set will be highlighted with red circles
locus.zoom(data = GWAS, # GWAS data
           region = c(GWAS$CHR, min(GWAS$BP), max(GWAS$BP)), # genome region
           offset_bp = 0, # offset
           ld.file = LD, # LD data
           genes.data = UCSC_GRCh37_Genes_UniqueList.txt, # a list of genes with coordinates
           plot.title = "Association of CTSB with PD in Europeans", # plot title
           file.name = "GWAS_FINEMAP_locus_plot.pdf", # file name
           secondary.snp = finemap_result$cred1, # FINEMAP SNPs
           secondary.label = TRUE) # SNPs in 95% credible set will be highlighted with red

