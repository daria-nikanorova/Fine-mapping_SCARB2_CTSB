#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Load packages 
library(coloc)
library(data.table)
library(colochelpR)
library(dplyr)
library(dbplyr)
library(tibble)
library(purrr)
library(locuscomparer)
library(ggplot2)
library(stringr)

# Read arguments
path_to_gwas <- args[1]
path_to_eqtl <- args[2]
output_name <- args[3]

# Upload gwas data
gwas_data <- read.table(path_to_gwas, sep = '\t', header = T)

# Check for complete cases
gwas_data <- gwas_data[complete.cases(gwas_data), ]

# Get rid of a SNP with beta == 0
gwas_data <- gwas_data[gwas_data$b != 0, ]

# Make varBetas
gwas_data$varbeta = gwas_data$se*gwas_data$se

# Find MA beta and MAF
gwas_data$b_MA <- gwas_data$b
gwas_data[gwas_data$freq > 0.50, ]$b_MA <- gwas_data[gwas_data$freq > 0.50, ]$b_MA *-1

gwas_data$MAF <- gwas_data$freq
gwas_data[gwas_data$freq > 0.50, ]$MAF <- 1 - gwas_data[gwas_data$freq > 0.50, ]$MAF

# Select required columns
gwas_data <- gwas_data %>% 
  mutate(s = N_cases / (N_cases + N_controls)) %>% 
  dplyr::select(rsID, A1, A2, freq, b_MA, MAF, varbeta, se, p, bp, s) %>%  
  dplyr::rename(SNP = rsID)

# Read QTL data
eqtl_data <- read.table(path_to_eqtl, sep = '\t', header = T, stringsAsFactors = TRUE)

# Check for complete cases
eqtl_data <- eqtl_data[complete.cases(eqtl_data), ]

# Make varBetas
eqtl_data$varbeta = eqtl_data$SE*eqtl_data$SE

# Find MA beta and MAF
eqtl_data$b_MA <- eqtl_data$b
eqtl_data[eqtl_data$Freq > 0.50, ]$b_MA <- eqtl_data[eqtl_data$Freq > 0.50, ]$b_MA *-1

eqtl_data$MAF <- eqtl_data$Freq
eqtl_data[eqtl_data$Freq > 0.50, ]$MAF <- 1 - eqtl_data[eqtl_data$Freq > 0.50, ]$MAF

# Select required columns
eqtl_data <- eqtl_data %>% 
  dplyr::select(SNP, A1, BP, A2, Freq, Gene, Probe, SE, p, MAF, varbeta, b_MA)

# join GWAS + eQTLs
df <- dplyr::inner_join(gwas_data, eqtl_data, by = "SNP", suffix = c(".gwas", ".eqtl"))

# loop over genes and apply coloc for each gene in a locus
genes <- unique(df$Gene)
my.res <- lapply(genes, function(x) {
  df_sub <- df[df$Gene == x, ]
  
  # run coloc.abf() for each gene 
  res <- coloc.abf(dataset1=list(beta=df_sub$b_MA.gwas, 
                                 varbeta=df_sub$varbeta.gwas, 
                                 s =df_sub$s, 
                                 type="cc"),
                   dataset2=list(beta=df_sub$b_MA.eqtl, 
                                 varbeta=df_sub$varbeta.eqtl, 
                                 N=nrow(df_sub), 
                                 MAF=df_sub$MAF.eqtl, 
                                 type="quant"))
  # if result is significant - save SNPs in credible set...
  if (res$summary[6] > 0.8){
    o <- order(res$results$SNP.PP.H4,decreasing=TRUE)
    cs <- cumsum(res$results$SNP.PP.H4[o])
    w <- which(cs > 0.95)[1]
    top_snps_index <- res$results[o,][1:w,]$snp
    top_snps_index_num <- as.numeric(str_replace(top_snps_index, "SNP.", ""))
    top_snps <- df_sub[top_snps_index_num, ]
    
    # save a table with SNPs
    top_snps_output_name <- paste0(output_name, "_top_snps_", x, ".tsv")
    write.table(top_snps, top_snps_output_name, sep = "\t", quote = F, row.names = F)
    
    # prepare datasets for locuscompare
    df_sub_eqtl <- df_sub %>% 
      dplyr::select(SNP, p.eqtl) %>%
      dplyr::mutate(SNP = as.character(SNP)) %>% 
      dplyr::rename(rsid = SNP, pval = p.eqtl)
    df_sub_gwas <- df_sub %>% 
      dplyr::select(SNP, p.gwas) %>%
      dplyr::mutate(SNP = as.character(SNP)) %>%
      dplyr::rename(rsid = SNP, pval = p.gwas)
    
    # plot significant colocalization
    plot <- locuscompare(in_fn1=df_sub_gwas,
                         in_fn2=df_sub_eqtl,
                         title1="GWAS", 
                         title2="eQTL", 
                         genome="hg19", 
                         population = "EUR")
    
    # save plot
    fig_output_name <- paste0(output_name, "_top_snps_", x, ".pdf")
    ggsave(filename = fig_output_name, 
           plot = plot, 
           device = "pdf",
           width = 7, 
           height = 4)
  }
  res <- as.data.frame(t(res$summary))
  return (data.frame(gene = x, res))
})

# combine results
results <- do.call("rbind", my.res)

# save all
all_res_output_name <- paste0(output_name, ".tsv")
write.table(results, all_res_output_name, sep = "\t", quote = F, row.names = F)
