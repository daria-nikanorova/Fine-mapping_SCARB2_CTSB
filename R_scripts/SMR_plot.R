#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Add libraries
library(stringr)
library(dplyr)
library(ggplot2)

# Set path to working directory with .smr files from SMR analysis with different tissues
input_path <- args[1]
output_path <- args[2]

# Create a list of files from SMR analysis
SMR_files <- list.files(path = input_path, pattern="*.smr\\b", ignore.case = T, full.names = T)
print(SMR_files)
# Read them all and create a list
list_of_SMR_files <- lapply(SMR_files, function(x) read.table(x, header = TRUE, sep = "\t", ))

# Add tissue column from a name of a file
for (i in 1:length(SMR_files)){
  list_of_SMR_files[[i]] <- list_of_SMR_files[[i]] %>% 
    dplyr::mutate(tissue = basename(SMR_files[i]), .before = 1,
                  log_pSMR = -log10(p_SMR),
                  log_pHEIDI = -log10(p_HEIDI))
}

# Create a histogram for each tissue (for each .smr file you will get a histogram with the same name)
for (i in 1:length(SMR_files)){
  my_plot <- ggplot(list_of_SMR_files[[i]], aes(x = reorder(Gene, desc(log_pSMR)), y = log_pSMR, fill = Gene))+
    geom_col() +
    theme_classic() +
    geom_text(aes(label = round(log_pSMR, 2)), size=4) +
    scale_x_discrete(position = "bottom") +
    theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1, size=10, color = "black"),
          axis.text.y = element_text(size=10, color = "black"),
          axis.title.x = element_blank(),
          axis.ticks.y=element_line(color = "black"),
          strip.text.x = element_text(size=10, color = "black"),
          legend.position = "none") +
    geom_hline(yintercept = 4, linetype = 2, color = "#d73027", size = 0.7)
  
  # Save plot to the output directory
  plot_name <- str_replace(basename(SMR_files[i]), ".smr\\b", ".pdf")
  plot_path <- str_c(output_path, '/', plot_name)
  ggsave(filename = plot_path, plot = my_plot, device = "pdf", width = 6, height = 3)
}

# If you want to visualize the results of HEIDI test as well, just change the variable for y axis to log_pHEIDI