#!/usr/bin/env Rscript

library(stringr)
library(Biostrings)
library(tidyr)
library(ampir)

filename <- "/path/to/project_folder/raw_data/FASTA/input-TFO-library-rep02_S2_L001_R1_001_trimmed.fasta"
data <- readDNAStringSet(filename)

pattern <- ("CCTCCT(.....................)")

data_CCT <- as.data.frame(data)
data_CCT$name <- data@ranges@NAMES
data_CCT$seq <- str_extract(data,pattern)
data_CCT$x <- NULL
rownames(data_CCT) <- NULL
nrow(data_CCT)

data_CCT <- data_CCT %>% drop_na()
nrow(data_CCT)

output_file <- paste0(tools::file_path_sans_ext(filename),"_filt.fasta")
df_to_faa(data_CCT, output_file)
