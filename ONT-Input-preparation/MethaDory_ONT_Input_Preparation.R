#!/usr/bin/env Rscript
library(tidyverse)

suppressWarnings(suppressMessages(library(tidyverse)))

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: MethaDory_ONT_Input_Preparation.R <pseudoepic_directory> <output_file>
    Arguments:
       pseudoepic_directory   Path to the pseudoepic directory output of MethaDory-modkit.sh
       output_file            Path to the output file", call. = FALSE)
}
input_dir <- args[1]
output_file <- args[2]

res = list()

for (i in list.files(input_dir, full.names = TRUE)) {
  # Read data
  print(paste0("Processing ", i, "..."))
  
  df <- read.table(i, header = FALSE)
  
  df$TotalDepth = rowSums(df[, paste0("V", 12:18)])
  df$UsableFraction = df$V5 / df$TotalDepth
  
  df = as.data.frame(df[, c("V1", "V2", "V5", "V6", "V11", "TotalDepth", "UsableFraction", "V19")])
  
  # Filter for sites with more than 10 usable reads
  df = df[df$TotalDepth >= 13 & df$UsableFraction > 0.75, ]

  df = df[, names(df) %in% c("V19", "V11")]
  names(df) = c(basename(gsub(".pseudoepic.cpgID.bed.gz", "", i)), "IlmnID")
  
  res[[i]] = df
}

gc()

print(paste("Merging", length(list.files(input_dir)), " files"))

res2 = res %>%
  purrr::reduce(full_join, by = 'IlmnID')

res2 = res2 %>%
  relocate("IlmnID")

head(res2)
rm(res)
gc()

res2[, 2:ncol(res2)] = res2[, 2:ncol(res2)] / 100

write.table(res2, file = output_file, row.names = FALSE, sep = "\t", quote = FALSE)
