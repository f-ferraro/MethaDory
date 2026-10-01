#!/usr/bin/env Rscript
# Merge per-sample pseudoepic BEDs (produced by
# MethaDory-extract-methylation-pacbio.sh) into a single beta-value TSV
# keyed by Illumina probe ID, ready for MethaDory input.
#
# Expected per-sample columns (tab-separated, no header):
#   1 chrom  2 begin  3 end  4 mod_score  5 type  6 cov  7 IlmnID

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: MethaDory_PB_Input_Preparation.R <pseudoepic_directory> <output_file> [min_cov]
    Arguments:
       pseudoepic_directory   Directory of *.pseudoepic.cpgID.bed files
       output_file            Path to the output TSV
       min_cov                Minimum pb-CpG-tools coverage to keep a probe (default 4)",
       call. = FALSE)
}
input_dir   <- args[1]
output_file <- args[2]
min_cov     <- if (length(args) >= 3) as.integer(args[3]) else 4

files <- list.files(input_dir, full.names = TRUE,
                    pattern = "\\.pseudoepic\\.cpgID\\.bed$")
if (length(files) == 0) {
  stop("No *.pseudoepic.cpgID.bed files in: ", input_dir)
}

res <- list()
for (f in files) {
  message("Processing ", f)
  df <- read.table(f, header = FALSE, sep = "\t",
                   col.names = c("chr", "begin", "end", "mod_score",
                                 "type", "cov", "IlmnID"),
                   stringsAsFactors = FALSE)
  # Keep only combined-strand "Total" rows that pass coverage.
  df <- df[df$type == "Total" & df$cov >= min_cov, c("IlmnID", "mod_score")]
  sample_id <- basename(gsub("\\.pseudoepic\\.cpgID\\.bed$", "", f))
  names(df) <- c("IlmnID", sample_id)
  res[[f]] <- df
}

gc()

message("Merging ", length(res), " files")
out <- res %>%
  purrr::reduce(full_join, by = "IlmnID") %>%
  dplyr::relocate("IlmnID")

# pb-CpG-tools mod_score is 0-100 (percent methylated). Convert to beta.
out[, 2:ncol(out)] <- out[, 2:ncol(out)] / 100

write.table(out, file = output_file, row.names = FALSE,
            sep = "\t", quote = FALSE)
message("Wrote ", output_file)
