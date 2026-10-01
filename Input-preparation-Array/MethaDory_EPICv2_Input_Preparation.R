#!/usr/bin/env Rscript

.libPaths(c("~/R-champ-dev", .libPaths()))

suppressPackageStartupMessages({
  library(minfi)
  library(ChAMP)
  library(vroom)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: MethaDory_EPICv2_Input_Preparation.R <idat_directory> <output_file>
    Arguments:
       idat_directory   Directory of paired *_Grn.idat(.gz) and *_Red.idat(.gz) files
       output_file      Path to the output TSV", call. = FALSE)
}
input_dir   <- args[1]
output_file <- args[2]

idats <- list.files(input_dir, pattern = "_(Grn|Red)\\.idat(\\.gz)?$",
                    full.names = TRUE, recursive = TRUE)
if (length(idats) == 0) {
  stop("No idat files in: ", input_dir)
}

Basename <- unique(gsub("_(Grn|Red)\\.idat(\\.gz)?$", "", idats))
targets <- data.frame(Sample_Name = basename(Basename),
                      Basename = Basename)

rgSet <- read.metharray.exp(targets = targets, verbose = TRUE,
                            extended = FALSE, force = TRUE)

beta <- getBeta(mapToGenome(preprocessNoob(rgSet)))

detP <- detectionP(rgSet)
detP <- detP[rownames(beta), colnames(beta), drop = FALSE]

filtered <- champ.filter(beta = beta,
                         detP = detP,
                         M = NULL, pd = NULL, Meth = NULL, intensity = NULL,
                         UnMeth = NULL, beadcount = NULL,
                         filterBeads = FALSE,
                         arraytype = "EPICv2",
                         SampleCutoff = 0.1,
                         ProbeCutoff = 0.5,
                         detPcut = 0.01,
                         filterDetP = TRUE,
                         autoimpute = TRUE,
                         filterXY = FALSE,
                         filterNoCG = TRUE,
                         filterSNPs = TRUE,
                         population = NULL,
                         fixOutlier = TRUE,
                         filterMultiHit = TRUE)$beta

out <- data.frame(IlmnID = sub("_.*$", "", rownames(filtered)),
                  filtered, check.names = FALSE)
out <- out[!duplicated(out$IlmnID), ]

vroom_write(out, output_file)
message("Wrote ", output_file)
