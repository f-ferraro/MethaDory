# Fraction of the most variable CpGs kept for the sample QC PCA.
QC_PCA_TOP_FRACTION <- 0.01
QC_PCA_MANIFEST_FILE <- "../../data/support_files/manifest.qc_filtered.rds"
# Number of equal-width beta bins (0-1) for the QC beta-value density curves.
QC_DENSITY_BINS <- 100

#' Bin index (1..QC_DENSITY_BINS) of beta values; NA stays NA. Betas are rounded
#' to 3 decimals first so values sitting on a bin edge are not split by
#' floating point noise, and clamped so out-of-range values land in the end
#' bins.
qc_density_bin <- function(beta) {
  pmax(pmin(floor(round(beta * 1000) * QC_DENSITY_BINS / 1000) + 1, QC_DENSITY_BINS), 1)
}

#' Sample QC: PCA of each proband together with the control background
#'
#' Runs on the proband betas BEFORE imputation, against the genome-wide control
#' samples shipped for imputation (controls only). Each proband gets its own
#' PCA (proband + all controls), so a report never depends on the other samples
#' of the input file. Only CpGs measured in the proband are used (no imputed
#' values), sex chromosomes are dropped so the PCs are not driven by sex, and
#' the top `top_fraction` most variable of the remaining CpGs go into the PCA.
#'
#' The same scan also yields the beta-value density curves (histogram over
#' QC_DENSITY_BINS bins) shown under the PCA: each control over all autosomal
#' CpGs shared with the input, each proband over the ones it has measured.
#'
#' @param test_data Pre-imputation test data frame (IlmnID + one column per sample)
#' @param imputation_background Control beta data frame (IlmnID + one column per control)
#' @param test_sample_ids Sample IDs to process
#' @param top_fraction Fraction of most variable CpGs to keep (default 1%)
#' @param n_pcs Number of principal components to return
#' @param chunk_size Rows per chunk when scanning the control matrix
#' @return Named list (one entry per sample) with `scores` (SampleID, Group,
#'   PC1..PCn), `var_explained` (%), `n_cpgs_available`, `n_cpgs_used` and
#'   `density` (SampleID, Group, beta, density).
#'   Samples with too few measured CpGs are skipped.
compute_qc_pca <- function(test_data, imputation_background, test_sample_ids,
                           top_fraction = QC_PCA_TOP_FRACTION, n_pcs = 4,
                           chunk_size = 1e5) {
  control_ids <- setdiff(names(imputation_background), "IlmnID")

  cpgs <- intersect(test_data$IlmnID, imputation_background$IlmnID)
  if (file.exists(QC_PCA_MANIFEST_FILE)) {
    manifest <- readRDS(QC_PCA_MANIFEST_FILE)
    sex_cpgs <- manifest$ProbeID[manifest$chr %in% c("chrX", "chrY")]
    cpgs <- setdiff(cpgs, sex_cpgs)
  }
  bg_idx <- match(cpgs, imputation_background$IlmnID)
  test_idx <- match(cpgs, test_data$IlmnID)

  # Per-CpG sum and sum of squares over the controls, scanned in chunks so the
  # ~1M x 100 background is never held as a second full matrix. With these the
  # variance of (controls + proband) is obtained per proband without touching
  # the control matrix again.
  ctrl_sum <- numeric(length(cpgs))
  ctrl_sumsq <- numeric(length(cpgs))
  ctrl_bin_counts <- numeric(QC_DENSITY_BINS * length(control_ids))
  bin_offset <- rep((seq_along(control_ids) - 1) * QC_DENSITY_BINS, each = chunk_size)
  for (start in seq(1, length(cpgs), by = chunk_size)) {
    idx <- start:min(start + chunk_size - 1, length(cpgs))
    chunk <- as.matrix(imputation_background[bg_idx[idx], control_ids, drop = FALSE])
    ctrl_sum[idx] <- rowSums(chunk)
    ctrl_sumsq[idx] <- rowSums(chunk^2)
    # chunk is column-major, so bin_offset shifts each control to its own bins
    offsets <- if (length(idx) == chunk_size) bin_offset else
      rep((seq_along(control_ids) - 1) * QC_DENSITY_BINS, each = length(idx))
    ctrl_bin_counts <- ctrl_bin_counts +
      tabulate(qc_density_bin(chunk) + offsets, QC_DENSITY_BINS * length(control_ids))
  }

  bin_mids <- (seq_len(QC_DENSITY_BINS) - 0.5) / QC_DENSITY_BINS
  ctrl_bin_counts <- matrix(ctrl_bin_counts, nrow = QC_DENSITY_BINS)
  ctrl_density <- data.frame(
    SampleID = rep(control_ids, each = QC_DENSITY_BINS),
    Group = "Control",
    beta = rep(bin_mids, length(control_ids)),
    density = as.vector(sweep(ctrl_bin_counts, 2, colSums(ctrl_bin_counts), "/")) * QC_DENSITY_BINS
  )

  n <- length(control_ids) + 1
  qc_pca <- list()
  for (id in test_sample_ids) {
    x <- test_data[[id]][test_idx]
    # is.finite also drops CpGs where a control is NA (NA control sum)
    measured <- which(is.finite(x) & is.finite(ctrl_sum))
    n_top <- ceiling(top_fraction * length(measured))
    if (n_top < n_pcs) {
      cat("QC PCA skipped for", id, "- only", length(measured), "measured CpGs\n")
      next
    }

    total <- ctrl_sum[measured] + x[measured]
    variance <- (ctrl_sumsq[measured] + x[measured]^2 - total^2 / n) / (n - 1)
    top <- measured[order(variance, decreasing = TRUE)[seq_len(n_top)]]

    beta <- cbind(as.matrix(imputation_background[bg_idx[top], control_ids, drop = FALSE]),
                  x[top])
    pca <- prcomp(t(beta), center = TRUE, scale. = FALSE, rank. = n_pcs)

    scores <- data.frame(SampleID = c(control_ids, id),
                         Group = c(rep("Control", length(control_ids)), "Proband"),
                         pca$x, row.names = NULL)
    qc_pca[[id]] <- list(
      scores = scores,
      var_explained = (100 * pca$sdev^2 / sum(pca$sdev^2))[seq_len(ncol(pca$x))],
      n_cpgs_available = length(measured),
      n_cpgs_used = n_top,
      density = rbind(ctrl_density, data.frame(
        SampleID = id, Group = "Proband", beta = bin_mids,
        density = tabulate(qc_density_bin(x[measured]), QC_DENSITY_BINS) /
          length(measured) * QC_DENSITY_BINS))
    )
  }

  qc_pca
}

# Thresholds (% of CpGs missing before imputation) for the QC status.
QC_MISSING_WARN_PCT <- 5
QC_MISSING_FAIL_PCT <- 15

#' Sample QC: share of missing values before imputation
#'
#' Computed over the CpGs MethaDory actually uses (the imputation panel: SVM
#' predictors, NNET probes and all signature CpGs), so it is the share of the
#' model input that had to be imputed. CpGs absent from the input file count
#' as missing.
#'
#' @param test_data_pre Pre-imputation beta matrix (probes x samples), i.e.
#'   `prepare_imputation_data()$test_data`
#' @param sample_ids Test sample IDs
#' @return data.frame with SampleID, n_cpgs, n_missing, pct_missing and status
#'   (PASS < 5%, WARNING 5-15%, FAIL > 15%)
compute_pre_imputation_missing <- function(test_data_pre, sample_ids) {
  sample_ids <- intersect(sample_ids, names(test_data_pre))
  n_missing <- vapply(sample_ids, function(id) sum(is.na(test_data_pre[[id]])), numeric(1))
  pct_missing <- 100 * n_missing / nrow(test_data_pre)
  data.frame(
    SampleID = sample_ids,
    n_cpgs = nrow(test_data_pre),
    n_missing = n_missing,
    pct_missing = round(pct_missing, 1),
    status = ifelse(pct_missing > QC_MISSING_FAIL_PCT, "FAIL",
                    ifelse(pct_missing >= QC_MISSING_WARN_PCT, "WARNING", "PASS")),
    row.names = NULL
  )
}
