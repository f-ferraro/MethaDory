#' Create cell deconvolution table
#'
#' @param samples_of_interest_beta Input to calculate cell proportions
#' @return table for cell proportions
#'
create_cell_deconv_table <- function(samples_of_interest_beta) {

  rownames(samples_of_interest_beta) = samples_of_interest_beta$IlmnID
  samples_of_interest_beta$IlmnID = NULL
  BloodFrac.m <- epidish(beta.m = samples_of_interest_beta,
                         ref.m = centDHSbloodDMC.m,
                         method = "RPC")$estF

  bf = as.data.frame(t(BloodFrac.m))
  bf$CellType = rownames(bf)

  bf = pivot_longer(bf, -"CellType",
                    names_to="Proband",
                    values_to="CellProp")

  # EpiDISH returns "B" for B cells; keep that name so the proband points
  # land on the same x-axis category as the training background, which is
  # also produced from EpiDISH (no rename).
  return(bf)
}

#' Predict methylation age
#'
#' @param test_beta Input to calculate age
#' @return table with age predictions
predict_age <- function(test_beta) {
  print("Predicting age")
  rownames(test_beta) = test_beta$IlmnID
  test_beta$IlmnID = NULL

  ages = agep(test_beta, method='all')

  ages = tibble::rownames_to_column(ages, "Proband")
  ages = as.data.frame(ages)
  
  ages = ages[, names(ages) %like% "Proband|skin"]
  names(ages) = c("Proband", "Predicted age", "Missing CpGs for prediction")

  # Show the missing CpGs out of the clock's total (skin & blood clock, 391
  # CpGs), otherwise the count alone says little.
  n_clock_cpgs <- tryCatch({
    data("age_coefficients", package = "wateRmelon", envir = environment())
    sum(names(ageCoefs$SkinBlood) != "(Intercept)")
  }, error = function(e) 391)
  ages[["Missing CpGs for prediction"]] = paste(ages[["Missing CpGs for prediction"]], "/", n_clock_cpgs)
  
  return(ages)
}

#' Predict chromosomal sex table
#'
#' @param test_beta Input to calculate sex
#' @return table with sex predictions
predict_chr_sex_table <- function(test_beta) {
  print("Predicting chromosomal sex")
  rownames(test_beta) = test_beta$IlmnID
  test_beta$IlmnID = NULL

  # estimateSex misparses a tibble/grouped input (matrix coercion silently
  # drops sample names); force a plain data.frame so probe rownames + sample
  # colnames survive.
  x = estimateSex(as.data.frame(test_beta), do_plot=F)
  x = tibble::rownames_to_column(x, "Proband")
  return(x)
}