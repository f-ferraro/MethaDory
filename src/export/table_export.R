#' Generate HTML template for MethaDory export
#'
#' @param ... All the data components for the HTML report
#' @return Complete HTML string
generate_html_template <- function(welcome_html, references_html, prediction_plot_base64,
                                  prediction_table_data,
                                  age_table_data, dim_plots_base64, signatures,
                                  include_dim_plots = TRUE,
                                  signature_version = NA, qc_missing_html = "",
                                  qc_cell_sex_base64 = "", qc_pca_base64 = "",
                                  qc_density_base64 = "") {

  # Create JSON data for tables. na = "null" keeps every key present even when
  # a value is NA (e.g. pNNET_average for SVM-only signatures); otherwise
  # jsonlite drops the key for that row and DataTables warns about a
  # "Requested unknown parameter" for columns inferred from other rows.
  prediction_table_json <- jsonlite::toJSON(prediction_table_data, dataframe = "rows", na = "null")
  age_table_json <- jsonlite::toJSON(age_table_data, dataframe = "rows", na = "null")

  # Signature names contain characters that are illegal in a CSS/jQuery ID
  # selector - most notably "." (e.g. VACTERL.combined_Haghshenas2024,
  # chr17p11.2del_vanderLaan2025, RNU2-2.n35_Leitao2025), which jQuery reads as
  # a class separator, so $("#dimplot-VACTERL.combined_...") matched nothing and
  # those plots never became visible even though their <img> was in the file.
  # Slugify the DOM id (and the <option> value that keys it) while keeping the
  # original signature name as the visible label.
  dim_plot_slugs <- gsub("[^A-Za-z0-9_-]", "-", names(dim_plots_base64))
  # Guard against two signatures slugifying to the same id (make.unique's ".1"
  # suffix would smuggle a dot back in, so re-slugify afterwards).
  dim_plot_slugs <- gsub("[^A-Za-z0-9_-]", "-", make.unique(dim_plot_slugs))
  names(dim_plot_slugs) <- names(dim_plots_base64)

  # Create dimension plots navigation
  dim_nav_items <- paste0(
    sapply(names(dim_plots_base64), function(s) {
      paste0('<option value="', dim_plot_slugs[[s]], '">', s, '</option>')
    }),
    collapse = ""
  )

  # Create dimension plots content
  dim_plots_content <- paste0(
    sapply(names(dim_plots_base64), function(s) {
      paste0(
        '<div id="dimplot-', dim_plot_slugs[[s]], '" class="dimension-plot" style="display: none; text-align: center;">',
        '<h4>', s, '</h4>',
        '<img src="', dim_plots_base64[[s]], '" style="max-width: 100%; height: auto;">',
        '</div>'
      )
    }),
    collapse = "\n"
  )

  html_template <- paste0('
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MethaDory Analysis Report</title>

    <!-- Bootstrap CSS (minified) -->
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <!-- DataTables CSS (minified) -->
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.11.5/css/dataTables.bootstrap5.min.css">

    <style>
        .main-header {
            background-color: #3c8dbc;
            color: white;
            padding: 15px;
            margin-bottom: 20px;
        }

        .content-wrapper {
            padding: 20px;
        }

        .dimension-plot {
            margin: 20px 0;
            padding: 20px;
            border: 1px solid #ddd;
            border-radius: 5px;
        }

        .plot-navigation {
            margin: 20px 0;
            padding: 15px;
            background-color: #f8f9fa;
            border-radius: 5px;
        }

        .nav-tabs .nav-link.active {
            background-color: #3c8dbc;
            border-color: #3c8dbc;
            color: white;
        }

        /* Color coding for combined mean(SVM, NNET) scores */
        .score-low { background-color: white !important; }
        .score-medium { background-color: #FAD302 !important; }
        .score-high { background-color: #9E1E05 !important; color: white !important; font-weight: bold !important; }

        /* Fix tab content overflow */
        .tab-content {
            position: relative;
            z-index: 1;
        }

        .tab-pane {
            overflow: auto;
            max-height: 85vh;
        }

        /* Sticky tabs */
        .nav-tabs {
            position: sticky;
            top: 0;
            z-index: 1000;
            background-color: white;
            border-bottom: 1px solid #dee2e6;
            margin-bottom: 0;
        }

        /* Prediction plot container */
        .prediction-plot img {
            max-width: 100%;
            height: auto;
            max-height: 70vh;
            object-fit: contain;
        }

        /* Dimension plot container: the per-signature figure is much taller
           than the viewport, so no max-height - the tab pane scrolls. */
        .dimension-plot img {
            max-width: 100%;
            height: auto;
        }
    </style>
</head>
<body>
    <div class="main-header">
        <div class="container-fluid">
            <h1>MethaDory Analysis Report</h1>
            <p class="mb-0">Generated on ', Sys.Date(), '</p>
            <p class="mb-0">Signature version: ', if (is.na(signature_version)) "unknown" else signature_version, '</p>
        </div>
    </div>

    <div class="container-fluid content-wrapper">
        <ul class="nav nav-tabs" id="mainTabs" role="tablist">
            <li class="nav-item" role="presentation">
                <button class="nav-link active" id="welcome-tab" data-bs-toggle="tab" data-bs-target="#welcome" type="button" role="tab">Welcome</button>
            </li>
            <li class="nav-item" role="presentation">
                <button class="nav-link" id="qc-tab" data-bs-toggle="tab" data-bs-target="#qc" type="button" role="tab">QC</button>
            </li>
            <li class="nav-item" role="presentation">
                <button class="nav-link" id="prediction-plot-tab" data-bs-toggle="tab" data-bs-target="#prediction-plot" type="button" role="tab">Prediction Results Plot</button>
            </li>
            <li class="nav-item" role="presentation">
                <button class="nav-link" id="prediction-table-tab" data-bs-toggle="tab" data-bs-target="#prediction-table" type="button" role="tab">Prediction Table</button>
            </li>
            ', if(include_dim_plots && length(dim_plots_base64) > 0) '<li class="nav-item" role="presentation">
                <button class="nav-link" id="dimension-tab" data-bs-toggle="tab" data-bs-target="#dimension" type="button" role="tab">Visualizations for interpretation</button>
            </li>' else '', '
            <li class="nav-item" role="presentation">
                <button class="nav-link" id="references-tab" data-bs-toggle="tab" data-bs-target="#references" type="button" role="tab">References</button>
            </li>
        </ul>

        <div class="tab-content" id="mainTabContent">
            <!-- Welcome Tab -->
            <div class="tab-pane fade show active" id="welcome" role="tabpanel">
                <div class="mt-3">
                    ', welcome_html, '
                </div>
            </div>

            <!-- QC Tab: missing values, methylation age, cell proportions | chromosomal
                 sex, PCA vs controls, beta-value distribution -->
            <div class="tab-pane fade" id="qc" role="tabpanel">
                <div class="mt-3">
                    <h4>Missing values before imputation</h4>
                    ', qc_missing_html, '

                    <h4 class="mt-4">Methylation age prediction</h4>
                    <table id="ageTable" class="table table-striped table-bordered" style="width:100%">
                    </table>

                    ', if(nchar(qc_cell_sex_base64) > 0) paste0('<h4 class="mt-4">Cell proportions and chromosomal sex</h4>
                    <div class="text-center">
                        <img src="', qc_cell_sex_base64, '" style="max-width: 100%; height: auto;">
                    </div>') else '', '

                    ', if(nchar(qc_pca_base64) > 0) paste0('<h4 class="mt-4">PCA against controls</h4>
                    <p>Principal component analysis of the sample together with the control samples, computed on the beta values before imputation using the top 1% most variable autosomal CpGs measured in the sample. Controls are shown in grey and the sample in red. A sample lying far from the controls may be of poor quality (or come from a different tissue or platform) and its predictions should be interpreted with caution.</p>
                    <div class="text-center">
                        <img src="', qc_pca_base64, '" style="max-width: 100%; height: auto;">
                    </div>') else '', '

                    ', if(nchar(qc_density_base64) > 0) paste0('<h4 class="mt-4">Beta-value distribution against controls</h4>
                    <p>Distribution of the beta values before imputation over all autosomal CpGs measured in the sample. Controls are shown in grey and the sample in red. A distribution that departs from the two peaks (near 0 and 1) of the controls may indicate poor quality, or data from a different platform that is not on the array scale.</p>
                    <div class="text-center">
                        <img src="', qc_density_base64, '" style="max-width: 100%; height: auto;">
                    </div>') else '', '
                </div>
            </div>

            <!-- Prediction Plot Tab -->
            <div class="tab-pane fade" id="prediction-plot" role="tabpanel">
                <div class="mt-3 text-center prediction-plot">
                    <img src="', prediction_plot_base64, '" alt="Prediction Results Plot">
                </div>
            </div>

            <!-- Prediction Table Tab (combined mean(SVM, NNET) metapredictor) -->
            <div class="tab-pane fade" id="prediction-table" role="tabpanel">
                <div class="mt-3">
                    <table id="predictionTable" class="table table-striped table-bordered" style="width:100%">
                    </table>
                </div>
            </div>
            ', if(include_dim_plots && length(dim_plots_base64) > 0) paste0('<!-- Visualizations for interpretation Tab -->
            <div class="tab-pane fade" id="dimension" role="tabpanel">
                <div class="mt-3">
                    <div class="plot-navigation">
                        <label for="dimPlotSelect"><strong>Jump to plot:</strong></label>
                        <select id="dimPlotSelect" class="form-select" style="width: 300px; display: inline-block; margin-left: 10px;">
                            <option value="">Select a signature...</option>
                            ', dim_nav_items, '
                        </select>
                    </div>
                    <div id="dimensionPlots">
                        ', dim_plots_content, '
                    </div>
                </div>
            </div>') else '', '

            <!-- References Tab -->
            <div class="tab-pane fade" id="references" role="tabpanel">
                <div class="mt-3">
                    ', references_html, '
                </div>
            </div>
        </div>
    </div>

    <!-- Bootstrap JS -->
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
    <!-- jQuery -->
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <!-- DataTables JS -->
    <script src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/1.11.5/js/dataTables.bootstrap5.min.js"></script>

    <script>
        // Data for tables
        const predictionData = ', prediction_table_json, ';
        const ageData = ', age_table_json, ';

        // Function to apply color coding to combined mean(SVM, NNET) scores
        function formatScoreCell(data, type, row) {
            if (data === null || data === undefined || data === "") {
                return type === "display" ? "" : data;
            }
            if (type === "display") {
                const value = parseFloat(data);
                if (isNaN(value)) return "";
                const rounded = value.toFixed(2);
                let className = "score-low";

                if (value >= 0.5) {
                    className = "score-high";
                } else if (value >= 0.25) {
                    className = "score-medium";
                }

                return `<span class="${className}">${rounded}</span>`;
            }
            return data;
        }

        function formatRounded(data, type, row) {
            if (data === null || data === undefined || data === "") {
                return type === "display" ? "" : data;
            }
            if (type === "display") {
                const value = parseFloat(data);
                return isNaN(value) ? "" : value.toFixed(2);
            }
            return data;
        }

        $(document).ready(function() {
            // Initialize prediction table. Columns are inferred from the data
            // so it works whether NNET inference ran or only SVM is present.
            const baseCols = [
                { data: "SampleID", title: "Sample ID" },
                { data: "SVM", title: "Signature" }
            ];
            const numericCols = [
                ["mean_case",     "pCombined",        true],
                ["pSVM_average",  "pSVM Average",     false],
                ["pSVM_sd",       "pSVM SD",          false],
                ["pNNET_average", "pNNET Average",    false],
                ["pNNET_sd",      "pNNET SD",         false],
                ["n_svm",         "n SVM",            false],
                ["n_nnet",        "n NNET",           false],
                ["pct_na_pre",    "% NA pre-imputation", false]
            ];
            const sample = predictionData[0] || {};
            const dynamicCols = numericCols
                .filter(c => c[0] in sample)
                .map(c => ({
                    data: c[0],
                    title: c[1],
                    render: c[2] ? formatScoreCell : formatRounded
                }));
            const predTable = $("#predictionTable").DataTable({
                data: predictionData,
                columns: baseCols.concat(dynamicCols),
                pageLength: 25,
                order: [[2, "desc"]]
            });

            // Initialize age table
            $("#ageTable").DataTable({
                data: ageData,
                columns: Object.keys(ageData[0] || {}).map(key => ({
                    data: key,
                    title: key.replace(/_/g, " ").replace(/\\b\\w/g, l => l.toUpperCase()) // Clean up titles
                })),
                pageLength: 25
            });


            // Dimension plots navigation
            $("#dimPlotSelect").on("change", function() {
                const selectedPlot = this.value;

                // Hide all plots
                $(".dimension-plot").hide();

                if (selectedPlot) {
                    // Attribute selector rather than "#id" so ids are matched
                    // literally even if a signature name ever yields a
                    // character jQuery would treat as a CSS combinator.
                    const $target = $("[id=\'dimplot-" + selectedPlot + "\']");

                    // Show selected plot
                    $target.show();

                    // Scroll to the plot
                    if ($target.length > 0) {
                        $target[0].scrollIntoView({
                            behavior: "smooth",
                            block: "start"
                        });
                    }
                }
            });

            // Show first dimension plot by default if any exist
            const firstPlot = $(".dimension-plot").first();
            if (firstPlot.length > 0) {
                firstPlot.show();
                const firstId = firstPlot.attr("id").replace("dimplot-", "");
                $("#dimPlotSelect").val(firstId);
            }
        });
    </script>
</body>
</html>')

  return(html_template)
}

#' Create Excel tables export
#'
#' One workbook for every sample in `data_list$results`: Analysis_Info,
#' Predictions_Wide, Predictions_Long, Signatures_Summary and, when the sample
#' statistics were computed, Cell_Proportions, Methylation_Age, Chromosomal_Sex.
#' Written by MethaDory_html (generate_methadory_html_report.R). openxlsx is
#' called by namespace so a caller does not have to attach it.
#'
#' @param data_list list with `results` and, optionally, `cell_props`,
#'   `age_table`, `chr_sex_table`.
#' @param output_path the .xlsx file to write.
create_excel_export <- function(data_list, output_path) {
  cat("Generating Excel tables...\n")

  # Prepare prediction results table
  prediction_results <- data_list$results

  # Display: rename mean_case -> pCombined and drop whisker/abs_diff cols.
  drop_cols <- c("whisker_low", "whisker_high", "abs_diff")
  prediction_results <- prediction_results[, setdiff(names(prediction_results), drop_cols), drop = FALSE]
  names(prediction_results)[names(prediction_results) == "mean_case"] <- "pCombined"

  # Create wide format for easier reading
  value_cols <- intersect(c("pSVM_average", "pSVM_sd",
                            "pNNET_average", "pNNET_sd",
                            "pCombined",
                            "n_svm", "n_nnet", "pct_na_pre"),
                          names(prediction_results))
  prediction_wide <- prediction_results %>%
    select(SampleID, SVM, all_of(value_cols)) %>%
    pivot_wider(
      id_cols = SampleID,
      names_from = SVM,
      values_from = all_of(value_cols),
      names_sep = "_"
    )

  # Sample stats are optional: a caller may leave them NULL. Only build the
  # corresponding tables when the data is present.
  have_cell_props <- !is.null(data_list$cell_props)
  have_age <- !is.null(data_list$age_table)
  have_chr_sex <- !is.null(data_list$chr_sex_table)

  # Prepare cell proportions table
  cell_props_wide <- if (have_cell_props) {
    data_list$cell_props %>%
      pivot_wider(
        id_cols = Proband,
        names_from = CellType,
        values_from = CellProp
      )
  } else NULL

  # Prepare age predictions
  age_table_clean <- if (have_age) {
    tbl <- data_list$age_table
    names(tbl) <- gsub("\\.", "_", names(tbl))
    tbl
  } else NULL

  # Prepare chromosomal sex predictions
  chr_sex_clean <- if (have_chr_sex) data_list$chr_sex_table else NULL

  # Signatures summary uses the metapredictor mean(SVM, NNET). It is called
  # pCombined by now: mean_case was renamed above, so looking for "mean_case"
  # here always fell through to pSVM_average.
  score_col <- if ("pCombined" %in% names(prediction_results)) "pCombined" else "pSVM_average"
  signatures_summary <- prediction_results %>%
    group_by(SVM) %>%
    summarise(
      Max_score = round(max(.data[[score_col]], na.rm = TRUE), 3),
      Min_score = round(min(.data[[score_col]], na.rm = TRUE), 3),
      Avg_score = round(mean(.data[[score_col]], na.rm = TRUE), 3),
      N_Samples = n(),
      High_Confidence = sum(.data[[score_col]] >= 0.5, na.rm = TRUE),
      Medium_Confidence = sum(.data[[score_col]] >= 0.25 & .data[[score_col]] < 0.5, na.rm = TRUE),
      Low_Confidence = sum(.data[[score_col]] >= 0.05 & .data[[score_col]] < 0.25, na.rm = TRUE),
      Very_Low = sum(.data[[score_col]] < 0.05, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Avg_score))

  # Run metadata: record the signature version so results are traceable to the
  # exact signature set used.
  analysis_info <- data.frame(
    Field = c("signature_version", "signature_file", "generated"),
    Value = c(signature_version(), SIGNATURE_FILE, as.character(Sys.time())),
    stringsAsFactors = FALSE
  )

  # Create workbook with separate sheets
  wb <- openxlsx::createWorkbook()

  # Create separate sheets
  openxlsx::addWorksheet(wb, "Analysis_Info")
  openxlsx::addWorksheet(wb, "Predictions_Wide")
  openxlsx::addWorksheet(wb, "Predictions_Long")
  openxlsx::addWorksheet(wb, "Signatures_Summary")

  # Write data to sheets
  openxlsx::writeData(wb, "Analysis_Info", analysis_info, rowNames = FALSE)
  openxlsx::writeData(wb, "Predictions_Wide", prediction_wide, rowNames = FALSE)
  openxlsx::writeData(wb, "Predictions_Long", prediction_results, rowNames = FALSE)
  openxlsx::writeData(wb, "Signatures_Summary", signatures_summary, rowNames = FALSE)

  # Optional sample-stats sheets: only added when the stats were computed.
  if (have_cell_props) {
    openxlsx::addWorksheet(wb, "Cell_Proportions")
    openxlsx::writeData(wb, "Cell_Proportions", cell_props_wide, rowNames = FALSE)
  }
  if (have_age) {
    openxlsx::addWorksheet(wb, "Methylation_Age")
    openxlsx::writeData(wb, "Methylation_Age", age_table_clean, rowNames = FALSE)
  }
  if (have_chr_sex) {
    openxlsx::addWorksheet(wb, "Chromosomal_Sex")
    openxlsx::writeData(wb, "Chromosomal_Sex", chr_sex_clean, rowNames = FALSE)
  }

  # Add conditional formatting for pSVM values in summary
  openxlsx::conditionalFormatting(wb, "Signatures_Summary", cols = 2:4, rows = 2:(nrow(signatures_summary) + 1),
                       rule = ">=0.5", style = openxlsx::createStyle(bgFill = "#9E1E05", fontColour = "white"))
  openxlsx::conditionalFormatting(wb, "Signatures_Summary", cols = 2:4, rows = 2:(nrow(signatures_summary) + 1),
                       rule = ">=0.25", style = openxlsx::createStyle(bgFill = "#FAD302"))

  # Save workbook
  openxlsx::saveWorkbook(wb, output_path, overwrite = TRUE)
  cat(paste("Excel file saved:", output_path, "\n"))
}

#' Export the imputed methylation data: user samples + controls + real cases.
#'
#' One IlmnID x samples table of the beta values the predictions were made on
#' (the user samples after imputation), next to the array reference they are
#' compared with: the background controls and the real affected individuals.
#' Probes are the union across the three; a sample has NA where it has no value.
#'
#' @param imputed_data  data.frame, IlmnID + one column per user sample.
#' @param insilico_beta background controls (IlmnID column or CpG row names), or NULL.
#' @param output_path   the .tsv file to write.
export_imputed_data <- function(imputed_data, insilico_beta, output_path) {
  cat("\nExporting imputed methylation data (user samples + controls + cases)...\n")
  combined <- imputed_data

  if (!is.null(insilico_beta)) {
    controls <- insilico_beta
    if (!"IlmnID" %in% names(controls)) controls$IlmnID <- rownames(controls)
    combined <- merge(combined, controls, by = "IlmnID", all = TRUE)
    cat("  Added", ncol(controls) - 1, "control samples\n")
  }

  # The real cases come from the full affected-individuals matrix rather than
  # from the per-signature subsets used for plotting, which hold only each
  # signature's own probes.
  tryCatch({
    cases_beta <- readRDS("../../data/affectedindividuals/affectedindividuals_methadory.beta.rds")
    cases_meta <- readRDS("../../data/affectedindividuals/affectedindividuals_methadory.meta.rds")
    case_ids <- cases_meta$geo_accession[cases_meta$RealLabel != "control"]
    if (!"IlmnID" %in% names(cases_beta)) cases_beta$IlmnID <- rownames(cases_beta)
    cases_only <- cases_beta[, c("IlmnID", intersect(names(cases_beta), case_ids)),
                             drop = FALSE]
    combined <- merge(combined, cases_only, by = "IlmnID", all = TRUE)
    cat("  Added", ncol(cases_only) - 1, "real case samples\n")
  }, error = function(e) {
    cat("  Warning: Could not load real cases data:", e$message, "\n")
  })

  write.table(combined, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste(" Imputed data saved:", output_path, "\n"))
  cat(paste("  Total samples:", ncol(combined) - 1, "\n"))
  cat(paste("  Total CpGs:", nrow(combined), "\n"))
  invisible(output_path)
}
