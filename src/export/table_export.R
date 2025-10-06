#' Generate HTML template for MethaDory export
#'
#' @param ... All the data components for the HTML report
#' @return Complete HTML string
generate_html_template <- function(welcome_html, references_html, prediction_plot_base64,
                                  prediction_table_data, cell_prop_base64, chr_sex_base64,
                                  age_table_data, dim_plots_base64, signatures,
                                  include_cell_plot = TRUE, include_chr_sex_plot = TRUE, include_dim_plots = TRUE) {

  # Create JSON data for tables
  prediction_table_json <- jsonlite::toJSON(prediction_table_data, dataframe = "rows")
  age_table_json <- jsonlite::toJSON(age_table_data, dataframe = "rows")

  # Create dimension plots navigation
  dim_nav_items <- paste0(
    sapply(names(dim_plots_base64), function(s) {
      paste0('<option value="', s, '">', s, '</option>')
    }),
    collapse = ""
  )

  # Create dimension plots content
  dim_plots_content <- paste0(
    sapply(names(dim_plots_base64), function(s) {
      paste0(
        '<div id="dimplot-', s, '" class="dimension-plot" style="display: none; text-align: center;">',
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

        .table-psvm-filter {
            margin: 15px 0;
        }

        /* Color coding for pSVM values */
        .psvm-low { background-color: white !important; }
        .psvm-medium { background-color: #FAD302 !important; }
        .psvm-high { background-color: #9E1E05 !important; color: white !important; font-weight: bold !important; }

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

        /* Dimension plot container */
        .dimension-plot img {
            max-width: 100%;
            height: auto;
            max-height: 80vh;
            object-fit: contain;
        }
    </style>
</head>
<body>
    <div class="main-header">
        <div class="container-fluid">
            <h1>MethaDory Analysis Report</h1>
            <p class="mb-0">Generated on ', Sys.Date(), '</p>
        </div>
    </div>

    <div class="container-fluid content-wrapper">
        <ul class="nav nav-tabs" id="mainTabs" role="tablist">
            <li class="nav-item" role="presentation">
                <button class="nav-link active" id="welcome-tab" data-bs-toggle="tab" data-bs-target="#welcome" type="button" role="tab">Welcome</button>
            </li>
            <li class="nav-item" role="presentation">
                <button class="nav-link" id="prediction-plot-tab" data-bs-toggle="tab" data-bs-target="#prediction-plot" type="button" role="tab">Prediction Results Plot</button>
            </li>
            <li class="nav-item" role="presentation">
                <button class="nav-link" id="prediction-table-tab" data-bs-toggle="tab" data-bs-target="#prediction-table" type="button" role="tab">SVM Prediction Table</button>
            </li>
            ', if(include_cell_plot) '<li class="nav-item" role="presentation">
                <button class="nav-link" id="cell-prop-tab" data-bs-toggle="tab" data-bs-target="#cell-prop" type="button" role="tab">Cell Proportion Deconvolution</button>
            </li>' else '', '
            ', if(include_chr_sex_plot) '<li class="nav-item" role="presentation">
                <button class="nav-link" id="chr-sex-tab" data-bs-toggle="tab" data-bs-target="#chr-sex" type="button" role="tab">Chromosomal Sex Prediction</button>
            </li>' else '', '
            <li class="nav-item" role="presentation">
                <button class="nav-link" id="age-tab" data-bs-toggle="tab" data-bs-target="#age" type="button" role="tab">Methylation Age Prediction</button>
            </li>
            ', if(include_dim_plots && length(dim_plots_base64) > 0) '<li class="nav-item" role="presentation">
                <button class="nav-link" id="dimension-tab" data-bs-toggle="tab" data-bs-target="#dimension" type="button" role="tab">Dimension Reduction Plots</button>
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

            <!-- Prediction Plot Tab -->
            <div class="tab-pane fade" id="prediction-plot" role="tabpanel">
                <div class="mt-3 text-center prediction-plot">
                    <img src="', prediction_plot_base64, '" alt="Prediction Results Plot">
                </div>
            </div>

            <!-- SVM Prediction Table Tab -->
            <div class="tab-pane fade" id="prediction-table" role="tabpanel">
                <div class="mt-3">
                    <table id="predictionTable" class="table table-striped table-bordered" style="width:100%">
                    </table>
                </div>
            </div>

            ', if(include_cell_plot) paste0('<!-- Cell Proportion Tab -->
            <div class="tab-pane fade" id="cell-prop" role="tabpanel">
                <div class="mt-3 text-center">
                    <img src="', cell_prop_base64, '" style="max-width: 100%; height: auto;">
                </div>
            </div>') else '', '

            ', if(include_chr_sex_plot) paste0('<!-- Chromosomal Sex Tab -->
            <div class="tab-pane fade" id="chr-sex" role="tabpanel">
                <div class="mt-3 text-center">
                    <img src="', chr_sex_base64, '" style="max-width: 100%; height: auto;">
                </div>
            </div>') else '', '

            <!-- Age Prediction Tab -->
            <div class="tab-pane fade" id="age" role="tabpanel">
                <div class="mt-3">
                    <table id="ageTable" class="table table-striped table-bordered" style="width:100%">
                    </table>
                </div>
            </div>
            ', if(include_dim_plots && length(dim_plots_base64) > 0) paste0('<!-- Dimension Reduction Plots Tab -->
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

        // Function to apply color coding to pSVM values
        function formatPSVMCell(data, type, row) {
            if (type === "display") {
                const value = parseFloat(data);
                const rounded = value.toFixed(2);
                let className = "psvm-low";

                if (value >= 0.5) {
                    className = "psvm-high";
                } else if (value >= 0.25) {
                    className = "psvm-medium";
                }

                return `<span class="${className}">${rounded}</span>`;
            }
            return data;
        }

        function formatRounded(data, type, row) {
            if (type === "display") {
                return parseFloat(data).toFixed(2);
            }
            return data;
        }

        $(document).ready(function() {
            // Initialize prediction table
            const predTable = $("#predictionTable").DataTable({
                data: predictionData,
                columns: [
                    { data: "SampleID", title: "Sample ID" },
                    { data: "SVM", title: "SVM" },
                    {
                        data: "pSVM_average",
                        title: "pSVM Average",
                        render: formatPSVMCell
                    },
                    {
                        data: "pSVM_sd",
                        title: "pSVM SD",
                        render: formatRounded
                    }
                ],
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
                    // Show selected plot
                    $("#dimplot-" + selectedPlot).show();

                    // Scroll to the plot
                    $("#dimplot-" + selectedPlot)[0].scrollIntoView({
                        behavior: "smooth",
                        block: "start"
                    });
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