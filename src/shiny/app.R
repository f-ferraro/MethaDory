options(timeout = 10000)
Sys.setenv(R_DEFAULT_INTERNET_TIMEOUT = "10000")

# Add packages that are broken in pixi
# packages <- c("FDb.InfiniumMethylation.hg19", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
#               "ChAMPdata",  "GenomeInfoDb")
# 
# for (pkg in packages) {
#   if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
#     BiocManager::install(pkg)
#   }
# }

# Load required libraries
suppressPackageStartupMessages({
  library(BiocParallel)
  library(caret)
  library(circlize)
  library(ComplexHeatmap)
  library(data.table)
  library(DT)
  library(butcher)
  library(EpiDISH)
  library(wateRmelon)
  library(fs)
  library(ggplotify)
  library(methyLImp2)
  library(patchwork)
  library(PCAtools)
  library(plotly)
  library(shiny)
  library(shinybusy)
  library(shinydashboard)
  library(shinyFiles)
  # tidyverse components actually used (avoid pulling the meta-package).
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(tibble)
  library(scales)
  library(ggrepel)
  library(htmlwidgets)
  library(base64enc)
  library(jsonlite)
  library(markdown)
  library(reticulate)
})

# Set Shiny options
options(shiny.maxRequestSize = 500 * 1024^2)

# Source UI and Server components
source("ui.R")
source("server.R")

# Launch the Shiny application
shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
