Welcome to MethaDory
=======================

<img style="float: right;" src="methadory.png" width="200">

##### Prerequisites

To run Methadory

  * Place the `data` folder in the same directory where the MethaDory code resides. This folder must contain:
    * `affectedindividuals`, the folder containing beta data and meta data of example affected individuals (used for plotting).
    * `imputation`, the folder containing the data used for imputing missing values in the test samples.
    * `manifest.qc_filtered.rds`, of Illumina Epic V2 array manifest after filtering for QC as described in XXX. 
    * `merged_signatures_90DMRs.tsv`, text file containing information about the sites used for building the DNAm signatures and classifiers.
    * `background_training.cellprops.rds`, containing deconvoluted cell proportions from the samples used for the model training.
 
  * Prepare input data. MethaDory assumes you're providing a table with first column called `IlmnID` (containing the probe names from the array), and then a column for each sample. 
  * Select `svm-models` folder and upload your test data
  * Wait for processing to be completed
  * Select a signature or multiple to plot results
  * Export results as a table and / or as pdf


Source code can be found in [GitHub](https://github.com/f-ferraro/MethaDory).