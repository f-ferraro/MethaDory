Welcome to MethaDory!
=======================

<img style="float: right;" src="methadory.png" width="200">

##### Prerequisites

To run Methadory

  * Place the `data` folder in the same directory where the MethaDory code resides. This folder must contain:
    * `affectedindividuals`, folder containging beta data and meta data of example affected individuals
    * `imputation`, folder containg data used for imputing missing values in the test samples
    * `manifest.qc_filtered.rds`, Illumina Epic V2 array after filtering for QC as described in XXX. 
    * `merged_signatures_90DMRs.tsv`, containing information about the sites used for building the DNAm signatures and classifiers
 
  * Prepare test sample beta table. MethaDory assumes you're providing a table with first column called "IlmnID" (containg the probe names from the array), and then a column for each sample. 
  * Select classifier models and upload the test data
  * Wait for processing to be completed
  * Select a signature or multiple to plot results
  * Export results as a table and / or as pdf


Source code can be found in [GitHub](https://github.com/f-ferraro/MethaDory).