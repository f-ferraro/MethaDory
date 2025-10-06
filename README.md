# MethaDory

<img style="float: right;" src="src/shiny/html_imports/methadory.png" width="200">

`MethaDory` is a shiny app for testing of DNAm signatures as described in ['***Training with synthetic data provides accurate and openly-available DNA methylation classifiers for developmental disorders and congenital anomalies via MethaDory***' (Ferraro et al., 2025)](https://www.medrxiv.org/content/10.1101/2025.03.28.25324859v1). 



If you use MethaDory please cite our work and consider starring this repository to follow updates :) 

`MethaDory` is currently in beta testing and is in active development to provide more signatures and optimizations for long-read sequencing, so stay tuned for the latest version! 


## Installation

`MethaDory` is distributed with Pixi (v ≥ 0.55.0), which can be installed the instructions described [here](https://pixi.sh/latest/).

When pixi is available on your system, clone the `MethaDory` git page:

```
git clone git@github.com:f-ferraro/MethaDory.git
```

Then launch the app or one of the cli, installation will be performed automatically. 

## Running MethaDory
To ensure these command are executed correctly from anywhere, specify the full path to the `MethaDory/pixi.toml` included in the main `MethaDory` repository. You can omit this if you're in the main `MethaDory` directory.

You can launch the interactive app with:

```bash
pixi run --manifest-path MethaDory/pixi.toml MethaDory
```

MethaDory relies on a number of files provided in the `data` folder. This folder should be in the same directory where the MethaDory code resides. In the folder you will find:

- `affectedindividuals`, the folder containing filtered and anonymized beta data and meta data of example affected individuals (used for plotting).
- `imputationsamples`, the folder containing samples used for missing value imputation.
- `support_files`, a folder containing
  - `manifest.qc_filtered.rds`, of methylation array manifest after filtering for QC as described in [Ferraro et al., 2025](https://www.medrxiv.org/content/10.1101/2025.03.28.25324859v1).
  - `merged_signatures_90DMRs.tsv`, text file containing information about the sites used for building the DNAm signatures and classifiers.
  - `background_training.cellprops.rds`, containing deconvoluted cell proportions from the samples used for the model training.
 
**Trained classifiers are currently undergoing external validation before public release and are available upon request. Please reach out to test the `data` folder.**


Otherwise a number of CLI are also available. Please specificy full paths to all required inputs. 


### Table and PDF exports
```bash
pixi run --manifest-path MethaDory/pixi.toml MethaDory_cli <model_folder> <sample_file> <output_prefix> [options]


Arguments:
   model_folder: Full path to folder containing SVM model .rds files
   sample_file:  Full path to .tsv file containing sample data
   output_prefix: Prefix for output files (will create .xlsx and .pdf files)

 Options:
   --include-dim-plots      Include dimension reduction plots in PDF (default: TRUE)
   --include-cell-plots     Include cell deconvolution plots in PDF (default: TRUE)
   --include-chr-sex        Include chromosomal sex prediction plots in PDF (default: TRUE)
   --min-psvm               Minimum pSVM threshold for dimension plots (default: 0.05)
   --n-imputation-samples   Number of closest samples for imputation (default: 20)
   --n-samples-plots        Number of additional samples for visualization (default: 20)
   --export-imputed         Export imputed methylation data (default: FALSE)
   --help                   Show this help message
```


### Self-contained HTML (supported for single sample only)
```bash
pixi run --manifest-path MethaDory/pixi.toml MethaDory_html <model_folder> <sample_file> <output.html> [options]


 Arguments:
   model_folder: Full path to folder containing SVM model .rds files
   sample_file:  Full path to .tsv file containing sample data
   output_path:  Full name to the HTML that should be saved

 Options:
   --include-dim-plots      Include dimension reduction plots (default: TRUE)
   --include-cell-plots     Include cell deconvolution plots (default: TRUE)
   --include-chr-sex        Include chromosomal sex prediction plots (default: TRUE)
   --min-psvm               Minimum pSVM threshold for plots (default: 0.05)
   --n-imputation-samples   Number of closest samples for imputation (default: 20)
   --n-samples-plots        Number of additional samples for visualization (default: 20)
   --help                   Show this help message
```




The first time you run `MethaDory`, allow some time (~15') to install the necessary prerequisites.
If the local web browser doesn't open automatically, double click on the link shown on the terminal or paste it in a web browser of choice. 

`MethaDory` relies on web browser being installed and set as default in your system. If you get the error:

>Listening on http://127.0.0.1:... Error in utils::browseURL(appUrl) :  'browser' must be a non-empty character string

R is failing to find the browser or there is none set as default.
You can manually set your browser in `MethaDory` by adding to `src/shiny/app.R` e.g.:

```
options(browser="firefox")
```

On some systems, it might be required to enable the pixi postlinks to successfully complete the installation. In that case follow the terminal instructions from pixi.




## Input Data Format

`MethaDory` expects a tab-separated file with:
- **Column 1**: `IlmnID` (CpG probe names, e.g., `cg12345678`)
- **Remaining columns**: Sample names with beta values (0-1)

An example is provided in `data/demo.input.txt`:
```
IlmnID	GSM3173324	GSM3173369	GSM3173402	Missing50Sites
cg09499020	0.328	0.529	0.384
cg16535257	0.541	0.594	0.371	0.371
cg06325811		0.303	0.212
cg16619049	0.371	0.449		0
cg13938959	0.742	0.634	0.773
cg12445832	0.408	0.46	0.614	0.614
cg11527153	0.93	0.917	0.878
cg04195702	0.88	0.937	0.835	0.835
cg08128007	0.801	0.84	0.869
```

`MethaDory` will perform missing value imputation to ensure all necessary probes are present, however high missing probe rate increases computation time and can reduce the performance of the classifiers.

For Oxford Nanopore Data support see the `ONT-Input-Preparation` folder.