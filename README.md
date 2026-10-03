# MethaDory

<p align="center">
  <img src="src/shiny/html_imports/methadory.png" width="200">
</p>

`MethaDory` is a shiny app for testing of DNAm signatures as described in ['***Training with synthetic data provides accurate and openly-available DNA methylation classifiers for developmental disorders and congenital anomalies via MethaDory***' (Ferraro et al., 2025)](https://www.medrxiv.org/content/10.1101/2025.03.28.25324859v1). 



If you use MethaDory please cite our work and consider starring this repository to follow updates. 

`MethaDory` is currently in beta testing and is in active development to provide more signatures and optimizations for long-read sequencing, so stay tuned for the latest version! 


## Installation

`MethaDory` is distributed with Pixi (v ≥ 0.55.0), which can be installed the instructions described [here](https://pixi.sh/latest/).

When pixi is available on your system, clone the `MethaDory` git page:

```
git clone git@github.com:f-ferraro/MethaDory.git
```

Then launch the app or one of the cli, installation will be performed automatically.

The first time you run `MethaDory`, allow some time (~15') to install the necessary prerequisites. You might have to launch the command a couple of times to install all the dependencies. It might be required to enable the pixi postlinks to successfully complete the installation; in that case follow the terminal instructions from pixi.

## Running MethaDory
To ensure these command are executed correctly from anywhere, specify the full path to the `MethaDory/pixi.toml` included in the main `MethaDory` repository. You can omit this if you're in the main `MethaDory` directory. Please specify full paths to all required inputs.

MethaDory can be run in two ways, both built on the same analysis pipeline and producing the same predictions:

1. **Self-contained HTML report** - the primary output: one shareable file per sample, with every plot and table embedded, plus the result tables as one `.xlsx` workbook for the whole run.
2. **Interactive app** - explore the results, changing selections and thresholds interactively.

For help with plot interpretation see `Manual.md` document in this repo. 

### Data folder

MethaDory relies on a number of files provided in the `data` folder. This folder should be in the same directory where the MethaDory code resides. In the folder you will find:

- `affectedindividuals`, the folder containing filtered and anonymized beta data and meta data of example affected individuals (used for plotting).
- `imputationsamples`, the folder containing samples used for missing value imputation.
- `models`, the trained classifiers, organised per input platform (`arrays`, `ont`, `pb`). Each platform folder holds one model set containing an `SVM/` and an `NNET/` subfolder. The folder you point MethaDory at is the one containing `SVM/` and `NNET/` (e.g. `data/models/arrays/SVMa10b1_NNETa15b1_sig<...>`).
- `support_files`, a folder containing
  - `manifest.qc_filtered.rds`, of methylation array manifest after filtering for QC as described in [Ferraro et al., 2025](https://www.medrxiv.org/content/10.1101/2025.03.28.25324859v1).
  - `merged_signatures_90DMRs.tsv`, text file containing information about the sites used for building the DNAm signatures and classifiers.
  - `background_training.cellprops.rds`, containing deconvoluted cell proportions from the samples used for the model training.



### Predictions

Each signature is scored by two independent classifier families, a support vector machine (SVM) and a neural network (NNET), each trained over several checkpoints. Reported results include the per-family averages (`pSVM_average`, `pNNET_average`) with their standard deviations, and the metapredictor score `pCombined`, the mean of the two. `pCombined` is the score used for ranking, for the confidence bins in the summary tables, and for the `--min-p` threshold that decides which signatures get dimension reduction plots.

### Running the Container

```
docker run -it -v ./:/data ghcr.io/f-ferraro/MethaDory:latest MethaDory_cli /app/data/models/<platform>/<model>/ /data/demo.input.txt /data/output
```



### 1. Self-contained HTML report (primary output)

The recommended way to run MethaDory is the HTML report generator. It runs the full pipeline and writes a self-contained `.html` file for each sample with the prediction plot and table, the dimension reduction figures and the sample-level QC panels all embedded in it.

Multi-sample input files are supported: the pipeline runs once over the whole file and the samples are then rendered one at a time. In that case the third argument is the **output directory**, and one report per sample is written to it as `<SampleName>.MethaDory-output.html`. With a single-sample input the third argument is the `.html` file to write, as before. Running in this way is faster than running multiple samples at the same time. We reccomend splitting the samples per technology, i.e. test at once only arrays, or ONT, or PacBio data.

```bash
pixi run --manifest-path MethaDory/pixi.toml MethaDory_html <model_folder> <sample_file> <output.html|output_dir> [options]


 Arguments:
   model_folder: Full path to the model folder containing the SVM/ and NNET/ subfolders
   sample_file:  Full path to .tsv file containing sample data
   output_path:  Single sample: full name of the HTML file to write.
                 Several samples: the output directory, which will receive one
                 <SampleName>.MethaDory-output.html per sample.

 Options:
   --include-dim-plots      Include dimension reduction plots (default: TRUE)
   --include-cell-plots     Include cell deconvolution plots (default: TRUE)
   --include-chr-sex        Include chromosomal sex prediction plots (default: TRUE)
   --min-p                  Keep a signature for the per-signature dimension plots when its combined
                            score (pCombined, the mean of the SVM and NNET scores) is at or above
                            this value, between 0 and 1 (default: 0.20). Tables and the prediction
                            plot always show every signature.
   --n-imputation-samples   Number of closest samples for imputation (default: 20)
   --n-samples-plots        Number of additional samples for visualization (default: 20)
   --export-xlsx            Also write the result tables as an .xlsx workbook (default: TRUE)
   --export-imputed         Also write the imputed methylation data (user samples + controls +
                            real cases) as <report>_imputed.tsv (default: FALSE)
   --help                   Show this help message
```

Files written next to the report(s):

| File | When | Content |
|---|---|---|
| `<report>.html`, or `<SampleName>.MethaDory-output.html` per sample | always | the self-contained report, one per sample |
| `<report>.xlsx`, or `<input file name>.MethaDory-output.xlsx` | unless `--export-xlsx FALSE` | the result tables, **one workbook for the whole run** with all samples together: predictions (long and wide), per-signature summary, cell proportions, methylation age, chromosomal sex |
| `<report>_imputed.tsv`, or `<input file name>.MethaDory-output_imputed.tsv` | with `--export-imputed TRUE` | `IlmnID` x samples: your samples after imputation, plus the background controls and the real cases they are compared with |

The tables are written before the reports are rendered, so a report that fails does not cost the numbers.

Example, for a single Nanopore proband:

```bash
pixi run --manifest-path MethaDory/pixi.toml MethaDory_html \
  /path/to/data/models/ont/SVMa10b1_NNETa15b1_sig<...> \
  /path/to/sample.tsv /path/to/report.html
```

Example, for a file holding several samples, drawing figures only for signatures scoring 0.25 or more and exporting the imputed data:

```bash
pixi run --manifest-path MethaDory/pixi.toml MethaDory_html \
  /path/to/data/models/arrays/SVMa10b1_NNETa15b1_sig<...> \
  /path/to/samples.tsv /path/to/reports/ \
  --min-p 0.25 --export-imputed TRUE
```

This writes `reports/<SampleName>.MethaDory-output.html` for each sample, plus `reports/samples.MethaDory-output.xlsx` and `reports/samples.MethaDory-output_imputed.tsv` for the run.

### 2. Interactive app

The app is the exploratory front-end: it runs the same pipeline, but lets you re-select probands, signatures and thresholds and redraw the figures without recomputing. Launch it with:

```bash
pixi run --manifest-path MethaDory/pixi.toml MethaDory
```

#### App options

Beyond selecting the model folder and uploading the input file, the sidebar exposes:

- **Number of closest samples for imputation**: how many nearest background samples are used to impute missing probes (default 20).
- **Select Proband(s) / Select Signature(s)**: restrict which samples and signatures are plotted.
- **Minimum mean(SVM,NNET) for dimension plots**: score threshold above which a signature gets a PCA/heatmap panel (default 0.05).
- **Number of additional samples to use for PCA and heatmap**: reference samples added per group to the dimension plots (default 20).
- **Export Options**: *Download Results* writes the tables, *Download Plots* writes a self-contained HTML report reflecting the current selections.

#### If the app doesn't open

If the local web browser doesn't open automatically, double click on the link shown on the terminal or paste it in a web browser of choice.

`MethaDory` relies on web browser being installed and set as default in your system. If you get the error:

>Listening on http://127.0.0.1:... Error in utils::browseURL(appUrl) :  'browser' must be a non-empty character string

R is failing to find the browser or there is none set as default.
You can manually set your browser in `MethaDory` by adding to `src/shiny/app.R` e.g.:

```
options(browser="firefox")
```


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


`MethaDory` uses models trained on Noob-normalized arrays. Please use this normalization on your data before testing new samples for best performance. See for an example the `Input-preparation-Array` folder.
For Oxford Nanopore data support see the `Input-preparation-ONT` folder.
For PacBio data support see the `Input-preparation-PB` folder.
