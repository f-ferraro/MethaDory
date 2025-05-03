# MethaDory

<img style="float: right;" src="methadory.png" width="200">

MethaDory is a shiny app for testing of DNAm signatures as described in ['***Training with synthetic data provides accurate and openly-available DNA methylation classifiers for developmental disorders and congenital anomalies via MethaDory***' (Ferraro et al., 2025)](https://www.medrxiv.org/content/10.1101/2025.03.28.25324859v1). 



If you use MethaDory please cite our work and consider starring this repository to follow updates :) 

MethaDory is currently in beta testing and is in active development to provide more signatures and optimizations for long-read sequencing, so stay tuned for the latest version! 

**Trained classifiers are currently undergoing external validation before public release and are available upon request. Please reach out to test them.**

## Installation

MethaDory is distributed with Pixi, which can be installed the instructions described [here](https://pixi.sh/latest/).

For Linux and MacOS this boils down to:

```
curl -fsSL https://pixi.sh/install.sh | bash
```
You might need to restart your terminal for the installation to take effect

Then clone the MethaDory git page:

```
git clone git@github.com:f-ferraro/MethaDory.git
```

## Running MethaDory 

From the directory where you cloned MethaDory, run:

```
pixi run MethaDory
```
The first time you run this, allow some time (~15') to install the necessary prerequisites and lanch the app. After installation, you can launch MethaDory with the same command.
If the local web browser doesn't open automatically, double click on the link shown on the terminal. 

MethaDory relies on web browser being installed and set as default in your system. If you get the error:

>Listening on http://127.0.0.1:... Error in utils::browseURL(appUrl) :  'browser' must be a non-empty character string

R is failing to find the browser or there is none set as default. 
You can manually set your browser in MethaDory by adding to `ShinyApp_MethylClassifier.r` e.g.:

```
options(browser="firefox")
```


## MethaDory input data
MethaDory assumes you're providing an input table with the first column called `IlmnID` (containing the probe names, e.g. `cp123456` from an Illumina methylation array), followed by a number of columns each representing a sample. MethaDory will perform missing values imputation for you to ensure all the necessary probes are present.

For ONT data, you can find the scripts needed to tabulate methylation from a bam file and format them for MethaDory in the [ONT-Input-preparation](https://github.com/f-ferraro/MethaDory/tree/main/ONT-Input-preparation) folder.
