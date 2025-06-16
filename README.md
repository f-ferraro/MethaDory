# MethaDory

<img style="float: right;" src="methadory.png" width="200">

MethaDory is a shiny app for testing of DNAm signatures as described in ['***Training with synthetic data provides accurate and openly-available DNA methylation classifiers for developmental disorders and congenital anomalies via MethaDory***' (Ferraro et al., 2025)](https://www.medrxiv.org/content/10.1101/2025.03.28.25324859v1). 



If you use MethaDory please cite our work and consider starring this repository to follow updates :) 

MethaDory is currently in beta testing and is in active development to provide more signatures and optimizations for long-read sequencing, so stay tuned for the latest version! 

Trained classifiers are currently undergoing external validation before public release and are available upon request. Please reach out to test them.

## Installation

MethaDory environment is distributed with Pixi, which can be installed with following the instructions described [here](https://pixi.sh/latest/).

For Linux and MacOS this boils down to:

```
curl -fsSL https://pixi.sh/install.sh | bash
```
You might need to restart your terminal for the installation to take effect

### Clone the repository 

To obtain MethaDory clone its git page

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


## MethaDory input data
MethaDory assumes you're providing an input table with the first column called `IlmnID` (containing the probe names, e.g. `cp123456` from an Illumina methylation array), followed by a number of columns each representing a sample. MethaDory will perform missing values imputation for you to ensure all the necessary probes are present.

You will also need the trained classifiers, which are currently undergoing external validation before public release and are available upon request. Please reach out to test them.
