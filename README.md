# MethaDory

MethaDory is a shiny app for testing of DNAm signatures as described in XXX. 

If you use MethaDory please cite our work and consider starring this repository :) 


## Installation

MethaDory environment is distributed as Pixi that can be easily installed with following the instructions described [here](https://pixi.sh/latest/).

For Linux and MacOS this boils down to:

```
curl -fsSL https://pixi.sh/install.sh | bash
```
You might need to restart your terminal for the installation to take effect

### Clone the repository 

Now clone the MethaDory git page and download ChAMPdata v2.38.0

```
git clone git@github.com:f-ferraro/MethaDory.git #or download the repository locally, unzip it
cd MethaDory #open the folder of the repository
wget https://bioconductor.org/packages/release/data/experiment/src/contrib/ChAMPdata_2.38.0.tar.gz
```

## Running MethaDory 

```
pixi run rscript ShinyApp-MethylClassifier.r
```
The command will take care of installing the necessary prerequisites and launching the app.
If the local web browser doesn't open automatically, double click on the link shown on the terminal. 
