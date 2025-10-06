# MethaDory Analysis Report

## About MethaDory

Changes in DNA methylation (DNAm) at characteristic genomic loci - **DNAm signatures** - have been associated with genomic variants and environmental exposures. These DNAm signatures enable molecular classification or diagnosis for numerous disorders and have shown potential as biomarker to estimate disease prognosis and to assess therapeutic responses.

**MethaDory** is a tool to recognize DNAm signatures associated with various diseases or exposures in blood samples from patients. If you use MethaDory, please cite our work: [Ferraro et al., 2025](https://www.medrxiv.org/content/10.1101/2025.03.28.25324859v1)

## Interpreting Results

The **Prediction Results Plot** tab shows a graphical representation of any DNAm signatures that scored at least as high as the threshold set at runtime. You can find the raw scores in tabular form in **SVM Prediction Table**. A positive score alone may not sufficient to classify a sample, this is possible only with matching phenotype and genetic variant (for genetic disorders).

Be also aware that similar signatures are associated with related genes and overlapping pathways.

The **Dimension Reduction Plots** tab shows the sample of interest together with control samples from healthy individuals and samples from individuals affected by the disorder of the specific signature (or synthetic cases), to further corroborate the results. A genuine positive sample should cluster together with samples from the specific disorder and separately from controls.

For quality control purposes **Cell Proportion Deconvolution** and **Methylation Age prediction** can be inspected.
