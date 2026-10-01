## Methylation data preprocessing from Illumina arrays for MethaDory analysis

`MethaDory` classifiers were trained on Noob-normalized Illumina methylation
arrays. We provide scripts to normalize raw IDAT files, filter low-quality
probes and write the MethaDory input table, one for EPICv1 and one for EPICv2.

### Installation

The scripts require R with [minfi](https://bioconductor.org/packages/minfi/), [vroom](https://cran.r-project.org/package=vroom) and the GitHub development build of [ChAMP](https://github.com/YuanTian1991/ChAMP), which is the only version with EPICv2 support. The scripts load ChAMP from `~/R-champ-dev`, so install it there to keep any Bioconductor ChAMP installation intact:

```r
remotes::install_github("YuanTian1991/ChAMPdata", lib = "~/R-champ-dev")
remotes::install_github("YuanTian1991/ChAMP",     lib = "~/R-champ-dev")
```

EPICv2 data also need the minfi annotation packages:

```r
BiocManager::install(c("IlluminaHumanMethylationEPICv2manifest",
                       "IlluminaHumanMethylationEPICv2anno.20a1.hg38"))
```

### Array data

Place the paired `<sample>_Grn.idat(.gz)` and `<sample>_Red.idat(.gz)` files
of all samples in one directory (subdirectories are searched too). The sample
name is taken from the file name without the `_Grn`/`_Red` suffix.

For EPICv1 arrays:

```bash
Rscript MethaDory_EPICv1_Input_Preparation.R \
          <idat_directory> \
          <output_file>
```

For EPICv2 arrays:

```bash
Rscript MethaDory_EPICv2_Input_Preparation.R \
          <idat_directory> \
          <output_file>
```

```
Usage: MethaDory_EPICv1_Input_Preparation.R <idat_directory> <output_file>
       MethaDory_EPICv2_Input_Preparation.R <idat_directory> <output_file>
  Arguments:
       idat_directory   Directory of paired *_Grn.idat(.gz) and *_Red.idat(.gz) files
       output_file      Path to the output TSV
```

