## Methylation data extraction from ONT data for MethaDory analysis

### Installation

Similarly to MethaDory, also the dependencies for these scripts can be obtained with Pixi and are included in the `pixi.toml` in this directory.

The script assumes that the required `ilmn.epic.v2.manifest.annot.1pd.bed` is uncompressed and in the same directory. 

### Usage

To extract methylation data and annotate it with the CpG id from Illumina EPIC v2 you can use the wrapper:

```bash
pixi run Methadory-methylation.sh

Usage: MethaDory-methylation.sh --bam_path FILE --sampleID ID --ref_genome DIR --threads N
 Arguments:
 --bam_path   Path to the ONT BAM file
 --sampleID   Sample ID for output files
 --ref_genome Path to the reference genome directory
 --threads    Number of threads to use
 
```

This will produce two folders: `modkitoutput` and `pseudoepic`. The latter is the input for the next step.

To reformat and merge multiple files together in a single MethaDory input file, you can run R script:

```bash

pixi run Rscript MethaDory_ONT_Input_Preparation.R

Usage: MethaDory_ONT_Input_Preparation.R <pseudoepic_directory> <output_file>
  Arguments:
       pseudoepic_directory   Path to the `pseudoepic` directory output of MethaDory-modkit.sh
       output_file            Path to the output file
```

*remember that MethaDory accepts files up to 300MB*
