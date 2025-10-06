## Methylation data extraction from ONT data for MethaDory analysis
`MethaDory` can also be used with Oxford Nanopore Technologies (ONT) data. We provide the scripts to extract and reformat the data use the preprocessing tool.

### Installation

The dependencies for these scripts are included in the `pixi.toml` of MethaDory.

We provide the required `ilmn.epic.s.manifest.annot.1pd.bed.gz` containing the coordinates we use.

### ONT data 
To ensure these command are executed correctly, specify the full path to the `MethaDory/pixi.toml` from the main `MethaDory` repository.

For each sample of interest:

```bash
pixi run --manifest-path MethaDory/pixi.toml MethaDory_ont_extract \
  --bam_path FILE \
  --sampleID ID \
  --ref_genome DIR \
  --threads N \
  --manifest MANIFEST

Usage: MethaDory-extract-methylation.sh --bam_path FILE --sampleID ID --ref_genome DIR --manifest MANIFEST --threads N
 Arguments:
 --bam_path   Full path to the BAM file
 --sampleID   Sample ID for output files
 --ref_genome Full path to the reference genome directory
 --manifest   Full path to the manifest file
 --threads    Number of threads to use
 
```

This will produce two folders: `modkitoutput` and `pseudoepic`. If this command will be used for multiple samples and these folder already exist, the output will be written in them preserving the existing data from other samples. The `pseudoepic` folder is the input for the next step.

To reformat and merge multiple files together in a single `MethaDory` input file, you can run R script:

```bash
pixi run --manifest-path MethaDory/pixi.toml MethaDory_ont_input \
  <pseudoepic_directory> \
  <output_file>

Usage: MethaDory_ONT_Input_Preparation.R <pseudoepic_directory> <output_file>
  Arguments:
       pseudoepic_directory   Full path to the `pseudoepic` directory output of MethaDory-modkit.sh
       output_file            Full path to the output file
```

*remember that MethaDory accepts files up to 300MB*
