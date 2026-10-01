## Methylation data extraction from PacBio data for MethaDory analysis

`MethaDory` can also be used with PacBio HiFi data preprocessed with
[pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools). We provide
scripts to reformat those per-CpG BEDs into the MethaDory input table.

### Installation

The dependencies (`bedtools`, R + dplyr + purrr) are included in the
`pixi.toml` of MethaDory.

The same EPIC manifest used for the ONT pipeline works here too:
`../ONT-Input-preparation/ilmn.epic.s.manifest.annot.1pd.bed.gz`. Uncompress
it after downloading if needed.

### PacBio data

pb-CpG-tools emits two BED variants per sample:

- `<sample>.pb.model.bed.combined.bed.gz` &mdash; `--pileup-mode model`
- `<sample>.pb.count.bed.combined.bed.gz` &mdash; `--pileup-mode count`

pb-CpG-tools defaults to `model`, so that's the default mode of this script.
Both share `chrom / begin / end / mod_score / type / cov` as columns 1&ndash;6;
only the downstream columns differ. `mod_score` is on a 0&ndash;100 scale in
both cases and is converted to a beta value (0&ndash;1) by the R step.

For each sample of interest:

```bash
pixi run --manifest-path MethaDory/pixi.toml MethaDory_pb_extract \
  --pb_file FILE \
  --sampleID ID \
  --manifest MANIFEST \
  [--mode model|count]
```

```
Usage: MethaDory-extract-methylation-pacbio.sh \
       --pb_file FILE --sampleID ID --manifest MANIFEST [--mode model|count]
 Arguments:
 --pb_file    Path to the pb-CpG-tools BED (gzipped or plain).
 --sampleID   Sample ID for output filenames.
 --manifest   Path to the EPIC manifest BED.
 --mode       'model' (default) or 'count'.
```

This writes one `pseudoepic/<sampleID>.pseudoepic.cpgID.bed` per call. If you
process multiple samples into the same working directory, all of them
accumulate in `pseudoepic/` &mdash; existing files from other samples are
preserved.

To reformat and merge multiple samples into a single MethaDory input table,
run the R script:

```bash
pixi run --manifest-path MethaDory/pixi.toml MethaDory_pb_input \
  <pseudoepic_directory> \
  <output_file> \
  [min_cov]
```

```
Usage: MethaDory_PB_Input_Preparation.R \
       <pseudoepic_directory> <output_file> [min_cov]
  Arguments:
       pseudoepic_directory   Directory of *.pseudoepic.cpgID.bed files
       output_file            Output TSV path
       min_cov                Minimum coverage to keep a probe (default 4)
```

The output is a tab-separated table with `IlmnID` as the first column and one
column per sample, with beta values in [0, 1].

*Remember that MethaDory accepts files up to 300 MB.*
