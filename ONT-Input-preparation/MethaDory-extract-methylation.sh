#!/bin/bash

usage() {
  echo "Usage: $0 --bam_path FILE --sampleID ID --ref_genome DIR --manifest MANIFEST --threads N"
  echo " Arguments:"
  echo " --bam_path   Path to the BAM file"
  echo " --sampleID   Sample ID for output files"
  echo " --ref_genome Path to the reference genome directory"
  echo " --manifest   Path to the manifest file"
   echo " --threads    Number of threads to use"
  exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --bam_path)
      bam_path="$2"
      shift 2
      ;;
    --sampleID)
      sampleID="$2"
      shift 2
      ;;
    --ref_genome)
      ref_genome="$2"
      shift 2
      ;;
    --threads)
      threads="$2"
      shift 2
      ;;
    --manifest)
      manifest="$2"
      shift 2
      ;;
    *)
      usage
      ;;
  esac
done

# Check if all required args are provided
if [[ -z "$bam_path" || -z "$sampleID" || -z "$ref_genome" || -z "$threads" || -z "$manifest" ]]; then
  usage
fi

# Create output directories
mkdir -p modkitoutput
mkdir -p pseudopedic

# Run modkit
modkit pileup "$bam_path" \
  --threads "$threads" \
  --combine-mods --cpg --combine-strands \
  --ref "${ref_genome}" \
  --filter-threshold 0.5 \
  "modkitoutput/${sampleID}.modkit.unstranded.combined.bed" \
  --include-bed $manifest

# Intersect with EPIC manifest and reformat
bedtools intersect \
  -a "modkitoutput/${sampleID}.modkit.unstranded.combined.bed" \
  -b $manifest \
  -wb -wa | \
  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$22}' \
  > "pseudopedic/${sampleID}.pseudoepic.cpgID.bed"
