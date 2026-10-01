#!/bin/bash
# Reformat a pb-CpG-tools per-CpG BED into a "pseudoepic" BED keyed by
# Illumina EPIC probe IDs, ready for MethaDory_PB_Input_Preparation.R.

usage() {
  cat <<EOF
Usage: $0 --pb_file FILE --sampleID ID --manifest MANIFEST [--mode model|count]
 Arguments:
 --pb_file    Path to the pb-CpG-tools output BED (gzipped or plain).
              Expects the standard pb-CpG-tools "combined" BED, in either
              pileup-mode=model (default for pb-CpG-tools) or pileup-mode=count.
 --sampleID   Sample ID used for output filenames.
 --manifest   Path to the EPIC manifest BED (e.g.
              ../ONT-Input-preparation/ilmn.epic.s.manifest.annot.1pd.bed[.gz]).
 --mode       'model' (default) or 'count'. Selects which column carries the
              Illumina probe ID after bedtools intersect; the rest of the
              downstream pipeline is mode-agnostic.
EOF
  exit 1
}

mode="model"
while [[ $# -gt 0 ]]; do
  case $1 in
    --pb_file)  pb_file="$2";  shift 2 ;;
    --sampleID) sampleID="$2"; shift 2 ;;
    --manifest) manifest="$2"; shift 2 ;;
    --mode)     mode="$2";     shift 2 ;;
    *) usage ;;
  esac
done

if [[ -z "$pb_file" || -z "$sampleID" || -z "$manifest" ]]; then usage; fi


case "$mode" in
  model) ilmn_col=13 ;;
  count) ilmn_col=14 ;;
  *) echo "Bad --mode: $mode"; usage ;;
esac

mkdir -p pseudoepic

zcat -f "$pb_file" \
  | grep -v '^#' \
  | bedtools intersect -a - -b "$manifest" -wb -wa \
  | awk -v OFS="\t" -v c="$ilmn_col" '{print $1,$2,$3,$4,$5,$6,$c}' \
  > "pseudoepic/${sampleID}.pseudoepic.cpgID.bed"

echo "Wrote pseudoepic/${sampleID}.pseudoepic.cpgID.bed" >&2
