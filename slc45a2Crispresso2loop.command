#!/bin/bash

shopt -s nullglob

for i in *.fastq # loop thru files in directory that finish by .fastq (we unzipped all before)
do

  FILE="$i"

  ### get the well
  # will assume always first in filename, before first underscore _
  # e.g. A01_xxx
  WELL="$(echo "$i" | cut -d'_' -f 1)"

  ### in the folder, find files that start with WELL
  fastqpath=$(dirname "$FILE")
  READS=$(find "$fastqpath" -type f -name "$WELL*")

  # check that exactly two fastq files were found
  readscount=$(echo "$READS" | wc -l)

  if [ "$readscount" -ne 2 ]; then
    echo "Error: there should be exactly TWO fastq files starting with '$WELL', but found $readscount."
    exit 1
  fi

  # will assume one is R1, one is R2
  # but will not assume which is which to be safe
  first=$(echo "$READS" | head -n 1)
  second=$(echo "$READS" | head -n 2 | tail -n 1)
  # keep only filename
  firstnm=$(basename "$first")
  # if first file has R1 in its name
  if [[ "$firstnm" == *"R1"* ]]; then
    FWD="$first"
    RVS="$second"
  else
    RVS="$first"
    FWD="$second"
  fi

  ### in summary:
  echo
  echo
  echo "---- [ CRISPResso2 on "$"$FWD & $RVS ] ----"
  echo
  echo

  # now run CRISPResso2

  # previously, was doing actual prime-editing mode
  # (see README.md in 23117_sequencing/231117_miseq/plasmidsaurus for details on how I built the command)
  # command was:

  # CRISPResso --fastq_r1 "$FWD" --fastq_r2 "$RVS" --amplicon_seq GTACAGTCTGGTGTGGCTCATAAGCCCCATTTTGGGTTTTATCCTACAGCCCGTCATCGGCTCGGCGAGCGACTACTGTAGGTCGTCATAAGGCCGAAGGAGACCGTACATACTCTTACTGGGGATTCTGATGTTAGTGGGCATGACTTTATTTCTAAATGGAGATGCAGTCACAACAGGTGGGTGA --amplicon_name slc45a2TAA --prime_editing_pegRNA_spacer_seq gactactgtaggtcgtcata --prime_editing_pegRNA_extension_seq tctccttcggccccatgacgacctacagt --prime_editing_pegRNA_scaffold_seq gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtgggaccgagtcggtcc
  
  # now, I prefer to simply use CRISPResso for alignment,
  # and call all the mutations myself
  CRISPResso --fastq_r1 "$FWD" --fastq_r2 "$RVS" --amplicon_seq GTACAGTCTGGTGTGGCTCATAAGCCCCATTTTGGGTTTTATCCTACAGCCCGTCATCGGCTCGGCGAGCGACTACTGTAGGTCGTCATAAGGCCGAAGGAGACCGTACATACTCTTACTGGGGATTCTGATGTTAGTGGGCATGACTTTATTTCTAAATGGAGATGCAGTCACAACAGGTGGGTGA
done

shopt -u nullglob