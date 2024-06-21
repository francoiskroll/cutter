#!/bin/bash

shopt -s nullglob

for i in *.fastq # loop thru files in directory that finish by .fastq (we unzipped all before)
do

  FILE="$i"
  DIR="$(echo "$i" | cut -d'_' -f 4)" # DIR (direction) should be R1 = Forward reads or R2 = Reverse reads

  if [ "$DIR" = "R1" ]; then # do I have the Forward reads?
  # if yes: need to find the Reverse reads
    FWD="$FILE"

    HALF1="$(echo "$FWD" | cut -d'R' -f 1)" # everything before 'R'
    HALF2="$(echo "$FWD" | cut -d'_' -f 5)" # everything after '_'
    RVS="$(echo $HALF1$"R2_"$HALF2)" # based on the name of the Forward file, this should be the name of the Reverse file
  fi

  WELL="$(echo "$FWD" | cut -d'_' -f 1 | cut -d'-' -f 2)" # first cut gets everything before first '_'; second cut gets everything after first '-'

  echo
  echo
  echo "---- [ CRISPResso2 on "$"$WELL ] ----"
  echo
  echo

  # now run CRISPResso2
  # see README.md in 23117_sequencing/231117_miseq/plasmidsaurus for details on how I built the command
  CRISPResso --fastq_r1 "$FWD" --fastq_r2 "$RVS" --amplicon_seq GTACAGTCTGGTGTGGCTCATAAGCCCCATTTTGGGTTTTATCCTACAGCCCGTCATCGGCTCGGCGAGCGACTACTGTAGGTCGTCATAAGGCCGAAGGAGACCGTACATACTCTTACTGGGGATTCTGATGTTAGTGGGCATGACTTTATTTCTAAATGGAGATGCAGTCACAACAGGTGGGTGA --amplicon_name slc45a2TAA --prime_editing_pegRNA_spacer_seq gactactgtaggtcgtcata --prime_editing_pegRNA_extension_seq tctccttcggccccatgacgacctacagt --prime_editing_pegRNA_scaffold_seq gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtgggaccgagtcggtcc

done

shopt -u nullglob
