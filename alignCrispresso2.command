#!/bin/bash

# 31/01/2025
# was making a new script for each new amplicon
# this script alignCrispresso2.command will

shopt -s nullglob

#####################

### function checkPath to check that a path exists
checkPath() {
    local file_path="$1"

    # Check if the file exists
    if [ -e "$file_path" ]; then
        :
    else
        echo "Error. File does not exist: $file_path"
        exit 1
    fi
}

### function checkDir to check path leads to an existing directory
checkDir() {
    local folder_path="$1"

    # Check if the file exists
    if [ -d "$folder_path" ]; then
        :
    else
        echo "Error. Folder does not exist: $folder_path"
        exit 1
    fi
}

#####################
# read the arguments
# only one now is path to folder with fasta references, typically /refseqs/
while getopts a: flag
do
    case "${flag}" in
        a) refpath=${OPTARG};;
    esac
done

# check directory exists
checkDir "$refpath"

#####################
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
  # if first file has _R1 in its name
  if [[ "$firstnm" == *"_R1"* ]]; then
    FWD="$first"
    RVS="$second"
  else
    RVS="$first"
    FWD="$second"
  fi

  echo "Forward reads = $FWD"
  echo "Reverse reads = $RVS"


  ### find the fasta reference
  # get the amplicon name
  # will assume filename looks like well_amplicon_R1/R2.fastq
  # where 'amplicon' could also have underscores
  # so remove anything before first underscore and anything after last underscore
  # > remove everything before the first underscore
  AMP="${FILE#*_}"
  # > remove everything after the last underscore
  AMP="${AMP%_*}"
  # add .fa to the name
  REF="$(echo "$AMP"$".fa")"

  echo "Reference should be $REF"

  # need to find it in /refseqs/ folder
  # path should be
  refp="$(echo "$refpath"$"/""$REF")"

  # check that it exists
  checkPath "$refp"

  ### get the actual sequence
  refseq=$(tail -n +2 "$refp") # +2 skips first row which is fasta header


  #####################
  ### in summary:
  echo
  echo
  echo "---- [ CRISPResso2 on "$"$FWD & $RVS ] ----"
  echo " amplicon sequence is "$refseq" "
  echo
  echo

  # now run CRISPResso2
  # create output folder, in parent folder of folder containing the reads
  # get full path of FWD file
  fwdp=$(realpath "$FWD")
  # get the folder in which it is
  readdir="$(dirname "$fwdp")"
  # get the parent folder
  pardir="$(dirname "$readdir")"
  # this is where we create output folder
  outdir="$(echo "$pardir"$"/crispresso")"
  mkdir "$outdir"

  # previously, was doing actual prime-editing mode
  # (see README.md in 23117_sequencing/231117_miseq/plasmidsaurus for details on how I built the command)
  # command was:
  # CRISPResso --fastq_r1 "$FWD" --fastq_r2 "$RVS" --amplicon_seq GTACAGTCTGGTGTGGCTCATAAGCCCCATTTTGGGTTTTATCCTACAGCCCGTCATCGGCTCGGCGAGCGACTACTGTAGGTCGTCATAAGGCCGAAGGAGACCGTACATACTCTTACTGGGGATTCTGATGTTAGTGGGCATGACTTTATTTCTAAATGGAGATGCAGTCACAACAGGTGGGTGA --amplicon_name slc45a2TAA --prime_editing_pegRNA_spacer_seq gactactgtaggtcgtcata --prime_editing_pegRNA_extension_seq tctccttcggccccatgacgacctacagt --prime_editing_pegRNA_scaffold_seq gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtgggaccgagtcggtcc
  
  # now, I prefer to simply use CRISPResso for alignment,
  # and call all the mutations myself
  CRISPResso --fastq_r1 "$FWD" --fastq_r2 "$RVS" --amplicon_seq "$refseq" --output_folder "$outdir"
done

shopt -u nullglob