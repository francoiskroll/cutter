#!/bin/bash

# requires picard
# conda version had a conflict with java
# homebrew version worked
# brew install picard-tools

# filter alignments in a BAM file

# currently handles
  # minimum Phred score
  # minimum reference span (previously read length)
    # Note, reference span is how much of the reference the read covers
    # it includes deletions/insertions but *not* soft-clipping, as this is not officially part of the alignment
    # e.g. read has 150 bp aligned + 10 bp deletion + 20 bp soft-clipped, reference span is 150 + 10 = 160 bp
    # Note, this is a better solution than looking at read length (i.e. number of characters, without considering the alignment)
    # Indeed we would have a risk of removing reads with big deletions & we would be counting soft-clipping even though it is not part of the alignment
  # ~~~ maximum read length (not currently in this version, see in v1 if need to add it back)
  # maximum proportion soft-clipped
  # for CRISPR experiments: read covering the likely double-strand break site (4 bp next to N of the NGG/CCN PAM) ± padding, e.g. 20 bp means at least 20 bp on each side of double-strand break site
  # only primary alignments (or all)

# explanation
  # filterBam.command -i XXX.bam -a ./refseqs/ -e int -f int -s float -d int -p (yes or no) -o XXX.bam
    # -i = input = bam file to process
    # -a = path to folder with reference sequences in fasta
    # -e = PhrEd score = minimum Phred score
    # -f = floor = minimum read span
    # -s = maximum proportion soft-clipped
    # -d = padding around double-strand break site
    # -p = keep only primary alignments? yes or no
    # -o = output = name of bam file to output


# COMMENTS
# to turn OFF a filter, simply do not provide the corresponding parameter
# ! remove unmapped reads by default – otherwise why filter?

# examples

# filterBAM.command -i C01.bam -e 40 -f 140 -s 0.2 -d 20 -p no -o CO1_filtered.bam

# filterBAM.command -i C01.bam -e 40 -s 0.2 -d 20 -p no -o CO1_filtered.bam
# would turn off minimum reference span filter

### VERSIONS ###

# v2: relies more heavily on R script readsToThrow.R to name the reads to remove
# easier as can edit there if want to add filters

# v3: added filter to exclude reads not covering likely DSB site ± padding, see readsToThrow.R

######################################

### function checkPath to check paths given
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

######################################

# read the flags/arguments
while getopts i:a:e:f:s:d:p:o: flag
do
    case "${flag}" in
        i) bam=${OPTARG};;
        a) refpath=${OPTARG};;
        e) min_phred=${OPTARG};;
        f) min_readspan=${OPTARG};;
        s) max_softprop=${OPTARG};;
        d) dsbpad=${OPTARG};;
        p) primary=${OPTARG};;
        o) out=${OPTARG};;
    esac
done

### check paths given by user
checkPath "$bam"
checkDir "$refpath"

# Note, because of the way command line arguments are read into R, we need every parameter to exist
# so simply turn any parameter NOT given by user into 'off' so we know we should skip this filter in R
###
if [ -z "$min_phred" ]
then
  min_phred=$"off"
fi

###
if [ -z "$min_readspan" ]
then
  min_readspan=$"off"
fi

###
if [ -z "$max_softprop" ]
then
  max_softprop=$"off"
fi

###
if [ -z "$dsbpad" ]
then
  dsbpad=$"off"
fi

# slightly different here, if not given we should *not* remove secondary alignments, so turn to NO
if [ -z "$primary" ]
then
  primary=$"no"
fi


echo
echo "Input BAM = $bam";
echo "Minimum Phred score = $min_phred";
echo "Minimum reference span (bp) = $min_readspan";
echo "Maximum proportion soft-clipped = $max_softprop";
echo "Padding around DSB site (bp) = $dsbpad";
echo "Keep only primary alignments = $primary";
echo "Output BAM = $out";
echo

######################################

nreads=$(samtools view -c $bam)
echo "Number of reads: $nreads"

######################################

# remove unmapped reads
# prepare file, should be in folder with input bam file
outdir=$(dirname "$bam")
bamm="$(echo "$outdir"$"/tmp0.bam")" # bamm for BAM Mapped
samtools view -b -F 4 $bam > $bamm

nreads=$(samtools view -c $bamm)
echo "Number of alignments: $nreads"

######################################

# will write a temporary SAM to remove secondary alignments and/or filter with readsToThrow.R script
# again, should in the folder with the input BAM
tmp1="$(echo "$outdir"$"/tmp1.sam")"
samtools view -h -o $tmp1 $bamm # convert to sam for below

######################################

# remove secondary alignments if needed
  # write temporary SAM = without secondary alignments
  # (or just a copy of tmp2.sam if -primary no)
# prepare file name
tmp2="$(echo "$outdir"$"/tmp2.sam")"

if [ $primary = "yes" ]
then

  # counting number of alignments before
  before=$(samtools view -c $tmp1)

  # removing secondary alignments
  samtools view -h -F 256 -F 4 -F 2048 $tmp1 > $tmp2
    # See https://broadinstitute.github.io/picard/explain-flags.html
    # Note I think - F 4 (unmapped) -F 256 (not primary) are either not used or were removed at some point before
      # but safer to keep
    # -F 2048 = supplementary alignment, seems to include -2064 as well
    # but if I explicitly add -F 2064, it removes all Reverse (-F 16) for some reason

  # counting number of alignments after
  after=$(samtools view -c $tmp2)

  # how many did we remove?
  ndel=$( expr $before - $after )
  echo "Number of secondary alignments removed: $ndel"

else # if No (or anything else)
  # just store copy of tmp1 as tmp2 so we continue below
  cp $tmp1 $tmp2

fi

# make a BAM copy for below (after readsToThrow.R, for FilterSamReads)
tmp2b="$(echo "$outdir"$"/tmp2.bam")"
samtools view -bS $tmp2 -o $tmp2b

######################################

### now give tmp2.sam to readsToThrow.R
# removereads.txt = list of read names to be removed, called by readsToThrow.R
removereadstxt="$(echo "$outdir"$"/removereads.txt")" # write an empty file
touch "$removereadstxt"

### get complete paths
# for tmp2.sam
tmp2full=$(realpath "$tmp2")

# for folder with fasta references
refpathfull=$(realpath "$refpath")

# for removereads.txt
removereadstxtfull=$(realpath "$removereadstxt")

# now run the R script
### directory where the present script is
# this will allow us to locate R script, as cannot make it run from anywhere
utilsdir=$(dirname "$0")
pathRscript="$(echo "$utilsdir"$"/readsToThrow.R")"

Rscript "$pathRscript" \
  $tmp2full \
  $refpathfull \
  $min_phred \
  $min_readspan \
  $max_softprop \
  $dsbpad \
  $removereadstxtfull

# tell user number of reads which we will remove after readsToThrow.R
ndel=$(wc -l < "$removereadstxt")
echo "Number of alignments not passing filters:$ndel"

# now exclude these reads from tmp2b.bam (the reads listed in removereads.txt)
tmp3=$"tmp3.bam"

# previous solution; but takes ages: samtools view -h $tmp3b | grep -vf removereads.txt | samtools view -bS -o $tmp4
# v2
echo # FilterSamReads still talks even though quiet, so leave some space
echo

# ! skip that step if 0 reads to remove
if [ $ndel -ne 0 ] # if not 0 reads to remove
then

picard FilterSamReads \
    I=$tmp2b \
    O=$tmp3 \
    READ_LIST_FILE="$removereadstxt" FILTER=excludeReadList \
    use_jdk_deflater=true \
    use_jdk_inflater=true \
    quiet=TRUE

else
  # just store copy of tmp2b as tmp3 so we can continue below
  cp $tmp2b $tmp3

fi

# Note: R script readsToThrow.R will write read names in removereads.txt
# but one read can have multiple alignments in the SAM/BAM file
# in that case picard's FilterSamReads will then remove multiple rows (alignments)
# this is by design, as we are interested here in filtering reads, not alignments
# as they will then be converted back into fastq in alignFilterBackSingle.command

# sort tmp2.bam and write final bam file
# index the final bam file
samtools sort $tmp3 > $out
samtools index $out

# tell final number of alignments
fin=$(samtools view -c $out)
echo
echo
echo "Final number of alignments: $fin"
echo
echo
echo

######################################

# clean after ourselves
rm $bamm
rm $tmp1
rm $tmp2
rm $tmp2b
rm $tmp3
rm $removereadstxt