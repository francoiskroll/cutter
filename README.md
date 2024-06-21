# miseqUtils

Main script is **alignFilterBack.command**, which aligns MiSeq fastq reads to a fasta reference (uses `bwa mem`), then optionally filters the resulting BAM file and converts it back to fastq files; one R1 containing Forward reads, one R2 containing Reverse reads. To know which fasta reference to use for each pair of fastq files, it reads a config.xlsx file with two columns: _well_ & _ref_.

Assuming filtering is ON, running on one pair of fastq files creates:
* first alignment: `well_reference.bam`, e.g. _D01_psen1_1.bam_, in directory _bam_.
* alignment after filtering: `well_reference_filt.bam`, e.g. _D01_psen1_1_filt.bam_, in directory _bamfilt_.
* fastq files after filtering, back-converted from the filtered bam file: `well_reference_R1.fastq.gz` & `well_reference_R2.fastq.gz`, e.g. _D01_psen1_1_R1.fastq.gz_ & _D01_psen1_1_R2.fastq.gz_, in directory _filterfastq_.

Directories _bam_, _bamfilt_, _filterfastq_ are created in the folder where the config file is.

For filtering, alignFilterBack uses **filterBam.command**, which can also be used alone.

Reference should be in fasta format, all in lowercase, except PAM (if any) in uppercase. Multiple PAMs are supported.  
e.g. file _cog1Ex12.fa_, which is:
```
> cog1Ex12
tgtcaataaccgcacgttcaacagtgtgaaaggtaagatggagctgtattctcaccggatctgactgctggccaaAGGaaggatgttatatggttCCTgagagttgatgctgttgctcctagacgtgaactgcttctctgtgcccgtcaacagccctaacaacac
```
There are two PAMs: `AGG` and `CCT`.

## Dependencies
* `samtools`
* `picard`  
I had to install with homebrew, conda version had a conflict with Java.
* `bwa-mem`

### alignFilterBack.command

* `-c` path to .xlsx config file. It must have two columns (and only two): _well_ and _ref_. For example, it could look like:

    | well | ref |
    |:---|:---|
    |A01|geneX_amplicon2.fa|
    |B11|geneY_amplicon1.fa|

* `-r` path to directory containing fastq reads.

* `-a` path to directory containing the reference sequences (fasta format) listed in the config file.

* `-l` whether to filter the alignment or no. To turn filtering OFF, do not include the flag; to turn filtering OFF, include the flag.

If filtering is ON (flag `-l` is present), you may also include those parameters: 

* `-e` minimum PhrEd score for an alignment to be kept, e.g. `-e 40`.

* `-f` minimum read span for an alignment to be kept, e.g. `-f 100` means the alignment should span at least 100 bp of the reference. Looking at read span, rather than read length, makes it agnostic to possible insertions/deletions. For example, a 100-bp deletion in a 150-bp wild-type amplicon would create a 50-bp read, but it would still span the full 150 bp of the reference.

* `-s` maximum proportion Soft-clipped for the alignment to be kept, e.g. `-s 0.2` would remove any alignment which has more than 20% of its length soft-clipped.

* `-d` padding around double-strand break site in bp, e.g. `-d 20` would remove any alignment which does not cover at least double-strand break site ± 20 bp. For this filter to work, the PAM should be in capital letters in the fasta reference, e.g. `...agggattAGGacct...`. Multiple PAMs are supported.

* `-p` whether (`-p yes`) or not (`-p no`) to keep only primary alignments.

#### Some examples:

To align a bunch of fastq files without filtering:
```
alignFilterBack.command -c ./config.xlsx -r ./reads/ -a ./refseqs/
```

With every possible filtering:
```
alignFilterBack.command -c ./config.xlsx -r ./reads/ -a ./refseqs/ -l -e 40 -f 100 -s 0.2 -d 20 -p yes
```

### filterBam.command

* `-i` path to input BAM file.

* `-a` path to directory containing the reference sequences (fasta format). Only needs the reference which was used to create the BAM file.

* `-e` minimum PhrEd score for an alignment to be kept, e.g. `-e 40`.

* `-f` minimum read span for an alignment to be kept, e.g. `-f 100` means the alignment should span at least 100 bp of the reference. Looking at read span, rather than read length, makes it agnostic to possible insertions/deletions. For example, a 100-bp deletion in a 150-bp wild-type amplicon would create a 50-bp read, but it would still span the full 150 bp of the reference.

* `-s` maximum proportion Soft-clipped for the alignment to be kept, e.g. `-s 0.2` would remove any alignment which has more than 20% of its length soft-clipped.

* `-d` padding around double-strand break site in bp, e.g. `-d 20` would remove any alignment which does not cover at least double-strand break site ± 20 bp. For this filter to work, the PAM should be in capital letters in the fasta reference, e.g. `...agggattAGGacct...`. Multiple PAMs are supported.

* `-p` whether (`-p yes`) or not (`-p no`) to keep only primary alignments.

* `-o` path to BAM file to output.

#### Example

```
filterBam.command -i ./bam/D01_cog1Ex12.bam -a ./refseqs/ -e 40 -f 100 -s 0.2 -d 20 -p yes -o D01_cog1Ex12_filt.bam
```

## Other scripts

### slc45a2Crispresso2loop.command
For ZPRI project, this script is to run CRISPResso2 analysis on _slc45a2_ TAA>TGG prime editing MiSeq samples, with standard pegRNA. It simply runs on a loop in a folder of fastq, there are no flags.

CRISPResso2 command is:
```
CRISPResso --fastq_r1 "$FWD" --fastq_r2 "$RVS" --amplicon_seq GTACAGTCTGGTGTGGCTCATAAGCCCCATTTTGGGTTTTATCCTACAGCCCGTCATCGGCTCGGCGAGCGACTACTGTAGGTCGTCATAAGGCCGAAGGAGACCGTACATACTCTTACTGGGGATTCTGATGTTAGTGGGCATGACTTTATTTCTAAATGGAGATGCAGTCACAACAGGTGGGTGA --amplicon_name slc45a2TAA --prime_editing_pegRNA_spacer_seq gactactgtaggtcgtcata --prime_editing_pegRNA_extension_seq tctccttcggccccatgacgacctacagt --prime_editing_pegRNA_scaffold_seq gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtgggaccgagtcggtcc
```