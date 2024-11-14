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


### miscellaneous notes

Read with only a mutation at edge is counted as reference, see alleleToMutation.

classifyReads, expedit = NA will never call edit (perhaps obvious), only two possible categories become reference and mutated.


### about scaffold detection
v2: read class "scaffold" trumps all others, i.e. agnostic to edit or other mutations, if scaffold is there, it will be label read as "scaffold".

Potentially important for scaffold direction: always align to 5'–3' (genome direction) fasta reference sequence, so we can assume reads always represent the Fw genome.

After RTT, scaffold (standard) is CCTGGC... (reading in same direction as RT, i.e. from PBS to RTT).  
* If PE occurs in Forward strand: incorporates scaffold in Forward strand as GGACCG... (reading from the RTT-templated nucleotides)  
* If PE occurs in Reverse strand: incorporates scaffold in Reverse strand as GGACCG... (reading from the RTT-templated nucleotides), which appears in the Forward strand as CCTGGC... (reading from the RTT templated nucleotides) or ...CGGTCC reading 5'–3' genome.

Defined RHA position (`rhapos`) as last nucleotide templated from the RTT (maybe "ligation position" would be a better name).
* PE on Forward strand: expect scaffold to be `rhapos` + 1, substitution or insertion that _starts_ with G.
* PE on Reverse strand: expect scaffold to be `rhapos` - 1, substitution or insertion that _ends_ with C (as we are looking at Forward strand reads).  

! also impacts if want to look at `start` or `stop` for scaffold and whether sequence _starts_ or _ends_ with a given nucleotide. I set start to be always lower number (so more 5' on reference/genome) & stop to be always higher number.
* PE on Forward strand: scaffold mutation's `start` is closer to `rhapos`, and mutation _starts_ with G. 
* PE on Forward strand: scaffold mutation's `stop` is closer to `rhapos`, and mutation _ends_ with C.  

TODO: need to check again scaffdetectwin. Probably did not take this into account, especially if asymmetrical.
TODO: for filtering, did I take this into account? Looking at all positions between start and stop should be OK.
TODO: I am not sure about decision to add filtered-out reads as reference reads. Probably they should just be thrown out? i.e. they will be removed from mutated reads & total reads. Currently, by switching their labels to "ref", we removing them from mutated reads, but keeping them in total reads. ! those "manually edited to ref" reads will have ref & ali sequence as NA. While sequences directly called as ref before filterMutations (including those which were called as ref because the mutation was at the edge) have the actual ref & ali sequences.


### detectMHdel

MH stands for microhomology.

Runs on complete mutation table but currently only detects MH for deletions. If type is ref or sub or ins, returns all NA.

* `mut` mutation table, created by `callMutations`.

* `minMHlen` minimum length allowed for the MH, in bp. Default is 2 bp.

* `cutpos` position of cut site in reference sequence used for alignment (first nucleotide is #1). This position should point to exactly the nucleotide before the cut. For example, for reference sequence (PAM in uppercase) `attctagactNGGcattca`, we expect `cutpos=7` (position of `g`).

Typical case looks like:  
`MH - inner sequence - MH`  
and deleted sequence is
`MH - inner sequence` or `inner sequence - MH`

Columns it adds to mutation table:

* `MHside` side of the deleted MH. The MH left in the sequence after deletion is at the opposite side. Presumably, this is an arbitrary decision from the alignment algorithm as it could align the deletion shifted to the left or to the right. Value `none` with all other columns (below) being `NA` indicates a deletion that was analysed but that did not return any MH.

* `MHbp` length of the MH in bp.

* `leftMHstart` start position of the left MH.

* `leftMHstop` stop position of the left MH.

* `rightMHstart` start position of the right MH.

* `rightMHstop` stop position of the right MH.

* `innerSeq` inner sequence. This is the sequence flanked by the two MHs (cf. schematic above).

* `innerbp` length of inner sequence in bp.

* `leftFlapSeq` sequence of the 3'-flap on the left of the cut. I do not know for sure this is what the flap looked like as there have could been resection of the end before annealing of the MH. More precisely: sequence after the left MH and before the cut.

* `leftFlapbp` length of the left flap in bp.

* `rightFlapSeq` sequence of the 3'-flap on the right of the cut. I do not know for sure this is what the flap looked like as there have could been resection of the end before annealing of the MH. More precisely: sequence after the cut and before the right MH.

* `rightFlapbp` length of the left flap in bp.

All positions (`leftMHstart`, `leftMHstop`, `rightMHstart`, `rightMHstop`) refer to the reference sequence used for alignment (first nucleotide is #1). This is not always the same as the `ref` stored in the mutation table. In the case of an insertion, hyphens are added in the `ref` sequence for alignment, which shifts the positions.

`detectMHdel` runs on complete mutation table row by row and adds columns to it.

There are two, not mutually exclusive, exceptions to the typical case represented above:  

* **exception #1**: left MH is directly next to right MH, for example (real): 
    ```
    TAGGTCGTCATGGGGCCG
    TAGGTC---ATGGGGCCG
    ```
    (MH is `GTC`)  
    In this case, there is no inner sequence as the two MHs are directly one after the other. `innerSeq` is thus an empty string `""` of length (`innerbp`) 0. There are no left and right 3'-flap, so `leftFlapSeq`, `leftFlapbp`, `rightFlapSeq`, `rightFlapbp` are all NA.

* **exception #2**: cut is not between the two MHs, i.e. left MH stops after the cut or right MH starts before the cut. For example (real):
    ```
    GGCGAGCGACTACTGTAGGTCGTCATG*GGG*CC
    GGCG------------------TCATG*GGG*CC
                           ^
    ```
    (MH is `CG`, PAM is `*GGG*`, `cutpos` nucleotide is indicated with `^`)

    I cannot easily explain such case with MMEJ (how would MH anneal if both MHs are on the same side after cut?), but the fact that an MH is present seems to be a telltale sign. An option is that Cas9 did not cut at the expected site, but somewhere between the two MHs. I decided to still record those MH as they are potentially interesting. `innerSeq` and `innerbp` are still recorded, but we cannot tell anything about the 3'-flaps as it is unclear where the cut occurred, so `leftFlapSeq`, `leftFlapbp`, `rightFlapSeq`, `rightFlapbp` are all NA.

`detectMHdel` tries both left MH and right MH for each deletion. For each side, it looks for the longest MH, starting with `minMHlen`. The maximum MH size is set as the length of the deletion. If the deletion has MH on both sides (I think this is extremely rare), it compares both and keeps the longest one. If both are the same size (I do not know if this would ever happen), it will all record the left one.

### version history

* v1

* v2
Detection of scaffold incorporation, currently just looks at substitutions or insertions that starts with G around end of RTT (end of RHA). Can control detection with `rhapos` and `scaffdetectwin`.

* v3  
Detection of microhomology for deletions, see `detectMHdel`.