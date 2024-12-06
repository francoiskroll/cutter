refff <- 'GTACAGTCTGGTGTGGCTCATAAGCCCCATTTTGGGTTTTATCCTACAGCCCGTCATCGGCTCGGCGAGCGACTACTGTAGGTCGTCATGGGGCCGAAGGAGACCGTACATACTCTTACTGGGGATTCTGATGTTAGTGGGCATGACTTTATTTCTAAATGGAGATGCAGTCACAACAGGTGGGTGA'
refff <- strsplit(refff, '')[[1]]
refff[93]


# create composite column srid which is sample_rid so we have a unique identifier for each read across samples
mut <- mut %>%
  mutate(srid=paste(sample, rid, sep='_'), .after='rid')


stringdist(newseq, searchSeq, method='lcs')

ta <- drop(attr(adist(newseq, searchSeq, counts=TRUE), "trafos"))
stri_sub(searchSeq, stri_locate_all_regex(ta, "M+")[[1]])

## get all forward substrings of 'b'
sb <- stri_sub(newseq, 1, 1:nchar(newseq))
## extract them from 'a' if they exist
sstr <- na.omit(stri_extract_all_coll(a, sb, simplify=TRUE))
## match the longest one
sstr[which.max(nchar(sstr))]
# [1] "hel"


library(Biostrings)

out1 <- pmatchPattern(newseq, searchSeq, maxlength.out=1L)

out2 <- pmatchPattern(searchSeq, newseq, maxlength.out=1L)

# ! output is different depending on order
# newseq, searchSeq: GCCGAAGG
# searchSeq, newseq: CGAAGGAG
# so both outputs are very good (8 nt each) but I do not understand logic behind

as.character(out1[[1]])
as.character(out2[[1]])



# -------------------------------------------------------------------------

###### stuff to sort


mut %>%
  filter(type=='ins') %>%
  filter(mutid=='ins_89_89_9_NA_AGCCGAAGG') %>%
  filter(sample=='240902_F05_slc45a2WT')

which(mut$sample=='240902_F05_slc45a2WT' & mut$mutid=='ins_89_89_9_NA_AGCCGAAGG')

mut[9890,]

insro <- mut[9890,]

cutpos <- 87
cutdist <- 12
allowStretchMatches <- 5
extendNewlySeq <- 3
searchWindowStarts <- 20


# TODO
# ! will need to take into account extendNewlySeq when deciding on longest match
# here, the minimum is 3, even if newly synthesised was completely random

# -------------------------------------------------------------------------

### test
mut <- read.csv('~/Dropbox/miseqUtils/mutcallsTest.csv')

library(tictoc)
tic()
# mut <- mut[sample(nrow(mut), 1000), ]
mutin <- detectTemplatedIns(mut=mut,
                            cutpos=87,
                            cutdist=12, 
                            minLCSbp=1,
                            allowStretchMatches=5,
                            extendNewlySeq=3,
                            searchWindowStarts=20)
toc()
# 1k rows (226 unique alignments): 127 sec
# so 1.8 sec per alignment
# full dataset has 7k alignments, so 3.5 hours
# no, solution is to take only ins rows,
# then do mutuni etc.
# only insertions: 33 alignments, 15 sec
# so still 2.2 per alignment, but fewer


### TODO
# I think issue with minLCSbp setting; hardcoded?
# no actually looks fine now, it returns many which are 3 bp

# is this the right way to output the results? Because we have row for each mut, including del & sub
# I think correction should be: do not check mutation type initially, we will just run on every alignment
# then add all the columns to mut
# then overwrite in rows which are ref or del or sub, all NA

# check search sequence of 164 bp?

mutin %>%
  filter(!is.na(lcbp)) %>%
  View(.)

mutin %>%
  filter(lcbp<5) %>%
  View(.)



tmp <- mut[sample(nrow(mut), 1000), ]
tic()
mutins <- detectTemplatedIns(mut=tmp,
                             cutpos=87,
                             cutdist=12, 
                             minLCSbp=1,
                             allowStretchMatches=5,
                             extendNewlySeq=3,
                             searchWindowStarts=20)
toc()
# 1k rows: from 36 sec to 6-7 sec
# so divided by 6
# complete dataset should take 14 min, which is not fast but not too bad either



tmp <- mut[sample(nrow(mut), 500), ]
tmp <- mut[9890:11000,]
tic()
mutins <- detectTemplatedIns(mut=tmp,
                   cutpos=87,
                   cutdist=12, 
                   minLCSbp=1,
                   allowStretchMatches=5,
                   extendNewlySeq=3,
                   searchWindowStarts=20)
toc()
# very slow!
# 1k: 36 sec
# 2k: 64 sec
# all dataset (151k rows) would take 90 min

### now, how to add to mut

# to take an example
mut1 <- mut %>%
  filter(ref==mut[9890, 'ref'] & ali==mut[9890, 'ali'])
nrow(mut1) # 1754 rows

mutins1 <- mutins %>%
  filter(ref==mut[9890, 'ref'] & ali==mut[9890, 'ali'])
nrow(mutins1) # 1 row

# from mutins, only keep ref, ali, and the new columns
# then we can join by matching ref & ali
# this will automatically copy all the new columns when alignment is the same


# TODO I think issue with minLCSbp because only returns quite long matches but it should not

### approach to speed up would be to only look at unique alignments
# complete dataset is 151k rows
mutuni <- mut %>%
  distinct(ref, ali, .keep_all=TRUE)
# 7,7k
# which brings it to 4 min, I guess that is fine

mut %>%
  filter(type=='ins') %>%
  filter(mutid=='ins_89_89_9_NA_AGCCGAAGG') %>%
  filter(sample=='240902_F05_slc45a2WT')

which(mut$sample=='240902_F05_slc45a2WT' & mut$mutid=='ins_89_89_9_NA_AGCCGAAGG')

mut[9890,]

insro <- mut[9890,]

cutpos <- 87
cutdist <- 12
allowStretchMatches <- 5
extendNewlySeq <- 3
searchWindowStarts <- 20



tmp[tmp$type %in% c('ref', 'sub', 'del'),
       c('lcbp', 'lcseq', 'lcdir', 'lcStart', 'lcStop')] <- list(999, NA, NA, NA, 888)




# make cartoon, some checks -----------------------------------------------

mut <- read.csv('~/Dropbox/cutter/mutcallsTest.csv')

# find example reverse match
# top 1 insertion at slc45a2, see studyRepair.ai
# found one #7127

oref <- gsub('-', '', mut$ref[7127])

detectTemplatedIns_one(alrow=mut[7127,],
                       oref=oref,
                       cutpos=87,
                       cutdist=12,
                       minLCSbp=1,
                       allowStretchMatches=5,
                       extendNewlySeq=3,
                       searchWindowStarts=20)