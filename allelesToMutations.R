#####################################################
# ~ ZPRI: CRISPResso2 alleles to mutations ~
#
# utility to convert Allele frequencies from CRISPResso2 to table of unique mutations
# 
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# v0

# v1
# records reference & mutated sequence aligned
# this will make it easier to analyse repair


# packages ----------------------------------------------------------------

library(data.table)


# function allelesToMutations() -------------------------------------------

# overarching function
# loops through table of alleles
# and convert it into table of read-by-read mutations

# alpath: path to one Alleles_frequency_table, created by CRISPResso2
allelesToMutations <- function(alpath) {
  
  ### import
  # als for alleles table
  als <- as.data.frame(fread(alpath))
  
  ### loop through it row by row
  # each row, call alleleToMutation()
  # mutl is a list where each element is all mutations for one read (alignment)
  mutl <- lapply(1:nrow(als), function(rowi) {
    
    # for that row,
    # read REFERENCE sequence
    ref <- als[rowi, 'Reference_Sequence']
    # split it directly
    ref <- strsplit(ref, split='')[[1]]
    
    # read ALIGNED sequence
    ali <- als[rowi, 'Aligned_Sequence']
    # split it directly
    ali <- strsplit(ali, split='')[[1]]
    
    # read number of reads
    nreads <- als[rowi, '#Reads']
    
    return( alleleToMutation(ref=ref,
                             ali=ali,
                             nreads=nreads,
                             rowi=rowi) )
    
  })
  
  # put everything in one big dataframe
  mutdf <- data.frame(do.call(rbind, mutl))
  
  # fix some of the types
  mutdf$start <- as.integer(mutdf$start)
  mutdf$stop <- as.integer(mutdf$stop)
  mutdf$bp <- as.integer(mutdf$bp)
  
  ### calculate total number of reads for this sample, i.e. coverage
  # and add it as column
  splcov <- sum(als$`#Reads`)
  cat('\t \t \t \t >>> Total number of reads:', splcov, '\n')
  mutdf <- mutdf %>%
    mutate(splcov=splcov)
  
  # return
  return(mutdf)
  
}


# function alleleToMutation() ---------------------------------------------

# for one alignment (one row of allele table), record mutations

# ref = one aligned reference, split
# ali = one aligned sequence, split
alleleToMutation <- function(ref,
                             ali,
                             nreads,
                             rowi) {
  
  ### REFERENCE
  refdf <- recordRef(ref=ref,
                     ali=ali)
  
  ### DELETIONS
  deldf <- recordDel(ref=ref,
                     ali=ali)
  
  ### INSERTIONS
  insdf <- recordIns(ref=ref,
                     ali=ali)

  ### SUBSTITUTIONS
  subdf <- recordSub(ref=ref,
                     ali=ali)
  
  ### EXCEPTION: once in a while, there is a mutation at the edge of the read
  # so recordDel/Ins/Sub skips it
  # ! we will count the read as reference
  if(length(refdf)==0 & length(deldf)==0 & length(insdf)==0 & length(subdf)==0) {
    # cat('\t \t \t \t \t >>> Only a mutation at the edge of the read found, counting this read as reference.\n')
    # 18/09/2024; skipping some verbose, could make a total count for each sample, but for each read it too much
    refdf <- c(type='ref',
               start=NA,
               stop=NA,
               bp=NA,
               refseq=NA,
               altseq=NA,
               ref=paste0(ref, collapse=''),
               ali=paste0(ref, collapse=''))
  }
  
  ### pool & return
  # pool
  mutdf <- as.data.frame(rbind(refdf, deldf, insdf, subdf))
  # sometimes messes up the column names, add them back
  colnames(mutdf) <- c('type', 'start', 'stop', 'bp', 'refseq', 'altseq', 'ref', 'ali')
  # delete any row names
  row.names(mutdf) <- NULL
  # if more than one read like this, multiply the rows
  mutdfm <- do.call('rbind', replicate(nreads, mutdf, simplify=FALSE))
  
  # add read IDs
  # which is original row.1, 2, 3, ...
  mutdfm$rid <- sprintf('r%i.%i', rowi, rep(1:nreads, each=nrow(mutdf)))
  # 16/10/2024: originally was 1.1, 1.2, etc.
  # but read.csv imports it as numeric, not character
  # which then creates issues because e.g. 1.100 becomes 1.1, so then multiple 1.1 reads
  # write.csv correctly puts "" around rid, the issue is from read.csv
  # so only solution without every time thinking about it is to add a letter here
  
  return( mutdfm )
  
}


# function recordRef() ----------------------------------------------------

# ref = reference sequence (aligned), split
# ali = aligned sequence, split

# if aligned sequence is the same as reference,
# we record special "ref" mutation

recordRef <- function(ref,
                      ali) {
  
  # put back the sequences together
  refc <- paste0(ref, collapse='') # c for collapsed
  alic <- paste0(ali, collapse='')
  
  if( identical(refc , alic) ) {
    # cat('\t \t \t \t >>> Detected reference allele.\n')
    # 18/09/2024; skipping some verbose, could make a total count for each sample, but for each read it too much
    return(c(type='ref',
             start=NA,
             stop=NA,
             bp=NA,
             refseq=NA,
             altseq=NA,
             ref=refc,
             ali=alic))
  }
  # else, do nothing
}


# function recordDel() ----------------------------------------------------

# ref = reference sequence (aligned), split
# ali = aligned sequence, split

recordDel <- function(ref,
                      ali) {
  
  # find all deleted positions
  # these are positions where aligned sequence is a gap -
  delpos <- which(ali=='-')
  
  ### if no deletion, we return an empty dataframe in the correct format
  if(length(delpos)==0) {
    return(c(type=character(),
             start=integer(),
             stop=integer(),
             bp=integer(),
             refseq=character(),
             altseq=character()))
  } # return would stop function here, so consider below is `else`
  
  # each deletion starts when there is a break in the sequence
  # e.g. 9, 10, 11, 22, 23, 24
  # means one deletion starts at #9 and one deletion starts at #22
  # we can test this with diff (which is = n+1 minus n)
  # e.g. 4, 5, 6, 10, 11, 12 will give 1, 1, 4, 1, 1
  # any position that is not 1 is a break in the sequence
  # as can see, need to do +1 to the diff positions
  # also record 1 as first gap is certain to be start of first deletion
  # _i for index of each deletion start
  # ! indices of delpos, not sequence positions
  delpos_i <- c(1, (which(diff(delpos)!=1)+1))
  # cat('\t \t \t \t >>> Detected', length(delpos_i), 'deletion(s). \n')
  # 18/09/2024; skipping some verbose, could make a total count for each sample, but for each read it too much

  # now we split delpos in a list, where each element is a continuous sequence of positions
  # e.g. deletion1 is positions 1, 2, 3; deletion2 is positions 11, 12; deletion3 is position 22; etc.
  # _s for sequences of deletions
  delpos_s <- splitAt(delpos, delpos_i)
  
  # now we loop through the sequences and record each deletion
  # seqpos is sequence of positions
  # delrs is a list of "recorded deletions"
  delrs <- lapply(delpos_s, function(seqpos) {
    
    # ! is this deletion that at the edge of the read?
    # e.g. aligned sequence is AGGGACTAA--- or ---AGGGACTAA
    # in this case we probably do not have the entire deletion
    # so if we counted naively, we would count 3 bp deletion above
    # while it is likely that the deletion starts upstream or stops downstream so we do not see it entirely
    # we can tell if position #1 or position #`length of alignment` is in the deleted positions
    # in this case, we should not record this deletion
    if(1 %in% seqpos | length(ali) %in% seqpos) {
      # cat('\t \t \t \t >>> Deletion at the edge of the read, skipping record as it may be truncated. \n')
      # 18/09/2024; skipping some verbose, could make a total count for each sample, but for each read it too much
      return()
    }
    
    # prepare reference sequences
    # if length bp is 1, just that nucleotide
    # if length bp is more than 1, paste & collapse
    # (paste & collapse when length is 1 creates a glitch)
    if(length(seqpos)==1) {
      refseq <- ref[seqpos]
    } else if(length(seqpos)>1) {
      refseq <- paste0(ref[seqpos], collapse='')
    }
    
    # ! we need to correct the start/stop positions
    # see fixpos() for details
    start <- fixpos(ref=ref, pos=seqpos[1])
    stop <- fixpos(ref=ref, pos=seqpos[length(seqpos)])
    
    return(c(type='del',
             start=start,
             stop=stop,
             bp=length(seqpos),
             refseq=refseq,
             altseq=NA,
             ref=paste0(ref, collapse=''), # return reference & aligned sequence, pasted back together
             ali=paste0(ali, collapse='')))
  })
  # note, altseq is undetermined as the deleted sequence is not present anymore
  
  deldf <- data.frame(do.call(rbind, delrs))
  
  return(deldf)

}


# function recordIns() ----------------------------------------------------

# will be fairly similar to recordDel

# ref = reference sequence (aligned), split
# ali = aligned sequence, split

recordIns <- function(ref,
                      ali) {
  
  # find all inserted positions
  # these are positions where reference sequence is a gap -
  inspos <- which(ref=='-')
  
  ### if no insertion, we return an empty dataframe in the correct format
  if(length(inspos)==0) {
    return(c(type=character(),
             start=integer(),
             stop=integer(),
             bp=integer(),
             refseq=character(),
             altseq=character()))
  } # return would stop function here, so consider below is `else`
  
  # each insertion starts when there is a break in the sequence
  # e.g. 9, 10, 11, 22, 23, 24
  # means one insertion starts at #9 and one insertion starts at #22
  # we can test this with diff (which is = n+1 minus n)
  # e.g. 4, 5, 6, 10, 11, 12 will give 1, 1, 4, 1, 1
  # any position that is not 1 is a break in the sequence
  # as can see, need to do +1 to the diff positions
  # also record 1 as first gap is certain to be start of first insertion
  # _i for index of each insertion start
  # ! indices of inspos, not sequence positions
  inspos_i <- c(1, (which(diff(inspos)!=1)+1))
  # cat('\t \t \t \t >>> Detected', length(inspos_i), 'insertion(s). \n')
  # 18/09/2024; skipping some verbose, could make a total count for each sample, but for each read it too much
  
  # now we split inspos in a list, where each element is a continuous sequence of positions
  # e.g. insertion1 is positions 1, 2, 3; insertion2 is positions 11, 12; insertion3 is position 22; etc.
  # _s for sequences of insertions
  inspos_s <- splitAt(inspos, inspos_i)
  # now we loop through the sequences and record each insertion
  # seqpos is sequence of positions
  # insrs is a list "recorded deletions"
  insrs <- lapply(inspos_s, function(seqpos) {
    
    # ! is this insertion at the edge of the read?
    # e.g. reference sequence is AGGGACTAA--- or ---AGGGACTAA
    # in this case we probably do not have the entire deletion
    # so if we counted naively, we would count 3 bp insertion above
    # while it is likely that the insertion starts upstream or stops downstream so we do not see it entirely
    # we can tell if position #1 or position #`length of alignment` is in the inserted positions
    # in this case, we should not record this insertion
    if(1 %in% seqpos | length(ref) %in% seqpos) {
      # cat('\t \t \t \t >>> Insertion at the edge of the read, skipping record as it may be truncated. \n')
      # 18/09/2024; skipping some verbose, could make a total count for each sample, but for each read it too much
      return()
    }
    
    # prepare alternative sequences
    # if length bp is 1, just that nucleotide
    # if length bp is more than 1, paste & collapse
    # (paste & collapse when length is 1 creates a glitch)
    if(length(seqpos)==1) {
      altseq <- ali[seqpos]
    } else if(length(seqpos)>1) {
      altseq <- paste0(ali[seqpos], collapse='')
    }
    
    # ! we need to correct the start/stop positions
    # see fixpos() for details
    start <- fixpos(ref=ref, pos=seqpos[1])
    stop <- fixpos(ref=ref, pos=seqpos[length(seqpos)])
    
    return(c(type='ins',
             start=start,
             stop=stop,
             bp=length(seqpos),
             refseq=NA,
             altseq=altseq,
             ref=paste0(ref, collapse=''), # return reference & aligned sequence, pasted back together
             ali=paste0(ali, collapse='')))
  })
  # note, refseq is undetermined as the inserted sequence was not present in the reference
  
  insdf <- data.frame(do.call(rbind, insrs))
  
  return(insdf)
  
}



# function recordSub() ----------------------------------------------------

# ref = aligned reference, split
# ali = aligned sequence, split

recordSub <- function(ref,
                      ali) {
  
  # check again ref & ali have same length bp
  if(length(ref) != length(ali))
    stop('\t \t \t \t >>> Error: reference and aligned sequence do not have the same length (bp).\n')
  
  # for substitutions,
  # we simply compare ref vs ali
  # but ignore where gaps
  # i.e. substitution is only when two nucleotides do not match
  subpos <- which(ref!=ali) # these are all positions where aligned does not match reference
  # from this, remove where there are gaps in either the aligned or the reference sequence
  subpos <- subpos[! subpos %in% sort(unique(c(which(ref=='-'), which(ali=='-'))))] # i.e. keep only subpos which are *not* positions with gap
  
  ### if no substitution, we return an empty dataframe in the correct format
  if(length(subpos)==0) {
    return(c(type=character(),
             start=integer(),
             stop=integer(),
             bp=integer(),
             refseq=character(),
             altseq=character()))
  } # return would stop function here, so consider below is `else`
  
  # ! make substitutions as long as possible,
  # i.e. we should not record one substitution 1bp and another one next to it 1bp, but rather one substitution of 2bp
  # we can tell where the transitions happen (i.e. each start of a new substitution) using diff
  # see recordDel for logic
  # _i for index of each substitution start
  # ! indices of subpos, not sequence positions
  subpos_i <- c(1, (which(diff(subpos)!=1)+1))
  # cat('\t \t \t \t >>> Detected', length(subpos_i), 'substitution(s). \n')
  # 18/09/2024; skipping some verbose, could make a total count for each sample, but for each read it too much
  # now we split subpos in a list, where each element is a continuous sequence of positions
  # e.g. substitution1 is positions 1, 2, 3; substitution2 is positions 11, 12; substitution3 is position 22; etc.
  # _s for sequences of substitutions
  subpos_s <- splitAt(subpos, subpos_i)
  # now we loop through the sequences and record each substitution
  # seqpos is sequence of positions
  subrs <- lapply(subpos_s, function(seqpos) {
    
    # ! is this substitution at the edge of the read?
    # e.g. reference AGGGACTAACCC vs aligned AGGGACTAAGGG
    # in this case we probably do not have the entire substitution
    # so if we counted naively, we would count 3 bp insertion above
    # while it is likely that the insertion starts upstream or stops downstream so we do not see it entirely
    # we can tell if position #1 or position #`length of alignment` is in the inserted positions
    # in this case, we should not record this insertion
    if(1 %in% seqpos | length(ref) %in% seqpos) {
      # cat('\t \t \t \t >>> Substitution at the edge of the read, skipping record as it may be truncated. \n')
      # 18/09/2024; skipping some verbose, could make a total count for each sample, but for each read it too much
      return()
    }
    
    # prepare reference/alternative sequences
    # if length bp is 1, just that nucleotide
    # if length bp is more than 1, paste & collapse
    # (paste & collapse when length is 1 creates a glitch)
    if(length(seqpos)==1) {
      refseq <- ref[seqpos]
      altseq <- ali[seqpos]
    } else if(length(seqpos)>1) {
      refseq <- paste0(ref[seqpos], collapse='')
      altseq <- paste0(ali[seqpos], collapse='')
    }
    
    # ! we need to correct the start/stop positions
    # see fixpos() for details
    start <- fixpos(ref=ref, pos=seqpos[1])
    stop <- fixpos(ref=ref, pos=seqpos[length(seqpos)])
    
    return(c(type='sub',
             start=start,
             stop=stop,
             bp=length(seqpos),
             refseq=refseq,
             altseq=altseq,
             ref=paste0(ref, collapse=''), # return reference & aligned sequence, pasted back together
             ali=paste0(ali, collapse='')))
  })
  
  subdf <- data.frame(do.call(rbind, subrs))
  
  return(subdf)
  
}


# function fixpos() -------------------------------------------------------

# we need final positions to always be in reference to absolute reference (i.e. genome)
# but if there are insertions, aligned reference will get gaps ---, increasing its length bp vs. absolute reference
# so when we record a position, we should account for this
# which is simply subtracting the number of gaps (---) in aligned reference found before that position
fixpos <- function(ref,
                   pos) {
  # note, works when no gap because it gives integer(0) which is length 0
  # note, e.g. C C A - T and pos is -, we *should* correct it so that insertion start marks nucleotide just before
  # it currently does correct it, in e.g. it would give 4 -- corrected to -- 3
  return( pos - length(which(ref[1:pos]=='-')) )
}



# function splitAt() ------------------------------------------------------
# from https://stackoverflow.com/questions/16357962/r-split-numeric-vector-at-position

splitAt <- function(x,
                    pos) {
  return( unname(split(x, cumsum(seq_along(x) %in% pos))) )
}
  
