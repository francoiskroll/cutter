#####################################################
# ~ cutter: simulate random deletions ~
#
# see lots of notes in fitDelLengths_notes.R
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################
library(ggplot2)
library(dplyr)
library(MASS)
# v0

# expects mut, mutation table


# simulateDel -------------------------------------------------------------

simulateDel <- function(mut,
                        nreads=1000,
                        mincov=100,
                        cutpos,
                        cutDelbp=3,
                        awayfromCut=4) {
  
  ### fit distribution of deletion lengths with a log-normal function
  logparams <- fitDelLengths(mut,
                             mincov=mincov)
  
  ### now simulate `nreads` with each a deletion
  # step below is surprisingly slow...
  simdelL <- lapply(1:nreads, function(i) {
    cat('\t \t \t \t simulating read', i, 'of', nreads, '\n')
    simulateDel_one(logparams=logparams,
                    cutpos=cutpos,
                    cutDelbp=3,
                    awayfromCut=4)
  })
  simdel <- do.call(rbind, simdelL)
  
  # correct some column types
  simdel$start <- as.integer(simdel$start)
  simdel$stop <- as.integer(simdel$stop)
  simdel$bp <- as.integer(simdel$bp)
  
  ### add rid
  # normally done by alleleToMutation
  # here, will simply be r1.1, r2.1, r3.1, etc.
  # until e.g. r1000.1 if nreads=1000
  simdel$rid <- sprintf('r%i.1', 1:nrow(simdel))
  
  ### add coverage
  # copying from allelesToMutations so we get the same columns
  # calculate total number of reads for this sample, i.e. coverage
  # and add it as column
  # here, we know this is exactly `nreads`
  cat('\t \t \t \t >>> Total number of reads:', nreads, '\n')
  simdel <- simdel %>%
    mutate(splcov=nreads)
  
  ### run mutation table through filterMutations
  # with all filters off
  # just to get the same columns
  simdel <- filterMutations(muttb=simdel,
                            minnreads=NA,
                            cutpos=NA,
                            rhapos=NA,
                            rhadist=NA,
                            controltb=NA,
                            callSubs=TRUE)
  
  ### add rundate, well, locus
  # these we know are present in mut
  simdel <- simdel %>%
    mutate(well='Z99', .before=1) %>%
    mutate(locus='simulated', .before=1) %>%
    mutate(rundate=999999, .before=1)
  
  ### also create sample column
  # which is rundate_well_locus
  simdel <- simdel %>%
    mutate(sample=paste(rundate, well, locus, sep='_'), .before=1)
  
  ### add last columns as NA
  # find the ones we are missing from mut
  # these were added from meta
  cols2add <- colnames(mut) [! colnames(mut) %in% colnames(simdel)]
  
  for(colnm in cols2add) {
    simdel <- simdel %>%
      mutate(tmp=NA, .before='well')
    # now put the correct column name
    colnames(simdel)[which(colnames(simdel)=='tmp')] <- colnm
  }
  
  ### add simulated reads to mut
  # first put the columns of simdel in the same order as in mut
  
  # check columns are identical (before looking at their positions)
  if(!identical( sort(colnames(mut)) , sort(colnames(simdel))))
    stop('\t \t \t \t >>> Error simulateDel: all columns of mut are not in simdel, or vice-versa.\n')
  
  # now put the columns of simdel in the same order
  simdel <- simdel[, match(colnames(simdel), colnames(mut)) ]
  # check correct
  if(!identical(colnames(simdel), colnames(mut)))
    stop('\t \t \t \t >>> Error simulateDel: columns of simdel are not in the same order as columns of mut.\n')
  
  # can now rbind
  mutsim <- rbind(mut, simdel)
  # and we return result
  return(mutsim)
  
}



# simulateDel_one ---------------------------------------------------------

simulateDel_one <- function(logparams,
                            cutpos,
                            cutDelbp,
                            awayfromCut) {
  
  # if delOK did not get switch to TRUE below
  # then simulate a new deletion
  delOK <- FALSE
  
  while(!delOK) {
    
    ### from the fitted log-normal distribution, we draw one deletion length
    # and round it to closest integer
    # ! if this integer is 0, we draw again
    len <- 0
    while(len==0) {
      len <- round(rlnorm(n=1,
                          meanlog=logparams[1],
                          sdlog=logparams[2]))
    }
    
    ### recover original reference sequence
    # i.e. the reference sequence used for alignment
    # (in alignment, the reference sequence can get --- to indicate insertions)
    # get all the 'true' original reference sequences:
    # which we do by just removing the --- from the aligned reference seauences
    orefs <- gsub('-', '', mut$ref) # "original references"
    # check they are all the same, after we remove NA
    orefs <- orefs[!is.na(orefs)]
    
    if(length(unique(orefs))!=1)
      stop('\t \t \t \t >>> Error simulateDel_one: issue when computing overall reference sequence.\n')
    # we have the original reference sequence
    oref <- unique(orefs)
    
    ### 
    # deletion must remove the 'cutpos' nucleotide
    # which is PAM minus 4
    # cf. https://elifesciences.org/articles/59683#fig2s1, it is indeed the most frequently deleted nucleotide
    # except if 1 or 2 nt (flexible), then we can leave a bit of margin
    # but we still want the deletion to be very close, say ± 4 bp (flexible too)
    
    if(len < cutDelbp) {
      # if deletion is very small, 
      # we do not necessarily have to remove `cutpos` nucleotide
      # deletion most shifted to the left that is acceptable is
      # e.g. if cutpos is #87 and length is 2 bp
      # then most shifted to the left is #83, #84
      # because 4 bp (`awayfromCut`) from #84 includes #87, which is cutpos
      delstart <- cutpos - awayfromCut - len + 1 + 1 # a bit by trial and error to get what I want...
      delposleft <- delstart : (delstart+len-1)
      
      # then we just shift the deletion to the right 1 bp at a time
      # we do it until the left-most deleted position is `awayfromCut` away from cutpos
      # note, window of possible deletions is not symmetrical around `cutpos`
      # this is on purpose, because cutpos labels the nucleotide left of the cut
      # it *is* symmetrical if you consider the exact cut position (i.e. where the DNA chain is broken)
      shiftNtimes <- (cutpos + awayfromCut) - delstart
      
      # starting at 0 allows to get delposleft as first element of the list
      # then we shift one by one to the right
      delposL <- lapply(0:shiftNtimes, function(le) {
        delposleft + le
      })
      
    } else {
      # if deletion is fairly big, it must delete `cutpos` position
      # possible deletion most shifted to the left (so when last deleted nucleotide is cutpos) is
      delposleft <- (cutpos - len + 1) : cutpos # vector of all the deleted positions
      
      # we just shift the deletion to the right 1 bp at a time
      # `len` times
      
      # starting at 0 allows to get delposleft as first element of the list
      # then we shift one by one to the right
      delposL <- lapply(0:len, function(le) {
        delposleft + le
      })
      # so this is our set of possible deletions of length `len`
      # all delete the `cutpos` nucleotide
    }
    
    ### draw one deletion from the set
    drawi <- sample(1:length(delposL), 1)
    delpos <- delposL[[drawi]]
    
    # we now have our simulated deletion
    
    ### ! check the positions !
    # ! the deleted positions may be impossible
    # that is, we may be deleting negative positions (below the start of the read)
    # or positions beyond the end of the read
    # if any of these occur, repeat the procedure
    
    # if smallest deleted position is negative or 0
    # or if largest deleted position is past the end of the read
    # this is not OK and we should repeat the procedure
    if( min(delpos)<1 | max(delpos)>nchar(oref) ) {
      cat('\t \t \t \t impossible deletion, repeating procedure.\n')
      delOK <- FALSE
    } else {
      delOK <- TRUE
    }
  }
  
  # if deletion is OK, finish
  if(delOK) {
    ### create the aligned sequence
    # we split the reference sequence
    refsp <- strsplit(oref, split='')[[1]]
    # then swap the deleted nucleotides to -
    alisp <- refsp
    alisp[delpos] <- '-'
    # the reference sequence stays the same
    
    print('ref')
    print(refsp)
    
    print('ali')
    print(alisp)
    
    print(delpos)
    
    # check they are the same size just to be safe
    if(length(refsp)!=length(alisp))
      stop('\t \t \t \t >>> Error simulateDel_one: reference sequence and simulated aligned sequence do not have the same length.\n')
    
    # now we can give the reference & aligned sequences directly to recordDel from allelesToMutations.R
    # this will give the row of mut
    # this is what we return
    # overarching function can run simulateDel_one many times to create data from the simulated sample
    return( recordDel(ref=refsp,
                      ali=alisp) )
  }
  
}



# fitDelLengths -----------------------------------------------------------

# takes a big dataset as input
# gives as many samples & loci as possible
# will fit a log-normal distribution of deletion lengths
# reproduces key steps of fitDelLengths_notes.R

# will hard-code a default fit at some point
# so this function is only if want to fit it again, e.g. different model or experimental condition

# e.g.
# fitDelLengths(mut=read.csv('~/Dropbox/cutter/mutcallsTest.csv'),
#               mincov=100)

fitDelLengths <- function(mut,
                          mincov=100) {
  
  #### TODO
  # I am not confident about the approach below anymore
  # shouldn't we be taking frequencies as out of total deleted reads in each sample?
  # take a sample which is 10% mutated, and 90% of mutations are a 11-bp deletion
  # then, frequency of 11-bp deletion is 9%, so the 11-bp deletion is rare
  # but really it should mean the 11-bp deletion is very frequent, as it occured 9 times out of 10
  
  
  ### get all the deletions
  del <- mut %>%
    filter(type=='del')
  
  ### normalise all coverage to total 1000 reads
  # co1k for counts out of 1000x coverage
  del <- del %>%
    mutate(co1k=round(freq*1000))
  
  # now, we want unique deletions within each sample
  # and we count it "co1k" times
  delu <- del %>%
    distinct(sample, mutid, .keep_all=TRUE)
  # we loop through these deletions,
  # and each time repeat their lengths "co1k" times
  # essentially we re-create the dataset as if every sample had exactly coverage 1000x
  lens <- unlist(sapply(1:nrow(delu), function(ro) {
    # from this row, get the deletion length & the normalised count
    # and return length, `normalised count` times
    # e.g. we get 8 bp deletion, and normalised count is 4
    # we return 8 8 8 8
    rep( delu[ro, 'bp'] , delu[ro, 'co1k'])
  }))
  # unlist: we can throw all together, no need to keep track of samples etc.
  # we now have all the deletion lengths, as if every sample had the same coverage
  
  ### fit log-normal distribution
  flog <- fitdistr(lens, 'lognormal', lower=1)
  # lower: lower bound, which is 1 (bp) here
  
  ### plot fit on density plot
  gglog <- ggplot(del, aes(x=bp)) + 
    geom_density() +
    stat_function(
      fun=dlnorm,
      args=list(meanlog=flog$estimate[1],
                sdlog=flog$estimate[2]),
      colour='red'
    )
  print(gglog)
  
  ### return the log-normal parameter
  # vector meanlog, sdlog
  return(c(flog$estimate[1],
           flog$estimate[2]))
}


# freqBins ----------------------------------------------------------------

# small utilities function
# to turn some data (e.g. deletion lengths) into frequency within bins
# v: a vector of numerics
# every: new bin every x, i.e. controls the breaks

freqBins <- function(v,
                     every,
                     lower=0,
                     last=NA) {
  
  if(is.na(last)) {
    # %% is modulo: e.g. 57 %% 5 gives 2 because 55 can be divided by 5, then remains 2
    # so we can tell how much we need to add to make it divisible based on the remainder
    remainder <- max(v, na.rm=TRUE) %% every
    
    if(remainder != 0) {
      toadd <- every - (max(v, na.rm=TRUE) %% every) + lower
    } else {
      toadd <- 0
    }
    # say every is 5
    # if max(v) was
    # 51, gives 4 > correct
    # 52, gives 3 > correct
    # 53, gives 2 > correct
    # 54, gives 1 > correct
    # 55, gives 5 > incorrect, we should add nothing
    # just catch this case
    lastbreak <- max(v, na.rm=TRUE) + toadd
  } else {
    lastbreak <- last
  }
  
  breaks <- seq(from=lower, to=lastbreak, by=every)
  
  his <- hist(v, breaks=breaks, plot=FALSE)
  hisd <- data.frame(upbound=his$breaks[2:length(his$breaks)],
                     counts=his$counts)
  # add frequencies
  hisd$freq <- hisd$counts / sum(hisd$counts)
  # dataframe
  # breaks: last value of each bin
  # counts: number of observations in this bin
  
  return(hisd)
  
}