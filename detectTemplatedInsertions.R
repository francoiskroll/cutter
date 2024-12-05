#####################################################
# ~ miseqUtils/CUTTER: detect templated insertions ~
#
# (which is characteristic of polq activity)
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# v0

# expects mutation table written by callMutations
library(Biostrings)
library(dplyr)
library(stringi)

# example -----------------------------------------------------------------

# mut <- read.csv('~/Dropbox/miseqUtils/mutcallsTest.csv')


# detectTemplatedIns ------------------------------------------------------

# overarching function which loops through rows of mut
# it is mainly a wrapper for detectTemplatedIns_one

# in detectMHdel, we could just test every row because detection is fast
# but detection of templated insertions is slow
# approach here is to only look at *unique* alignments where there is at least one insertion
# then copy the results of the detection in mut wherever the alignment is the same

# complete test dataset 151k rows: 8 min
# it is 1074 unique alignments with insertion

# will leave for now; but for the record: a further improvement would be to trim the alignment around cut position
# then keep only unique trimmed
# currently, we would repeat the detection if e.g. there is 1-bp substitution 150 bp from the cut position
# which is a waste of time as it does not change the outcome of the detection

# for 1k rows taken as random from test dataset: improvement is ~ 36 sec to ~ 6-7 sec
detectTemplatedIns <- function(mut,
                               cutpos,
                               cutdist,
                               allowStretchMatches=5,
                               extendNewlySeq=3,
                               searchWindowStarts=20,
                               minLCSbp=5) {
  
  ### "summarise" mut
  # only keep rows recording an insertion
  # then of those rows, just keep one example of each alignment
  # or in other words, each alignment (combination of "ref" and "ali") should be present only once
  mutuni <- mut %>%
    filter(type=='ins') %>%
    distinct(ref, ali, .keep_all=TRUE)
  # we work on mutuni below
  # then we will copy results into mut so that each row of mut gets annotated
  
  cat('\t \t \t \t >>>', nrow(mutuni), 'unique alignments, i.e. pairs of aligned reference & read.\n')
  
  ### recover original reference sequence
  # first, recover the original reference seauence,
  # i.e. the reference sequence used for alignment
  # (in alignment, the reference sequence can get --- to indicate insertions)
  # get all the 'true' original reference sequences:
  # which we do by just removing the --- from the aligned reference seauences
  orefs <- gsub('-', '', mutuni$ref) # "original references"
  # check they are all the same, after we remove NA
  orefs <- orefs[!is.na(orefs)]
  
  if(length(unique(orefs))!=1)
    stop('\t \t \t \t >>> Error detectTemplatedIns: issue when computing overall reference sequence.\n')
  # we have the original reference sequence
  oref <- unique(orefs)
  
  ### loop through every row of mutuni (so every unique alignment)
  # if type is ref or sub or del, return all NA
  # if type is ins, run detection of templated insertions
  tempInsL <- lapply(1:nrow(mutuni), function(rowi) {
    
    # here, ignore whether row we were given is labelled del or sub or ins
    # we just test every unique alignment
    cat('\t \t \t \t unique aligment', rowi, '/', nrow(mutuni), '\n')
    
    ### run detection of templated insertion on this alignment
    # alrow is alignment row
    alrow <- mutuni[rowi,]
    
    # detectTemplatedIns_one simply returns LCS (Longest Common Sub-string/sequence)
    # i.e. just a string
    detectedLCS <- detectTemplatedIns_one(alrow=alrow,
                                          oref=oref,
                                          cutpos=cutpos,
                                          cutdist=cutdist,
                                          minLCSbp=minLCSbp,
                                          allowStretchMatches=allowStretchMatches,
                                          extendNewlySeq=extendNewlySeq,
                                          searchWindowStarts=searchWindowStarts)
    
  })
  # from lapply, we get a big list where each element is detection (if any) of templated insertions from one row of mut
  # note detectTemplatedIns_one will always return a vector, even if detection inconclusive
  # make it into a dataframe
  tempIns <- as.data.frame(do.call(rbind, tempInsL))
  
  # fix some column types
  tempIns$lcbp <- as.integer(tempIns$lcbp)
  tempIns$lcStart <- as.integer(tempIns$lcStart)
  tempIns$lcStop <- as.integer(tempIns$lcStop)

  # add templated insertions info to mutuni
  # check same number of rows
  if(nrow(mutuni) != nrow(tempIns))
    stop('\t \t \t \t >>> Error detectTemplatedIns: not same number of rows when adding information about templated insertions to mutuni.\n')

  mutunins <- as.data.frame(cbind(mutuni, tempIns))
  
  # finally, add all the detections to mut
  # we only looked at unique alignments to save time
  # now can just copy the detection when the alignment is the same
  
  # from mutins, only keep ref, ali, and the new columns
  # then we can join by matching ref & ali
  # this will automatically copy all the new columns when alignment is the same
  mutins <- mutunins %>%
    select(ref, ali, lcbp, lcseq, lcdir, lcStart, lcStop) %>%
    left_join(x=mut, y=., by=c('ref', 'ali'))
  
  # one alignment will often have multiple rows, as each row is a mutation
  # e.g. could be an insertion & two substitutions, so three rows for one alignment
  # we copied the detection wherever the alignment is the same,
  # so in the example we will have same columns lcbp, lcseq, etc. for all three rows
  # framework is to detect templated *insertions*, so this may be confusing
  # I think better to only keep detections on rows which are insertion mutations
  # so now, wherever not ins, overwrite all NA
  
  mutins[mutins$type %in% c('ref', 'sub', 'del'),
         c('lcbp', 'lcseq', 'lcdir', 'lcStart', 'lcStop')] <- list(NA, NA, NA, NA, NA)
  # so now only rows where type is 'ins' can have a detection of templated insertions
  # (makes more sense)
  
  return(mutins)
  
}


# detectTemplatedIns_one --------------------------------------------------

# detection of templated insertion on one 'insertion' row of mut

detectTemplatedIns_one <- function(alrow,
                                   oref,
                                   cutpos,
                                   cutdist,
                                   minLCSbp,
                                   allowStretchMatches,
                                   extendNewlySeq,
                                   searchWindowStarts) {
  
  # split reference & aligned sequences
  refsp <- strsplit(alrow$ref, '')[[1]]
  alisp <- strsplit(alrow$ali, '')[[1]]
  # check they are the same length
  if(length(refsp) != length(alisp))
    stop('\t \t \t \t >>> Error detectTemplatedIns_one: something unexpected about the alignment; length of reference sequence is the not the same as the length of the aligned sequence.\n')
  
  ### step1: get the newly synthesised sequence
  # insertion is often part of a longer newly-synthesised sequence
  # which appears as substitutions - insertion
  # or insertion - substitutions
  
  # ! pitfall: if e.g. there is a sequencing error 100 bp from the cut, we risk taking this as part of the newly synthesised sequence
  # we need to limit ourselves to some region around the cut, this is why we ask the user for cutpos & cutdist
  # current 'cutpos' given by the user is on original reference used for alignment,
  # but we will need cut position on aligned reference
  # (position may be shifted due to insertions which appear as --- in aligned reference)
  # small function cutpos_aligned converts cutpos (on original reference) to cutpos on aligned reference
  cutpos_aref <- cutpos_aligned(cutpos=cutpos,
                                refsp=refsp)
  
  # previously, was approaching as: take position of modified nucleotide most to the left - up to - position of modified nucleotide most to the right
  # but does not work well because positions always refer to reference used for alignment (without any -) so does not give directly the newly synthesised sequence
  # better to go back to the alignment and look at first modified nucleotide / last modified nucleotide
  
  # finding the most likely "newly synthesised sequence" is the job of function newlySynthesised
  # newly synthesised sequence is a substring (sequence) of the aligned read
  # essentially from the first mismatch to the last mismatch (but there are a bit more details to it)
  newseq <- newlySynthesised(refsp=refsp,
                             alisp=alisp,
                             cutpos_aref=cutpos_aref,
                             cutdist=cutdist,
                             allowStretchMatches=allowStretchMatches,
                             extendNewlySeq=extendNewlySeq)
  
  # if the mutation is too far from the cut position (threshold controlled with cutdist)
  # newseq returns "MUT_TOO_FAR"
  # in this case, return NA
  if(newseq=='MUT_TOO_FAR') {
    return( c(lcbp=0,
              lcseq=NA,
              lcdir=NA,
              lcStart=NA,
              lcStop=NA) )
    # here, returning length of longest common substring (lcbp) as 0
    # obviously not actually 0 but helps to mark that we did analyse this mutation
    # but search for templated synthesis was inconclusive
    # (here because mutation was too far from cut position)
  }
  
  ### step2: get the search window
  # i.e. window of the reference sequence that we think could have been template for the newly synthesised sequence (will say newseq)
  # setting searchWindowStarts is the distance from cut we start the search window at
  # we will then add the length of the newly synthesised sequence
  # from looking at a few slc45a2 insertions: of nucleotides that served as template, closest is 1--9 bp from cut
  # so 20 bp is probably a good start
  # precise definition of searchWindowStarts:
  # of all the nucleotides that could have been a template for the newly synthesised sequence, how far (from the cut) do we allow the closest one to be
  # "when searching for the template, how far should we look"
  # "closest one", i.e. do not worry about accounting for the length of the insertion/newly synthesised sequence
  # e.g. if 10 bp, it means that we do not consider a potential template that is found at 11 bp (or further) from the cut
  
  # so, for search window, we go from cut up to + searchWindowStarts + length of newly synthesised sequence
  # in both ways (leftward & rightward)
  
  # search window:
  # start = cutpos - searchWindowStarts - length of newly synthesised sequence
  searchStart <- cutpos - (searchWindowStarts + nchar(newseq) - 1)
  # here minus 1, to be fully precise, we should account for the fact that cutpos is nucleotide before the cut
  # so without minus 1 we would take one nucleotide extra
  # and window would not be symmetrical (not same number of bp on either side of cut position)
  
  # ! if below 1, edit to 1
  if(searchStart < 1) {
    cat('\t \t \t \t >>> Start of search window was negative so edited to 1, i.e. the start of the reference.\n')
    searchStart <- 1
  }
  
  # stop = cutpos + searchWindowStarts + length of newly synthesised sequence
  searchStop <- cutpos + searchWindowStarts + nchar(newseq)
  # ! if above length of the reference, edit to length of the reference
  if(searchStop > nchar(oref)) {
    cat('\t \t \t \t >>> End of search window would go past the end of the reference so edited to last position of reference.\n')
    searchStop <- nchar(oref)
  }
  
  # the search window sequence is
  searchSeqFw <- substr(oref, start=searchStart, stop=searchStop)
  cat('\t \t \t \t >>> Length of search sequence is', nchar(searchSeqFw), 'bp, centered on the cut position.\n')
  cat('\t \t \t \t >>> search window:', searchSeqFw, '\n')
  
  # ! we also want to try search window from right to left
  # (not reverse-complement, just reverse)
  # I have seen examples where it is fairly convincing that it is the template
  searchSeqRv <- stringi::stri_reverse(searchSeqFw)
  
  # we now have searchSeqFw & searchSeqRv, our search windows
  
  ### step3: find if a significant chunk of newly synthesised sequence is from the search window (i.e. the flanking regions)
  # well known algorithmic problem called "longest common substring", or LCS
  # I have tried a few options found on SO
  # from Biostrings, pmatchPattern gives good outputs
  # no real explanation of the difference between pattern & subject though
  # on test sequence, swapping newseq / searchSeq as pattern or subject give different outputs, but both good
  # will try both and decide later which one we report
  
  # so first searching in Fw sequence, try both options
  lcs_fw1 <- as.character( pmatchPattern(newseq, searchSeqFw, maxlength.out=1L)[[1]] )
  lcs_fw2 <- as.character( pmatchPattern(searchSeqFw, newseq, maxlength.out=1L)[[1]] )
  
  # second searching in Rv sequence, try both options
  lcs_rv1 <- as.character( pmatchPattern(newseq, searchSeqRv, maxlength.out=1L)[[1]] )
  lcs_rv2 <- as.character( pmatchPattern(searchSeqRv, newseq, maxlength.out=1L)[[1]] )
  
  # get all the results together
  lcs <- c(lcs_fw1, lcs_fw2, lcs_rv1, lcs_rv2)
  names(lcs) <- c('fw1', 'fw2', 'rv1', 'rv2')
  
  # which is the longest?
  maxi <- which.max(nchar(lcs))
  lcseq <- as.character(lcs[maxi])
  lcnm <- names(lcs)[maxi]
  # looks like fw1 & fw2, and rv1 & rv2 will often give same length
  # here, which.max will just take the first one, so fw1 or rv1
  # ! this is arbitrary !
  # or in other words: I think guaranteed that we found the longest match
  # but there may be an equally good (long) solution
  
  # here, apply length threshold
  # below a certain length could be just chance, and not templated synthesis from flanking sequences
  # so we want to report only if it is above a certain length, controlled by minLCSbp
  if(nchar(lcseq) < minLCSbp) {
    cat('\t \t \t \t common substring is only', nchar(lcseq), 'bp long, which is below minLCSbp threshold, returning 0 for lcbp & NA for the rest.\n')
    return( c(lcbp=0,
              lcseq=NA,
              lcdir=NA,
              lcStart=NA,
              lcStop=NA) )
    # here, returning length of longest common substring (lcbp) as 0
    # obviously not actually 0 but helps to mark that we did analyse this mutation
    # but search for templated synthesis was inconclusive
    # (here because length of common substring was too short)
  }
  
  # locate the substring in original reference sequence
  # ! depends if Forward or Reverse match
  if(startsWith(lcnm, 'fw')) {
    lcdir <- 'fw'
    lcseqFw <- lcseq
  } else if(startsWith(lcnm, 'rv')) {
    lcdir <- 'rv'
    lcseqFw <- stringi::stri_reverse(lcseq)
  }
  
  # now locate it
  lcpos <- stri_locate(oref, regex=lcseqFw, mode='first')
  lcStart <- lcpos[1]
  lcStop <- lcpos[2]
  
  # check this really gives the longest common sequence from above
  if (!identical(substr(oref, lcStart, lcStop), lcseqFw))
    stop('\t \t \t \t >>> Error detectTemplatedIns_one: found longest common substring (sequence) to be ',
         lcseq, ' but the positions in the reference sequence do not give the same substring (sequence); they give ', substr(oref, lcStart, lcStop))
  
  # we are ready to return
  # lcseq: longest common substring (sequence);
  # string is always find in aligned read in forward direction (left to right)
  # but may need to be flipped to be found in reference
  # lcdir: fw if lcseq is found as it is in reference sequence; rv if lcseq has to be reversed to be found in reference sequence
  # lcStart & lcStop: positions in *original* reference sequence which give the longest common substring;
  # stop will always be after start but start:stop in reference may give the longest common substring in reverse, if lcdir is "rv"
  return(c(lcbp=nchar(lcseq),
           lcseq=lcseq,
           lcdir=lcdir,
           lcStart=lcStart,
           lcStop=lcStop))
  
}



# newlySynthesised --------------------------------------------------------

# refsp = aligned reference, split
# alisp = aligned read, split
# cutpos_aref = cut position on aligned reference
# cutdist
# allowStretchMatches
# extendNewlySeq

newlySynthesised <- function(refsp,
                             alisp,
                             cutpos_aref,
                             cutdist,
                             allowStretchMatches,
                             extendNewlySeq) {
  
  # walk through the alignment and compare each position
  seqcompare <- sapply(1:length(refsp), function(pos) {
    return(refsp[pos]==alisp[pos])
  })
  # series of TRUE & FALSE
  # TRUE if the sequences match at this position
  # FALSE if they do not
  
  # which are the modified positions?
  # any FALSE, so this can be mismatch, insertion, or deletion
  modpos <- which(!seqcompare)
  
  # use modified position that is closest to the cut position as 'seed'
  dist_fromcut <- modpos-cutpos_aref # gives distance from cut; signed so negative means before (left) / positive means after (right)
  mutseed <- cutpos_aref + dist_fromcut[which.min(abs(dist_fromcut))]
  # we have position (in alignment) of the modified nucleotide that is closest to the cut
  
  # if modified nucleotide that is closest to the cut
  # is in fact quite far from the cut position
  # we should let it go, as it is unlikely to be a Cas9-generated mutation
  
  # below: abs difference gives distance between mutseed and cut
  if( abs(mutseed-cutpos_aref) > cutdist ) {
    print(mutseed)
    print(cutdist)
    cat('\t \t \t \t mutation seed is returning 0 for lcbp & NA for the rest.\n')
    return('MUT_TOO_FAR')
    # then in detectTemplatedIns_one
    # we will understand what this means
  }
  
  # now we will try to extend the newly synthesised sequence in either directions
  # general logic is: if it continues being mismatches (could be insertion or substitution), it is probably part of the newly synthesised sequence
  # but we could get "lucky matches", which are positions where the newly synthesised sequence matches the original sequence,
  # even though the newly synthesised sequence was untemplated or copied from a region elsewhere
  # so we should not stop as soon as we have a single match
  # we allow a few matches; "a few" being allowStretchMatches
  
  ### extend LEFTward the furthest we can
  # we start at seed position
  leftmostMod <- mutseed
  extend <- TRUE
  # and try to extend more to the left
  # for each new position to the left, we look a few nucleotides ahead
  # if any modified nucleotide ahead, we continue extending
  while(extend) {
    # below: counting number of mismatches (cf. sum of !) from here to 5 nucleotides ahead
    # 5 is allowStretchMatches
    # e.g. if were allowing 0 match (i.e. we say "newly synthesised sequence cannot have any match with original sequence"),
    # then we would have to look at next position each time
    # if it matches, we would stop at the current position and call it "start of newly synthesised sequence"
    
    # more generally, whenever number of mismatches is 0 in the few nucleotides ahead, we are done
    # i.e. we found the start of the newly synthesised sequence
    if(sum(!seqcompare[(leftmostMod - allowStretchMatches) : leftmostMod]) > 0) {
      # if there are one or more mismatches, we extend more
      leftmostMod <- leftmostMod - 1
      
    } else {
      # no more mismatches ahead!
      # so we stop extending
      extend <- FALSE
      # precisely, we should stop at the nucleotide before the one we are trying now
      leftmostMod <- leftmostMod + 1
    }
  }
  
  
  ### extend RIGHTward the furthest we can
  # we start at seed position
  rightmostMod <- mutseed
  extend <- TRUE
  # and try to extend more to the right
  # for each new position to the right, we look  a few nucleotides ahead
  # if any modified nucleotide ahead, we continue extending
  while(extend) {
    # below: counting number of mismatches (cf. sum of !) from here to 5 nucleotides ahead
    # 5 is allowStretchMatches
    # e.g. if were allowing 0 match (i.e. we say "newly synthesised sequence cannot have any match with original sequence"),
    # then we would have to look at next position each time
    # if it matches, we would stop at the current position and call it "end of newly synthesised sequence"
    
    # more generally, whenever number of mismatches is 0 in the few nucleotides ahead, we are done
    # i.e. we found the end of the newly synthesised sequence
    if(sum(!seqcompare[(rightmostMod + allowStretchMatches) : rightmostMod]) > 0) {
      # if there are one or more mismatches, we extend more
      rightmostMod <- rightmostMod + 1
      
    } else {
      # no more mismatches ahead!
      # so we stop extending
      extend <- FALSE
      # precisely, we should stop at the nucleotide before the one we are trying now
      rightmostMod <- rightmostMod - 1
    }
  }
  
  # we now have first (left-most) mismatch and last (right-most) mismatch
  # guaranteed not to have other mismatches nearby (a bit before left-most or a bit after right-most)
  # precisely: there are no mismatches in the 'allowStretchMatches' nt before left-most
  # and there are no mismatches in the 'allowStretchMatches' nt after right-most
  
  # we can still extend them a little to allow "lucky matches" at the edges
  # i.e. nucleotides which were actually newly synthesised during repair but just happen to match original sequence
  # if "lucky matches" were random (should be approx. correct I think?),
  # probability of one should be 1/4, followed by another one is 1/4 * 1/4, etc
  # so probability that we are correct to extend by 2 nt is 6%
  # probability that we are correct to extend by 3 nt is 1.5%
  # (correct meaning the extra nucleotides we take were indeed newly synthesised during repair)
  # I think 3 nt extra on either side is reasonable
  # but will leave it flexible with extendNewlySeq
  
  # we can now finally extract our "newly synthesised sequence"
  newsp <- alisp[ (leftmostMod - extendNewlySeq) : (rightmostMod + extendNewlySeq) ]
  newseq <- paste0(newsp, collapse='')
  # note this newseq could include hyphens
  # we can leave them in
  
  cat('\t \t \t \t >>> newly synthesised sequence:', newseq, '\n')
  
  return(newseq)
  
}


# cutpos_aligned ----------------------------------------------------------

# small function to calculate cut position on aligned reference
# i.e. account for the --- which mark insertions
# not as simple as just adding number of ---,
# need to look at insertions starting before cutpos and add their lengths

# cutpos = cut position on original reference sequence
# refsp = aligned reference sequence, split

cutpos_aligned <- function(cutpos,
                           refsp,
                           alisp) {
  
  # we need to add the length of any insertion that starts before that position (cutpos)
  # find all "inserted" positions
  hypos <- which(refsp=='-') # hypos for hyphen positions
  
  # correct those positions so that they refer to the original reference
  # will make it a bit simpler later so all positions use the same reference
  # and we do not have to worry about correcting positions of later insertions to account for previous ones
  # we have function fixpos() in allelesToMutations.R for this
  ohypos <- sapply(hypos, fixpos, ref=refsp)
  # added o for 'original reference'
  
  # for example,
  # ...tgtaggtcgtcat---------gGGG
  # last t is #89
  # hypos would be 90, 91, ..., 98
  # correction by fixpos will give 89, 89, ..., 89
  # i.e. will mark, on the original reference, the last nucleotide before the insertion
  
  # ! there could be multiple insertions before cutpos
  # if so, find their boundaries
  # add NA at the start otherwise we get original length minus 1 
  # (which I find confusing to account for later)
  ohypos_dif <- c(NA, diff(ohypos))
  
  # diff gives, at each position n, (n+1) - n
  # so usually 0 (because insertion positions, after correction with fixpos, are e.g. 89, 89, 89, 89, ...)
  # but, if there are more than one insertion, then we expect e.g. 0 0 0 4 0 0 0 0
  # so wherever greater than 0 there is a jump in the positions, i.e. a new insertion starts
  # we should keep separate the hyphen positions of each insertion
  # we have small function splitAt from allelesToMutations.R which is exactly for that purpose
  ohypos_byins <- splitAt(ohypos, which(ohypos_dif>0))
  # hyphen positions, split for each insertion
  
  # small example:
  # positions with hyphens, after fixpos, are 89, 89, 89, 92, 92, 92, 101, 101, 101
  # so diff will give NA 0 0 3 0 0 9 0 0
  # so positions to split at are #4 & #7
  # splitAt will give a list of three elements:
  # 1: 89 89 89
  # 2: 92 92 92
  # 3: 101 101 101
  
  # we have cutpos & positions of each insertion all on original reference
  # now we look at the start position of each insertion
  # if this insertion started before cutpos, we need add the total length of this insertion to cutpos
  toadd <- lapply(ohypos_byins, function(ohyposi) {
    # about <= cutpos below, instead of < cutpos
    # example:
    # agtaaC*ttaNGGacc
    # cut is at *; cutpos points to C (#6)
    # if insertion like so:
    # agtaaC----ttaNGGacc
    # its positions (after fixpos) will be 6 6 6 6; same as cutpos
    # we cannot easily tell where the cut (or ligation, rather) was!
    # could be after the insertion: agtaaC----*ttaNGGacc
    # (left end was extended)
    # > in which case we need to correct cutpos
    # or before the insertion: agtaaC*----ttaNGGacc
    # (right end was extended)
    # > in which case we do not need to correct cutpos
    # no answer will be always correct! will (arbitrarily) assume first case,
    # i.e. we correct cutpos if an insertion starts at that position
    # so <= (instead of <)
    if(min(ohyposi) <= cutpos) {
      # if an insertion starts before cutpos, record its length
      return(length(ohyposi))
    } else {
      # if not, record 0 as the length
      return(0)
    }
  })
  # we have, for each insertion, number of positions to add to correct cutpos
  # sum them & add them to cutpos
  cutpos_aref <- cutpos + sum(unlist(toadd))
  # we now have cutpos on aligned reference
  
  # cat('\t \t \t \t >>> received cutpos =', cutpos, '; converted to cutpos on aligned reference =', cutpos_aref, '\n')
  
  return(cutpos_aref)
  
}