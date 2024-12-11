#####################################################
# ~ miseqUtils/CUTTER: detect probable repair on mutations ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# v0
# currently only microhomologies on deletions

# expects mutation table written by callMutations


# example -----------------------------------------------------------------

# mut <- read.csv('~/Dropbox/miseqUtils/mutcallsTest.csv')

# mutmh <- detectMHdel(mut=mut,
#                      minMHlen=2,
#                      cutpos=87)


# detectMHdel -------------------------------------------------------------

# overarching function for detection of microhomologies in deletions
# runs detectMHdel_one for each unique deletion in mut (written by callMutations)

detectMHdel <- function(mut,
                        minMHlen=2,
                        cutpos) {
  
  # ! in mut, mutation positions (start & stop) refer to original reference sequence
  # *not* to reference seauence aligned (column ref)
  # distinction is important because in the case of an insertion, ref gets additional "-"
  # to keep positions always referring to the same reference, we have a function "fixPos" in allelesToMutations to correct for this
  # but this means that in functions below we should always work with the original reference, not the reference listed in that row of mut
  # so that user does not have to give manually the reference sequence, we can recover it here ourselves
  # only case where reference is not actually the original reference sequence is if it has any "-"
  # get all the 'true' original reference sequences:
  orefs <- mut[!grepl('-', mut$ref, fixed=TRUE), 'ref'] # "overall references"
  # check they are all the same, after we removed NA
  orefs <- orefs[!is.na(orefs)]
  if(length(unique(orefs))!=1)
    stop('\t \t \t \t >>> Error detectMHdel: issue when computing overall reference sequence.\n')
  # we have our overall reference sequence
  oref <- unique(orefs)
  
  ### loop through every row of mut
  # if type is ref or sub or ins, return all NA
  # if type is del, run MH detection
  delmhL <- lapply(1:nrow(mut), function(rowi) {
    
    if(mut[rowi, 'type'] %in% c('ref', 'sub', 'ins')) {
      return(c(MHside=NA, # differentiating vs. none which is when a deletion was actually analysed but did not return any MH
               MHseq=NA,
               MHbp=NA,
               leftMHstart=NA,
               leftMHstop=NA,
               rightMHstart=NA,
               rightMHstop=NA,
               innerSeq=NA,
               innerbp=NA,
               leftFlapSeq=NA,
               leftFlapbp=NA,
               rightFlapSeq=NA,
               rightFlapbp=NA))
    }
    # rest is 'else' (because did not 'return' yet), which is = del
    
    # if not returned, it means we have a deletion
    ### run MH detection
    delrow <- mut[rowi,]
    
    detectedMH <- detectMHdel_one(delrow=delrow,
                                  oref=oref,
                                  minMHlen=minMHlen)
    
    return( recordMH(delrow=delrow,
                     oref=oref,
                     cutpos=cutpos,
                     detectedMH=detectedMH) )
    
  })
  delmh <- as.data.frame(do.call(rbind, delmhL))
  
  # fix some column types
  delmh$MHbp <- as.integer(delmh$MHbp)
  delmh$leftMHstart <- as.integer(delmh$leftMHstart)
  delmh$leftMHstop <- as.integer(delmh$leftMHstop)
  delmh$rightMHstart <- as.integer(delmh$rightMHstart)
  delmh$rightMHstop <- as.integer(delmh$rightMHstop)
  delmh$innerbp <- as.integer(delmh$innerbp)
  delmh$leftFlapbp <- as.integer(delmh$leftFlapbp)
  delmh$rightFlapbp <- as.integer(delmh$rightFlapbp)
  
  # add MH deletion info to mut
  # check same number of rows
  if(nrow(mut) != nrow(delmh))
    stop('\t \t \t \t >>> Error detectMHdel: not same number of rows when adding MH information to mut.\n')
  
  mutmh <- as.data.frame(cbind(mut, delmh))
  
  return(mutmh)

}


# recordMH ----------------------------------------------------------------
# separate function to record all the information we want about the microhomology
# probably a bit easier to edit than recording this information directly when we detect the microhomology and having to pass it up across multiple functions

# it receives the output of detectMH as input
# i.e. a vector with MHside and MHseq, will call it detectedMH

# cutpos: nucleotide before the cut

recordMH <- function(delrow,
                     oref,
                     cutpos,
                     detectedMH) {
  
  # ! if getting detectedMH as NA, NA
  # should just return everything as NA
  if( is.na(detectedMH[1]) & is.na(detectedMH[2]) ) {
    return(c(MHside='none', # will return as 'none' so can easily see those which were analysed but did not return any MH
             MHseq=NA,
             MHbp=NA,
             leftMHstart=NA,
             leftMHstop=NA,
             rightMHstart=NA,
             rightMHstop=NA,
             innerSeq=NA,
             innerbp=NA,
             leftFlapSeq=NA,
             leftFlapbp=NA,
             rightFlapSeq=NA,
             rightFlapbp=NA))
  }
  # below is "else", i.e. we did not return above yet
  
  # from detectedMH, we already have
  ### MHside
  # this is side of MH still in sequence, i.e. not deleted
  MHside <- as.character(detectedMH['MHside'])
  ### MHseq
  MHseq <- as.character(detectedMH['MHseq'])
  
  # also record
  ### length of MH
  MHbp <- nchar(MHseq)
  
  ### start & stop positions of leftMH
  # if MHside (which MH was deleted) is left, then we want the first few positions of deletion
  # if MHside (which MH was deleted) is right, then we want positions just before deletion
  if(MHside=='left') {
    leftMHstart <- delrow$start
    leftMHstop <- delrow$start + MHbp - 1
  } else if(MHside=='right') {
    leftMHstart <- delrow$start - 1 - MHbp + 1
    leftMHstop <- delrow$start - 1
  }
  # check that it is correct
  if(substr(oref, leftMHstart, leftMHstop) != MHseq)
    stop('\t \t \t \t >>> Error: positions leftMHstart to leftMHstop in ref do not give the MH sequence.\n')

  ### start & stop positions of rightMH
  # if MHside (which MH was deleted) is left, then we want positions just after deletion
  # if MHside (which MH was deleted) is right, then we want last few positions of deletion
  if(MHside=='left') {
    rightMHstart <- delrow$stop + 1
    rightMHstop <- delrow$stop + MHbp
  } else if(MHside=='right') {
    rightMHstart <- delrow$stop - MHbp + 1
    rightMHstop <- delrow$stop
  }
  # check that it is correct
  if(substr(oref, rightMHstart, rightMHstop) != MHseq)
    stop('\t \t \t \t >>> Error: positions rightMHstart to rightMHstop in ref do not give the MH sequence.\n')
  
  ### "inner" sequence
  
  # take care of odd cases first
  # all are when cut is not between left & right MH
  # we could simply not consider those events
  # but I think they are interesting to record because some do have MH
  
  # #1
  # in some cases, we can have leftMH directly flanking rightMH
  # and deletion is exactly one MH
  # e.g.
  # AGGTCGTCAT
  # AGGTC---AG
  # (MH is GTC)
  
  # #2
  # (which can include some cases of #1)
  # we have cut site not between the MH
  # I do not know how to explain this event
  # a possibility is that Cas9 did not cut where expected
  
  ## case #1 ##
  # easier to take care of this situation separately
  # so if just after leftMHstop is directly rightMHstart
  if(leftMHstop+1==rightMHstart) {
    # innerSeq is an empty string
    innerSeq <- ''
    # of length 0
    innerbp <- 0
    # leftFlapSeq is NA
    leftFlapSeq <- NA
    # of length NA
    leftFlapbp <- NA
    # rightFlapSeq is NA
    rightFlapSeq <- NA
    # of length NA
    rightFlapbp <- NA
  
  ## case #2 ##
  # if leftMH stops *after* cut
  # or if rightMH starts *before* cut
  # both are including =, same position still means MH stops after cut or starts before cut
  # e.g. (real example) ...agta*tggGGG; MH is atg but cut is here *
  # this gives cutpos = 4 (should select the nucleotide just before) & rightMHstart as 4 (selects the a of atg)
  } else if( (leftMHstop>=cutpos) | (rightMHstart<=cutpos) ) {
    # we can record innerSeq as normal (cf. `else` below)
    innerSeq <- substr(oref, leftMHstop+1, rightMHstart-1)
    innerbp <- nchar(innerSeq)
    
    # but we record flaps as NA, as we do not know where the cut occurred
    leftFlapSeq <- NA
    leftFlapbp <- NA
    rightFlapSeq <- NA
    rightFlapbp <- NA
    
  } else {
    # if not, we proceed as normal
    # this is sequence between the two MH
    # we now calculated all the positions,
    # so we can just take from leftMHstop + 1 up to rightMHstart - 1
    innerSeq <- substr(oref, leftMHstop+1, rightMHstart-1)
    
    ### length of "inner" sequence
    innerbp <- nchar(innerSeq)
    
    ### sequence between leftMH and cut
    # I think safest is to get the actual sequence and count the number of nt
    leftFlapSeq <- substr(oref, leftMHstop+1, cutpos)
    # cutpos should be precisely the nucleotide before the cut
    # check sequence makes sense, should be first few nucleotides of innerSeq
    if (!identical( leftFlapSeq , substr(innerSeq, 1, nchar(leftFlapSeq)) ))
      stop('\t \t \t \t >>> leftFlapSeq does not give first n nucleotides of innerSeq.\n')
    
    ### distance between leftMH and cut
    leftFlapbp <- nchar(leftFlapSeq)
    
    ### sequence between cut and rightMH
    rightFlapSeq <- substr(oref, cutpos+1, rightMHstart-1)
    # check sequence makes sense, should be last few nucleotides of innerSeq
    if (!identical( rightFlapSeq , substr(innerSeq, nchar(innerSeq)-nchar(rightFlapSeq)+1, nchar(innerSeq)) ))
      stop('\t \t \t \t >>> rightFlapSeq does not give last n nucleotides of innerSeq.\n')
    # this should also be the same as innerSeq - leftFlapSeq
    # or in other words, innerSeq should be exactly leftFlapSeq + rightFlapSeq
    if(!identical(innerSeq ,  paste0(leftFlapSeq, rightFlapSeq, collapse='')))
      stop('\t \t \t \t >>> innerSeq is not exactly leftFlapSeq + rightFlapSeq.\n')
    
    ### distance between cut and rightMH
    rightFlapbp <- nchar(rightFlapSeq)
  }
  
  ### return all info in one vector
  return(c(MHside=MHside,
           MHseq=MHseq,
           MHbp=MHbp,
           leftMHstart=leftMHstart,
           leftMHstop=leftMHstop,
           rightMHstart=rightMHstart,
           rightMHstop=rightMHstop,
           innerSeq=innerSeq,
           innerbp=innerbp,
           leftFlapSeq=leftFlapSeq,
           leftFlapbp=leftFlapbp,
           rightFlapSeq=rightFlapSeq,
           rightFlapbp=rightFlapbp))
  
}



# detectMHdel_one ---------------------------------------------------------

# will find best possible (longest) MH
# returns a vector with two information:
# MHside: left or right, this tells whether MH *deleted* (i.e. not left in sequence) was the left or right MH
# remember, I think this is always (?) an arbitrary decision by the alignment algorithm, but could be interesting to tell whether it has a preference or not
# MHseq: sequence of the detected MH
detectMHdel_one <- function(delrow,
                            oref,
                            minMHlen) {
  
  ### try left MH
  # from minimum size until we fail
  # fail = 0 means we did not fail yet
  nMH <- minMHlen
  # preallocate to NA, i.e. no MH detected
  leftMH <- NA
  
  fail <- 0
  while(fail==0) {
    # leftMHnext is testing +1 longer MH
    leftMHnext <- tryLeftMH(delrow=delrow,
                            oref=oref,
                            nMH=nMH)
    
    # if NA, we failed
    # so we keep the previous leftMH
    # if not, we update leftMH and we try +1 longer MH
    if(is.na(leftMHnext)) {
      fail <- 1
    } else {
      leftMH <- leftMHnext
      nMH <- nMH + 1
      # and while loop will restart
    }
  }
  
  ### try right MH
  # from minimum size until we fail
  # fail = 0 means we did not fail yet
  nMH <- minMHlen
  # preallocate to NA, i.e. no MH detected
  rightMH <- NA
  
  fail <- 0
  while(fail==0) {
    rightMHnext <- tryRightMH(delrow=delrow,
                              oref=oref,
                              nMH=nMH)
    
    # if NA, we failed
    # so we keep the previous rightMH
    # if not, we update rightMH and we try +1 longer MH
    if(is.na(rightMHnext)) {
      fail <- 1
    } else {
      rightMH <- rightMHnext
      nMH <- nMH + 1
      # and while loop will restart
    }
  }
  
  ### compare left & right MH
  # we know have left & right MH
  # compare which one is best (i.e. exists & is longest)
  if(is.na(leftMH) & is.na(rightMH)) {
    # cat('\t \t \t \t no MH detected\n') # slows down to print all of them
    return(c(side=NA,
             MHseq=NA))
    
  } else if(!is.na(leftMH) & is.na(rightMH)) {
    cat('\t \t \t \t left MH:', leftMH, '\n')
    return(c(MHside='left',
             MHseq=leftMH))
    
  } else if(is.na(leftMH) & !is.na(rightMH)) {
    cat('\t \t \t \t right MH:', rightMH, '\n')
    return(c(MHside='right',
             MHseq=rightMH))
    
  } else if(!is.na(leftMH) & !is.na(rightMH)) {
    # I am guessing this is a rare occurence
    # if it happens, report it & take the longest one
    bothMH <- c(leftMH, rightMH)
    bothMHi <- which.max(nchar(bothMH))
    longestMH <- bothMH[bothMHi]
    
    # note, I guess it is possible to have left & right MH both the same size; example aligment:
    # GTAGTAAGGCAGTGTAGTA
    # GTA-------------GTA
    # (could even be different MH sequence)
    # here, will arbitrarily take left one
    # because e.g. which.max(c(9, 9)) says 1
    
    if(bothMHi==1) {
      MHside <- 'left'
      cat('\t \t \t \t MH on both sides! Longest is left MH:', longestMH, '\n')
    } else if(bothMHi==2) {
      MHside <- 'right'
      cat('\t \t \t \t MH on both sides! Longest is right MH:', longestMH, '\n')
    }
    return(c(MHside=MHside,
             MHseq=longestMH))
  }
  
}


# tryLeftMH ---------------------------------------------------------------

tryLeftMH <- function(delrow,
                      oref,
                      nMH) {
  
  # ! microhomology length (nMH) cannot be longer than the deletion, does not make sense
  # maximum size allowed for the MH is exactly the size of the deletion
  if(nMH > nchar(delrow$refseq)) {
    return(NA)
  }
  
  # preallocate to NA
  leftMH <- NA
  # take first (left) n nt of deleted sequence
  tryMH <- substr(delrow$refseq, 1, nMH)
  
  # positions of tryMH will be first n positions of deletion
  # Same is for same side, i.e. positions of the tryMH on the same side
  startposSame <- delrow$start
  stopposSame <- delrow$start + nchar(tryMH) - 1
  
  # in the reference sequence, we know already that at these positions are the same nucleotides
  # i.e. the "tryMH" sequence
  # can check this
  if( ! identical( substr(oref, startposSame, stopposSame) , tryMH ) )
    stop('\t \t \t \t >>> Error tryLeftMH: in the reference sequence, expected to find the same sequence as the deleted sequence, but that was not the case.\n')
  
  # check at the opposite side of this tryMH, just after the deletion (i.e. not in the deleted sequence)
  # if we find the 'tryMH' sequence, this is evidence of MMEJ
  
  # those positions will be end of deletion + 1 >> end of deletion + length of tryMH
  # Oppo is for opposide, i.e. positions of the tryMH on the opposite side
  startposOppo <- delrow$stop + 1
  stopposOppo <- delrow$stop + 1 + nchar(tryMH) - 1
  
  # if this sequence is the same as tryMH, this is evidence of MMEJ
  # we record leftMH as this sequence
  if( substr(oref, startposOppo, stopposOppo) == tryMH) {
    leftMH <- tryMH
  }
  
  # return
  return(leftMH)
}


# tryRightMH --------------------------------------------------------------

tryRightMH <- function(delrow,
                       oref,
                       nMH) {
  
  # ! microhomology length (nMH) cannot be longer than the deletion, does not make sense
  # maximum size allowed for the MH is exactly the size of the deletion
  if(nMH > nchar(delrow$refseq)) {
    return(NA)
  }
  
  ### try right MH
  # preallocate to NA
  rightMH <- NA
  
  # take last (right) n nt of deleted sequence
  tryMH <- substr(delrow$refseq, nchar(delrow$refseq)-nMH+1, nchar(delrow$refseq))
  
  # positions of tryMH will be last two positions of deletion
  # Same is for same side, i.e. positions of the tryMH on the same side
  startposSame <- delrow$stop - nchar(tryMH) + 1
  stopposSame <- delrow$stop
  
  # in the reference sequence, we know already that at these positions are the same nucleotides
  # i.e. the "tryMH" sequence
  # can check this
  if( ! identical( substr(oref, startposSame, stopposSame) , tryMH ) )
    stop('\t \t \t \t >>> Error tryRightMH: in the reference sequence, expected to find the same sequence as the deleted sequence, but that was not the case.\n')
  
  # check at the opposite side of this tryMH, just before the deletion (i.e. not in the deleted sequence)
  # if we find the 'tryMH' sequence, this is evidence of MMEJ
  
  # those positions will be start of deletion - length of MH
  # until start of deletion - 1
  # Oppo is for opposide, i.e. positions of the tryMH on the opposite side
  startposOppo <- delrow$start - nchar(tryMH)
  stopposOppo <- delrow$start - 1
  
  # if this sequence is the same as tryMH, this is evidence of MMEJ
  # we record rightMH as this sequence
  if( substr(oref, startposOppo, stopposOppo) == tryMH) {
    rightMH <- tryMH
  }
  
  # return
  return(rightMH)
}