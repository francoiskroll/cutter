#####################################################
# ~ cutter: simulateIns ~
#
# ... ... ...
# Francois Kroll 2025
# francois@kroll.be
#####################################################

# v0
# exponential fitted in fitIns_notes.R on Cas9MiSeqDB v0
# rate = 0.06643421

# TODO: note, we overestimate almost certainly
# because we will generate many random sequences
# while a real sample would only have a few possible alleles
# right, but then each random match would only contribute 1 out of 1000 in the simulated sample
# I guess the probabilities cancel each other?


# simulateIns -------------------------------------------------------------

# ! assumes that mut was already ran through detectTemplatedIns

simulateIns <- function(mut,
                        nreads=1000,
                        extendNewlySeq=3,
                        searchWindowStarts=20,
                        minLCSbp=5,
                        fit_exprate=0.06643421) {
  
  # ! can only simulate a sample for *one* locus
  # check that we only have one locus
  if(length(unique(mut$locus))>1) stop('\t \t \t \t Error simulateIns: can only simulate a sample for *one* locus at a time.\n')
  
  ### get overall reference sequence
  # we simply remove the - hyphens
  orefs <- gsub('-', '', mut$ref)
  
  # check that all the same
  oref <- unique(orefs)
  
  if(length(oref)>1)
    stop('\t \t \t \t >>> Error simulateIns: more than one reference sequence in this mutation table.\n')
  
  # get cutpos
  cutpos <- unique(mut$cutpos)
  # check only one
  if(length(cutpos)>1)
    stop('\t \t \t \t >>> Error simulateIns: more than one cutpos in this mutation table.\n')
  
  ### now simulate `nreads` newly synthesised sequences
  # step below is surprisingly slow...
  siminsL <- lapply(1:nreads, function(i) {
    cat('\t \t \t \t simulating newly synthesised sequence', i, 'of', nreads, '\n')
    simulateIns_one(mut=mut,
                    extendNewlySeq=extendNewlySeq,
                    searchWindowStarts=searchWindowStarts,
                    minLCSbp=minLCSbp,
                    cutpos=cutpos,
                    oref=oref,
                    fit_exprate=fit_exprate)
  })
  simins <- data.frame(do.call(rbind, siminsL))
  # >> we get a dataframe with columns newlyseq, newlyseq_bp, etc
  # including the LCS etc.
  # as in detectTemplatedIns
  
  ### add rid
  # normally done by alleleToMutation
  # here, will simply be r1.1, r2.1, r3.1, etc.
  # until e.g. r1000.1 if nreads=1000
  simins$rid <- sprintf('r%i.1', 1:nrow(simins))
  
  ### add coverage
  # copying from allelesToMutations so we get the same columns
  # calculate total number of reads for this sample, i.e. coverage
  # and add it as column
  # here, we know this is exactly `nreads`
  cat('\t \t \t \t >>> Total number of reads:', nreads, '\n')
  simins <- simins %>%
    mutate(splcov=nreads)
  
  ### add rundate, well, locus, type
  # these we know are present in mut
  simins <- simins %>%
    mutate(well='Z99', .before=1) %>%
    mutate(locus='simulated', .before=1) %>%
    mutate(rundate=999999, .before=1) %>%
    mutate(type='ins', .before=1)
  
  ### also create sample column
  # which is rundate_well_locus
  simins <- simins %>%
    mutate(sample=paste(rundate, well, locus, sep='_'), .before=1)
  
  ### add last columns as NA
  # find the ones we are missing from mut
  # these were added from meta
  cols2add <- colnames(mut) [! colnames(mut) %in% colnames(simins)]
  
  for(colnm in cols2add) {
    simins <- simins %>%
      mutate(tmp=NA, .before='well')
    # now put the correct column name
    colnames(simins)[which(colnames(simins)=='tmp')] <- colnm
  }
  
  # refill cutpos column
  simins$cutpos <- cutpos
  
  ### add simulated reads to mut
  # first put the columns of simins in the same order as in mut
  
  # check columns are identical (before looking at their positions)
  if(!identical( sort(colnames(mut)) , sort(colnames(simins))))
    stop('\t \t \t \t >>> Error simulateIns: all columns of mut are not in simins, or vice-versa.\n')
  
  # now put the columns of simins in the same order
  simins <- simins[, colnames(mut) ]
  # check correct
  if(!identical(colnames(mut), colnames(simins)))
    stop('\t \t \t \t >>> Error simulateIns: columns of simins are not in the same order as columns of mut.\n')
  
  # can now rbind
  mutsim <- rbind(mut, simins)
  # and we return result
  return(mutsim)
  
}
# we now added a simulate sample
# which has 1,000 simulated newly synthesised sequences
# includes detection of templated insertions



# simulateIns_one ---------------------------------------------------------

simulateIns_one <- function(mut,
                            extendNewlySeq,
                            searchWindowStarts,
                            minLCSbp,
                            cutpos,
                            oref,
                            fit_exprate) {
  
  ### find the "extendNewlySeq" extensions
  # "extendNewlySeq" extensions for this locus are
  # left extension: cutpos-2, cutpos-1, cutpos
  # right extension: cutpos+1, cutpos+2, cutpos+3
  
  ### left extension seq is
  extLpos <- cutpos - seq(extendNewlySeq-1, 0) # positions of right extension sequence
  extL <- substring(oref, extLpos[1], extLpos[length(extLpos)])
  
  # check its length
  if(nchar(extL)!=extendNewlySeq)
    stop('\t \t \t \t Error simulateIns: simulated left extension sequence is not the same length as extendNewlySeq.\n')
  
  ### right extension seq is
  extRpos <- cutpos + seq(1, extendNewlySeq) # positions of right extension sequence
  extR <- substring(oref, extRpos[1], extRpos[length(extRpos)])
  
  # check its length
  if(nchar(extR)!=extendNewlySeq)
    stop('\t \t \t \t Error simulateIns: simulated right extension sequence is not the same length as extendNewlySeq.\n')
  
  ### draw a length for the simulated newly synthesised sequence
  # we need the drawn length to be below extendNewlySeq * 2 + 1 (typically 7 bp)
  minlen <- extendNewlySeq*2 + 1
  
  ### from the fitted exponential distribution, we draw one newly synthesised sequence length
  # and round it to closest integer
  # ! if this integer is below , we draw again
  len <- 0
  while(len < minlen) {
    len <- round(rexp(n=1,
                      rate=fit_exprate))
  }
  # >> we have a length above 7 bp
  
  # then we should randomly generate a sequence this long:
  nran <- len - 2 * extendNewlySeq # number of bp to draw randomly
  # as we will then add the extensions
  
  ### simulate the newly seq without the extension sequences
  seqran <- paste(sample(c('A', 'C', 'G', 'T'),
                         size=nran,
                         replace=TRUE),
                  collapse='')
  
  # now add the extension sequences
  newlyseq_sim <- paste0(extL, seqran, extR, collapse='')
  # >> we have the simulated extension sequence
  
  ### now run searchLCS() from detectTemplatedIns.R
  # will prepare searchWindow and look for LCS
  # we can return directly its output
  return( searchLCS(newseq=newlyseq_sim,
                    oref=oref,
                    minLCSbp=minLCSbp,
                    searchWindowStarts=searchWindowStarts,
                    cutpos=cutpos) )
}
