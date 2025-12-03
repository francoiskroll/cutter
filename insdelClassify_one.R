#####################################################
# ~ ZPRI: mutations to read labels ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# expect mut, mutation table created by filterMutations

# v0
# insdelClassify gives four labels: 
# no indel
# insertion
# deletion
# both

insdelClassify_one <- function(mut,
                               unwantedSubs) {
  ### check mut is giving data for just one sample
  if(length(unique(mut$sample))>1)
    stop('\t \t \t \t >>> Error: attempting to give mutation data for more than one sample.\n')
  
  # record rundate, sid, locus, well, grp sample, sample's coverage
  rundate <- unique(mut$rundate)
  sid <- unique(mut$sid)
  locnm <- unique(mut$locus)
  wellnm <- unique(mut$well)
  grpnm <- unique(mut$grp)
  splnm <- unique(mut$sample)
  
  # check that there is only one sample's coverage
  splcov <- unique(mut$splcov)
  if(length(splcov)!=1)
    stop('\t \t \t \t >>> Error: 0 or multiple sample coverage.\n')
  
  # gather unique read IDs
  uids <- unique(mut$rid)
  # loop through unique read IDs
  rlabsL <- lapply(uids, function(ui) {
    
    ## preallocate vector which is 1) insertion?; 2) deletion?
    # preallocate as FALSE, FALSE so = no indel read
    rlab <- c(ins=FALSE, del=FALSE)
    
    ## gather all the detected mutations in this read
    muti <- mut %>%
      filter(rid==ui)
    
    # if mutation "ref" is present
    if('ref' %in% unique(muti$mutid)) {
      # there should not be anything else, check that
      if(nrow(muti)!=1)
        stop('\t \t \t \t >>> Error: for read', ui,', there is mutation "ref" and other mutations, which does not make sense.\n')
      # and we can return vector rlab as it is, i.e. FALSE, FALSE
      return(rlab)
    }
    # we are done here as we returned, below is `else`
    # so everything below: read has at least one mutation
    # and we know "ref" is not there
    
    ## is there an insertion?
    if('ins' %in% muti$type) {
      # then switch ins flag to TRUE
      rlab['ins'] <- TRUE
    
    ## is there a deletion?
    } else if ('del' %in% muti$type) {
      # then switch del flag to TRUE
      rlab['del'] <- TRUE
    }
    
    ## return labels for this read
    return(rlab)
  })
  
  # create a simpler dataframe
  rlabs <- data.frame(do.call(rbind, rlabsL))
  rlabs <- rlabs %>%
    mutate(rid=uids, .before=1)
  
  # now add columns to classify reads
  # read categories
  rcats <- sapply(1:nrow(rlabs), function(r) {
    # row is:
    rl <- rlabs[r,]
    
    # if no insertion, no deletion > noindel
    if(!rl['ins'] & !rl['del']) {
      return('noindel')
    # if insertion & no deletion > insertion
    } else if(rl['ins'] & !rl['del']) {
      return('insertion')
    # if no insertion & deletion > deletion
    } else if(!rl['ins'] & rl['del']) {
      return('deletion')
    # if insertion & deletion > both
    } else if(rl['ins'] & rl['del']) {
      return('both')
    }
  })
  
  rlabs <- rlabs %>%
    mutate(cat=rcats)
  
  ### add meta columns
  # add from right to left in final columns
  rlabs <- rlabs %>%
    mutate(grp=grpnm, .before=1) %>%
    mutate(well=wellnm, .before=1) %>%
    mutate(locus=locnm, .before=1) %>%
    mutate(sid=sid, .before=1) %>%
    mutate(rundate=rundate, .before=1) %>%
    mutate(sample=splnm, .before=1) %>%
    mutate(splcov=splcov, .after='grp')
  
  return(rlabs)
  
}
