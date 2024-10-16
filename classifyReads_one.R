#####################################################
# ~ ZPRI: mutations to read labels ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# expect mut, mutation table created by filterMutations

classifyReads_one <- function(mut,
                              expedit) {
  
  ### check one giving data for one sample
  if(length(unique(mut$sample))>1)
    stop('\t \t \t \t >>> Error: attempting to give mutation data for more than one sample.\n')
  # record rundate, sid, locus, well, grp sample, sample's coverage
  rundate <- unique(mut$rundate)
  sid <- unique(mut$sid)
  locnm <- unique(mut$locus)
  wellnm <- unique(mut$well)
  grpnm <- unique(mut$grp)
  splnm <- unique(mut$sample)
  
  # check only one sample's coverage
  splcov <- unique(mut$splcov)
  if(length(splcov)!=1)
    stop('\t \t \t \t >>> Error: 0 or multiple sample coverage.\n')
  
  # gather unique read IDs
  uids <- unique(mut$rid)
  # loop through unique read IDs
  rlabsL <- lapply(uids, function(ui) {
    
    ## preallocate vector which is 1) prime-edit present?; 2) other mutation present?
    # preallocate as FALSE, FALSE so = reference read
    rlab <- c(expedit=FALSE, mut=FALSE) # read labels
    
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
    } # we are done here, below is `else`
    # so everything below: read has at least one mutation
    # and we know "ref" is not there
    
    ## is prime-edit present?
    if(expedit %in% muti$mutid) {
      # then switch edit label to TRUE
      rlab['expedit'] <- TRUE
    }
  
    ## other mutation(s) present?
    # can simply count number of rows
    # if edit is present and > 1 rows, then other mutation
    # if edit is absent and > 0 rows, then other mutation
    if(rlab['expedit'] & nrow(muti)>1) {
      rlab['mut'] <- TRUE
    } else if(!rlab['expedit'] & nrow(muti)>0) {
      rlab['mut'] <- TRUE
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
    
    if(!rl['expedit'] & !rl['mut']) {
      return('reference')
    } else if(rl['expedit'] & !rl['mut']) {
      return('pure')
    } else if(rl['expedit'] & rl['mut']) {
      return('impure')
    } else if(!rl['expedit'] & rl['mut']) {
      return('mutated')
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
