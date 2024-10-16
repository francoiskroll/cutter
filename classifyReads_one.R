#####################################################
# ~ ZPRI: mutations to read labels ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# expect mut, mutation table created by filterMutations

# v1
# v2: scaffold detection
# scaffdetectWin: distance before/after end of RHA to detect scaffold
# for now, detect scaffold just means insertion or deletion that starts with G
# ! definitions of classes are slightly different:
# pure = edit / no mutation / no scaffold
# impure = edit / mutation / no scaffold
# mutated = no edit / mutation / no scaffold
# scaffold = edit or not / mutation / scaffold
# so scaffold class is agnostic to edit or not, i.e. whenever scaffold is detected that is the class the read gets

classifyReads_one <- function(mut,
                              expedit,
                              rhapos,
                              scaffdetectwin=c(-2,+1)) {
  
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
    # rlab <- c(expedit=FALSE, mut=FALSE) # read labels; v1
    rlab <- c(expedit=FALSE, mut=FALSE, scaffold=FALSE) # v2
    
    ## gather all the detected mutations in this read
    muti <- mut %>%
      filter(rid==ui)
    
    # if mutation "ref" is present
    if('ref' %in% unique(muti$mutid)) {
      # there should not be anything else, check that
      if(nrow(muti)!=1)
        stop('\t \t \t \t >>> Error: for read', ui,', there is mutation "ref" and other mutations, which does not make sense.\n')
      # and we can return vector rlab as it is, i.e. FALSE, FALSE, FALSE
      return(rlab)
    } # we are done here as we returned, below is `else`
    # so everything below: read has at least one mutation
    # and we know "ref" is not there
    
    ## is prime-edit present?
    if(expedit %in% muti$mutid) {
      # then switch edit flag to TRUE
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
    
    ## detect scaffold incorporation
    # current definition is insertion or substitution that starts within scaffdetectwin
    # and altseq starts with G
    # convert window given by user to absolute positions
    # make sure both positions are positive
    scaffdetectwin <- abs(scaffdetectwin)
    scaffdetectpos <- (rhapos-scaffdetectwin[1]) : (rhapos+scaffdetectwin[2])
    # e.g. if rhapos is 103 and scaffdetectwin is (-2, +1),
    # will give 101, 102, 103, 104
    # do not look at whether expedit is there or not
    # but most of the time should be there if scaffold incorporation is present
    scaff <- muti %>%
      filter(start %in% scaffdetectpos) %>%
      filter(type %in% c('sub', 'ins')) %>%
      filter(startsWith(altseq, 'G'))
    if(nrow(scaff)>0) {
      rlab['scaffold'] <- TRUE
      # if scaffold is present, mut flag is guaranteed to be TRUE
      # check
      if(!rlab['mut']) stop('\t \t \t \t Error classifyReads_one: detected scaffold but "mut" flag for this read is FALSE, which does not make sense.\n')
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
    } else if(rl['expedit'] & !rl['mut'] & !rl['scaffold']) { # if edit and no mutation, no scaffold
      return('pure')
    } else if(rl['expedit'] & rl['mut'] & !rl['scaffold']) { # if edit and mutation, but no scaffold
      return('impure')
    } else if(!rl['expedit'] & rl['mut'] & !rl['scaffold']) { # if no edit, but mutation (no scaffold)
      return('mutated')
    } else if(rl['mut'] & rl['scaffold']) { # if mutated & scaffold (agnostic re edit)
      return('scaffold')
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
