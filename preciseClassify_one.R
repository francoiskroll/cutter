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
# for now, detect scaffold just means insertion or substitution that starts with G
# ! definitions of classes are slightly different:
# pure = edit / no mutation / no scaffold
# impure = edit / mutation / no scaffold
# mutated = no edit / mutation / no scaffold
# scaffold = edit or not / mutation / scaffold
# so scaffold class is agnostic to edit or not, i.e. whenever scaffold is detected that is the class the read gets

# v3: correction of scaffold detection, taking into account PE strand

# v4: scaffDetect argument, to turn ON or OFF scaffold detection

# v5: changing name to preciseClassify_one
# will allow multiple modes in classifyReads, mode='precise' to detect classify reads when expected edit (e.g. prime editing)

# v6: added argument unwantedSubs
# changed slightly how detected unwanted mutations as a result

preciseClassify_one <- function(mut,
                                scaffDetect,
                                scaffdetectwin=c(-2,+1),
                                unwantedSubs=FALSE) {
  
  ### check mut is giving data for just one sample
  if(length(unique(mut$sample))>1)
    stop('\t \t \t \t >>> Error: attempting to give mutation data for more than one sample.\n')
  
  # scaffDetect: we checked in classifyReads that columns rhapos & pestrand look correct
  # we also turned expedit to NA if column "expedit" was not present
  
  # record rundate, sid, locus, well, grp sample, sample's coverage
  rundate <- unique(mut$rundate)
  sid <- unique(mut$sid)
  locnm <- unique(mut$locus)
  wellnm <- unique(mut$well)
  grpnm <- unique(mut$grp)
  splnm <- unique(mut$sample)
  # also record rhapos & pestrand (if scaffDetect is ON)
  if(scaffDetect) {
    rhapos <- unique(mut$rhapos)
    pestrand <- unique(mut$pestrand)
  }
  # for expedit, record it if present or turn to NA if not
  if('expedit' %in% colnames(mut)) {
    expedit <- unique(mut$expedit)
  } else {
    cat('\t \t \t \t >>> classifyReads: column "expedit" not in mut, will assume that there is no expected edit to detect.\n')
    expedit <- NAs
  }
  
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
    
    # note, read can still be called reference below
    # for example, a read that just has one substitution and no edit will be called reference (assuming unwantedSubs is FALSE)
    
    ## is prime-edit present?
    if(expedit %in% muti$mutid) {
      # then switch edit flag to TRUE
      rlab['expedit'] <- TRUE
    }
    
    ## other mutation(s) present?
    # can simply count number of rows,
    # after we remove expected edit (if it is present)
    # is expected edit present? if yes, remove it
    muti_unw <- muti # will store unwanted mutations
    if(rlab['expedit']) {
      muti_unw <- muti_unw %>%
        filter(mutid!=expedit)
    }
    
    # should we count substitutions as unwanted mutation or not?
    # if not, then exclude substitutions
    if(!unwantedSubs) {
      muti_unw <- muti_unw %>%
        filter(type!='sub')
    }
    
    # if any rows left in unwanted mutations, then we have some unwanted mutations
    # flip 'mut' flag to TRUE
    if(nrow(muti_unw)>0) {
      rlab['mut'] <- TRUE
    }
    # note, if read only has a substitution (no expected edit)
    # and unwantedSubs is FALSE
    # then all three flags are still FALSE and read will be called reference below
    
    
    ## detect scaffold incorporation
    if(scaffDetect) {
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
      
      ### notes in README
      # PE on Forward: insertion or substitution, whose `start` is close to rhapos and *starts* with G
      # PE on Reverse: insertion or substitution, whose `stop` is close to rhapos and *ends* with G
      if(pestrand=='forward') {
        scaff <- muti %>%
          filter(start %in% scaffdetectpos) %>%
          filter(type %in% c('sub', 'ins')) %>%
          filter(startsWith(altseq, 'G'))
        
      } else if(pestrand=='reverse') {
        scaff <- muti %>%
          filter(stop %in% scaffdetectpos) %>% # *** here stop, not start
          filter(type %in% c('sub', 'ins')) %>%
          filter(endsWith(altseq, 'C')) # *** here endsWith, not startsWith
      }
      
      if(nrow(scaff)>0) {
        rlab['scaffold'] <- TRUE
        # note v6: scaffold can be present and mut flag be FALSE if unwantedSubs is FALSE
        # in other words: unwantedSubs will exclude the scaffold incorporation (if it is only a substitution)
        # if nothing else in the read, then mut is FALSE
        # but scaffold detection can still detect the scaffold
      }
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
    # from row, only keep expedit / mut / scaffold flags
    # as a vector of booleans
    rl <- as.logical( rlabs[r, c('expedit', 'mut', 'scaffold')] )
    # add back the names
    names(rl) <- c('expedit', 'mut', 'scaffold')
    
    # no edit, no unwanted mutations, no scaffold = reference
    if(!rl['expedit'] & !rl['mut'] & !rl['scaffold']) {
      return('reference')
    # edit, no unwanted mutations, no scaffold = pure edit
    } else if(rl['expedit'] & !rl['mut'] & !rl['scaffold']) {
      return('pure')
    # edit, but also unwanted mutations (that is not scaffold) = impure edit
    } else if(rl['expedit'] & rl['mut'] & !rl['scaffold']) {
      return('impure')
    # no edit and unwanted mutations (that is not scaffold) = mutated
    # this category can also be called "indel" is unwantedSubs is FALSE
    } else if(!rl['expedit'] & rl['mut'] & !rl['scaffold']) { # if no edit, but mutation (no scaffold)
      return('mutated')
    # category "scaffold" trumps everything
    # if we see scaffold, we call that read scaffold, whatever the situation is regarding edit/other unwanted mutations
    } else if(rl['scaffold']) { # if scaffold (agnostic re edit or unwanted mutations)
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