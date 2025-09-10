#####################################################
# ~ ZPRI: mutations to read labels ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# expect mut, mutation table created by filterMutations

# v0
# based on preciseClassify_one.R v7
# but simple, just three categories: reference, edit, unwanted mutation (= indel if unwantedSubs=FALSE)

simpleClassify_one <- function(mut,
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
  
  ### get the original reference sequence
  # take all reference sequences that are not NA
  # remove any -
  orefs <- gsub('-', '', mut[!is.na(mut$ref), 'ref'])
  # check they are all the same
  if(length(unique(orefs))>1)
    stop('\t \t \t \t >>> Error preciseClassify_one: more than one original reference sequence in this sample\'s mutation table, which does not make sense.\n')
  oref <- unique(orefs)
  
  # for expedit, record it if present or turn to NA if not
  if('expedit' %in% colnames(mut)) {
    expedit <- unique(mut$expedit)
  } else {
    cat('\t \t \t \t >>> classifyReads: column "expedit" not in mut, will assume that there is no expected edit to detect.\n')
    expedit <- NA
  }
  
  # check only one sample's coverage
  splcov <- unique(mut$splcov)
  if(length(splcov)!=1)
    stop('\t \t \t \t >>> Error: 0 or multiple sample coverage.\n')
  
  
  ######################### below is lapply read by read within that sample
  
  # gather unique read IDs
  uids <- unique(mut$rid)
  # loop through unique read IDs
  rlabsL <- lapply(uids, function(ui) {
    
    # for debug
    # cat('ui is', ui, '\n')
    
    ## preallocate vector which is
    # 1) prime-edit present?; 2) other mutation present?
    # preallocate as FALSE, FALSE so = reference read
    rlab <- c(expedit=FALSE, mut=FALSE)
    
    ## gather all the detected mutations in this read
    muti <- mut %>%
      filter(rid==ui)
    
    ### if mutation "ref" is present ###
    if('ref' %in% unique(muti$mutid)) {
      # there should not be anything else, check that
      if(nrow(muti)!=1)
        stop('\t \t \t \t >>> Error: for read', ui,', there is mutation "ref" and other mutations, which does not make sense.\n')
      # and we can return vector rlab as it is, i.e. FALSE, FALSE
      return(rlab)
    } # we are done for this read as we called return, below is `else`
    # so everything below: read has at least one mutation
    # and we know "ref" is not there
    
    # note, read could still be called reference below
    # for example, a read that just has one substitution and no edit will be called reference (assuming unwantedSubs is FALSE)
    
    ### is prime-edit present? ###
    if(expedit %in% muti$mutid) {
      # then switch edit flag to TRUE
      rlab['expedit'] <- TRUE
    }
    
    ### other mutation(s) present? ###
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
    # note, if read only has substitution(s) (no expected edit)
    # and unwantedSubs is FALSE
    # then both flags are still FALSE and read will be called reference below
    
    ## return labels for this read
    return(rlab)
  }) # closes lapply to check read by read
  
  # create a simpler dataframe
  rlabs <- data.frame(do.call(rbind, rlabsL))
  
  # add back the read ids
  rlabs <- rlabs %>%
    mutate(rid=uids, .before=1)
  
  # we now have the read-by-read labels
  # assign categories based on those labels
  rcats <- sapply(1:nrow(rlabs), function(r) {
    # from row, only keep the flags, i.e. expedit / mut
    # just a vector of three booleans
    rl <- as.logical( rlabs[r, c('expedit', 'mut')] )
    # add back the names
    names(rl) <- c('expedit', 'mut')
    
    ### no edit, no unwanted mutations = REFERENCE
    if(!rl['expedit'] & !rl['mut']) {
      return('reference')
      
      #### *edit*, no unwanted mutations = EDIT
    } else if(rl['expedit'] & !rl['mut']) {
      return('edit')
      # in most cases this category will be = PURE EDIT from preciseClassify_one
      # as we decided below to call MUTATED even if edit is present
      
      ### *edit*, but also *unwanted mutations* = MUTATED
      # ! this is an arbitrary choice for simpleClassify_one
      # we decide that category MUTATED is over EDIT
      # i.e. if we have edit & a mutation, we ignore the edit and call the read mutated
    } else if(rl['expedit'] & rl['mut']) {
      return('mutated')
      
      ### no edit and *unwanted mutations* = MUTATED
      # if unwantedSubs is FALSE, this category can also be called INDEL
    } else if(!rl['expedit'] & rl['mut']) { # if no edit, but mutation
      return('mutated')
    }
  })
  
  # add the read categories to rlabs
  rlabs <- rlabs %>%
    mutate(cat=rcats)
  
  ### add meta columns to rlabs
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