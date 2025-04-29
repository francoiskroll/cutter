#####################################################
# ~ ZPRI: mutations to read labels ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# expect mut, mutation table created by filterMutations

# copy of preciseClassify_one on 29/04/2025
# here, frameshiftClassify_one to classify indels when doing standard CRISPR

frameshiftClassify_one <- function(mut) {
  
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
    
    ## preallocate vector which is 1) indel?; 2) frameshift?
    # preallocate as FALSE, FALSE so = reference read
    rlab <- c(indel=FALSE, frameshift=FALSE)
    
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
    
    ## is there an indel?
    if( 'ins' %in% muti$type | 'del' %in% muti$type ) {
      # then switch indel flag to TRUE
      rlab['indel'] <- TRUE
      
      # and count whether it is a frameshift
      # make sure to look only at indels (not substitutions)
      indels <- muti %>%
        filter(type %in% c('ins', 'del'))
      
      # for deletions, set bp as negative lengths
      indels[which(indels$type=='del'), 'bp'] <- - indels[which(indels$type=='del'), 'bp']
      
      # now we can simply sum the lengths
      totlen <- sum(indels$bp)
      # and check whether it is a frameshift
      # i.e. should not be a multiple of three
      if(totlen%%3 != 0) {
        rlab['frameshift'] <- TRUE
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
    # row is:
    rl <- rlabs[r,]
    
    if(!rl['indel'] & !rl['frameshift']) {
      return('reference')
    } else if(rl['indel'] & !rl['frameshift']) { # if indel but no frameshift, indel_inframe
      return('indel_inframe')
    } else if(rl['indel'] & rl['frameshift']) { # if indel and frameshift, indel_frameshift
      return('indel_frameshift')
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
