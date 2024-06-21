#####################################################
# ~ ZPRI: mutations to read labels ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# expect mut, mutation table created by filterMutations

classifyReads <- function(mut) {
  # gather unique read IDs
  uids <- unique(mut$rid)
  # loop through unique read IDs
  rlabsL <- lapply(uids, function(ui) {
    
    ## preallocate vector which is 1) prime-edit present?; 2) other mutation present?
    # preallocate as FALSE, FALSE so reference read
    rlab <- c(edit=FALSE, mut=FALSE) # read labels
    
    ## gather all the detected mutations in this read
    muti <- mut %>%
      filter(rid==ui)
    
    ## is prime-edit present?
    if('sub_90_91_2_AA_GG' %in% muti$mutid) {
      # then switch edit label to TRUE
      rlab['edit'] <- TRUE
    }
    
    ## other mutation(s) present?
    # can simply count number of rows
    # if edit is present and > 1 rows, then other mutation
    # if edit is absent and > 0 rows, then other mutation
    if(rlab['edit'] & nrow(muti)>1) {
      rlab['mut'] <- TRUE
    } else if(!rlab['edit'] & nrow(muti)>0) {
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
    
    if(!rl['edit'] & !rl['mut']) {
      return('reference')
    } else if(rl['edit'] & !rl['mut']) {
      return('pure')
    } else if(rl['edit'] & rl['mut']) {
      return('impure')
    } else if(!rl['edit'] & rl['mut']) {
      return('mutated')
    }
  })
  
  rlabs <- rlabs %>%
    mutate(cat=rcats)
  
  return(rlabs)
  
}
