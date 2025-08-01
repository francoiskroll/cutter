#####################################################
# ~ miseqUtils: function classifyReads ~
#
#
# Francois Kroll 2024
# francois@kroll.be
#####################################################

source('~/Dropbox/cutter/reversecomplement.R') # preciseClassify_one needs it

# v1

# v2: with classifyReads_one_v2.R which includes scaffold detection

# v3: added 'pestrand' argument ('forward' or 'reverse') to adjust scaffold detection.
# in v2, would expect starting with G, which is only correct if PE is happening on Forward strand
# see more notes in README.md

# v4: can turn ON or OFF scaffold detection with argument scaffDetect

# v5: takes rhapos from meta file

# v6: added argument unwantedSubs

classifyReads <- function(mut,
                          mode='precise',
                          scaffDetect=FALSE,
                          whichScaff='std',
                          scaffWin=NA,
                          unwantedSubs=FALSE,
                          exportpath) {
  
  ### check export path ends with .csv
  if(!is.na(exportpath)) {
    frmt <- substr(exportpath, start=nchar(exportpath)-2, stop=nchar(exportpath))
    if(frmt!='csv')
      stop('\t \t \t \t >>> exportpath should end with .csv.\n')
  }
  
  if(mode=='precise') {
    ### if scaffold detection is ON,
    # check rhapos is given in the mut table so we know where to look
    # check pestrand is giving in the mut table
    if(scaffDetect) {
      ## check rhapos
      if(! 'rhapos' %in% colnames(mut) ) stop('\t \t \t \t >>> Error classifyReads: scaffold detection is ON (scaffDetect=TRUE) but column "rhapos" is not in mut table so do not know where to look.\n')
      
      ## check pestrand
      if(! 'pestrand' %in% colnames(mut) ) stop('\t \t \t \t >>> Error classifyReads: scaffold detection is ON (scaffDetect=TRUE) but column "pestrand" is not in mut table. Make sure it is given in the meta.xlsx file when running callMutations.\n')
      if(! all(mut$pestrand %in% c('forward', 'reverse'))) stop('\t \t \t \t >>> Error classifyReads: in column "pestrand", only acceptable values are "forward" and "reverse".\n')
    }
  }
  
  ### import mut
  # if is a character, assume we are given a path,
  # if not, assume we are given R object directly
  if(is.character(mut)) {
    # check it is .csv
    frmt <- substr(mut, start=nchar(mut)-2, stop=nchar(mut))
    if(frmt!='csv')
      stop('\t \t \t \t >>> If giving a path for mut, it should end with .csv.\n')
    mut <- read.csv(mut)
  }
  
  ### split into mutations by sample
  # to get a list where each slot is one sample
  mutL <- split(mut, mut$sample)
  # then we pass each mutation table to classifyReads_one
  rlabL <- lapply(1:length(mutL), function(muti) {
    cat('\t \t \t \t >>> Sample', muti, 'out of', length(mutL), '\n')
    mut <- mutL[[muti]]
    
    ### mode PRECISE
    if(mode=='precise') {
      
      # we should only have one 'pestrand'
      if(length(unique(mut$pestrand))>1)
        stop('\t \t \t \t >>> Error classifyReads: more than one "pestrand" value within one sample, which does not make sense.\n')
      
      return( preciseClassify_one(mut,
                                  scaffDetect=scaffDetect,
                                  pestrand=unique(mut$pestrand),
                                  whichScaff=whichScaff,
                                  scaffWin=scaffWin,
                                  unwantedSubs=unwantedSubs) )
      
      ### mode FRAMESHIFT
    } else if(mode=='frameshift') {
      return( frameshiftClassify_one(mut=mut) )
    }
    
  })
  
  ### gather in one dataframe
  rlab <- do.call(rbind, rlabL)

  # make sure it does not add row names
  row.names(rlab) <- NULL
  # export
  if(!is.na(exportpath)) {
    write.csv(rlab, exportpath, row.names=FALSE)
    cat('\t \t \t \t >>> Wrote', exportpath, '\n')
  } else {
    cat('\t \t \t \t >>> Export is OFF.')
  }

  # return
  invisible(rlab)
  
}
