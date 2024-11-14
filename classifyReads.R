#####################################################
# ~ miseqUtils: function classifyReads ~
#
#
# Francois Kroll 2024
# francois@kroll.be
#####################################################

# v1

# v2: with classifyReads_one_v2.R which includes scaffold detection

# v3: added 'pestrand' argument ('forward' or 'reverse') to adjust scaffold detection.
# in v2, would expect starting with G, which is only correct if PE is happening on Forward strand
# see more notes in README.md

# v4: can turn ON or OFF scaffold detection with argument scaffDetect
# TODO (maybe): could list events to detect

classifyReads <- function(mut,
                          expedit,
                          scaffDetect=FALSE,
                          rhapos,
                          pestrand='forward',
                          scaffdetectwin=c(-2,+1),
                          exportpath) {
  
  ### check export path ends with .csv
  if(!is.na(exportpath)) {
    frmt <- substr(exportpath, start=nchar(exportpath)-2, stop=nchar(exportpath))
    if(frmt!='csv')
      stop('\t \t \t \t >>> exportpath should end with .csv.\n')
  }
  
  ### if scaffold detection is ON, check rhapos is given so we know where to look
  if(scaffDetect) {
    if(is.na(rhapos)) stop('\t \t \t \t >>> Error classifyReads: when scaffold detection is ON (scaffDetect=TRUE), need to give rhapos.\n')
    # pestrand & scaffdetectwin have defaults so do not have to check
  }
  
  ### check pestrand is forward or reverse
  if(!pestrand %in% c('forward', 'reverse'))
    stop('\t \t \t \t >>> Error classifyReads: pestrand can only be "forward" or "reverse".\n')
  
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
    return(classifyReads_one(mut=mut,
                             expedit=expedit,
                             scaffDetect=scaffDetect,
                             rhapos=rhapos,
                             pestrand=pestrand,
                             scaffdetectwin=scaffdetectwin))
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
