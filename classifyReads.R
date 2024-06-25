#####################################################
# ~ miseqUtils: function classifyReads ~
#
#
# Francois Kroll 2024
# francois@kroll.be
#####################################################

classifyReads <- function(mut,
                          expedit,
                          exportpath) {
  
  ### check export path ends with .csv
  frmt <- substr(exportpath, start=nchar(exportpath)-2, stop=nchar(exportpath))
  if(frmt!='csv')
    stop('\t \t \t \t >>> exportpath should end with .csv.\n')
  
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
  # we run function classifyReads_one on each mutation table one by one
  rlabL <- lapply(1:length(mutL), function(muti) {
    cat('\t \t \t \t >>> Sample', muti, 'out of', length(mutL), '\n')
    mut <- mutL[[muti]]
    return(classifyReads_one(mut=mut,
                             expedit=expedit))
  })
  
  ### gather in one dataframe
  rlab <- do.call(rbind, rlabL)
  # make sure it does not add row names
  row.names(rlab) <- NULL
  # export
  write.csv(rlab, exportpath, row.names=FALSE)
  # return
  invisible(rlab)
  
}
