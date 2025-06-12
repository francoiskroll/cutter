#####################################################
# ~ ZPRI: reverse-complement function ~
#
#
# Francois Kroll 2025
# francois@kroll.be
#####################################################

reversecomplement <- function(seq) {
  
  # make sure it is a string all uppercase
  seq <- as.character(seq)
  seq <- toupper(seq)
  
  # split it
  seqsplit <- strsplit(seq, split='')[[1]]
  
  # reverse it
  seqrev <- rev(seqsplit)
  
  # run through it and get complementary of each nt
  seqrc <- sapply(1:length(seqrev), function(i) {
    if(seqrev[i]=='A') {
      return('T')
    } else if(seqrev[i]=='C') {
      return('G')
    } else if(seqrev[i]=='G') {
      return('C')
    } else if(seqrev[i]=='T') {
      return('A')
    } else if(seqrev[i]=='U') {
      return('A')
    }
  })
  
  # bind it all together
  seqrc <- paste0(seqrc, collapse='')
  
  # return
  return(seqrc)
}