#####################################################
# ~ miseqUtils/CUTTER: detect probable repair on mutations ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# here, detect 1-bp insertion of nucleotide next to cut
# see e.g. Fig 1b https://www.nature.com/articles/s41467-024-54566-7
# these disappear in absence of lig4, i.e. are a signature of NHEJ


# example -----------------------------------------------------------------

# mut <- read.csv('~/Dropbox/miseqUtils/mutcallsTest.csv')

# mutmh <- detectMHdel(mut=mut,
#                      minMHlen=2,
#                      cutpos=87)



# detectIns1bp ------------------------------------------------------------

# overarching function for detection of 1-bp ins
# runs detectIns1bp_one for each unique 1-bp insertion in mut (written by callMutations)

# 17/03/2025, deleted cutpos as argument, should take it from mut table


#' detectNHEJins
#'
#' @param mut
#' @param fwNGG
#'
#' @returns
#' @export
#'
#' @examples
detectNHEJins <- function(mut,
                          fwNGG) {

  ### loop through every row of mut
  mutinsL <- lapply(1:nrow(mut), function(rowi) {

    # case we want to find is 1-bp or 2-bp insertion at cutpos
    # that is nucleotide `cutpos` (or cutpos & cutpos - 1)
    # do we have to change something when reverse strand is cleaved?
    # alignment is always forward so I think no?
    # if insertion, 1 bp, at `cutpos`
    if(mut[rowi, 'type']=='ins' & mut[rowi, 'bp'] %in% c(1,2) & mut[rowi, 'start']==mut[rowi, 'cutpos']) {
      # get original reference sequence
      oref <- gsub('-', '', mut[rowi, 'ref'])
      # split it
      orefsi <- strsplit(oref, split='')[[1]]

      # get length of insertion (can only be 1 or 2 bp)
      bp <- mut[rowi, 'bp']

      # get cutpos
      cutpos <- mut[rowi, 'cutpos']

      # if forward strand has NGG PAM, we expect insertion to be the `cutpos` nucleotide
      # i.e. adjacent the cut on the left
      if(fwNGG) {

        if(bp==1) {
          expectnt <- orefsi[cutpos]
        } else if(bp==2) {
          expectnt <- orefsi[(cutpos-1):cutpos]
        }
        # if reverse, we expect insertion to be the `cutpos` nucleotide + 1
        # i.e. adjacent the cut on the right
      } else {
        if(bp==1) {
          expectnt <- orefsi[ cutpos+1 ]
        } else if(bp==2) {
          expectnt <- orefsi[ (cutpos+1):(cutpos+2) ]
        }
      }

      # make one string
      expectnt <- paste(expectnt, collapse='')

      # is the inserted nucleotide the expected nucleotide?
      if( mut[rowi, 'altseq'] == expectnt ) {
        return(NHEJins=TRUE)
      } else {
        return(NHEJins=FALSE)
      }
    } else {
      return(NHEJins=FALSE)
    }

  })
  mutins <- as.data.frame(do.call(rbind, mutinsL))
  colnames(mutins) <- 'NHEJins'

  # add NHEJ ins info to mut
  # check same number of rows
  if(nrow(mut) != nrow(mutins))
    stop('\t \t \t \t >>> Error detectNHEJ: not same number of rows when adding 1-bp ins information to mut.\n')

  mutins <- as.data.frame(cbind(mut, mutins))

  return(mutins)

}
