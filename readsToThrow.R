# function to loop through reads in a SAM file
# and for each read, check some filtering parameters

# !!! hard-coded path path2fastarefs !!!


# packages & small functions ----------------------------------------------

library(stringr)

substrEnding <- function(x, n){ # x = string (or vector of strings); n = last n characters
  substr(x, nchar(x)-n+1, nchar(x))
}



# function refSpan(...) ---------------------------------------------------

# from a row of a SAM file, gives the reference span
# e.g. this read spans 100 bp of the reference

refSpan <- function(samrow) {
  
  # get the read's CIGAR
  cigar <- as.character(samrow[6]) # ! assumes it is column 6
  
  cisp <- str_extract_all(cigar, '[0-9]+[A-Z]')[[1]] # splits the CIGAR into its components; e.g. 59S129M32S >> 59S 129M 32S
  
  # remove any hard or soft-clipping, we do not want to count this as spanning the reference
  # for example, 10H59S129M12I10D10H, split would be 10H 59S 129M 12I 10D 10H
  # we want to keep only 129M 12I 10D
  toremove <- which(sapply(cisp, function(sp) {substrEnding(sp, 1)}) %in% c('H', 'S'))
  
  # if any CIGAR element to remove, do so
  if (length(toremove) > 0) { # ! important to put in If condition, otherwise if nothing to remove it will delete everything
    cisp <- cisp[-toremove]
  }
  
  # now sum components to get read's span
  span <- sum(as.numeric(unlist(strsplit(cisp, '[A-Z]')))) # sum the number parts, e.g. 59 + 129 + 32 = 220; this gives the reference span
  span <- as.numeric(span) # just to be safe
  
  # return
  return(span)
}


# function refSpanPositions(...) ------------------------------------------

# from a row of a SAM file, gives the positions of the reference spanned by this read
# e.g. this read spans from 4th to 130th position of the reference

refSpanPositions <- function(samrow) {
  
  # get the left-most position of the alignment
  leftpos <- as.numeric(samrow[4]) # which should be the 4th column
  # this is the first position, i.e. where the alignment starts
  
  # second position is earliest position + reference span of this read
  rightpos <- leftpos + refSpan(samrow)
  
  # return both positions as a vector
  spanpos <- c(leftpos, rightpos)
  return(spanpos)
  
}


# function whichUpper(...) ------------------------------------------------
# given a string, returns positions of letters in uppercase

# smaller function isUpper(...)
# given a character, returns TRUE if uppercase, FALSE is lowercase
whichUpper <- function(string) {
  
  # split the string
  splits <- strsplit(string, split='')[[1]]
  
  # loop through each character and check if it is in uppercase
  isupper <- sapply(splits, function(char) {
    grepl("^[[:upper:]]+$", char)
  })
  
  # convert booleans into positions of uppercase characters
  upperpos <- as.numeric(which(isupper))
  
  # return these positions
  return(upperpos)
  
}


# function DSBpos(...) ----------------------------------------------------

# from a row of a SAM file,
# -- reads the fasta sequence used to align these reads
# -- gets the position of the PAM sequence, ! assumes the reference sequence is in lowercase, except the PAM in UPPERcase
# -- returns PAM position minus 4 as likely double-strand break position on the reference sequence

DSBpos <- function(refname, path2fastarefs) {
  
  # file should be same with .fa at the end
  refname <- paste0(refname, '.fa')
  
  # look for this sequence in folder with fasta references
  path2ref <- list.files(path2fastarefs, full.names=TRUE)[which(list.files(path2fastarefs)==refname)]
  
  if (length(path2ref)==0)
    stop('\t \t \t \t >>> Error: could not find fasta reference ', refname, '\n')
  
  # import the reference
  # this seems to be the most versatile import method, read.table would fail for formatting details
  refseq <- as.character(data.table::fread(path2ref)[1,])
  
  # get positions of uppercase characters
  pamposs <- whichUpper(refseq) # will return the 3 positions of the PAM
  
  # number of positions has to be a multiple of 3 (allows multiple PAMs)
  if (length(pamposs) %% 3 != 0) stop('\t \t \t \t >>> Number of nucleotides in uppercase is not a multiple of three, but each PAM should be three nucleotides? \n')
  
  # split positions of each PAM
  # how many PAMs do we have?
  npam <- length(pamposs)/3
  # split positions into PAM1, PAM2, ...
  pamposl <- split(pamposs, rep(1:npam, each=3))
  
  # for each PAM,
  # two possibilities
  # 1-- gRNA is +
  # then PAM is NGG and DSB site is at -4 on the left (position of N - 4)
  
  # 2-- gRNA is -
  # then PAM is CCN and DSB site is at +4 on the right (position of N + 4)
  
  dsbposs <- unlist(lapply(1:length(pamposl), function(pami) {
    
    cat('\n')
    cat('\t \t \t \t Looking at PAM', pami, 'of', length(pamposl), '\n')
    pamposs <- pamposl[[pami]]
    
    if (substr(refseq, start=min(pamposs)+1, stop=max(pamposs)) == 'GG') {
      cat('\t \t \t \t PAM is NGG so assuming gRNA is + \n')
      pampos <- min(pamposs) # get start of PAM, i.e. the N (leftmost position)
      dsbpos <- pampos - 4
      
    } else if (substr(refseq, start=min(pamposs), stop=max(pamposs)-1) == 'CC') {
      cat('\t \t \t \t PAM is CCN so assuming gRNA is - \n')
      pampos <- max(pamposs) # get start of PAM, i.e. the N (rightmost position)
      dsbpos <- pampos + 4
    }
    
    cat('\t \t \t \t PAM position on the reference = ', pampos, '\n')
    cat('\t \t \t \t Likely DSB site = ', dsbpos, '\n')
    
    return(dsbpos)
    
  }))
  
  return(dsbposs)
}



# function windowCoverageFilter(...) --------------------------------------
# given some alignments (rows of a SAM file), check the window coverage filter
# 'promiscuous' filter
# i.e. returns TRUE if *any* alignment in input is good, FALSE if *all* are bad

# input
# rows of SAM file to check as a dataframe (actual rows, not indices)
# window of interest

anyCoverageWindowValid <- function(samrows,
                                   woi) {
  
  # preallocate vector of TRUE/FALSE
  # TRUE = alignment covers window / FALSE = alignment does not cover window
  covorno <- rep(NA, nrow(samrows)) # covers or not
  
  # go through rows given one by one
  sapply(1:nrow(samrows), function(ri) { # ri for row index
    
    samrow <- samrows[ri,]
    
    spanstart <- refSpanPositions(samrow)[1]
    spanstop <- refSpanPositions(samrow)[2]
    
    if (woi[1] > spanstart & woi[2] < spanstop) { # if alignment *does* cover window
      # add TRUE to covorno
      covorno[ri] <<- TRUE
    } else { # if alignment does *not* cover window
      covorno[ri] <<- FALSE
    }
  })
  
  # check covorno is as expected
  if(length(which(is.na(covorno))) > 0) stop('\t \t \t \t >>> Error: some NA when double-checking window coverage of reads \n')
  
  # now look if any TRUE, if yes then just return TRUE
  if(length(which(covorno)) > 0) {
    return (TRUE)
  } else (
    return (FALSE) # if not, return FALSE, i.e. no valid alignment
  )
}

# function readNametoThrow() ----------------------------------------------

# given one row from a SAM file, checks whether we need to throw that read or not based on filters
# if no need to throw the read; does not do anything
# if need to throw the read; outputs its name

# Note; probably useful in the future: how to split a CIGAR in its component:
# str_extract_all('29S45M633H59S', '[0-9]+[A-Z]')

readNametoThrow <- function(samrow, # samrow = one row of a correctly-imported SAM file
                            min_phred,
                            min_readspan,
                            max_softprop,
                            woi) {
  
  ###
  
  # 1- check Phred score
  
  if (min_phred != 'off') {
    # get the alignment's Phred score
    phred <- as.numeric(samrow[5]) # ! assumes it is in column 5
    
    if (phred < min_phred) { # if it is below minimum Phred score, name that read
      cat('\t \t \t \t >>> Read', samrow[1], 'excluded because its Phred score is', phred, '\n')
      return(as.character(samrow[1]))
    } # Note; if Phred is below threshold, read name is output and function stops here due to return() call
    # i.e. below is only read if Phred is okay;
  }
  
  ###
  
  # 2- check the read's reference span
  
  if (min_readspan != 'off') {
    if (refSpan(samrow) < min_readspan) { # if read's span is below minimum length, name that read
      cat('\t \t \t \t >>> Read', samrow[1], 'excluded because its reference span is', refSpan(samrow), 'bp \n')
      return(as.character(samrow[1]))
    }
  }
  
  ###
  
  # 3- check the read's soft-clipping
  
  if (max_softprop != 'off') {
    
    # get the read's CIGAR
    cigar <- as.character(samrow[6]) # ! assumes it is column 6
    
    if(length(grep(pattern='S', cigar)) != 0) { # if read is soft-clipped; i.e. there is S in CIGAR
      tmp <- str_extract_all(cigar, '[0-9]+S')[[1]] # extracts the S components; eg. 59S129M32S >> 59S 32S
      softbp <- sum(as.numeric(unlist(strsplit(tmp, '[A-Z]'))))
      # split at any letter (here will be only S), so example above: 59 32
      # convert to numeric & sum
      
      # divide by read length to get proportion soft-clipped
      softp <- as.numeric(softbp/nchar(samrow[10]))
      
      if (softp > max_softprop) { # if proportion soft-clipped is above unwanted threshold; get the read name
        cat('\t \t \t \t >>> Read', samrow[1], 'excluded because', round(softp * 100), '% of it is soft-clipped \n')
        return(as.character(samrow[1]))
      }
    }
  }
  
  ###
  
  # 4- check if the read covers the window of interest for Cas9 experiments, which is DSB position ± padding
  
  if (woi[1] != 'off') {
    # window of interest has to be covered by the read
    # or in other words, left position of WOI has to be after alignment starts and right position has to be before alignment stops
    
    spanstart <- refSpanPositions(samrow)[1]
    spanstop <- refSpanPositions(samrow)[2]
    
    if (spanstart==0 & spanstop==0) { # if 0 to 0, just name that read to be removed
      
      cat('\t \t \t \t >>> Read', samrow[1],
          'excluded because only covers ', spanstart,
          'to',  spanstop, '\n')
      
      return(as.character(samrow[1]))
        
    } else if (! (woi[1] > spanstart & woi[2] < spanstop)) { # if read does not cover completely window
      # flag it for second check
      return(paste0('CHECK_', as.character(samrow[1])))
      }
    }
  
  # to get until here where nothing happens, read has to pass all filters above
  
}




# main function -----------------------------------------------------------


namesOfReadsToThrow <- function(sampath,
                                min_phred,
                                min_readspan,
                                max_softprop,
                                dsbpad,
                                outpath) {
  
  # check the path looks ok
  if(substrEnding(sampath, 4) != '.sam') stop('\t \t \t \t >>> That does not look like the path to a SAM file')
  
  # check the max_softprop argument looks ok
  if (max_softprop != 'off') {
    if(!(max_softprop > 0 & max_softprop < 1)) stop('\t \t \t \t >>> Maximum proportion soft-clipped must be 0~1, e.g. 0.2')
  }
  
  # if all good -- import the SAM file
  numcols <- max(count.fields(sampath, sep='\t'), na.rm=TRUE)
  sam <- read.table(file=sampath, fill=TRUE, comment.char='@', sep='\t',
                    col.names=sprintf('col%i', 1:numcols))
  
  # check there are alignments
  if (nrow(sam)==0) stop('\t \t \t \t >>> STOP: there is no alignments in this SAM file \n')
  
  # if we are filtering on coverage of window, get the positions now
  # so we do not compute the window of interest at each alignment
  # ! assumes that same reference sequence throughout the entire SAM file
  if (dsbpad != 'off') {
    
    # get the names of all reference sequences
    refs <- as.character(sam[,3])
    
    # should be either * if read unmapped or the name of the reference
    # exclude any *
    if(length(which(refs=='*'))>0) {
      refs <- refs[-which(refs=='*')]
    }
    
    # check there is only one reference name
    if(length(unique(refs)) == 0) stop('\t \t \t \t Error: there is no reference sequence in this SAM file \n')
    if(length(unique(refs)) > 1) stop('\t \t \t \t Error: there is more than one reference sequence in this SAM file \n')
    
    # if ok, take the name of that reference
    refname <- unique(refs)
    
    cat('\t \t \t \t >>> Reference is', refname, '\n')
    
    # get DSB position for this reference
    dsbpos <- DSBpos(refname=refname, path2fastarefs=path2fastarefs)
    # compute window of interest (WOI)
    # there may be multiple DSB positions (if multiple PAMs)
    # the solution is simply to take min & max, will work if one or multiple DSB positions
    woi <- c(min(dsbpos) - dsbpad, max(dsbpos) + dsbpad)
    cat('\n \t \t \t \t Window of interest = ', woi[1], ':', woi[2], '\n')
    
  } else { # if dsbpad is OFF
    woi <- 'off'
  }
  
  
  # apply readNametoThrow() function to all the reads
  # i.e. readNametoThrow runs on each row of the SAM file
  # if any filter is trigged, readNametoThrow returns the name of the read and apply() makes it run on the next row
  excludereads <- as.character(unlist(apply(sam[1:10], 1,
                                            readNametoThrow, min_phred, min_readspan, max_softprop, woi)))
  
  # if we are filtering based on window coverage, there is a chance some reads were flagged for second check
  # below: if window coverage is ON
  if (dsbpad != 'off') {
    
    # first, keep only unique readnames
    # if some reads were labelled CHECK by window coverage filter, this keeps only one
    # e.g. if both forward and reverse read were labelled CHECK,
    # it would be a waste of time to double-check both below, so keep only one time this readname
    excludereads <- unique(excludereads)
    
    # are there any reads flagged for second check?
    r2ch <- which(substr(excludereads, start=1, stop=6)=='CHECK_')
    # if need to check read, first six characters of name are CHECK_
    # r2ch will be indices of excludereads where we need to double-check
    
    if (length(r2ch) > 0) {
      
      cat('\n \n \t \t \t \t >>>', length(r2ch),' reads flagged for second check by window coverage filter \n \n')
      
      # now go through reads to check one by one
      sapply(r2ch, function(r2) {
        # r2 is an index of excludereads
        
        # note, if this read was already excluded for another reason, it is also a waste of time to double-check it
        # gets its name without CHECK_ in front (the original readname)
        rnam <- substr(excludereads[r2], start=7, stop=999)
        # is it already in excludereads? i.e. was it excluded because of another filter
        # if yes -- just exclude this read by re-writing its full name and move on
        if (rnam %in% excludereads) {excludereads[r2] <<- rnam}
        
        # if still needs to check, get all the alignments of this read
        rali <- which(sam[,1] == rnam) # rali for read's alignments
        # rali will give row indices of SAM file
        
        # check window coverage filter of these alignments
        anyok <- anyCoverageWindowValid(sam[rali,], woi=woi) # will return TRUE if any ok, FALSE if no
        
        if (anyok) { # if OK, switch to NA so we do not exclude this read
          excludereads[r2] <<- NA # do not remove the read here so we do not affect the indexing
          cat('\t \t \t \t \t >>> Read ', rnam, 'rescued as other alignment covers window \n')
        } else { # i not, switch to normal readname so we exclude this read
          excludereads[r2] <<- rnam
          cat('\t \t \t \t \t >>> Read ', rnam, 'removed, other alignments did not help \n')
        }
      })
      
    }

  }

  # remove any reads switched to NA above (if any created)
  if (length(which(is.na(excludereads))) > 0) {
    excludereads <- excludereads[-which(is.na(excludereads))]
  }

  # only keep unique names to be safe
  excludereads <- unique(excludereads)
  
  # say how many we are removing
  cat('\n \n \t \t \t \t >>> Excluding total', length(excludereads), 'reads \n \n \n')
  
  # write list of reads to throw as text file
  write.table(excludereads, file=outpath,
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  
}



# run function namesOfReadsToThrow ----------------------------------------

# small function to import each numerical argument
# if argument says 'off', keep it this way
# if argument is something else, assume it was correctly provided and turn it into a numeric
importNumArg <- function(arg) {
  if (arg!='off') {
    arg <- as.numeric(arg) # ! as.numeric() important, otherwise arguments are imported as character which creates odd outcomes
  }
  return(arg)
}

# first read the Terminal arguments
args <- commandArgs(trailingOnly=TRUE)

cat('\n \n readsToThrow.R got arguments: \n')

### path to SAM file ###
sampath <- as.character(args[1]) # read 1st argument from Terminal = path to SAM file
cat('\t Path to SAM file = ', sampath, '\n')

### path to fasta references ###
path2fastarefs <- as.character(args[2])
cat('\t Path to fasta references = ', path2fastarefs, '\n')

### Phred filter ###
min_phred <- importNumArg(args[3])
cat('\t Minimum Phred score = ', min_phred, '\n')

### Reference span filter ###
min_readspan <- importNumArg(args[4]) 
cat('\t Minimum reference span (bp) = ', min_readspan, ' \n')

### Soft-clipping filter ###
max_softprop <- importNumArg(args[5])
cat('\t Maximum proportion soft-clipped = ', max_softprop, '\n')

### Coverage of DSB site ###
dsbpad <- importNumArg(args[6])
cat('\t Minimum padding around DSB site = ', dsbpad, '\n')

### Output path ###
outpath <- as.character(args[7])
cat('\t Output =', outpath, '\n')


# run the main function namesOfReadsToThrow
namesOfReadsToThrow(sampath=sampath,
                    min_phred=min_phred,
                    min_readspan=min_readspan,
                    max_softprop=max_softprop,
                    dsbpad=dsbpad,
                    outpath=outpath)

###

# sampath <- '~/Dropbox/phd/220524_miseq/testfilter/G10_clu_3_example.sam'
# min_phred <- 'off'
# min_readspan <- 'off'
# max_softprop <- 0.2
# dsbpad <- 20
# outpath <- '~/Dropbox/phd/220524_miseq/testfilter/listreads.txt'
