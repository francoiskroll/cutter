#####################################################
# ~ miseqUtils: function callMutations ~
#
#
# Francois Kroll 2024
# francois@kroll.be
#####################################################

library(openxlsx)

# loops through CRISPResso output folders

# copath = folder that contains the CRISPResso results directories

# v1
# v2: expects columns rundate, well, locus to be in meta file, but will add any other column
# note, callMutations did not change when added scaffold detection

callMutations <- function(copath,
                          metapath,
                          minnreads=NA,
                          cutpos=NA,
                          cutdist=NA,
                          rhapos=NA,
                          rhadist=NA,
                          controltb=NA,
                          callSubs=TRUE,
                          exportpath) {
  
  if(!is.na(exportpath)) {
    ### check export path ends with .csv
    frmt <- substr(exportpath, start=nchar(exportpath)-2, stop=nchar(exportpath))
    if(frmt!='csv')
      stop('\t \t \t \t >>> exportpath should end with .csv.\n')
  }

  ### import meta file
  # check file exists
  if(!file.exists(metapath))
    stop('\t \t \t \t >>> Error: cannot find ', metapath, '.\n')
  # check it finishes by .xlsx
  frmt <- substr(metapath, start=nchar(metapath)-3, stop=nchar(metapath))
  if(frmt!='xlsx')
    stop('\t \t \t \t >>> Expecting meta file to be .xlsx.\n')
  # import it
  meta <- read.xlsx(metapath)
  
  # it has to have at least columns rundate & well & locus
  if(!'rundate' %in% colnames(meta))
    stop('\t \t \t \t >>> Error: expecting column "rundate" in meta file.\n')
  if(!'well' %in% colnames(meta))
    stop('\t \t \t \t >>> Error: expecting column "well" in meta file.\n')
  if(!'locus' %in% colnames(meta))
    stop('\t \t \t \t >>> Error: expecting column "locus" in meta file.\n')
  
  ### find CRISPResso output directories
  dirs <- list.dirs(copath)
  # first one is current directory, skip it
  dirs <- dirs[2:length(dirs)]
  # which ones are CRISPResso result directories?
  dirs <- dirs[startsWith(basename(dirs), 'CRISPResso')]
  
  ### find alleles table
  # loop through the directories,
  mutL <- lapply(1:length(dirs), function(di) {
    cat('\n')
    # path to Alleles_frequency_table.zip should be:
    alzip <- paste(dirs[di], 'Alleles_frequency_table.zip', sep='/')
    # check we found it
    if(!file.exists(alzip)) {
      cat('\t \t \t \t >>> Warning: no Alleles_frequency_table.zip in folder', dirs[di], '. Skipping this sample.\n')
      return()
    }
      
    # unzip it in same folder
    unzip(alzip, exdir=  dirname(alzip))
    # check unzipped file exists
    # should be same as zip, but with .txt
    altxt <- paste0(substr(alzip, start=1, stop=nchar(alzip)-3), 'txt')
    if(!file.exists(altxt))
      stop('\t \t \t \t >>> Error: after unzipping, expecting Alleles_frequency_table.txt in folder', dirs[di], '.\n')
    
    ### convert alleles table to mutation table
    cat('\t \t \t \t >>> Calling mutations from', altxt,'\n')
    muttb <- allelesToMutations(alpath=altxt)
    
    ### filter the detected mutations
    cat('\t \t \t \t >>> Filtering mutation calls.\n')
    mutf <- filterMutations(muttb=muttb,
                            minnreads=minnreads,
                            cutpos=cutpos,
                            cutdist=cutdist,
                            rhapos=rhapos,
                            rhadist=rhadist,
                            controltb=controltb,
                            callSubs=callSubs)
    
    ### add column well & column locus
    # from meta,
    # unique well names are:
    wells <- unique(meta$well)
    # unique locus names are:
    loci <- unique(meta$locus)
    # from the folder name,
    ## try to find well name
    # assumption: well name is in the folder name between some _
    dirsplit <- unlist(strsplit(basename(dirs[di]), '_'))
    wellnm <- dirsplit[which(dirsplit %in% wells)][1] # we take first occurence of that looks like well name
    if(length(wellnm)==0)
      stop('\t \t \t \t >>> Error: did not find well name in directory name', dirs[di], '\n')
    ## try to find locus name
    locnm <- dirsplit[which(dirsplit %in% loci)][1] # we take first occurence of that looks like locus name
    if(length(locnm)==0)
      stop('\t \t \t \t >>> Error: did not find locus name in directory name', dirs[di], '\n')
    ## check that it makes sense in comparison with meta
    # i.e. that this well coordinate has this locus
    metarow <- intersect(which(meta$well==wellnm), which(meta$locus==locnm))
    if(length(metarow)==0)
      stop('\t \t \t \t >>> Error: there is no row in meta file that has well ', wellnm, ' and locus ', locnm,'.\n')
    if(length(metarow)>1)
      stop('\t \t \t \t >>> Error: there are multiple rows in meta file that has well ', wellnm, ' and locus ', locnm,'. Please make well/locus unique.\n')
    
    ## add minimal meta information to mutation table
    # add from right to left
    mutf <- mutf %>%
      mutate(well=meta[metarow,'well'], .before=1) %>%
      mutate(locus=meta[metarow, 'locus'], .before=1) %>%
      mutate(rundate=meta[metarow,'rundate'], .before=1)
    # also add rundate_well_locus to use as unique sample ID
    # I think this should always be unique, even if pooling multiple runs
    mutf <- mutf %>%
      mutate(sample=paste(rundate, well, locus, sep='_'), .before=1)
    
    # now add other meta information found in meta file
    # this allows it to be any number of other columns so free to add new ones for specific experiments
    metacols <- colnames(meta)
    # the extra meta columns are:
    metacols <- metacols[which(!metacols %in% c('rundate', 'locus', 'well'))]
    
    # add extra columns
    # cannot make mutate work within an sapply loop, despite help from ChatGPT
    # but below works
    for(colnm in metacols) {
      mutf <- mutf %>%
        mutate(tmp=meta[metarow, colnm], .before='well')
      # now put the correct column name
      colnames(mutf)[which(colnames(mutf)=='tmp')] <- colnm
    }
    
    ### return filtered mutations
    return(mutf)
    
  })
  # we get list mutL which is table of filtered mutations for each sample
  # gather everything in one dataframe
  # we have locus & well to keep track of which mutation is from which sample
  mut <- do.call(rbind, mutL)
  
  # save this dataframe
  if(!is.na(exportpath)) {
    write.csv(mut, exportpath, row.names=FALSE)
    cat('\t \t \t \t >>> Wrote', exportpath, '\n')
  } else {
    cat('\t \t \t \t >>> Export is off.')
  }

  # we also return mut
  invisible(mut)
}
