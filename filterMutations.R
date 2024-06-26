#####################################################
# ~ ZPRI: filter mutations aligned by CRISPResso2 then listed with allelesToMutations.R ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################


# source & packages -------------------------------------------------------

library(tidyr)
library(dplyr)
library(tibble)


# function filterMutations() ----------------------------------------------

# runs on table created by function allelesToMutations()

# muttb = table of mutations created by function allelesToMutations()

# minnreads = minimum number of reads calling one mutation for it to be kept
  # e.g. minnreads = 3, and 5 reads have a 2-bp deletion at this position, so we keep that deletion

# cut_window = number of bp before/after cut site where mutation has to start for it to be kept
  # e.g. cut_window is 15 bp, so we consider a 30-bp window centered on cut site; if 2-bp deletion starts at cut site + 2 bp, we keep it
  # note, one might think this defeats the purpose of having long reads, but that is not quite right because the mutation itself can be long
  # e.g. if we imagine a 200-bp deletion that starts at cut site + 2 bp, it would have been impossible to detect it with short reads

# controltb = table of mutations for control sample(s)
  # if one = give as data frame
  # if multiple = give as list of data frames

# callsSubs = whether or not to call substitutions (i.e. keep them in output, assuming they pass the filters)

# (will do one for now, but want to allow multiple)

# defaults for filter is NA, i.e. filter is OFF

filterMutations <- function(muttb,
                            minnreads=NA,
                            cutpos=NA,
                            cutdist=NA,
                            rhapos=NA,
                            rhadist=NA,
                            controltb=NA,
                            callSubs=TRUE) {
  
  ### check settings
  # if given cut, need to give cutdist too
  # and vice-versa
  if(is.na(cutpos) & !is.na(cutdist))
    stop('\t \t \t \t >>> Error filterMutations: please give cutpos parameter, or set both cutpos and cutdist to NA to turn filter OFF.\n')
  
  if(!is.na(cutpos) & is.na(cutdist))
    stop('\t \t \t \t >>> Error filterMutations: please give cutdist parameter, or set both cutpos and cutdist to NA to turn filter OFF. \n')
  
  ### record total number of reads from muttb
  # this is calculated as sum of #Reads column in CRISPResso's Alleles
  # so will include reference reads & mutation we discarded here
  # sanity check: check it is just one total
  if(length(unique(muttb$splcov))!=1)
    stop('\t \t \t \t filterMutations: More than one value in splcov column, which does not make sense.\n')
  splcov <- unique(muttb$splcov)
  
  ### create unique mutation ID
  muttb <- muttb %>%
    mutate(mutid=paste(type, start, stop, bp, refseq, altseq, sep='_'), .before=1)
  # for "ref" mutation, better to just have mutid as "ref"
  muttb[which(muttb$mutid=='ref_NA_NA_NA_NA_NA'), 'mutid'] <- 'ref'
  
  ### keep ref mutations aside
  # we should never remove them
  # so keep them aside then we will put them back after
  tbref <- muttb %>%
    filter(mutid=='ref')
  
  # store the real mutations
  muttb <- muttb %>%
    filter(mutid!='ref')
  
  ### now we can start filtering
  cat('\t \t \t \t >>> Before filtering:', length(unique(muttb$mutid)),'unique mutations.\n')
  
  ### store the names of the mutated reads
  rmutids <- unique(muttb$rid)
  cat('\t \t \t \t >>> Before filtering:', length(rmutids),'mutated reads.\n')
  
  ### filter1: at least n reads calling the mutation
  # talr for tally of number of reads
  if(!is.na(minnreads) & nrow(muttb)!=0) { # also need to check we have any mutations to filter
    talr <- muttb %>%
      group_by(mutid) %>%
      tally()
    
    # keep these mutations
    mut2keep <- as.character(unlist(talr[which(talr$n >= minnreads), 'mutid']))
    muttb <- muttb[muttb$mutid %in% mut2keep , ]
    
    # left with
    cat('\t \t \t \t >>> After filter minimum number of reads:', length(unique(muttb$mutid)) ,'unique mutations left. \n')
  }
  
  
  ### filter2: cut_buffer
  # i.e. how close to the cut the indel has to start to be considered as potentially Cas9-made
  # ampliCan default is 5 bp
  # which I extended at some point to 12 bp
  # comment at the time (from oct2021MiSeq_main.R)
  # "I extended it after noticing a sample with a 1-bp deletion at 8 bp of the PAM (clearly made by Cas9 as not present in 3 SCR larvae)"
  
  # note, situation may be different here?
  # indels generated from nick repair, I think we do expect them to be all close to cut
  # but e.g. scaffold insertion, could be further?
  
  # e.g. cutbuffer 15 bp
  # we want *start* at or after 35 & before or at 65
  # stop can be anything; which gives us some of the advantages of long reads
  # in practice, the mutations are not precisely centered
  # see Fig 2.-suppl. 1 of my eLife paper
  # there are more of them towards spacer, upstream
  # so taking this into account could be a tiny bit more accurate
  
  if(!is.na(cutpos) & !is.na(cutdist) & nrow(muttb)!=0) { # also need to check we have any mutations left to filter
    cat('\t \t \t \t >>> Filtering by cut window ON. \n')
    # previous version, filter was looking whether start position was within window cut - cutdist until cut + cutdist
    # the issue with this solution is:
    # imagine a 60-bp deletion centered on the cut site
    # it effectively starts at cut - 30 bp, so if cutdist is 15 bp (a reasonable setting), it will get removed
    # while it is exactly centered on the cut!
    # >> need to consider whichever mutated position is the closest to the cut
    # i.e. consider every position from start to stop
    # note, for insertion, start is always = stop, so we just look where it starts
    
    # we loop through muttb
    dsb_keepornot <- sapply(1:nrow(muttb), function(nro) {
      
      # for this mutation, what are all the mutated positions
      mutposs <- muttb[nro, 'start'] : muttb[nro, 'stop']
      
      # substract cut position, so we get absolute distance from cut
      # then look if the closest mutated position is below (or equal to) cutdist
      # TRUE = some mutated position is less than `cutdist` away from cut
      # FALSE = some mutation position is less than `cutdist` away from cut
      return( min(abs(mutposs - cutpos)) <= cutdist )
      
    })
  }
  # we get vector dsb_keepornot which is TRUE/FALSE, i.e. keep that row or not
  

  ### if user wants to filter by end of RHA position
  # exactly same logic as for DSB position
  if(!is.na(rhapos) & !is.na(rhadist) & nrow(muttb)!=0) { # also need to check we have any mutations left to filter
    cat('\t \t \t \t >>> Filtering by RHA window ON.\n')
    # we loop through muttb
    rha_keepornot <- sapply(1:nrow(muttb), function(nro) {
      
      mutposs <- muttb[nro, 'start'] : muttb[nro, 'stop']
      
      return( min(abs(mutposs - rhapos)) <= rhadist )
      
    })
  }
  
  
  ### at this stage, we may or not have vectors dsb_keepornot & rha_keepornot
  # merge in one vector keepornot
  # if both DSB & RHA filters are on:
  if( exists('dsb_keepornot') & exists('rha_keepornot') ) {
    # we do OR operation, i.e. return TRUE if TRUE for DSB or RHA
    keepornot <- dsb_keepornot | rha_keepornot
  } else if( exists('dsb_keepornot') & !exists('rha_keepornot') ) {
    keepornot <- dsb_keepornot
  } else if( !exists('dsb_keepornot') & exists('rha_keepornot') ) {
    keepornot <- rha_keepornot
  }
  
  
  ### now ready to filter
  if(exists('keepornot')) {
    # check we have one boolean for each row
    if(length(keepornot)!=nrow(muttb)) stop('\t \t \t \t >>> Error: filtering by cut and/or RHA window: we expect one boolean for each row.\n')
    # remove those rows
    muttb <- muttb[which(keepornot),]
    # left with
    cat('\t \t \t \t >>> After filter cut window:', length(unique(muttb$mutid)) ,'unique mutations left. \n')
  }
  
  
  ### filter3: remove mutations also found in control
  # as these are frequent nanopore error at this position (e.g. TTTT, nanopore is likely to add/remove a T)
  # or mutations already in the genome before, so not caused by Cas9
  
  # simplest option would be to remove completely a mutation that is present in uninjected
  # in theory that is the best solution, as it is unlikely that a given mutation can be artefactual AND Cas9-mediated
  # but the possibility of contamination raises a concern: there may be some templates from injected in the uninjected sample
  # in this case, we would end up deleting all this mutation from the injected sample
  # e.g. injected sample has 10-bp deletion next to cut site, present in 15% of reads
  # if uninjected sample has 2 reads with this 10-bp deletion, this is very likely to be contamination
  # but solution above would delete that mutation from the injected sample entirely, so not good in practice
  
  # I think criterion should be based on ratio of frequencies
  # where ratio = frequency mutation injected / frequency mutation uninjected
  # i.e. how many times more frequently is the mutation found in injected compared to uninjected
  # options:
  # freq_ratio is = 1: means *same* frequency of this mutation in uninjected
    # I think certain this mutation is an artefact, we should exclude it
  # freq_ratio is < 1: means *more* of this mutation in uninjected
    # I think certain this mutation is an artefact, I do not see any scenario where uninjected sample would have *more* of a mutation
    # ! note if the goal is to repair a mutation, this may delete the mutation
    # e.g. reference is AA, mutated fish is GG, and our goal is to do GG>AA (back to wild-type)
    # then we expect mutation to be less present in injected, hence this filter will delete it
    # but here we are focusing on indel quantification
  # freq_ratio is > 1: means *less* of this mutation in uninjected;
    # there are 3 options
  # 1/ it is contamination from injected sample
      # in this case, we should ideally *not* edit the counts from the injected sample
      # how to recognise contamination?
      # if e.g. 10% of the uninjected reads are from the injected sample, this would give ratio > 10
      # e.g. if it is a real mutation and the frequency is 50% in injected, and 10% of the reads of the uninjected sample are from the injected sample
      # then this mutation's frequency in uninjected should be 50% * 0.1 = 5%, so ratio should be 50% / 5% = 10
      # I think 10% contamination is already high, but we can say > 5 to be safe
      # (which represents 20% of reads are from the other sample, i.e. very high contamination which should be obvious in IGV)
      # in summary: if ratio > 5, the much lower
  # 2/ it is artefact in both, just unlucky that mutation has less
    # here, we expect ratio to be fairly similar, so error rate 5x from above seems appropriate
  
  # 3/ same exact mutation can be generated by artefact and by Cas9
    # I think this is unlikely, but could be problematic:
    # in this case, injected should be higher but not by a lot
    # e.g. 2x would mean that in injected, 50% of these mutations are artefact, 50% are Cas9-generated
    # according to criterion above (we exclude mutation if ratio below 5x), we exclude this mutation entirely
    # I think this will remain a blindspot, unfortunately
  
  # if filter is ON, i.e. we are given a/some control sample mutation table(s)
  # if given one: it should be a data frame, jumps below
  # if given multiple: it should be a list of data frames, deal with this case now
  if( inherits(controltb, 'list') & nrow(muttb)!=0) {  # also need to check we have any mutations left to filter
    # note, data frame is also a list in R, so need to test if *actually* a list with inherits
    
    # ! next, we will count number of unique reads
    # so we need to make reads unique
    # currently, read 1.1 will be in control sample 1, and control sample 2, etc.
    # we will simply add a digit to the rid
    # e.g. 1.4.67 will mean control sample 1 . line 4 from allele table . read 73
    controltb <- lapply(1:length(controltb), function(i) {
      controltb[[i]] %>%
        mutate(rid=paste0(i, '.', rid))
    })
    # now pool
    controltb <- data.frame(do.call(rbind, controltb))
  }

  # so now in theory controltb can only be a single data frame or NA if filter is OFF
  if(is.data.frame(controltb)) {
    
    # record total number of reads in control sample
    covcon <- length(unique(controltb$rid))
    cat('\t \t \t \t >>> Total coverage control sample(s):', covcon, 'reads \n')
    
    # create unique mutation ID for control sample
    controltb <- controltb %>%
      mutate(mutid=paste(type, start, stop, bp, refseq, altseq, sep='_'), .before=1)
    
    # I think makes sense not to filter control,
    # but rather use it as a big database, and each mutation found in injected is like a query
    
    cat('\t \t \t \t >>> Control sample:', length(unique(controltb$mutid)),'unique mutations. \n')
    
    # tally mutations control
    talcon <- controltb %>%
      group_by(mutid) %>%
      tally(name='ncon')
    
    # add frequency control
    talcon <- talcon %>%
      mutate(freqcon=ncon/covcon, .after=ncon)
    
    # tally mutations injected
    talr <- muttb %>%
      group_by(mutid) %>%
      tally()
    
    # add frequency to injected counts
    talr <- talr %>%
      mutate(freq=n/cov, .after=n)
    
    # we can merge the control counts/frequencies to the sample counts/frequencies
    talr <- left_join(talr, talcon, by='mutid')
    
    # now calculate ratio of frequencies
    talr <- talr %>%
      mutate(freq_ratio=freq/freqcon)
    # freq_ratio is how many times more of this mutation in injected vs uninjected
    # e.g. 2x is "this mutation is found twice as often in injected vs uninjected"
    
    # now filter the sample's mutations
    # in summary from above, only mutations we want to keep:
    # NA, i.e. this mutation is not found in injected, which makes us more confident it is Cas9-generated
    # freq_ratio > 5, this mutation is at least 5x more frequent in injected than uninjected
    # which is consistent with real mutation from injected that contaminated in uninjected
    talr <- talr %>%
      filter(is.na(freq_ratio) | freq_ratio>5)
    
    # keep these mutations
    muttb <- muttb[muttb$mutid %in% talr$mutid , ]
    
    # left with
    cat('\t \t \t \t >>> After filter comparison to control sample:', length(unique(muttb$mutid)) ,'unique mutations left. \n')
    
  }
  
  
  ### should we keep substitutions?
  if(!callSubs) {
    
    muttb <- muttb %>%
      filter(type!='sub')
    
    # left with
    cat('\t \t \t \t >>> After removing substitutions:', length(unique(muttb$mutid)) ,'unique mutations left. \n')
  }
  
  ###### filtering is done
  
  ### record the mutated read IDs now
  rmutidsFIN <- unique(muttb$rid)
  cat('\t \t \t \t >>> After filtering:', length(rmutidsFIN),'mutated reads.\n')
  # so we filtered out those read IDs:
  filtids <- rmutids[which(!rmutids %in% rmutidsFIN)]
  if(length(filtids)>0) {
    cat('\t \t \t \t \t >>> so we filtered out', length(filtids),'reads; assuming those reads should be counted as reference.\n')
    # preparing a data frame with those read IDs
    ref2add <- data.frame(mutid='ref',
                          type='ref',
                          start=NA,
                          stop=NA,
                          bp=NA,
                          refseq=NA,
                          altseq=NA,
                          rid=filtids,
                          splcov=splcov)
    # adding those reference reads to the reference reads we already kept aside above
    tbref <- rbind(tbref, ref2add)
  }
  
  ### filtering is done, add back ref mutations
  muttb <- rbind(tbref, muttb)
  
  ### calculate frequencies
  # tally number of reads that has each mutation
  talr <- muttb %>%
    group_by(mutid) %>%
    tally(name='nreads') %>%
    mutate(freq=nreads/splcov)
  
  # "wish-list": could add information from control sample too
  # e.g. frequency in control sample
  
  # add that information to final mutation table
  muttb <- left_join(muttb, talr, by='mutid')
  
  
  ### return
  return(muttb)
  
}
