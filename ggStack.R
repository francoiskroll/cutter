#####################################################
# ~ ZPRI: read labels to stack plot ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# heavily based on function ggFrameshift created for previous projects

# ggFrameshift(...)

# given ampliCan results (config_summary.csv), plot frameshift plots per sample
# in the style of Figure 2A https://elifesciences.org/articles/59683/figures#content

# ! expects sample ID to be exactly: gene.locus_scr/ko_num
# e.g. psen2.1_ko_8
# e.g. psen2.3_scr_1

# arguments

# amplican = full path to main results from ampliCan (typically config_summary.csv)
# or R dataframe of it
# should guess if given the path or the dataframe directly

# onlygene = plot all samples of that gene
# default 'all' = plot all samples
# ! currently only supports all genes or one

# onlysource = only plot samples of a certain source
# 'source' column added to amplican results file
# default 'all' = plot all samples
# if no source column, just leave default

# covFilter = whether or not to filter low-coverage samples
# default = TRUE, i.e. yes filter out low coverage samples (see below)

# mincov = minimum coverage (above or equal to), in number of pairs of reads
# default = 30 x (i.e. 30x and above are okay)
# pairs of reads, so 1 pair = 1x, not 2x (i.e. what ampliCan reports, but not IGV)
# exclude the sample from the plot if coverage is lower

# mutated_col = colour to use for mutated proportion
# default = #5a6974, grey

# frameshift_col = colour to use for frameshifted proportion
# default = #f1876b, orange

# exportOrNo = whether to export the plot to pdf or no
# default = TRUE = yes, export to pdf

# width = width of the pdf plot (mm)
# by default = 159.1 mm = A4 minus margins

# height = height of the pdf plot (mm)
# by default = 82.03 mm = 1/3rd of A4 minus margins

# exportfull = full output path for plot export


# functions & packages ----------------------------------------------------

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(stringr)

sem <- function(x) sd(x)/sqrt(length(x))


# function ggFrameshift(...) ----------------------------------------------

ggStack <- function(rcats,
                    onlygene='all',
                    onlylocus='all',
                    onlysource='all',
                    locusorder=NA,
                    covFilter=TRUE,
                    mincov=30,
                    mutated_col='#5a6974',
                    frameshift_col='#f1876b',
                    xtext=TRUE,
                    ytextOrNo=TRUE,
                    ynameOrNo=TRUE,
                    annotateOrNo=FALSE,
                    annotateSize=2,
                    exportOrNo=TRUE,
                    width=159.1,
                    height=82.03,
                    exportfull) {
  
  
  # import ampliCan results -------------------------------------------------
  if(is.character(amplican)) { # if string then assume we are being given the path
    amp <- read.csv(amplican)
  } else {
    amp <- amplican
  }
  # if not, we are given the dataframe so we do not have to import it
  
  # split the sample ID into its components ---------------------------------
  
  # make sure sample ID is character and not factor
  amp$ID <- as.character(amp$ID)
  
  col_locus <- as.character(lapply(strsplit(amp$ID, '_'), function(x) x[1])) # e.g. psen2.1_ko_8 >> psen2.1
  col_gene <- as.character(lapply(strsplit(col_locus, '\\.'), function(x) x[1])) # psen2.1 >> psen2
  col_locusnum <- as.character(lapply(strsplit(col_locus, '\\.'), function(x) x[2])) # psen2.1 >> 1
  col_grp <- as.character(lapply(strsplit(amp$ID, '_'), function(x) x[2])) # e.g. psen2.1_ko_8 >> ko
  
  # for sample_num, I have not been consistent
  # either psen2.1_ko_8 or psen2.1_ko8
  
  grp_spl <- sapply(amp$ID, function(id) {
    # which(...) finds the first _ (+1 for position just after)
    # then we take everything after
    return ( substr( id, start=which(strsplit(id, '')[[1]] == '_')[1]+1 , stop=nchar(id) ) )
    # psen2.1_ko_8 >> ko_8
    # psen2.1_ko8 >> ko8
  })
  
  # now get the grp
  col_grp <- str_extract(grp_spl, pattern=regex('[A-Z a-z]+')) # regex, any uppercase or lowercase letter, one time or more (+)
  # now get the sample number
  col_spl <- str_extract(grp_spl, pattern=regex('\\d+')) # regex, any digit, one time or more (+)
  
  # do another version of grp_spl, we want e.g. ko8 for all
  col_grpspl <- paste0(col_grp, col_spl)
  
  # add them to data
  amp <- amp %>%
    mutate(samplenum=col_spl, .after='ID') %>%
    mutate(grp=col_grp, .after='ID') %>%
    mutate(grpsamplenum=col_grpspl, .after='ID') %>%
    mutate(locusnum=col_locusnum, .after='ID') %>%
    mutate(locus=col_locus, .after='ID') %>%
    mutate(gene=col_gene, .after='ID')
  
  
  # keep only source we want to plot ----------------------------------------
  if (onlysource != 'all') {
    amp <- amp %>%
      subset(source==onlysource)
  }
  
  
  # keep only gene we want to plot ------------------------------------------
  if (onlygene[1] != 'all') {
    amp <- amp %>%
      subset(gene %in% onlygene)
  }
  
  # keep only locus we want to plot ------------------------------------------
  if (onlylocus[1] != 'all') {
    amp <- amp %>%
      subset(locus %in% onlylocus)
  }
  
  
  # filter low coverage samples ---------------------------------------------
  if (covFilter) {
    
    # find their positions
    lcis <- which(amp$Reads_Filtered < mincov) # low coverage indices
    
    # tell user we are throwing them out
    if (length(lcis) != 0) { # are there any samples to throw?
      for (i in 1:length(lcis)) {
        cat('\t \t \t \t >>> Excluding sample', as.character(amp[lcis[i], 'ID']), 'because its coverage is', as.numeric(amp[lcis[i], 'Reads_Filtered']), 'x \n')
      }
    }
    
    # exclude them
    amp <- amp %>%
      subset(Reads_Filtered >= mincov)
    
  }
  
  
  # calculate frameshift/nonframeshift --------------------------------------
  
  # total number of reads = Reads_Filtered
  # stacked barplot is in fact two barplots on top of each other
  # so one read cannot be in two categories at once (cannot be in both edited & frameshifted)
  # >> need to compute reads *not* frameshifted, and add that category on top of reads frameshifted
  # so total height will be reads frameshifted + reads mutated but not frameshifted = all reads mutated
  
  amp$edit <- amp$Reads_Edited/amp$Reads_Filtered # proportion of filtered reads that have edits
  amp$frameshift <- amp$Reads_Frameshifted/amp$Reads_Filtered # proportion of filtered reads that have frameshift
  # ! this is not the same as proportion of edited reads that have frameshift
  # this is because total height of the plot (100% of reads) should be = total number of reads (Reads_Filtered)
  # not total number of edited reads
  
  # proportion of reads edited but not frameshifted
  amp$nonframeshift <- (amp$Reads_Edited - amp$Reads_Frameshifted) / amp$Reads_Filtered
  # i.e. reads edited but not frameshifted (so left = reads mutated but indel is a multiple of 3), as a proportion of total number of reads
  
  
  # control order of samples in plot ----------------------------------------
  
  sampleorder <- as.character(amp[order(amp$locus, -rank(amp$grp)), 'ID'])
  # i.e. order first by locus (alphabetical order), then by grp (inversed alphabetical order, so scr then ko)
  amp$ID <- factor(amp$ID, levels=sampleorder)
  
  # set order of loci, if given by the user
  if(!is.na(locusorder[1])) {
    amp$locus <- factor(amp$locus, levels=locusorder)
  }
  
  
  # pivot to long format ----------------------------------------------------
  # for plotting
  
  ampl <- amp %>%
    select(ID, gene, locus, locusnum, grpsamplenum, grp, samplenum, frameshift, nonframeshift) %>%
    pivot_longer(cols=-c(ID, gene, locus, locusnum, grpsamplenum, grp, samplenum),
                 names_to='type',
                 values_to='pro')
  
  
  # control order of type in plot -------------------------------------------
  
  ampl$type <- factor(ampl$type, levels=c('nonframeshift', 'frameshift'))
  
  
  
  # prepare annotations -----------------------------------------------------
  
  # we want to write on the bars the column grpsamplenum, i.e. ko1, ko2, etc.
  # turn half of them (the non-frameshift) to NA, otherwise it will write them twice
  ampl[which(ampl$type=='nonframeshift'), 'grpsamplenum'] <- NA
  
  
  
  # plot --------------------------------------------------------------------
  
  framestack <- ggplot (data=ampl, aes (x=ID, y=pro, fill=type)) +
    facet_grid(~locus, scales='free_x', space='free') + # space='free' lets subplots being different width,
    # which results in bar width staying equal regardless of number of samples for each locus
    geom_col(width=0.8) +
    
    # annotate grp_samplenum on bar
    {if(annotateOrNo) (geom_text(aes(x=ID, y=0.03, label=grpsamplenum, colour=grp), angle=90, hjust=0, size=2))} +
    # should write KO in white (because usually in front of orange/grey), SCR in black (because usually in front of white background)
    {if(annotateOrNo) (scale_colour_manual(values=c('white', 'black')))} +
    
    scale_fill_manual(drop=FALSE, values=c(mutated_col, frameshift_col)) +
    theme_minimal() +
    theme(
      strip.text.x=element_blank(),
      panel.grid.major.x=element_blank(),
      panel.grid.minor=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=9, margin=margin(t=0, r=-1, b=0, l=0)),
      axis.text.x=element_text(size=7, angle=90, hjust=1, vjust=0.5, margin=margin(t=-10, r=0, b=0, l=0)),
      axis.text.y=element_text(size=7, margin=margin(t=0, r=-1, b=0, l=0)),
      legend.position='none') +
    coord_cartesian(ylim=c(0,1.0)) +
    scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                       labels=c(0, 25, 50, 75, 100)) +
    
    {if(ytextOrNo) ylab('% reads')} +
    {if(!ytextOrNo) ylab('')} +
    
    {if(!xtext) theme(axis.text.x=element_blank())} +
    
    {if(!ytextOrNo) theme(axis.text.y=element_blank())}
  
  print(framestack)
  
  
  
  # export plot -------------------------------------------------------------
  
  if (exportOrNo) {
    ggsave(exportfull, plot=framestack, width=width, height=height, unit='mm', device=cairo_pdf)
  }
  
  
  # summary statistics ------------------------------------------------------
  
  # will simply do summary statistics scr vs ko on the data we have at this stage
  # that means either user asked for a single gene and we will calculate summary statistics just for this gene;
  # or user asked for all samples to be plotted and we will calculate summary statistics for all samples
  
  # tell user summary stats
  cat('\n \n \t \t \t \t >>> SUMMARY STATS... \n')
  
  cat('\n \n \t \t \t by group \n')
  
  ampsum <- amp %>%
    select(ID, gene, locus, locusnum, grp, samplenum, edit, frameshift, nonframeshift) %>%
    pivot_longer(cols=-c(ID, gene, locus, locusnum, grp, samplenum),
                 names_to='type',
                 values_to='pro') %>%
    group_by(grp, type) %>%
    summarise_at(vars(pro),
                 list(
                   mean= ~ mean(.),
                   sd= ~ sd(.),
                   sem= ~ sem(.),
                   nspl= ~ length(.)
                 ))
  
  ampsum %>%
    print()
  
  # give also complete mean / sd (sometimes easier to round properly to report in text)
  cat('\t \t \t \t full means: ')
  print(ampsum$mean)
  cat('\n \t \t \t \t full sds: ')
  print(ampsum$sd)
  
  
  cat('\n \n \t \t \t by locus \n')
  
  ampsum <- amp %>%
    select(ID, gene, locus, locusnum, grp, samplenum, edit, frameshift, nonframeshift) %>%
    pivot_longer(cols=-c(ID, gene, locus, locusnum, grp, samplenum),
                 names_to='type',
                 values_to='pro') %>%
    group_by(locus, grp, type) %>%
    summarise_at(vars(pro),
                 list(
                   mean= ~ mean(.),
                   sd= ~ sd(.),
                   sem= ~ sem(.),
                   nspl= ~ length(.)
                 ))
  
  ampsum %>%
    print()
  
  # give also complete mean / sd (sometimes easier to round properly to report in text)
  cat('\t \t \t \t full means: ')
  print(ampsum$mean)
  cat('\n \t \t \t \t full sds: ')
  print(ampsum$sd)
  
  
  
  # return processed amplican results ---------------------------------------
  return(amp)
  
} # closes the function