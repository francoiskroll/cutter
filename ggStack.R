#####################################################
# ~ ZPRI: reads labels to stack plot ~
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


ggStack <- function(rlab,
                    onlylocus='all',
                    locusorder=NA,
                    catorder=c('reference', 'mutated', 'impure', 'pure'),
                    catcols=c('#aeb3b4', '#fcb505', '#a3bbdb', '#417dcd'),
                    mincov=NA,
                    xtextOrNo=TRUE,
                    ytextOrNo=TRUE,
                    ynameOrNo=TRUE,
                    exportOrNo=TRUE,
                    width=159.1,
                    height=82.03,
                    exportpath) {
  
  ### check exportpath is a pdf
  frmt <- substr(exportpath, start=nchar(exportpath)-2, stop=nchar(exportpath))
  if(frmt!='pdf')
    stop('\t \t \t \t >>> exportpath should end with .pdf.\n')
  
  ### import reads labels, if necessary
  if(is.character(rlab)) { # if string then assume we are being given the path
    rlab <- read.csv(rlab)
  }
  # if not a string, we are given the dataframe so we do not have to import it
  
  ### tally by sample and read category
  rtal <- rlab %>%
    group_by(sample, locus, well, splcov, cat) %>%
    tally(name='nreads')
  
  ### as check, separately, count total number of reads for each sample
  rtot <- rlab %>%
    group_by(sample) %>%
    tally(name='ntot')
  # add the total to tally from above
  rtal <- right_join(rtal, rtot, by='sample')
  # check that splcov is always = ntot
  if(!identical(rtal$splcov, rtal$ntot))
    stop('\t \t \t \t >>> Sample\'s coverage calculated at the start is not equal to adding reads from the different categories.\n')
  # can delete ntot
  rtal$ntot <- NULL
  
  ### calculate proportions
  rtal <- rtal %>%
    mutate(catpro=nreads/splcov)
  
  ### keep only locus we want to plot
  if (onlylocus[1] != 'all') {
    rtal <- rtal %>%
      filter(locus %in% onlylocus)
  }
  
  ### filter low coverage samples
  if (!is.na(mincov)) {
    
    # find their positions
    lcis <- which(rtal$splcov < mincov) # low coverage indices
    
    # tell user we are throwing them out
    if (length(lcis) != 0) { # are there any samples to throw?
      for (i in 1:length(lcis)) {
        cat('\t \t \t \t >>> Excluding sample', as.character(rtal[i, 'sample']),
            'because its coverage is', as.numeric(rtal[i, 'splcov']), 'x.\n')
      }
    }
    
    # exclude them
    rtal <- rtal %>%
      filter(splcov >= mincov)
    
  }
  
  ### control order of loci, if given by the user
  if(!is.na(locusorder[1])) {
    rtal$locus <- factor(rtal$locus, levels=locusorder)
  }
  
  ### control order of type in plot
  # categories will in stacked barplot from top to bottom
  rtal$cat <- factor(rtal$cat, levels=catorder)
  
  ### plot
  ggstack <- ggplot(data=rtal, aes (x=sample, y=catpro, fill=cat)) +
    facet_grid(~locus, scales='free_x', space='free') + # space='free' lets subplots being different width,
    # which results in bar width staying equal regardless of number of samples for each locus
    geom_col(width=0.8) +
    scale_fill_manual(drop=FALSE, values=catcols) +
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
    
    {if(ynameOrNo) ylab('% reads')} +
    {if(!ynameOrNo) ylab('')} +
    
    {if(!xtextOrNo) theme(axis.text.x=element_blank())} +
    {if(!ytextOrNo) theme(axis.text.y=element_blank())}
  
  print(ggstack)
  
  
  ### export plot
  if (exportOrNo) {
    ggsave(exportpath, plot=ggstack, width=width, height=height, unit='mm', device=cairo_pdf)
  }
  
  
  ### return processed tallies
  invisible(rtal)
  
} # closes the function