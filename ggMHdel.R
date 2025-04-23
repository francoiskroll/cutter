#####################################################
# ~ miseqUtils/cutter: plot stacked barplot of deletions with microhomology ~
#
#
# Francois Kroll 2025
# francois@kroll.be
#####################################################

library(colorspace)
library(ggplot2)
library(dplyr)
library(tidyr)


# ggMHdel -----------------------------------------------------------------
# colours correspond to length of MH

ggMHdel <- function(mut,
                    min_del_nreads=50, # arbitrary
                    colourLight='#f1b9c7',
                    legendOrNo=TRUE,
                    titleOrNo=TRUE,
                    xtextOrNo=TRUE,
                    ytextOrNo=TRUE,
                    ynameOrNo=TRUE,
                    exportpath,
                    width=110,
                    height=65) {
  
  ### keep only reads with deletion
  del <- mut %>%
    filter(type=='del')
  
  ### tally, by sample & length of MH (MHbp), number of reads
  mhtal <- del %>%
    group_by(sample, rundate, locus, sid, grp, well, splcov, MHbp,
             .drop=FALSE) %>%
    tally(name='nreads')
  
  # replace NA in MHbp by 0, easier to deal with later
  mhtal[which(is.na(mhtal$MHbp)), 'MHbp'] <- 0
  
  # same tally but total number of reads with deletion
  deltal <- del %>%
    group_by(sample,
             .drop=FALSE) %>%
    tally(name='del_nreads')
  
  # merge both
  mhtal <- left_join(mhtal, deltal, by='sample')
  
  # calculate frequencies & plot
  mhtal <- mhtal %>%
    mutate(catpro=nreads/del_nreads)
  
  # exclude samples with less than min_del_nreads reads with deletion (arbitrary)
  mhtal <- mhtal %>%
    filter(del_nreads >= min_del_nreads)
  
  # set order of MHbp
  # we want 0 first, then decreasing
  # e.g. 0, 6, 4, 3, 2, 1
  MHbp_ord <- sort(unique(mhtal$MHbp), decreasing=TRUE)
  # remove 0 if present
  MHbp_ord <- MHbp_ord[-which(MHbp_ord==0)]
  # and add it as first
  MHbp_ord <- c(0, MHbp_ord)
  mhtal$MHbp <- factor(mhtal$MHbp, levels=MHbp_ord)
  
  # prepare the colours
  # user gave us the lightest colour
  # then we make it darker for each additional bp of MH
  MHbp_ord_no0 <- MHbp_ord[-which(MHbp_ord==0)]
  catcols <- sapply(MHbp_ord_no0, function(bp) {
    darken(col=colourLight, amount=0.08*bp)
  })
  # add colour for 0
  catcols <- c('#5a6974', catcols)
  
  #catcols <- c('#5a6974', '#c54867', '#db5072', '#e2738e', '#e996aa', '#f1b9c7')

  ggMhbp <- ggplot(mhtal, aes(x=sample, y=catpro, fill=MHbp)) +
    facet_grid(~grp, scales='free_x', space='free') +
    geom_col(width=0.8) +
    scale_fill_manual(drop=FALSE, values=catcols) +
    theme_minimal() +
    theme(
      strip.text=element_blank(),
      panel.grid.minor=element_blank(),
      axis.title.x=element_blank(),
      # axis.title.y=element_text(size=9, margin=margin(t=0, r=-1, b=0, l=0)),
      axis.title.y=element_blank(),
      # axis.text.x=element_text(size=7, angle=90, hjust=1, vjust=0.5, margin=margin(t=-10, r=0, b=0, l=0)),
      axis.text.y=element_text(size=7, margin=margin(t=0, r=-1, b=0, l=0)),
      axis.text.x=element_text(size=7, angle=90, hjust=1, vjust=0.5, margin=margin(t=-10, r=0, b=0, l=0))) +
    
    {if(!legendOrNo) theme(legend.position='none')} +
    {if(ynameOrNo) ylab('% of reads with deletion')} +
    {if(!ynameOrNo) ylab('')} +
    
    {if(!xtextOrNo) theme(axis.text.x=element_blank())} +
    {if(!ytextOrNo) theme(axis.text.y=element_blank())} +
    {if(!titleOrNo) theme(strip.text.x=element_blank())}
  
  print(ggMhbp)
  
  ggsave(exportpath, ggMhbp, width=width, height=height, units='mm')
  
}
