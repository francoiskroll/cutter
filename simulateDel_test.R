#####################################################
# ~ cutter: test simulation of deletion ~
#
#
# Francois Kroll 2024
# francois@kroll.be
#####################################################

source('~/Dropbox/cutter/detectMHdel.R')


library(ggplot2)

mut <- read.csv('~/Dropbox/cutter/mutcallsTest.csv')
nrow(mut)

mutsim <- simulateDel(mut=mut,
                      nreads=1000,
                      mincov=100,
                      cutpos=87,
                      cutDelbp=3,
                      awayfromCut=4)
nrow(mutsim)
# added 1k rows, as expected


# histogram of deletion lengths -------------------------------------------

# what we want precisely is within each sample, frequency of this deletion out of all deletions
# because: final barplot is 100% = all deleted reads
# and because simulateDel simulates a sample of 100% deleted reads

# only keep deletions
del <- mutsim %>%
  filter(type=='del')

# add column sample_rid
# which gives a unique read id for the entire dataset
del <- del %>%
  mutate(srid=paste(sample, rid, sep='_'), .after='rid')

# are there duplicated reads?
# this would be a single reads with two (or more) deletions called
which(duplicated(del$srid))
# there are two reads
# which each have two deletions, so present in two rows each

# will just keep distrinct srid
# which means we lose the second deletion call for each of the duplicated reads
# it is so few reads, does not matter much
del <- del %>%
  distinct(srid, .keep_all=TRUE)
# 84045 >> 84043 rows
# seems correct

# calculate del_nreads per sample
# which is, per sample, total number of deleted reads
# below works because each row of delo represents a single read
deltal <- del %>%
  group_by(sample, .drop=FALSE) %>%
  tally(name='del_nreads')

# merge back to all deletions
del <- left_join(del, deltal, by='sample')

# simulated sample should have del_nreads 1000
del %>%
  filter(locus=='simulated')
# yes

# now we can calculate within sample, for each unique deletion, frequency of reads with this deletion out of all of deleted reads
del <- del %>%
  mutate(delpro=nreads/del_nreads)

# from this, estimate how many reads we would have with this deletion
# if every sample had 1k coverage and every read were deleted
# which is simply delpro * 1k
del <- del %>%
  mutate(del_co1k=round(delpro*1000))

# get the unique deletions within each sample
delu <- del %>%
  distinct(sample, mutid, .keep_all=TRUE)

# ! now, split if simulated or simulated

### simulated
# below: sim for simulated
delusim <- delu %>%
  filter(locus=='simulated')

# we make one big vector with all the deletion lengths based on those counts
lensim <- unlist(sapply(1:nrow(delusim), function(ro) {
  # from this row, get the deletion length & the normalised count
  # and return length, `normalised count` times
  # e.g. we get 8 bp deletion, and normalised count is 4
  # we return 8 8 8 8
  rep( delusim[ro, 'bp'] , delusim[ro, 'del_co1k'])
}))

# prepare the histogram data
lensimh <- freqBins(v=lensim,
                    every=5,
                    lower=1)

gghsim <- ggplot(lensimh, aes(x=upbound, y=freq)) +
  geom_col() +
  scale_x_continuous(breaks=lenoh$upbound) +
  coord_cartesian(ylim=c(0, 0.5)) +
  theme_minimal()
gghsim


### observed
# below: o for observed
deluo <- delu %>%
  filter(locus!='simulated')

# we make one big vector with all the deletion lengths based on those counts
lenso <- unlist(sapply(1:nrow(deluo), function(ro) {
  # from this row, get the deletion length & the normalised count
  # and return length, `normalised count` times
  # e.g. we get 8 bp deletion, and normalised count is 4
  # we return 8 8 8 8
  rep( deluo[ro, 'bp'] , deluo[ro, 'del_co1k'])
}))

# prepare the histogram data

# make the maximum match what we have with simulated
lenoh <- freqBins(v=lenso,
                  every=5,
                  lower=1,
                  last=max(lensimh$upbound))

ggho <- ggplot(lenoh, aes(x=upbound, y=freq)) +
  geom_col() +
  scale_x_continuous(breaks=lenoh$upbound) +
  coord_cartesian(ylim=c(0, 0.5)) +
  theme_minimal()
ggho

ggsave('~/Dropbox/cutter/dev_plots/dellens_histSim.pdf', gghsim, width=100, height=100, units='mm')
ggsave('~/Dropbox/cutter/dev_plots/dellens_histObv.pdf', ggho, width=100, height=100, units='mm')



# can run detection of MHdel ----------------------------------------------
# on all dataset, including simulated sample

mutsimmh <- detectMHdel(mut=mutsim,
                        minMHlen=1,
                        cutpos=87)

# keep only deletions & tally
delmh <- mutsimmh %>%
  filter(type=='del')


# tally, by sample, number of reads with evidence of MMEJ
mhtal <- delmh %>%
  group_by(sample, rundate, locus, sid, grp, well, splcov, MHside,
           .drop=FALSE) %>%
  tally(name='nreads')

# same but total number of reads with deletion
deltal <- delmh %>%
  group_by(sample,
           .drop=FALSE) %>%
  tally(name='del_nreads')

# merge both
mhtal <- left_join(mhtal, deltal, by='sample')


# calculate frequencies & plot --------------------------------------------

mhtal <- mhtal %>%
  mutate(catpro=nreads/del_nreads)

# exclude samples with less than 50 reads with deletion (arbitrary)
mhtal <- mhtal %>%
  filter(del_nreads>=50)
# note, there is none!



# plot MHside -------------------------------------------------------------

unique(mhtal$MHside)
mhtal$MHside <- factor(mhtal$MHside, levels=c('none', 'right', 'left'))

catcols <- c('#5a6974', '#EE7163', '#fcb505')

ggMh <- ggplot(mhtal, aes(x=sample, y=catpro, fill=MHside)) +
  facet_grid(~grp, scales='free_x', space='free') +
  geom_col(width=0.8) +
  scale_fill_manual(drop=FALSE, values=catcols) +
  theme_minimal() +
  theme(
    #legend.position='none',
    panel.grid.minor=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9, margin=margin(t=0, r=-1, b=0, l=0)),
    axis.text.x=element_text(size=7, angle=90, hjust=1, vjust=0.5, margin=margin(t=-10, r=0, b=0, l=0)),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=-1, b=0, l=0)))
ggMh

ggsave(here('~/Dropbox/ZPRI/240902_sept2024MiSeq/slc45a2_stdCas9/plots/mhdel_wSim.pdf'),
       width=150, height=100, units='mm')


# plot MH length ----------------------------------------------------------

### tally, by sample, number of reads with evidence of MMEJ
mhtal <- delmh %>%
  group_by(sample, rundate, locus, sid, grp, well, splcov, MHbp,
           .drop=FALSE) %>%
  tally(name='nreads')

# replace NA in MHbp by 0, easier to deal with later
mhtal[which(is.na(mhtal$MHbp)), 'MHbp'] <- 0

# same but total number of reads with deletion
deltal <- delmh %>%
  group_by(sample,
           .drop=FALSE) %>%
  tally(name='del_nreads')

# merge both
mhtal <- left_join(mhtal, deltal, by='sample')

# calculate frequencies & plot
mhtal <- mhtal %>%
  mutate(catpro=nreads/del_nreads)

# exclude samples with less than 50 reads with deletion (arbitrary)
mhtal <- mhtal %>%
  filter(del_nreads>=50)
# note, there is none with current dataset

# ready to plot
# set order
unique(mhtal$MHbp)
mhtal$MHbp <- factor(mhtal$MHbp, levels=c(0, 1, 2, 3, 4, 5, 6))

catcols <- c('#5a6974', '#f8dce3', '#f2c2ce', '#eda8b9', '#be8694', '#77545d', '#473237')

ggMhbp <- ggplot(mhtal, aes(x=sample, y=catpro, fill=MHbp)) +
  facet_grid(~grp, scales='free_x', space='free') +
  geom_col(width=0.8) +
  scale_fill_manual(drop=FALSE, values=catcols) +
  theme_minimal() +
  theme(
    #legend.position='none',
    panel.grid.minor=element_blank(),
    axis.title.x=element_blank(),
    # axis.title.y=element_text(size=9, margin=margin(t=0, r=-1, b=0, l=0)),
    axis.title.y=element_blank(),
    axis.text.x=element_text(size=7, angle=90, hjust=1, vjust=0.5, margin=margin(t=-10, r=0, b=0, l=0)),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=-1, b=0, l=0)))
ggMhbp

ggsave(here('~/Dropbox/ZPRI/240902_sept2024MiSeq/slc45a2_stdCas9/plots/mhdel_wSim.pdf'),ggMhbp,
       width=150, height=100, units='mm')