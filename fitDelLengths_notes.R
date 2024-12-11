#####################################################
# ~ cutter: notes about fitting of deletion lengths ~
#
#
# Francois Kroll 2024
# francois@kroll.be
#####################################################

### 09/12/2024
# currently only using a set of deletions at slc45a2 for test
# may try on other loci later

library(ggplot2)
library(MASS)


# import & get deletion lengths -------------------------------------------

mut <- read.csv('~/Dropbox/cutter/mutcallsTest.csv')

del <- mut %>%
  filter(type=='del')
# note, was taking *unique* deletions before (using dplyr's distinct())
# either in the entire dataset or within samples
# I now think this was a mistake, as we want to model the probabilities as the dataset is
# not (currently) study more general properties of repair
# e.g. if 1-bp deletions are very frequent, when we simulate a sample we should generate many 1-bp deletions
# so I do not think we should lower the frequency of 1-bp deletions because many (within a sample) may be from the same repair event
# this "error" should be random anyways so corrected with sufficient samples
# "error" being for example
# "we overestimate the frequency of 1-bp deletions because some samples happened to have many of those
# because the edit occurred very early, so then many cells ended up with that mutation, or there were just many amplicons with this mutation from chance"
# both sources of error are random, sometimes we will underestimate because mutation occurred late in development or those amplicons were just not sampled often, etc"
# having said this, it could be a good idea to normalise for coverage (depth)
# as in this case we do have access to this (random) error
# (being some samples have higher coverage than others, so give inflated counts),
# so we can actually do something about it

# we have column "freq" in mut which is what we need
# it is number of reads with this mutation out of sample's coverage
# [to explore deeper here: there may be a small error from reads which have multiple deletions, probably rare...]
# as it is now it is basically = mutation's count if total coverage was 1
# a bit counterintuitive
# will * 1000 and round
# so we have read count if total coverage as 1000

# would be a good idea to delete samples with low coverage
# so we do not extrapolate probabilities of deletion lengths from low sample size (number of reads)
min(unique(del$splcov))
# min is 1043 which is large
# so no sample to exclude

# normalise all coverage to total 1000 reads
# co1k for counts out of 1000x coverage
del <- del %>%
  mutate(co1k=round(freq*1000))

# now, we want unique deletions within each sample
# and we count it "co1k" times
delu <- del %>%
  distinct(sample, mutid, .keep_all=TRUE)
# we loop through these deletions,
# and each time repeat their lengths "co1k" times
# essentially we re-create the dataset as if every sample had exactly coverage 1000x
lens <- unlist(sapply(1:nrow(delu), function(ro) {
  # from this row, get the deletion length & the normalised count
  # and return length, `normalised count` times
  # e.g. we get 8 bp deletion, and normalised count is 4
  # we return 8 8 8 8
  rep( delu[ro, 'bp'] , delu[ro, 'co1k'])
}))
# unlist: we can throw all together, no need to keep track of samples etc.
# we now have all the deletion lengths, as if every sample had the same coverage


# plot distribution -------------------------------------------------------

# as density plot
ggplot(del, aes(x=bp)) + 
  geom_density()

# or as histogram
# ! uses freqBins from simulateDel.R
lenh <- freqBins(v=lens,
                 every=5,
                 lower=1)
# ! important, smallest possible deletion is 1 bp, not 0
# so we should not waste one value of a bin on 0 bp as it does not exist
# so bins will go 1--5 , 6--11, etc.

ggplot(lenh, aes(x=upbound, y=counts)) +
  geom_col() +
  scale_x_continuous(breaks=lenh$upbound)


# fit distribution --------------------------------------------------------

# can try all options offered by fitdistr (package MASS)
# uses MLE fitting
# https://www.rdocumentation.org/packages/MASS/versions/7.3-61/topics/fitdistr

# "start" does not mean lower bound but some initial values for the parameters
# so that fitdistr can then iterate from them
# some distributions do accept a "lower" argument to set lower bound
# which we know here is 1 (bp), the smallest possible deletion

# fbet <- fitdistr(lens, 'beta') # >> give up, asks for start and is only defined for 0-1 range anyways, so would not work
fcau <- fitdistr(lens, 'cauchy', lower=1)
# fchi <- fitdistr(lens, 'chi-squared', lower=1) # >> asks for start, give up
fexp <- fitdistr(lens, 'exponential', lower=1) # does not require start, cf. documentation
fgam <- fitdistr(lens, 'gamma', lower=1)
fgeo <- fitdistr(lens, 'geometric', lower=1)
flog <- fitdistr(lens, 'lognormal', lower=1) # "log-normal" gives same result
floi <- fitdistr(lens, 'logistic', lower=1)
fbin <- fitdistr(lens, 'negative binomial', lower=1)
fnor <- fitdistr(lens, 'normal', lower=1) # does not require start, cf. documentation
fpoi <- fitdistr(lens, 'Poisson', lower=1) # does not require start, cf. documentation
ft <- fitdistr(lens, 't', lower=1)
fwei <- fitdistr(lens, 'weibull', lower=1)

# there are 13 distributions
# can proceed by elimination
# after fitting: 11 left


# assess fit with density on histogram ------------------------------------

### cauchy
ggplot(lenh, aes(x=upbound, y=freq)) +
  geom_col() +
  scale_x_continuous(breaks=lenh$upbound) +
  stat_function(fun=dcauchy,
                args=list(location=fcau$estimate[1],
                          scale=fcau$estimate[2]),
                colour='red')

ggplot(del, aes(x=bp)) + 
  geom_density() +
  stat_function(fun=dcauchy,
                args=list(location=fcau$estimate[1],
                          scale=fcau$estimate[2]),
                colour='red')
# not bad?


### exponential
ggplot(lenh, aes(x=upbound, y=freq)) +
  geom_col() +
  scale_x_continuous(breaks=lenh$upbound) +
  stat_function(fun=dexp,
                args=list(rate=fexp$estimate[1]),
                colour='red')

ggplot(del, aes(x=bp)) + 
  geom_density() +
  stat_function(fun=dexp,
                args=list(rate=fexp$estimate[1]),
                colour='red')
# not bad but does not catch the lower frequency of 1-bp
# which I think may be real
# will exclude, exponential can never catch this property


### gamma
ggplot(lenh, aes(x=upbound, y=freq)) +
  geom_col() +
  scale_x_continuous(breaks=lenh$upbound) +
  stat_function(
    fun=dgamma,
    args=list(shape=fgam$estimate[1],
              rate=fgam$estimate[2]),
    colour='red'
  )

ggplot(del, aes(x=bp)) + 
  geom_density() +
  stat_function(
    fun=dgamma,
    args=list(shape=fgam$estimate[1],
              rate=fgam$estimate[2]),
    colour='red'
  )
# interesting!


### geometric
ggplot(lenh, aes(x=upbound, y=freq)) +
  geom_col() +
  scale_x_continuous(breaks=lenh$upbound) +
  stat_function(
    fun=dgeom,
    args=list(prob=fgeo$estimate[1]),
    colour='red'
  )

ggplot(del, aes(x=bp)) + 
  geom_density() +
  stat_function(
    fun=dgeom,
    args=list(prob=fgeo$estimate[1]),
    colour='red'
  )

# no, something is off
# will exclude


### log-normal
ggplot(lenh, aes(x=upbound, y=freq)) +
  geom_col() +
  scale_x_continuous(breaks=lenh$upbound) +
  stat_function(
    fun=dlnorm,
    args=list(meanlog=flog$estimate[1],
              sdlog=flog$estimate[2]),
    colour='red'
  )

ggplot(del, aes(x=bp)) + 
  geom_density() +
  stat_function(
    fun=dlnorm,
    args=list(meanlog=flog$estimate[1],
              sdlog=flog$estimate[2]),
    colour='red'
  )
# not bad, although it misses the largest peak in density plot


### logistic
ggplot(lenh, aes(x=upbound, y=freq)) +
  geom_col() +
  scale_x_continuous(breaks=lenh$upbound) +
  stat_function(
    fun=dlogis,
    args=list(location=floi$estimate[1],
              scale=floi$estimate[2]),
    colour='red'
  )

ggplot(del, aes(x=bp)) + 
  geom_density() +
  stat_function(
    fun=dlogis,
    args=list(location=floi$estimate[1],
              scale=floi$estimate[2]),
    colour='red'
  )
# interesting!


### negative binomial
ggplot(lenh, aes(x=upbound, y=freq)) +
  geom_col() +
  scale_x_continuous(breaks=lenh$upbound) +
  stat_function(
    fun=dnbinom,
    args=list(size=fbin$estimate[1],
              mu=fbin$estimate[2]),
    colour='red'
  )

ggplot(del, aes(x=bp)) + 
  geom_density() +
  stat_function(
    fun=dnbinom,
    args=list(size=fbin$estimate[1],
              mu=fbin$estimate[2]),
    colour='red'
  )
# no, it is a discrete distribution, exclude


### normal
ggplot(lenh, aes(x=upbound, y=freq)) +
  geom_col() +
  scale_x_continuous(breaks=lenh$upbound) +
  stat_function(
    fun=dnorm,
    args=list(mean=fnor$estimate[1],
              sd=fnor$estimate[2]),
    colour='red'
  )

ggplot(del, aes(x=bp)) + 
  geom_density() +
  stat_function(
    fun=dnorm,
    args=list(mean=fnor$estimate[1],
              sd=fnor$estimate[2]),
    colour='red'
  )
# not bad


### poisson
ggplot(lenh, aes(x=upbound, y=freq)) +
  geom_col() +
  scale_x_continuous(breaks=lenh$upbound) +
  stat_function(
    fun=dpois,
    args=list(lambda=fpoi$estimate[1]),
    colour='red'
  )

ggplot(del, aes(x=bp)) + 
  geom_density() +
  stat_function(
    fun=dpois,
    args=list(lambda=fpoi$estimate[1]),
    colour='red'
  )
# again, discrete, exclude


### t
ggplot(lenh, aes(x=upbound, y=freq)) +
  geom_col() +
  scale_x_continuous(breaks=lenh$upbound) +
  stat_function(
    fun=dt,
    args=list(df=ft$estimate[3]),
    colour='red'
  )

ggplot(del, aes(x=bp)) + 
  geom_density() +
  stat_function(
    fun=dt,
    args=list(df=ft$estimate[3]),
    colour='red'
  )
# no


### weibull
ggplot(lenh, aes(x=upbound, y=freq)) +
  geom_col() +
  scale_x_continuous(breaks=lenh$upbound) +
  stat_function(
    fun=dweibull,
    args=list(shape=fwei$estimate[1],
              scale=fwei$estimate[2]),
    colour='red'
  )

ggplot(del, aes(x=bp)) + 
  geom_density() +
  stat_function(
    fun=dweibull,
    args=list(shape=fwei$estimate[1],
              scale=fwei$estimate[2]),
    colour='red'
  )
# not bad

### left in the race:
# cauchy
# gamma
# log-normal
# logistic
# normal
# weibull

# so 6 of them


# QQplot ------------------------------------------------------------------

### cauchy
draws <- rcauchy(n=10000, location=fcau$estimate[1], scale=fcau$estimate[2])

qqplot(draws, lens)
abline(0,1)
# weird...

### gamma
draws <- rgamma(n=10000, shape=fgam$estimate[1], rate=fgam$estimate[2])

qqplot(draws, lens)
abline(0,1)
# looks like a very bad fit

### log-normal
draws <- rlnorm(n=10000, meanlog=flog$estimate[1], sdlog=flog$estimate[2])

qqplot(draws, lens)
abline(0,1)
# looks OKish

### logistic
draws <- rlogis(n=10000, location=floi$estimate[1], scale=floi$estimate[2])

qqplot(draws, lens)
abline(0,1)
# looks OKish

### normal
draws <- rnorm(n=10000, mean=fnor$estimate[1], sd=fnor$estimate[2])

qqplot(draws, lens)
abline(0,1)
# not bad but it draws negative values...
# I think might exclude for this reason
# or should use half a truncated gaussian?
# truncate at 1 could make sense...

### weibull
draws <- rweibull(n=10000, shape=fwei$estimate[1], scale=fwei$estimate[2])

qqplot(draws, lens)
abline(0,1)
# good

# will exclude cauchy,
# cannot make sense of qqplot


# KS test -----------------------------------------------------------------

# KS test is not happy because of ties in the data though...
# will ignore but may be an issue

### gamma
ks.test(lens, 'pgamma', shape=fgam$estimate[1], rate=fgam$estimate[2]) # D = 0.31 / p very low

### log-normal
ks.test(lens, 'plnorm', meanlog=flog$estimate[1], sdlog=flog$estimate[2]) # D = 0.175 / p very low

### logistic
ks.test(lens, 'plogis', location=floi$estimate[1], scale=floi$estimate[2]) # D = 0.179 / p very low

### normal
ks.test(lens, 'pnorm', mean=fnor$estimate[1], sd=fnor$estimate[2]) # D = 0.26 / p very low

### weibull
ks.test(lens, 'pweibull', shape=fwei$estimate[1], scale=fwei$estimate[2]) # D = 0.172 / p very low

# based on D distance statistic, Weibull is the best
# but fairly similar between log-normal / logistic / weibull
# best density plot was logistic I think

# from quick reading wikipedia, all three look fairly convincing
# log-normal looks widely used in biology

# ! logistic goes to negative, easy reason to exclude
# weibull (0.02) goes a bit lower than log-normal (0.18)
# will pick log-normal


# where it the peak of log-normal? ----------------------------------------

ggplot(del, aes(x=bp)) + 
  geom_density() +
  stat_function(
    fun=dlnorm,
    args=list(meanlog=flog$estimate[1],
              sdlog=flog$estimate[2]),
    colour='red'
  ) +
  geom_vline(xintercept=3.5) # by trial and error
# mode of (rounded) values drawn is 3

draws <- round(rlnorm(n=10000, meanlog=flog$estimate[1], sdlog=flog$estimate[2]))
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(draws)

# I am sure someone better in statistics could do a better job but this seems a good approximation
# so approach will be to
# fit the log-normal distribution based on the dataset (for one locus)
# draw deletion lengths at random
# etc.

# TODO
# note: I think it would make sense to fit the log-normal of lengths distributions ONCE on a big dataset with many loci
# distribution would look like Fig. 2C https://elifesciences.org/articles/59683
# ideally should check if log-normal is still the best we can do
# rather than data from a specific locus
# will leave as it is for now
# as we want to approximate what kind of lengths deletions can do
# will do this later