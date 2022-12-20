
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multiAPCMC

<!-- badges: start -->
<!-- badges: end -->

multiAPCMC is a package that allows the user to fit multiple
age-period-cohort multiple classification (MultiAPCMC) models and test
them against a ‘leave-future-out’ set. In short, the package helps you
investigate a large range of possible APC model specifications and helps
you determine which is the best one for forecasting the future.

Age-period-cohort (APC) models are models that seek to decompose
outcomes into parts that are associated with time since birth (age),
calendar time (period), and time of birth (cohort). This is [tricky
business](https://maartenbijlsma.com/2017/11/25/a-hot-topic-and-a-futile-quest-the-recent-discussion-on-age-period-cohort-analysis-february-17-2014/)
because age, period and cohort are linearly related; if exact age of a
person is known at a particular point in time, that person’s birth date
can be calculated. This results in a ‘linear dependency problem’; linear
models cannot find a unique solution for this decomposition, and will
either give an error or they will remove terms from the model by fiat.
Nevertheless, APC scholars try to find solutions that result in
‘meaningful’ APC estimates. This is done by putting some constraint on
the model. Whether this is actually possible to do in a meaningful way
is up for debate (hint: I am rather skeptical). See the commentaries in
the [December 2013 issue of
Demography](https://link.springer.com/journal/13524/volumes-and-issues/50-6?page=1)
for a detailed discussion. This is especially relevant when we want to
give a substantive interpretation to, for example, cohort trends. Do
some birth cohorts have a higher risk of some outcome (e.g. lung cancer)
because they performed some activity (e.g. smoking) during their
lifetimes, relative to other birth cohorts (conditional on age and
calendar time).

However, APC models are not just used for interpretation. They are also
used for forecasting. This is where the multiAPCMC package comes in. If
some cohorts indeed have a higher risk of (e.g.) lung cancer, then this
would be useful information to know when making age-specific forecasts
of future cancer incidence. With this idea in mind (I presume) the
[Nordpred
software](https://www.kreftregisteret.no/en/Research/Projects/Nordpred/Nordpred-software/)
(Moeller et al. 2003) was developed by the Cancer Registry
(Kreftregisteret) of Norway. This package fits APC models where the age,
period, and cohort are entered into a generalized linear model as
categorical variables, and a continuous time trend referred to as
‘drift’ is also determined. In order for the model to find a unique
solution, a constraint has to be placed on the model. The Nordpred
software automatically does this by setting a second reference category
on the period and cohort dimensions. Specifically, it sets the
coefficient of both the first and last period in the time series to be
equal (i.e. have the same ‘effect’), and the coefficient of the oldest
and youngest birth cohort. However, this influences determines the drift
parameter, and thereby it has a potentially strong effect on
extrapolating the trend towards the future (forecasting). However, this
is not the only possible constraint of its type. We could instead set
the middle periods or middle cohorts to be equal (perhaps adjacent
periods or cohorts are more similar…?), or the first and the middle, and
so forth. The multiAPCMC package allows you to investigate this. Other
choices that can be investigated are the link function (log link or
power5 link) and the number of periods to use for forecasting. The
multiAPCMC package also allows any data resolution, as long as they are
in Lexis squares (1 year age by 1 year period, 2 year age by 2 year
period, … 5 year age by 5 year period, etc.).

Below I provide a demonstration of the use of this package, which comes
with mock data.

## Installation

You can install the development version of multiAPCMC like so:

``` r
# Will later add here how to install the package.
```

## Example

``` r
source('R/getmode.R')
source('R/multiAPCMC.multifit.R')
source('R/multiAPCMC.reclist.R')
source('R/multiAPCMC.datashape.R')
source('R/multiAPCMC.singlefit.R')
source('R/multiAPCMC.multipred.R')
source('R/multiAPCMC.singlepred.R')
source('R/multiAPCMC.multipredtest.R')
source('R/multiAPCMC.ranktable.R')
source('R/multiAPCMC.retrievemodel.R')
source('R/multiAPCMC.randpred.R')
source('R/multiAPCMC.predsummary.R')

# let's start by extracting one of the mock datasets and examining it:
data(multiAPCMC.example.data)

# this data has the following columns
head(multiAPCMC.example.data)
#>   Age Period  PY cases
#> 1   0   1989 119     1
#> 2   1   1989  93     1
#> 3   2   1989  77     1
#> 4   3   1989  80     2
#> 5   4   1989  78     2
#> 6   5   1989 105     3
# cases are the incidence of some outcome (e.g. cancer diagnosis or mortality)
# and PY stands for person-years at risk in that year and age category
# notice also that the data is in long format!
# this is different from Nordpred, where the data is in a matrix resembling
# a Lexis diagram

# it has the following age and period categories
unique(multiAPCMC.example.data$Age)
#>  [1]  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
#> [26] 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49
#> [51] 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74
#> [76] 75 76 77 78 79 80 81 82 83 84
unique(multiAPCMC.example.data$Period)
#>  [1] 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003
#> [16] 2004 2005 2006 2007 2008 2009

# that is, single year age categories
# and single year periods

# notice that Cohort has not yet been added
# we need the linear identity Cohort = Period - Age
# so to ensure that we have that, our functions will take care of that

# multiple classification models deal poorly with structural 0s.
# since incidences tend to be low for low ages, we could investigate
# at which age group we start having counts above 0.
multiAPCMC.example.data$cases[multiAPCMC.example.data$Age==0]
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
multiAPCMC.example.data$cases[multiAPCMC.example.data$Age==1]
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
multiAPCMC.example.data$cases[multiAPCMC.example.data$Age==2]
#>  [1] 1 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
multiAPCMC.example.data$cases[multiAPCMC.example.data$Age==4]
#>  [1] 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 1 2 2
multiAPCMC.example.data$cases[multiAPCMC.example.data$Age==6]
#>  [1] 3 3 3 3 3 3 3 2 2 3 2 2 3 2 2 2 2 2 2 2 2
multiAPCMC.example.data$cases[multiAPCMC.example.data$Age==8]
#>  [1] 4 4 4 5 4 3 3 4 3 4 4 3 3 3 3 3 3 3 3 3 3
multiAPCMC.example.data$cases[multiAPCMC.example.data$Age==20]
#>  [1]  7  8  7  8  9 11 10 11 14 14 17 15 18 18 20 21 18 18 21 18 18
multiAPCMC.example.data$cases[multiAPCMC.example.data$Age==30]
#>  [1]  4  4  5  6  7  8  8  9 11 13 18 20 21 25 31 36 33 42 44 48 48
# probably it is already fine at age 10, but let's say
# then that becomes the age at which we start estimation (we don't have to
# throw out younger ages; the package still uses information from younger ages.
# Namely: it takes the mean over those ages and projects it forward).
# we set this age as our startestage
startestage <- 30
# note that this is the 11th age group (since 0 is the first, i.e. 1).
# So if we have 5-year age by 5-year period
# categories, '11' would refer to age category '50-55' and '1' would refer to
# category '0 to 4'. Be mindful of this.

# now we have to decide which model parameters we want to investigate
# let's look at a large range
# we want to look at both link functions (log and power5)
vec.link <- c('power5','log')
# we want to look at taking into account a number of different period ranges:
vec.noperiod <- c(5,10,15,21)
# and we want to investigate all possible period and cohort constraints:
vec.refper <- c('extremes','outer','center', 'first middle','middle last',
                'first second', 'penultimate last')
vec.refcoh <- c('extremes','outer','center', 'first middle','middle last',
                'first second', 'penultimate last')

# let's fit all possible model combinations. We use the function 'multiAPCMC.multifit()'
# for this:
fitlist <- multiAPCMC.multifit(data=multiAPCMC.example.data,
                               startestage=startestage,
                               vec.link=vec.link,
                               vec.noperiod=vec.noperiod,
                               vec.refper=vec.refper,
                               vec.refcoh=vec.refcoh,
                               nomodelmax=1024)
#> [1] "fitting 392 APCMC models"
# nomodelmax ensures that we can't possibly fit more than 1024 models
# just in case we have too many combinations and we are too lazy
# to check first by hand how many combinations we have

# the object 'fitlist' now contains the output of 392 models!
# I wouldn't recommend opening it unless you have to.
# some other functions will help when interacting with this object.

# first, let's make predictions for each of those models.
# for that, we need to have an object that contains the personyears (PY)
# for future age and period categories
# the age categories will be the same in the future
# but of course the calendar years will, by definition, have increased
# future person-years could come from the (projected) estimates of population size
# from a national statistics office
data(multiAPCMC.example.futuredata)
# this is a dataset that contains the following information:
head(multiAPCMC.example.futuredata)
#>      Age Period  PY cases
#> 1786   0   2010 281     1
#> 1787   1   2010 338     1
#> 1788   2   2010 273     1
#> 1789   3   2010 315     1
#> 1790   4   2010 284     1
#> 1791   5   2010 284     2
# Age, Period and PY have the same meaning as before
# the column 'cases' is not needed, since we might not know future cases
# (if we did, why are we forecasting?)
# but those can come in handy later when we want to test our models!

# we can also set a 'cuttrend' parameter; this is an attenuation parameter
# that reduces the predictions by some factor
# it was introduced by Nordpred because predictions from a model with
# a log-link can sometimes become explosive (since the inverse of log is exp)
# and this holds that in check.
# I set it to 0 here, but it should be considered when the drift (=overall trend)
# is strongly upwards.

# for that, we will use the function 'multiAPCMC.multipred()'
predlist <- multiAPCMC.multipred(multiAPCMC.multifit.object=fitlist,
                                 futuredata=multiAPCMC.example.futuredata,
                                 cuttrend=rep(0,15))
#> [1] "Projecting 11 periods into the future"
# this function predicts the future following all of the 392 models!

# since predlist is also a very large object, I don't suggest interacting
# directly with it.

# so which model is best? For that, we will use the 'multiAPCMC.ranktable()' function
# this function compares the predicted incidence following each model to the true incidence
# in order to do that, we need to fit our models to only a part of the time trend that
# we have access to. For example, here we had data from 1989 to 2020. What we then did
# was take the first 21 years (1989 to 2009) and use it as 'training data'
# and then use the last 11 years (2010 to 2020) and use it as 'validation data.'
# we will get into the reasoning behind that a bit later.

# here we compare the predicted incidence to true incidence using
# root mean squared error (RMSE), but there are also other measures we can
# rank by (they will still get calculated, but the list will be ranked by RMSE)
# the object multiAPCMC.example.futuredata is placed as oos.data (out of sample data)
# (out of sample refers to the 'not part of the sample that the models were trained on)
# validate.what is a parameter that asks which part of the data to calculate RMSE for
# the out of sample data (data at the end of the time series), data that we did not
# fit on that was at the beginning of the time series (e.g. if we have 20 years of
# training data, but noperiod was 10, then we have 10 years BEFORE the data used for fitting
# that we can also use, since it wasn't trained on that data) -referred to as 'not periodfitted',
# or all the data including # the data we fitted on (not recommended).
# Generally, 'out of sample' (oos) is recommended.

# let's use the root-mean
ranktable <- multiAPCMC.ranktable(multiAPCMC.multifit.object=fitlist,
                                     multiAPCMC.multipred.object=predlist,
                                     rankhow='RMSE',
                                     oos.data=multiAPCMC.example.futuredata,
                                     validate.what='oos')
# let's look at the top 10
ranktable[1:10,]
#>     linkfunction noperiod   ref_period       ref_cohort   is.AIC oos.RMSE
#> 47        power5       21 first second         extremes 4973.594 10.06673
#> 383       power5       21 first second penultimate last 4973.594 10.10041
#> 103       power5       21 first second            outer 4973.594 10.12313
#> 159       power5       21 first second           center 4973.594 10.16307
#> 271       power5       21 first second      middle last 4973.594 10.17529
#> 215       power5       21 first second     first middle 4973.594 10.21666
#> 384          log       21 first second penultimate last 4972.366 10.47037
#> 272          log       21 first second      middle last 4972.366 10.50664
#> 48           log       21 first second         extremes 4972.366 10.51214
#> 200          log       21 first middle     first middle 4972.366 10.58938
#>      oos.MAE     driftSE is.pcAIC
#> 47  6.935561 0.009189754 236.8378
#> 383 6.972091 0.027821230 236.8378
#> 103 7.010530 0.009207170 236.8378
#> 159 7.066605 0.016378968 236.8378
#> 271 7.032294 0.009225716 236.8378
#> 215 7.140847 0.009268184 236.8378
#> 384 6.981997 0.207986718 236.7794
#> 272 7.014819 0.082801164 236.7794
#> 48  7.022835 0.082639430 236.7794
#> 200 6.130001 0.020423009 236.7794

# in this case, the top 10 models are all quite close together, as they have RMSE scores
# that are very similar. With real data, and depending on the constraints chosen, this
# does not have to be the case.
# we see that the model with a power5 link, 21 periods of fit, first-second for period references
# and extremes for cohort references, is the top model

# note that AIC, pAIC, and driftSE are straight from the models fitted on the training data
# they are not determined on the validation data.

ranktable[382:392,]
#>     linkfunction noperiod       ref_period   ref_cohort   is.AIC  oos.RMSE
#> 287       power5       21         extremes first second 4973.594  57.88050
#> 335       power5       21 penultimate last first second 4973.594  60.08507
#> 295       power5       21            outer first second 4973.594  63.12696
#> 319       power5       21      middle last first second 4973.594  66.85555
#> 303       power5       21           center first second 4973.594  67.28408
#> 312          log       21     first middle first second 4972.366 149.67118
#> 288          log       21         extremes first second 4972.366 188.38806
#> 336          log       21 penultimate last first second 4972.366 197.98628
#> 296          log       21            outer first second 4972.366 218.83234
#> 320          log       21      middle last first second 4972.366 237.20965
#> 304          log       21           center first second 4972.366 245.36393
#>      oos.MAE    driftSE is.pcAIC
#> 287 18.50168 0.08070200 236.8378
#> 335 19.45506 0.08088434 236.8378
#> 295 20.78033 0.08071851 236.8378
#> 319 22.40688 0.08071897 236.8378
#> 303 22.59416 0.08100834 236.8378
#> 312 30.33190 0.81908359 236.7794
#> 288 38.46044 0.81917222 236.7794
#> 336 40.46773 0.82044045 236.7794
#> 296 44.79456 0.81929798 236.7794
#> 320 48.57901 0.81930289 236.7794
#> 304 50.24993 0.82182224 236.7794
# just to demonstrate: the RMSE of the worst models is a lot worse than the RMSE
# of the best models

# let's look at the data and the best model, plus its forecast
# for this, we need to extract the predictions of the best model from the giant model object
# it would be a lot of work to find this object
# so we have the multiAPCMC.retrievemodel() function
# in this function, we can just enter the parameters of the best model from ranktable, i.e.:

# there are two ways to do this. Manually:
rank1mod <- multiAPCMC.retrievemodel(multiAPCMC.multipred.object=predlist,
                                     what="pred",
                                     link='log',
                                     noperiod=15,
                                     refper='first second',
                                     refcoh='extremes')

# alternatively, we could just take the elements of the first row of the ranktable
# and put them in multiAPCMC.retrievemodel()
rank1par <- ranktable[1,1:4]
rank1mod <- multiAPCMC.retrievemodel(multiAPCMC.multipred.object=predlist,
                                     what="pred",
                                     link=rank1par[1],
                                     noperiod=rank1par[2],
                                     refper=rank1par[3],
                                     refcoh=rank1par[4])
# either way, we now have the model predictions in rank1mod

head(rank1mod)
#>   Age Cohort Period cases observed periodfitted        rate   rate_high
#> 1   1      0      1     1        1            1 0.008403361 0.008403361
#> 2   2     -1      1     1        1            1 0.010752688 0.010752688
#> 3   3     -2      1     1        1            1 0.012987013 0.012987013
#> 4   4     -3      1     2        1            1 0.025000000 0.025000000
#> 5   5     -4      1     2        1            1 0.025641026 0.025641026
#> 6   6     -5      1     3        1            1 0.028571429 0.028571429
#>      rate_low  PY pred_cases pred_cases_low pred_cases_high or.Per or.Age
#> 1 0.008403361 119          1              1               1   1989      0
#> 2 0.010752688  93          1              1               1   1989      1
#> 3 0.012987013  77          1              1               1   1989      2
#> 4 0.025000000  80          2              2               2   1989      3
#> 5 0.025641026  78          2              2               2   1989      4
#> 6 0.028571429 105          3              3               3   1989      5

# note that Age, Cohort and Period columns in this object are recoded and have
# the earlier-mentioned constraints imposed on them
# so don't interact with those
# if you want to interact with the original Age or Period categories,
# use or.Age and or.Per

# let's make a plot with the original data first.
# let's just look by year, so we add all the incidences in a year together
# summing over the age groups:
# first, let's plot the real values
alldat <- rbind(multiAPCMC.example.data,multiAPCMC.example.futuredata)
inctot <- NULL
years <- sort(unique(alldat$Period))
for(k in 1:length(years)) {
  inctot[k] <- sum(alldat$cases[alldat$Period==years[k]])
}
plot(years,inctot,type='l', lwd=2)

# to get the same information from our best model's prediction, we can use
# the function multiAPCMC.predsummary().
# we just tell it which years we want a summary from, and which model
inctot.mod <- multiAPCMC.predsummary(years,rank1mod)
lines(years,inctot.mod$pred, col='red', lwd=2)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
# We see that our model fits really well for the fitted period
# this makes sense, our model was fitted on 21 years of the data (the training data)
# and then the model fits exactly right
# for the subsequent years (2010+) our model somewhat underpredicts initially
# but it does understand that a decline is happening despite an initial trend upwards
# note that the curve we see here is rather tricky: most trends will not look like this
# empirically speaking. So I think we model does decently.
# note that if a noperiod was chosen that was less than the total training data
# then the data from calendartime -before- the training data can also be mispredicted

# but this is a bit of an ugly plot, and it also doesn't include uncertainty
# from our estimates
# any self-respecting forecast should include uncertainty.
# It is an inherently probabilistic exercise.

# How do we get forecast uncertainty?
# since I will need to aggregate over categories, and I cannot just sum lower bounds
# or sum upper bounds from the 'rank1mod' object.
# instead, I need to get the original model first:
# for this, I use multiAPCMC.retrievemodel() again, but now with 'what' set to 'fit'
# instead of 'pred'

fitrank1mod <- multiAPCMC.retrievemodel(fitlist,
                                        what="fit",
                                        link=rank1par[1],
                                        noperiod=rank1par[2],
                                        refper=rank1par[3],
                                        refcoh=rank1par[4])
# with this object, I can do a parametric bootstrap
# this means I will create many predictions; each one will be a random draw
# from the so-called 'estimator distribution' of the model.
# This is very similar to taking draws from a posterior distribution in
# Bayesian estimation (see e.g. Gelman 'Regression and other stories')
# or McElreath ('Statistical Rethinking')

# what I will do is I will draw from this distribution
# then immediately summarize the incidences using multiAPCMC.predsummary()
# and then save that into an object called bsdat

# poissondraw makes the process even more probabilistic:
# not only do we predict the mean and its standard error
# but also we take into account that count variables
# follow a poisson distribution
# with poissondraw=FALSE this helps produce confidence intervals
# with poissondraw=TRUE it creates prediction intervals
# but with the large numbers that we have here, the prediction interval
# is only a fraction larger

# so let's draw 499 many times and save the quantiles
bssize <- 499
bsdat <- multiAPCMC.predsummary(years,multiAPCMC.randpred(multiAPCMC.singlefit.object=fitrank1mod,
                                                          futuredata=multiAPCMC.example.futuredata,
                                                          noproj=11,
                                                          cuttrend=rep(0,11)),
                                poissondraw=TRUE)
for(b in 1:bssize) {

  bsdat <- rbind(bsdat,multiAPCMC.predsummary(years,multiAPCMC.randpred(fitrank1mod,
                                                                        futuredata=multiAPCMC.example.futuredata,
                                                                        noproj=11,
                                                                        cuttrend=rep(0,11)),
                                                                        poissondraw=TRUE))
}

# from this, let's store some quantiles (this can probably be done in a
# more efficient way!)
predtable <- multiAPCMC.predsummary(years,rank1mod)
predtable$pred_025 <- NA
predtable$pred_975 <- NA
predtable$pred_100 <- NA
predtable$pred_900 <- NA
predtable$pred_800 <- NA
predtable$pred_200 <- NA
predtable$pred_700 <- NA
predtable$pred_300 <- NA
predtable$pred_400 <- NA
predtable$pred_600 <- NA
predtable$pred_450 <- NA
predtable$pred_550 <- NA

for(y in sort(unique(predtable$years))) {

  qs <- quantile(bsdat$pred[bsdat$years==y],probs=c(0.025,0.975,0.10,0.90,0.20,0.80,0.30,0.70,0.40,0.60,0.45,0.55))
  predtable$pred_025[predtable$years==y] <- qs[1]
  predtable$pred_975[predtable$years==y] <- qs[2]

  predtable$pred_100[predtable$years==y] <- qs[3]
  predtable$pred_900[predtable$years==y] <- qs[4]

  predtable$pred_800[predtable$years==y] <- qs[5]
  predtable$pred_200[predtable$years==y] <- qs[6]

  predtable$pred_700[predtable$years==y] <- qs[7]
  predtable$pred_300[predtable$years==y] <- qs[8]

  predtable$pred_400[predtable$years==y] <- qs[9]
  predtable$pred_600[predtable$years==y] <- qs[10]

  predtable$pred_450[predtable$years==y] <- qs[11]
  predtable$pred_550[predtable$years==y] <- qs[12]

}
predtable$obs <- inctot
predtable <- predtable[predtable$years >= 2000,]

# we need ggplot for this
library(ggplot2)

# lets make the transparency reflect the certainty level
# but in a cumulative way
prea <- c(0.05,0.20,0.40,0.60,0.80,0.90)
a <- prea - c(0,prea[-length(prea)])*0.75 
# the  0.75 is an inverse amplification value; lower it for stronger colors

# plot with confidence levels
p <- ggplot(predtable,aes(years))
p +
  geom_ribbon(aes(ymin=pred_025,ymax=pred_975), fill='darkorange2',alpha=a[1]) +
  geom_ribbon(aes(ymin=pred_100,ymax=pred_900), fill='darkorange2',alpha=a[2]) +
  geom_ribbon(aes(ymin=pred_200,ymax=pred_800), fill='darkorange2',alpha=a[3]) +
  geom_ribbon(aes(ymin=pred_300,ymax=pred_700), fill='darkorange2',alpha=a[4]) +
  geom_ribbon(aes(ymin=pred_400,ymax=pred_600), fill='darkorange2',alpha=a[5]) +
  geom_ribbon(aes(ymin=pred_450,ymax=pred_550), fill='darkorange2',alpha=a[6]) +
  geom_line(aes(y=obs),size=1.2)
```

<img src="man/figures/README-example2-1.png" width="100%" />

``` r
# using alpha allows for layered colours
```

There are various papers out there on why ‘leave-future-out’ validation
might be a good idea one reason is that it simply allows us to see how
our model behaves compared to the truth.

Another reason is that forecasting models have an implicit assumption
that some component of the current or past will continue in the future
(why else use data to predict the future?). Usually, that component is
the current trend. By using leave-future-out validation, we can actually
see which APC model creates the line that best fits the ‘current trend’
(the trend of the leave-future-out period). This is not normally clear
in APCMC models because the period factors cause a near-perfect fit to
the training data.

Hence, our assumption could be either: - (1) we continue extrapolating
this trend to the future that we don’t have data on yet. after all, we
started out believing that the future is best predicted by extrapolating
the current trend. So we are very transparent in our assumptions. - (2),
we believe that the parameters that were used to fit this model on the
training data are also the best parameters for fitting to the full data.
So we re-fit the APC model to the complete data using these parameters,
and then forecast into the unknown future.

When counts are low to medium, I believe these are decent models.
However, counts (incidence, etc.) are high, I find that these types of
models produce levels of uncertainty that are unrealistic. The truth
might fall outside the confidence intervals or prediction intervals. In
that case, a gamma-poisson age-cohort-drift model might be better (see
Bayesian Structural Time Series). In other cases, for example when there
is strong growth over time and we don’t want to us# e.g. the cuttrend
parameter, we see explosive incidence these are all drawbacks of these
types of models. In these, it may also be worth looking into other types
of forecasting methods such as Bayesian Structural Time Series.