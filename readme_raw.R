
# setwd('C:/MJB/github/apcexplore')
setwd('D:/GitHub/multiAPCMC/R')
source('multiAPCMC.multifit.R')
source('multiAPCMC.reclist.R')
source('multiAPCMC.datashape.R')
source('multiAPCMC.singlefit.R')
source('multiAPCMC.multipred.R')
source('multiAPCMC.singlepred.R')
source('getmode.R')
source('multiAPCMC.multipredtest.R')
source('multiAPCMC.ranktable.R')
source('multiAPCMC.retrievemodel.R')
source('multiAPCMC.randpred.R')
source('multiAPCMC.predsummary.R')

library(ggplot2)
library('haven')
library('MASS')

# let's start by extracting one of the mock datasets and examining it:
multiAPCMC.example.data <- multiAPCMC.example.data

# this data has the following columns
head(multiAPCMC.example.data)
# cases are the incidence of some outcome (e.g. cancer diagnosis or mortality)
# and PY stands for person-years at risk in that year and age category
# notice also that the data is in long format!
# this is different from Nordpred, where the data is in a matrix resembling
# a Lexis diagram

# it has the following age and period categories
unique(multiAPCMC.example.data$Age)
unique(multiAPCMC.example.data$Period)

# that is, single year age categories
# and single year periods

# notice that Cohort has not yet been added
# we need the linear identity Cohort = Period - Age
# so to ensure that we have that, our functions will take care of that

# multiple classification models deal poorly with structural 0s.
# since incidences tend to be low for low ages, we could investigate
# at which age group we start having counts above 0.
multiAPCMC.example.data$cases[multiAPCMC.example.data$Age==0]
multiAPCMC.example.data$cases[multiAPCMC.example.data$Age==1]
multiAPCMC.example.data$cases[multiAPCMC.example.data$Age==2]
multiAPCMC.example.data$cases[multiAPCMC.example.data$Age==4]
multiAPCMC.example.data$cases[multiAPCMC.example.data$Age==6]
multiAPCMC.example.data$cases[multiAPCMC.example.data$Age==8]
multiAPCMC.example.data$cases[multiAPCMC.example.data$Age==20]
multiAPCMC.example.data$cases[multiAPCMC.example.data$Age==30]
# probably it is already fine at really young ages because there
# are hardly any 0s. But let's say that we see that we don't get
# structural 0s anymore at e.g. age 20.
# Then that becomes the age at which we start estimation (we don't have to
# throw out younger ages; the package still uses information from younger ages.
# Namely: it takes the mean over those ages and projects it forward).
# we set this age as our startestage
startestage <- 30
# note that this is the 30th age group (since 0 is the first, i.e. 1).
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
multiAPCMC.example.futuredata <- multiAPCMC.example.futuredata
# this is a dataset that contains the following information:
head(multiAPCMC.example.futuredata)
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

# in this case, the top 10 models are all quite close together, as they have RMSE scores
# that are very similar. With real data, and depending on the constraints chosen, this
# does not have to be the case.
# we see that the model with a power5 link, 21 periods of fit, first-second for period references
# and extremes for cohort references, is the top model

# note that AIC, pAIC, and driftSE are straight from the models fitted on the training data
# they are not determined on the validation data.

ranktable[382:392,]
# just to demonstrate: the RMSE of the worst models is a lot worse than the RMSE
# of the best models

# let's look at the data and the best model, plus its forecast
# for this, we need to extract the predictions of the best model from the giant model object
# it would be a lot of work to find this object
# so we have the multiAPCMC.retrievemodel() function
# in this function, we can just enter the parameters of the best model from ranktable, i.e.:

# there are two ways to do this. They use the 'retrievemodel' function:
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
plot(years[1:21],inctot[1:21],type='l', lwd=2,
     xlim=c(1989,2020),
     ylim=c(400,1800),
     main='Fit on treacherous data')
lines(years[21:32],inctot[21:32],lty=2, lwd=2)

# to get the same information from our best model's prediction, we can use
# the function multiAPCMC.predsummary().
# we just tell it which years we want a summary from, and which model
inctot.mod <- multiAPCMC.predsummary(years,rank1mod)
lines(years,inctot.mod$pred, col='red', lty=2,lwd=3)

legend(1990,1700,
       legend=c('observed data',
                'unobserved future',
                'model fit and forecast'),
       lty=c(1,2,2),
       lwd=2,
       col=c('black','black','red'))

# We see that our model fits really well for the training period
# this makes sense, our model was fitted on 21 years of the data (the training data)
# and  since it has Period as a factor variable, the model fits period-time exactly right
# However, the subsequent years (2010+) our model somewhat underpredicts initially
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

# I can also do this a different way. If I specifically want to keep the
# same cohorts and periods as references

# The other two use datashape and specifically say which cohorts we want
# as the reference categories, and other parameters
# then uses that datashape to do a single fit with mutliAPCMC.singlefit()
# what reference categories did it have?
ranktable$refper1st[1]
ranktable$refper2nd[1]
ranktable$refcoh1st[1]
ranktable$refcoh2nd[1]

rank1ds <- multiAPCMC.datashape(multiAPCMC.example.data,
                                noperiod=21,
                                startestage=startestage,
                                refper1st="1989",
                                refper2nd="1990",
                                refcoh1st="1942",
                                refcoh2nd="1979")
fitrank1mod2 <- multiAPCMC.singlefit(rank1ds,link="power5")
# this effectively re-fits the model. But since it is just 1 GLM model
# it goes pretty fast.

# and of course, I can also this using the ranktable object:

# let's fit a single instance of that model
rank1ds <- multiAPCMC.datashape(multiAPCMC.example.data,
                                  noperiod=ranktable$noperiod[1],
                                  startestage=startestage,
                                  refper1st=ranktable$refper1st[1],
                                  refper2nd=ranktable$refper2nd[1],
                                  refcoh1st=ranktable$refcoh1st[1],
                                  refcoh2nd=ranktable$refcoh2nd[1])
fitrank1mod2 <- multiAPCMC.singlefit(rank1ds,link=ranktable$link[1])

# the advantage of explicitly referring to the reference cohorts and periods
# rather than to something like "extremes" or "outer" is that these references
# do not change when the time series changes somewhat (e.g. when you add)
# more years of data. This is handy for when you have a training and
# validation set.



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
predtable$obs <- c(inctot[1:21],rep(NA,length(22:32)))
predtable$future <- c(rep(NA,length(1:21)),inctot[22:32])
# predtable <- predtable[predtable$years >= 2000,]

# we need ggplot for this
library(ggplot2)

# lets make the transparency reflect the certainty level
# but in a cumulative way
prea <- c(0.05,0.20,0.40,0.60,0.80,0.90)
a <- prea - c(0,prea[-length(prea)])*0.75 # the  0.75 is an inverse amplification value; lower it for stronger colors

# using alpha allows for layered colours
# AA=cumsum(a)
# names(AA)<-c("a1","a2","a3","a4","a5","a6")
p <- ggplot(predtable,aes(years))
p +
  geom_ribbon(aes(ymin=pred_025,ymax=pred_975, fill='confidence density'),alpha=a[1]) +
  geom_ribbon(aes(ymin=pred_100,ymax=pred_900, fill='confidence density'), alpha=a[2]) +
  geom_ribbon(aes(ymin=pred_200,ymax=pred_800, fill='confidence density'), alpha=a[3]) +
  geom_ribbon(aes(ymin=pred_300,ymax=pred_700, fill='confidence density'), alpha=a[4]) +
  geom_ribbon(aes(ymin=pred_400,ymax=pred_600, fill='confidence density'), alpha=a[5]) +
  geom_ribbon(aes(ymin=pred_450,ymax=pred_550, fill='confidence density'), alpha=a[6]) +
  geom_line(aes(y=obs,colour='observed data', linetype='observed data'),size=1.2) +
  geom_line(aes(y=future,colour="unobserved future", linetype='unobserved future'),size=1.2) +
  ggtitle("Fit on treacherous data") +
  scale_color_manual(name = "data", values = c("observed data" = "black",
                                                   "unobserved future" = "black")) +
  scale_fill_manual(name='model', values = c("confidence density" = "darkorange2")) +
  scale_linetype_manual(name='data',values=c("observed data" = "solid",
                                             "unobserved future" = "dashed")) +
  theme(legend.position = c(0.15, 0.8),
        legend.background=element_blank())


# something with discrete_alpha
# https://ggplot2.tidyverse.org/reference/scale_alpha.html
# if you want to say exact levels in the legend

# there are various papers out there on why 'leave-future-out' validation might be a good idea
# one reason is that it simply allows us to see how our model behaves compared to the truth

# another reason is that forecasting models have an implicit assumption that some component
# of the current or past will continue in the future (why else use data to predict the future?)
# usually, that component is the current trend
# by using leave-future-out validation, we can actually see which APC model creates the line
# that best fits the 'current trend' (the trend of the leave-future-out period)
# this is not normally clear in APCMC models because the period factors cause a near-perfect fit
# to the training data.
# and then our assumption could be either:
# (1) we continue extrapolating this trend to the future that we don't have data on yet.
# after all, we started out believing that the future is best predicted by extrapolating
# the current trend. So we are very transparent in our assumptions.
# or (2), we believe that the parameters that were used to fit this model on the training data
# are also the best parameters for fitting to the full data. So we re-fit the APC model
# to the complete data using these parameters, and then forecast into the unknown future.

# when counts are low to medium, I believe these are decent models
# however, counts (incidence, etc.) are high, I find that these types of models
# produce levels of uncertainty that are unrealistic. The truth might fall outside
# the confidence intervals or prediction intervals
# in that case, a gamma-poisson age-cohort-drift model might be better (see Bayesian
# structural time series)
# in other cases, for example when there is strong growth over time and we don't want to use
# e.g. the cuttrend parameter, we see explosive incidence
# these are all drawbacks of these types of models.

# in these, it may be worth looking into other types of forecasting methods
# such as Bayesian Structural Time Series



