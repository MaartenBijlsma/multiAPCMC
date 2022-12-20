#'
#' @title multiAPCMC.singlepred
#'
#' @description makes a single set of predictions from a fitted multiAPCMC.singlefit object
#'
#' @param multiAPCMC.singlefit.object The list of fitted APCMC models that is produced by \code{multiAPCMC.multifit}.
#' @param futuredata a dataframe containing information on future periods and person-years (exposed population) in those years. Should have the same resolution as the data on which the models were also fitted (e.g. 1 year age group by 1 year periods, or 5-year age groups by 5-year periods, and so forth) and should have the columns Period, Age, and PY (person-years).
#' @param noproj number of periods to project to the future. Should have at most the length of the number of periods in futuredata.
#' @param cuttrend Cut trend in predictions? Should be a vector with at least a length equal to number of periods to predict, e.g. \code{c(0,0,0,0)} or \code{c(0,0.25,0.50,0.75)} or so if we project 4 periods into the future.
#'
#' @return returns a single set of predictions from an APCMC model.
#' @export
#'
#' @examples
#' multiAPCMC.example.data <- multiAPCMC.example.data
#' multifit.example <- multiAPCMC.multifit(multiAPCMC.example.data,
#'                                         startestage=5,
#'                                         vec.link=c('log','power5'),
#'                                         vec.noperiod=c(5,10),
#'                                         vec.refper=c('outer','center'),
#'                                         vec.refcoh=c('outer'),
#'                                         nomodelmax=1024)
#' # if multiAPCMC.ranktable shows that the model with a power5 link function
#' # a no-period of 10, period references of center and cohort references of 'outer'
#' # are best, then I could then retrieve this model as follows:
#' chosenmodel <- multiAPCMC.retrievemodel(multifit.example,
#'                          what="fit",
#'                          link='power5',
#'                          noperiod=10,
#'                          refper='center',
#'                          refcoh='outer')
#' multiAPCMC.singlepred(chosenmodel,
#'                     multiAPCMC.example.futuredata,
#'                     noproj=4,
#'                     cuttrend=rep(0,15))
multiAPCMC.singlepred <- function(multiAPCMC.singlefit.object,futuredata,noproj=4,cuttrend) {

  apcdata <- multiAPCMC.singlefit.object$apcdata
  startestage <- multiAPCMC.singlefit.object$startestage
  noperiod <- multiAPCMC.singlefit.object$noperiod
  refcoh1st <- multiAPCMC.singlefit.object$refcoh1st
  refcoh2nd <- multiAPCMC.singlefit.object$refcoh2nd
  refper1st <- multiAPCMC.singlefit.object$refper1st
  refper2nd <- multiAPCMC.singlefit.object$refper2nd

  dnoperiods <- length(unique(apcdata$PeriodTrue))

  # save some information to be able to later reconstruct
  # or.Per and or.Age
  ## reconstruct or.Age and or.Per
  # determine intercept for or.Per
  or.Per.int <- min(apcdata$or.Per,na.rm=TRUE)
  # determine intervalsize
  apcdata <- apcdata[order(apcdata$Age, apcdata$PeriodTrue),]
  or.Per.is <- getmode(diff(apcdata$or.Per))
  # determine intercept for or.Age
  or.Age.int <- min(apcdata$or.Age,na.rm=TRUE)
  # determine intervalsize
  apcdata <- apcdata[order(apcdata$PeriodTrue, apcdata$Age),]
  or.Age.is <- getmode(diff(apcdata$or.Age))

  apcdata$periodfitted <- ifelse(apcdata$pcode>(dnoperiods-noperiod),1,0)
  apcdata.fitted <- apcdata[apcdata$Age>=startestage,]
  apcdata.fitted <- apcdata.fitted[apcdata.fitted$pcode>(dnoperiods-noperiod),]

  # first we need to determine the number of age categories
  Age <- sort(unique(apcdata.fitted$Age))
  # then we need to know how many periods ahead and what name those periods should have:
  PeriodTrue <- max(apcdata$pcode)+1:noproj
  # then we need to determine how to attenuate those pcodes
  if(noproj > length(cuttrend)) {
    stop("cuttrend length is shorter than number of periods to project ahead")
  }
  pcode <- max(apcdata.fitted$pcode)+cumsum(1-cuttrend[1:noproj])

  # now that we have A and pcode, we can expand the grid
  # and then we can calculate cohort
  APtrue <- expand.grid(Age,PeriodTrue)
  AP <- expand.grid(Age,pcode)
  AP <- data.frame(cbind(AP,APtrue)[,-3])
  names(AP) <- c("Age","pcode","PeriodTrue")
  # pcode is now already attenuated, following the cuttrend we chose

  AP$CohortTrue <- AP$PeriodTrue - AP$Age # following the sane classification

  # which ones correspond to the ones we estimated and which ones do not?
  newcohorts <- unique(AP$CohortTrue)[(!unique(AP$CohortTrue) %in% unique(apcdata.fitted$CohortTrue))]
  # these ones should be set equal to the last (youngest) cohort that is estimated in the fitted model

  AP$Cohort <- AP$CohortTrue
  AP$Cohort[AP$Cohort %in% newcohorts] <- (newcohorts)[1]-1
  AP$Cohort <- as.factor(AP$Cohort)
  # the real values remain stored in CohortTrue
  # but this way the predict function knows what to do with them

  # unfortunately, the last observed cohort might also be a reference category
  # and we force two refs on one dimension to have the same designation, since
  # since having two refs on one dimension is not normal for R
  # so we need to make sure that we do so again just to be sure
  refcoh1st <- as.character(refcoh1st)
  refcoh2nd <- as.character(refcoh2nd)

  levels(AP$Cohort) <- c(levels(AP$Cohort),refcoh1st)
  AP$Cohort[AP$Cohort==refcoh2nd] <- refcoh1st

  # the Nordpred prediction function also automatically sets new periods to be equal to the last period
  # since the last period is a reference (technically, it doesn't do anything with Period
  # when it makes predictions, but that is equal to making it the ref, which is assumed to have an effect of 0
  # on top of the intercept)
  # so let's first set Period equal to last estimated Period
  # and then re-set the references, just like with Cohort
  AP$Period <- as.factor(as.character(max(apcdata.fitted$pcode)))
  levels(AP$Period) <- c(levels(AP$Period),refper1st)
  AP$Period[AP$Period==refper2nd] <- refper1st # and set it to be equal to first ref if it is a ref
  # this also needs to happen for apcdata if we want to make predictions for those
  # since some periods did not have estimated parameters due to the noperiod parameter
  apcdata$Period[apcdata$PeriodTrue <= (dnoperiods - noperiod)] <- refper1st
  # and the same for Cohort
  # Cohorts outside the range of fitted cohorts get set to recoh1st
  apcdata$Cohort[!apcdata$CohortTrue %in% unique(apcdata.fitted$CohortTrue)] <- refcoh1st

  # In order to join AP with apcdata I will remove or add some columns
  AP$periodfitted <- 0
  AP$cases <- NA
  apcdata$or.Age <- NULL # will be recalculated at the end
  apcdata$or.Per <- NULL # will be recalculated at the end
  apcdata$observed <- 1
  AP$PY <- 1 # will be assigned proper value later
  AP$observed <- 0
  AP <- rbind(apcdata,AP)
  AP$includepred <- ifelse(AP$Age < startestage,0,1)

  AP$y <- 1 # needed for the log-link part
  predobj <- predict(object=multiAPCMC.singlefit.object$glm,newdata=AP[AP$includepred==1,],type='link',se.fit=TRUE)
  AP$y <- NULL
  AP$inv_rate[AP$includepred==1] <- predobj$fit
  AP$inv_rate_low[AP$includepred==1] <- AP$inv_rate[AP$includepred==1] - predobj$se.fit*1.96
  AP$inv_rate_high[AP$includepred==1] <- AP$inv_rate[AP$includepred==1] + predobj$se.fit*1.96
  if(multiAPCMC.singlefit.object$link=='power5') {
    AP$rate <- AP$inv_rate^5
    AP$rate_low <- AP$inv_rate_low^5
    AP$rate_high <- AP$inv_rate_high^5
  } else {
    AP$rate <- exp(AP$inv_rate)
    AP$rate_low <- exp(AP$inv_rate_low)
    AP$rate_high <- exp(AP$inv_rate_high)
  }
  AP$inv_rate_low <- NULL
  AP$inv_rate_high <- NULL
  AP$inv_rate <- NULL
  AP$includepred <- NULL

  # let's also calculate the rates for the ages below startestage
  # for the training data
  AP$rate[AP$observed==1 & AP$Age < startestage] <- AP$cases[AP$observed==1 & AP$Age < startestage] / AP$PY[AP$observed==1 & AP$Age < startestage]
  AP$rate_low[AP$observed==1 & AP$Age < startestage] <- AP$rate_high[AP$observed==1 & AP$Age < startestage] <- AP$rate[AP$observed==1 & AP$Age < startestage]
  # if we want CIs for some of these, they would be calculable
  # but it would be pointless given how small these groups should be

  # now let's determine incidence in the lower age groups
  # following the Nordpred choice to use mean incidence in the last two periods

  # first we have to make a dataset with the lower ages
  underages <- unique(apcdata$Age)[unique(apcdata$Age) < startestage]
  PeriodTrue <- max(apcdata$pcode)+1:noproj
  APtrue_under <- data.frame(expand.grid(underages,PeriodTrue))
  names(APtrue_under) <- c("Age","pcode")
  # now I copy some irrelevant columns so that we can easily rbind later
  APtrue_under$PeriodTrue <- APtrue_under$pcode
  APtrue_under$Cohort <- APtrue_under$CohortTrue <- APtrue_under$PeriodTrue - APtrue_under$Age
  APtrue_under$Period <- APtrue_under$pcode
  APtrue_under$cases <- NA
  APtrue_under$PY <- 1
  APtrue_under$observed <- 0
  APtrue_under$periodfitted <- 0

  for(age in sort(unique(APtrue_under$Age))) {
    # For agegroups with little data, we use mean incidence for last two periods:
    APtrue_under$rate[APtrue_under$Age==age] <- mean((apcdata$cases/apcdata$PY)[apcdata$Age==age & apcdata$PeriodTrue>=max(apcdata$PeriodTrue-1)])
    APtrue_under$rate_low <- APtrue_under$rate_high <- APtrue_under$rate
    # we can safely assume 0 CI because we are not performing estimation here
    # even if we did assume so, this will have no relevant impact on the numbers
  }
  AP <- rbind(APtrue_under,AP)
  AP <- AP[order(AP$PeriodTrue,AP$Age),]

  # now in order to get predicted cases, we need to take the person-years into account
  # for the unobserved (out of sample: oos data) these were set to 1 so far
  futuredata$PeriodTrue <- futuredata$Period
  futuredata$Period <- NULL
  # sort both the same way
  futuredata <- futuredata[order(futuredata$PeriodTrue,futuredata$Age),]

  ## transform futuredata Period and Age so that the data can be merged with AP
  ap.pcodes <- sort(unique(AP$PeriodTrue[AP$observed==0]))
  py.pcodes <- sort(unique(futuredata$PeriodTrue))
  futuredata$PeriodTrue.temp <- futuredata$PeriodTrue
  for(pc in 1:length(ap.pcodes)) {
    futuredata$PeriodTrue[futuredata$PeriodTrue.temp==py.pcodes[pc]] <- ap.pcodes[pc]
  }
  ap.acodes <- sort(unique(AP$Age))
  py.acodes <- sort(unique(futuredata$Age))
  futuredata$Age.temp <- futuredata$Age
  for(ac in 1:length(ap.acodes)) {
    futuredata$Age[futuredata$Age.temp==py.acodes[ac]] <- ap.acodes[ac]
  }
  futuredata$PeriodTrue.temp <- NULL
  futuredata$Age.temp <- NULL

  # now I need to join on Age and PeriodTrue
  # because we already have PY and we want to keep it for the
  # observed data (training data)
  futuredata <- data.frame(Age=futuredata$Age,PeriodTrue=futuredata$PeriodTrue,PY=futuredata$PY)
  AP <- merge(AP,futuredata,by=c('PeriodTrue','Age'),all.x=TRUE)
  AP$PY <- AP$PY.x
  AP$PY[AP$observed==0] <- AP$PY.y[AP$observed==0]
  AP$PY.y <- AP$PY.x <- NULL

  # now multiply the rates with the person-years
  # the person-years should be entered by the user
  AP$pred_cases <- AP$rate * AP$PY
  AP$pred_cases_low <- AP$rate_low * AP$PY
  AP$pred_cases_high <- AP$rate_high * AP$PY

  ## reconstruct or.Age and or.Per
  AP$or.Per <- or.Per.int+(AP$PeriodTrue-1)*or.Per.is
  AP$or.Age <- or.Age.int+(AP$Age-1)*or.Age.is

  # then let's get rid of the columns that we don't need
  # we really just need A and P, but we can also keep C, and of course Rate, PY and cases
  # and or.Age and or.Per so that we know what the 'real' historical values are
  AP$Period <- AP$PeriodTrue
  AP$PeriodTrue <- AP$pcode <- NULL
  AP$Cohort <- AP$CohortTrue
  AP$CohortTrue <- NULL
  AP$y <- NULL

  AP <- AP[order(AP$Period,AP$Age),]

  return(AP)
}
