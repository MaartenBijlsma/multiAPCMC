#'
#' @title multiAPCMC.multifit
#'
#' @description Fits multiple specifications of an APCMC model at once.
#'
#' @param data a data frame containing the variables in the model.
#' @param startestage youngest age group to be included in the regression model. Note that this is based on how data is categorized. If you have have ages going from 0 to 75 with data on single year ages, then age 0 is the first age category and hence are age '1', age 1 is then '2', and so forth. If you have 5-year age categories and data starting from age 0 to 4, then '1' is age category 0 to 4, '2' is age category 5 to 9, and so forth.
#' @param vec.link a vector stating which link functions the user wants to explore; choices are log-link ("log") and power5-link ("power5")
#' @param vec.noperiod a vector containing the number of periods in the prediction base that the user wants to explore. For example, if you have 20 years of data, and data is categorized as 1 year period by 1 year age, then '10' means you only use the last 10 years of the data to extrapolate towards the future. Here you can specify multiple options, etc. c(3,5,7,10).
#' @param vec.refper which types of references to use for period. Options are 'extremes' (the first and last period), 'outer' (periods that represent roughly the 20\% and 80\% quantiles of all periods), 'center' (the middle categories), 'first middle' (the firsts period and the middle period), 'middle last' (the middle and the last category), 'first second' (the first and second period categories), 'penultimate last' (the penultimate and the last period category).
#' @param vec.refcoh which types of references to use for cohort. Options are 'extremes' (the first and last cohort), 'outer' (cohort that represent roughly the 20\% and 80\% quantiles of all cohorts), 'center' (the middle categories), 'first middle' (the firsts cohort and the middle cohort), 'middle last' (the middle and the last category), 'first second' (the first and second cohort categories), 'penultimate last' (the penultimate and the last cohort category).
#' @param nomodelmax outputs an error if more models would be fitted than this number.
#'
#' @return returns a nested list with all of the fitted models.
#' @importFrom stats quantile
#' @importFrom stats median
#' @export
#'
#' @examples
#' multiAPCMC.example.data <- multiAPCMC.example.data
#' multiAPCMC.multifit(multiAPCMC.example.data,
#'                     startestage=5,
#'                     vec.link=c('log'),
#'                     vec.noperiod=c(5,10),
#'                     vec.refper=c('outer','center'),
#'                     vec.refcoh=c('center'),
#'                     nomodelmax=1024)
#'
multiAPCMC.multifit_NEO <- function(data,startestage,vec.link,vec.noperiod,vec.refper,vec.refcoh,nomodelmax=1024) {
  depth = c(length(vec.link),
            length(vec.noperiod),
            length(vec.refper),
            length(vec.refcoh)
  )
  np.est <- multiAPCMC.reclist(depth)
  totalmodels <- length(vec.link)*length(vec.noperiod)*length(vec.refper)*length(vec.refcoh)
  print(paste(c("fitting",as.character(totalmodels),"APCMC models"),sep=" ",collapse=" "))
  if(totalmodels > nomodelmax) {
    stop("The total number of models that have to be estimated exceeds the total number of models that you currently allow to be estimated (nomodelmax parameter)")
  }
  ## run datashape2 here once, just so we get PeriodTrue and CohortTrue
  # which don't change even if we change the refs
  apcdata <- multiAPCMC.datashape_NEO(data,
                                noperiod=length(unique(data$Period)),
                                startestage=0,
                                refper1st=min(data$Period),
                                refper2nd=max(data$Period),
                                refcoh1st=min(data$Period-data$Age),
                                refcoh2nd=max(data$Period-data$Age))$apcdata

  # practically, those refs should always work
  # but it is not very elegant
  # now I need to fill those nested lists with estimated models
  agemode <- getmode(diff(apcdata$Age))
  permode <- permode <- getmode(diff(as.numeric(as.character(apcdata$Period[!(diff(data$Period)==0)]))))
  totalperiod <- length(unique(apcdata$PeriodTrue))
  totalcohort <- length(unique(apcdata$CohortTrue[apcdata$Age >= startestage]))
  refper1st <- NULL
  refper2nd <- NULL
  refcoh1st <- NULL
  refcoh2nd <- NULL
  for(l in 1:length(vec.link)) {
    link <- vec.link[l]
    for(n in 1:length(vec.noperiod)) {
      noperiod <- vec.noperiod[n]
      for(p in 1:length(vec.refper)) {
        if(noperiod > totalperiod) {
          stop("The total number of periods used to perform estimation (noperiod) is set to be higher than the number of periods detected in the dataframe")
        }
        if(vec.refper[p]=='extremes') {
          refper1st <- max(apcdata$PeriodTrue)-noperiod*permode+permode
          refper2nd <- max(apcdata$PeriodTrue)
        } else if(vec.refper[p]=='outer') {
          pseq <- seq(from=max(apcdata$PeriodTrue)-noperiod*permode+permode,to=max(apcdata$PeriodTrue),by=permode)
          refs <- round(quantile(pseq,c(0.2,0.8)))
          refper1st <- refs[1]
          refper2nd <- refs[2]
        } else if(vec.refper[p]=='center') {
          pseq <- seq(from=max(apcdata$PeriodTrue)-noperiod*permode+permode,to=max(apcdata$PeriodTrue),by=permode)
          refper1st <- floor(median(pseq))
          refper2nd <- refper1st + permode
        } else if(vec.refper[p]=='first middle') {
          pseq <- seq(from=max(apcdata$PeriodTrue)-noperiod*permode+permode,to=max(apcdata$PeriodTrue),by=permode)
          refper1st <- max(apcdata$PeriodTrue)-noperiod*permode+permode
          refper2nd <- floor(median(pseq))
        } else if(vec.refper[p]=='middle last') {
          pseq <- pseq <- seq(from=max(apcdata$PeriodTrue)-noperiod*permode+permode,to=max(apcdata$PeriodTrue),by=permode)
          refper1st <- floor(median(pseq))
          refper2nd <- max(apcdata$PeriodTrue)
        } else if(vec.refper[p]=='first second') {
          refper1st <- max(apcdata$PeriodTrue)-noperiod*permode+permode
          refper2nd <- refper1st + permode
        } else if(vec.refper[p]=='penultimate last') {
          refper1st <- max(apcdata$PeriodTrue) - permode
          refper2nd <- refper1st + permode
        } else {warning('refper not recognized')}
        for(c in 1:length(vec.refcoh)) {
          uniquecohorts <- sort(unique(apcdata$CohortTrue[apcdata$Age >= startestage &
                                                            apcdata$PeriodTrue > max(apcdata$PeriodTrue)-noperiod*permode]))
          if(vec.refcoh[c]=='extremes') {
            refcoh1st <- min(uniquecohorts)
            refcoh2nd <- max(uniquecohorts)
          } else if(vec.refcoh[c]=='outer') {
            refs <- round(quantile(uniquecohorts,c(0.2,0.8)))
            refcoh1st <- refs[1]
            refcoh2nd <- refs[2]
          } else if(vec.refcoh[c]=='center') {
            refcoh1st <- floor(median(uniquecohorts))
            refcoh2nd <- refcoh1st +1
          } else if(vec.refcoh[c]=='first middle') {
            refcoh1st <- min(uniquecohorts)
            refcoh2nd <- floor(median(uniquecohorts))
          } else if(vec.refcoh[c]=='middle last') {
            refcoh1st <- floor(median(uniquecohorts))
            refcoh2nd <- max(uniquecohorts)
          } else if(vec.refcoh[c]=='first second') {
            refcoh1st <- min(uniquecohorts)
            refcoh2nd <- refcoh1st +1
          } else if(vec.refcoh[c]=='penultimate last') {
            refcoh1st <- max(uniquecohorts)-1
            refcoh2nd <- refcoh1st+1
          } else {warning('refcoh not recognized')}
          apcobject <- multiAPCMC.datashape_NEO(data,
                                          noperiod=noperiod,
                                          startestage=startestage,
                                          refper1st=refper1st,
                                          refper2nd=refper2nd,
                                          refcoh1st=refcoh1st,
                                          refcoh2nd=refcoh2nd)
          np.est[[l]][[n]][[p]][[c]] <- multiAPCMC.singlefit(apcobject,link=link)
        }
      }
    }
  }
  return(list(
    startestage=startestage,
    vec.link=vec.link,
    vec.noperiod=vec.noperiod,
    vec.refper=vec.refper,
    vec.refcoh=vec.refcoh,
    np.est=np.est))
}
