#'
#' @title multiAPCMC.datashape
#'
#' @description The APCMC approach, following Clayton & Schifflers (Statistics in Medicine, 1987), works by setting additional reference categories on the time dimensions.
#' This is most easily accomplished by directly manipulating the factor coding for the time dimensions. Datashape does this.
#'
#' @param data a data frame containing the variables in the model.
#' @param noperiod the number of periods in the prediction base. For example, if you have 20 years of data,
#' and data is categorized as 1 year period by 1 year age, then '10' means you only use the last 10 years of the data to extrapolate towards the future.
#' @param startestage youngest age group to be included in the regression model. Note that this is based on how data is categorized.
#'  If you have have ages going from 0 to 75 with data on single year ages, then age 0 is the first age category and hence are age '1', age 1 is then '2', and so forth.
#'  If you have 5-year age categories and data starting from age 0 to 4, then '1' is age category 0 to 4, '2' is age category 5 to 9, and so forth.
#' @param refper1st first reference category of period, based on 'period codes'; earliest period is 1, second earliest is 2, etc.
#' @param refper2nd second reference category of period, based on 'period codes'; earliest period is 1, second earliest is 2, etc.
#' @param refcoh1st first reference category of cohort, based on 'cohort codes'; earliest cohort is 1, second earliest is 2, etc.
#' @param refcoh2nd second reference category of cohort, based on 'cohort codes'; earliest cohort is 1, second earliest is 2, etc.
#' @param cases the name (as a character) of the column containing the cluster identifiers.
#'
#' @return returns a data.frame with a structure that has the right constraints to subsequently be used by multiAPCMC.multifit and multiAPCMC.singlefit
#' @importFrom stats relevel
#' @export
#'
#' @examples
#' set.seed(100)
#' # datashape will generally not be called directly by a user.
#' # Instead, it will be called by multiAPCMC.multifit.
#' multiAPCMC.example.data <- multiAPCMC.example.data
#' apcdata <- multiAPCMC.datashape(multiAPCMC.example.data,
#'                                 noperiod=length(unique(multiAPCMC.example.data$Period)),
#'                                 startestage=1,
#'                                 refper1st=1,
#'                                 refper2nd=2,
#'                                 refcoh1st=1,
#'                                 refcoh2nd=2)
#'
multiAPCMC.datashape_NEO <- function(data,noperiod=4,startestage=5,refper1st,refper2nd,refcoh1st,refcoh2nd,cases=TRUE) {

  # check if columns are named as they should
  namesindf <- names(data)
  if(!"Age" %in% namesindf) {
    print("No column labelled 'Age' detected")
    stop()
  }
  if(!"Period" %in% namesindf) {
    print("No column labelled 'Period' detected")
    stop()
  }
  if(!"PY" %in% namesindf) {
    print("No column labelled 'PY' (person-years) detected")
    stop()
  }
  if(cases==TRUE) { # I do an additional check for this because we don't need it when making futureAPC objects with just PY
    if(!"cases" %in% namesindf) {
      print("No column labelled 'cases' detected")
      stop()
    } }

  # check if other columns are included that we don't want
  # and give warning if detected
  if(!all(namesindf %in% c("Age","Period","PY","cases"))) {
    warning("Columns other than Age, Period, cases and PY were detected, please
  note that this may result in errors when running later functions. Remove these
  additional columns from the data to be sure")
  }

  # determine number of age and period categories
  dnoagegr <- length(unique(data$Age))
  dnoperiods <- length(unique(data$Period))

  # Transform dataformat:
  data <- data[order(data$Period,data$Age),]
  data$or.Age <- data$Age # store untransformed age
  data$or.Per <- data$Period # store untransformed period
  agemode <- getmode(diff(data$Age))
  permode <- getmode(diff(data$Period[!(diff(data$Period)==0)]))
  data$acode <- (data$Age - min(data$Age))/agemode +1
  data$pcode <- (data$Period - min(data$Period))/permode +1
  data$Cohort <- data$Period - data$Age # this is a sane cohort classification

  # ! this could be even more improved because why make everything into codes?
  # we could just keep the ages, periods and cohorts by their real names
  # though since cohort is derived from age and period, they are born in a longer interval
  # so setting cohort ref becomes less 'natural'

  # ! set refs as character (if they already are, it does not matter)
  refper1st <- as.character(refper1st)
  refper2nd <- as.character(refper2nd)
  refcoh1st <- as.character(refcoh1st)
  refcoh2nd <- as.character(refcoh2nd)

  # I won't throw out data yet, since we need a WHOLE apcdataset for prediction later
  # but I do need to know which variables we can set as cohorts or periods
  # so I will just make a quick second dataset where we do throw out the data
  # based on startestage and noperiod
  data2 <- data[data$Age>=startestage,]
  data2 <- data2[data2$pcode>(dnoperiods-noperiod),]

  # ! check if the refs of period and cohort are actually in the list of period and cohorts
  # if not, quit and outputinfo

  if(refper1st %in% sort(unique(data2$Period)) == FALSE) {
    print("period refs need to be period categories")
    print("with this data, the following period categories can be chosen as refs")
    print(sort(unique(data2$Period))[-length(unique(data2$Period))])
    stop()
  }

  if(refper2nd %in% sort(unique(data2$Period)) == FALSE) {
    print("period refs need to be period categories")
    print("with this data, the following period categories can be chosen as refs")
    print(sort(unique(data2$Period))[-length(unique(data2$Period))])
    stop()
  }

  if(refcoh1st %in% sort(unique(data2$Cohort)) == FALSE) {
    print("cohort refs need to be cohort categories")
    print("with this data, the following cohort categories can be chosen as refs")
    print(sort(unique(data2$Cohort))[-length(unique(data2$Period))])
    stop()
  }

  if(refcoh2nd %in% sort(unique(data2$Cohort)) == FALSE) {
    print("cohort refs need to be cohort categories")
    print("with this data, the following cohort categories can be chosen as refs")
    print(sort(unique(data2$Cohort))[-length(unique(data2$Period))])
    stop()
  }
  rm(data2)

  # ! setting pcode as drift
  # ! setting Period and Cohort as factors
  # ! setting 2nd ref of Period and Cohort, rather than automatically last cat as ref

  # columns that store original information
  # this can be used when plotting or making preds
  data$PeriodTrue <- data$Period
  data$CohortTrue <- data$Cohort

  # the ones we use in the model should be factors (except pcode, which is drift)
  data$Period <- as.factor(data$Period)
  data$Cohort <- as.factor(data$Cohort)

  # set second refs equal to first refs
  data$Period[data$Period==refper2nd] <- refper1st
  data$Cohort[data$Cohort==refcoh2nd] <- refcoh1st

  # now that second and first refs have the same code, set that code as ref
  data$Period <- relevel(data$Period, ref=refper1st)
  data$Cohort <- relevel(data$Cohort, ref=refcoh1st)

  return(list(apcdata=data,
              noperiod=noperiod,
              startestage=startestage,
              refper1st=refper1st,
              refper2nd=refper2nd,
              refcoh1st=refcoh1st,
              refcoh2nd=refcoh2nd))
}



