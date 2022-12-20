#'
#' @title multiAPCMC.singlefit
#'
#' @description Fits multiple a single specification of an APCMC model. This function is called by multiAPCMC.multifit and hence it is not expected that users use the singlefit version directly, even when fitting just one model because multifit allows easy specification of references and such.
#'
#' @param multiAPCMC.datashape.object a data frame containing the variables in the model.
#' @param link the specified link function. Choices are log-link ("log") and power5-link ("power5")
#'
#' @return returns a single fitted APCMC model.
#' @importFrom stats poisson
#' @export
#'
#' @examples
#' multiAPCMC.example.data <- multiAPCMC.example.data
#' apcdata <- multiAPCMC.datashape(multiAPCMC.example.data,
#'                                 noperiod=length(unique(multiAPCMC.example.data$Period)),
#'                                 startestage=1,
#'                                 refper1st=1,
#'                                 refper2nd=2,
#'                                 refcoh1st=1,
#'                                 refcoh2nd=2)
#' multiAPCMC.singlefit(apcdata,link='power5')
#'
multiAPCMC.singlefit <- function(multiAPCMC.datashape.object,link="power5") {

  apcdata <- multiAPCMC.datashape.object$apcdata
  startestage <- multiAPCMC.datashape.object$startestage
  noperiod <- multiAPCMC.datashape.object$noperiod

  # Selecting data for regression:
  dnoperiods <- length(unique(apcdata$PeriodTrue))
  apcdata <- apcdata[apcdata$Age>=startestage,]
  apcdata <- apcdata[apcdata$pcode>(dnoperiods-noperiod),]

  # pcode is the drift parameter

  # Sett variable for use in estimation
  y <- apcdata$PY
  # Make link function:
  power5link <- poisson()
  power5link$link <- "0.2 root link Poisson family"
  power5link$linkfun <- function(mu)  { (mu/y)^0.2 }
  power5link$linkinv <- function(eta) { pmax(.Machine$double.eps, y*eta^5) }
  power5link$mu.eta <- function(eta)  { pmax(.Machine$double.eps, 5*y*eta^4) }

  # Estimation:
  if (link=="power5") {
    res.glm <- glm(cases~as.factor(Age)+pcode+as.factor(Period)+as.factor(Cohort) -1,family=power5link,data=apcdata)
  } else  if (link=="log") {
    res.glm <- glm(cases~as.factor(Age)+pcode+as.factor(Period)+as.factor(Cohort)+ offset(log(y)) -1,family=poisson(),data=apcdata)
  } else {
    stop("Unknown \"link\"")
  }

  # Set class and return results
  res <- list(glm=res.glm,
              apcdata=multiAPCMC.datashape.object$apcdata,
              noperiod=multiAPCMC.datashape.object$noperiod,
              startestage=multiAPCMC.datashape.object$startestage,
              refper1st=multiAPCMC.datashape.object$refper1st,
              refper2nd=multiAPCMC.datashape.object$refper2nd,
              refcoh1st=multiAPCMC.datashape.object$refcoh1st,
              refcoh2nd=multiAPCMC.datashape.object$refcoh2nd,
              link=link)
  return(res)
}
