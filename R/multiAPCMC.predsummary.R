#'
#' @title multiAPCMC.predsummary
#'
#' @description summarizes the outcome by year within some age range (i.e. adds all incidences together); use this together with randpred to get confidence intervals for an age range (see documentation).
#'
#' @param years a vector of the years that you want to have summaries for
#' @param modelpred the output of a retrievelmodel(multipred) or a randpred
#' @param agemin the lower bound of the age categories that you want to summarize the outcome over
#' @param agemax the upper bound o fthe age categories that you want to summarize the outcome over
#' @param poissondraw set to \code{TRUE} if you want to add a random draw from a Poisson distribution with lambda set to the summary value for that year (e.g. to build a prediction interval). \code{Default is FALSE}.
#'
#' @return returns a nested list with predictions from all of the models that were fitted with the \code{multiAPCMC.multifit.object}.
#' @importFrom stats rpois
#' @export
#'
#' @examples
#' multiAPCMC.example.data <- multiAPCMC.example.data
#' multifit.example <- multiAPCMC.multifit(multiAPCMC.example.data,
#'                                        startestage=5,
#'                                        vec.link=c('log','power5'),
#'                                        vec.noperiod=c(5,10),
#'                                        vec.refper=c('outer','center'),
#'                                        vec.refcoh=c('outer'),
#'                                        nomodelmax=1024)
#' chosenmodel <- multiAPCMC.retrievemodel(multifit.example,
#'                                        what="fit",
#'                                        link='power5',
#'                                        noperiod=10,
#'                                        refper='center',
#'                                        refcoh='outer')
#' years <- 1989:2020
#' multiAPCMC.predsummary(years,chosenmodel)
#'
multiAPCMC.predsummary <- function(years, modelpred, agemin=0, agemax=120, poissondraw=FALSE) {

  totcases <- totcases_low <- totcases_high <- NULL
  for(k in 1:length(years)) {
    totcases[k] <- sum(modelpred$pred_cases[modelpred$or.Per==years[k] & modelpred$or.Age >=agemin & modelpred$or.Age <=agemax])
  }

  if(poissondraw==TRUE) {
    totcases <- rpois(length(totcases),totcases)
  }

  return(data.frame(years=years,pred=totcases))
}
