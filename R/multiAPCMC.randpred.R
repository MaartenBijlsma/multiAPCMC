#'
#' @title multiAPCMC.randpred
#'
#' @description A function to simulate from a fitted APCMC model, used for parametric bootstrapping. APCMC models are built using generalized linear modelling machinery, and hence their coefficients are assumed to be multivariate normally distributed. This function uses the coefficients and variance-covariance matrix of a fitted APCMC model and uses them to make draw from a multivariate normal distribution. Do this many times to get a distribution of draws. This distribution can then be used to determine confidence intervals of predictions/forecasts.
#'
#' @param multiAPCMC.singlefit.object The list of fitted APCMC models that is produced by \code{multiAPCMC.multifit}.
#' @param futuredata a dataframe containing information on future periods and person-years (exposed population) in those years. Should have the same resolution as the data on which the models were also fitted (e.g. 1 year age group by 1 year periods, or 5-year age groups by 5-year periods, and so forth) and should have the columns Period, Age, and PY (person-years).
#' @param noproj number of periods to project to the future. Should have at most the length of the number of periods in futuredata.
#' @param cuttrend Cut trend in predictions? Should be a vector with at least a length equal to number of periods to predict, e.g. \code{c(0,0,0,0)} or \code{c(0,0.25,0.50,0.75)} or so if we project 4 periods into the future.
#'
#' @return returns a nested list with all of the fitted models.
#' @importFrom stats vcov
#' @importFrom stats glm
#' @importFrom stats predict
#' @export
#'
#' @examples
#' set.seed(100)
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
#' multiAPCMC.randpred(chosenmodel,
#'                     multiAPCMC.example.futuredata,
#'                     noproj=4,
#'                     cuttrend=c(0,0,0,0))
#'
#'#' @import MASS
multiAPCMC.randpred <- function(multiAPCMC.singlefit.object,
                     futuredata,
                     noproj=4,
                     cuttrend) {

  meancoefs <- multiAPCMC.singlefit.object$glm$coefficients
  vcovcoefs <- vcov(multiAPCMC.singlefit.object$glm)

  multiAPCMC.singlefit.object$glm$coefficients <- MASS::mvrnorm(1,meancoefs,vcovcoefs)

  multiAPCMC.singlepred(multiAPCMC.singlefit.object=multiAPCMC.singlefit.object,
                    futuredata=futuredata,
                    noproj=noproj,
                    cuttrend=cuttrend)

}
