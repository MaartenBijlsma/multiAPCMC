#'
#' @title multiAPCMC.retrievemodel
#'
#' @description retrieves a model, or its predictions, with these specifications from either a multiAPCMC.multifit object or from a multiAPCMC.multipred object.
#'
#' @param multiAPCMC.multipred.object either a multiAPCMC.multipred.object from the multiAPCMC.multipred function, or a multiAPCMC.multifit.object from the multiAPCMC.multifit function.
#' @param what should be \code{'pred'} if we want to retrieve predictions from a multiAPCMC.multipred object or it should be \code{'fit'} if we wish to retrieve a model from a multiAPCMC.multifit object.
#' @param link should be log-link ("log") or a power5-link ("power5").
#' @param noperiod the number of periods in the prediction base.
#' @param refper the period references. Options are 'extremes', 'outer', 'center', 'first middle', 'middle last', 'first second', 'penultimate last'.
#' @param refcoh the cohort references. Options are 'extremes', 'outer', 'center', 'first middle', 'middle last', 'first second', 'penultimate last'.
#'
#' @return returns a nested list with all of the fitted models.
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
#'                          
#'                                              
multiAPCMC.retrievemodel <- function(multiAPCMC.multipred.object,
                           what="pred",
                           link,
                           noperiod,
                           refper,
                           refcoh) {
  
  l <- which(multiAPCMC.multipred.object$vec.link==as.character(link))
  n <- which(multiAPCMC.multipred.object$vec.noperiod==as.numeric(noperiod))
  p <- which(multiAPCMC.multipred.object$vec.refper==as.character(refper))
  c <- which(multiAPCMC.multipred.object$vec.refcoh==as.character(refcoh))
  
  if(what=="pred") {
    return(multiAPCMC.multipred.object$np.pred[[l]][[n]][[p]][[c]])
  } else if(what=="fit") {
    return(multiAPCMC.multipred.object$np.est[[l]][[n]][[p]][[c]])
  }
  
}