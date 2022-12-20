#'
#' @title multiAPCMC.multipred
#'
#' @description Makes predictions for multiple specifications of an APCMC model at once.
#'
#' @param multiAPCMC.multifit.object The list of fitted APCMC models that is produced by \code{multiAPCMC.multifit}.
#' @param futuredata a dataframe containing information on future periods and person-years (exposed population) in those years. Should have the same resolution as the data on which the models were also fitted (e.g. 1 year age group by 1 year periods, or 5-year age groups by 5-year periods, and so forth) and should have the columns Period, Age, and PY (person-years).
#' @param cuttrend Cut trend in predictions? Should be a vector with at least a length equal to number of periods to predict, e.g. \code{c(0,0,0,0)} or \code{c(0,0.25,0.50,0.75)} or so if we project 4 periods into the future.
#'
#' @return returns a nested list with predictions from all of the models that were fitted with the \code{multiAPCMC.multifit.object}.
#' @export
#'
#' @examples
#'  multiAPCMC.example.data <- multiAPCMC.example.data
#' multiAPCMC.example.futuredata <- multiAPCMC.example.futuredata
#' multifit.example <- multiAPCMC.multifit(multiAPCMC.example.data,
#'                     startestage=5,
#'                     vec.link=c('log','power5'),
#'                     vec.noperiod=c(5),
#'                     vec.refper=c('outer'),
#'                     vec.refcoh=c('outer'),
#'                     nomodelmax=1024)
#' multiAPCMC.multipred(multifit.example,
#'                      multiAPCMC.example.futuredata,
#'                      cuttrend=rep(0,15))
#'
multiAPCMC.multipred <- function(multiAPCMC.multifit.object,futuredata,cuttrend) {

  # check if futuredata meets the noproj
  noproj <- length(unique(futuredata$Period))
  print(paste("Projecting",noproj,"periods into the future",sep=' '))

  # check if number of age groups in futuredata is equal to the number of age groups in the apcobject
  age_futuredata <- length(unique(futuredata$Age))
  per_futuredata <- length(unique(futuredata$Period))
  age_APCobject <- length(unique(multiAPCMC.multifit.object[[6]][[1]][[1]][[1]][[1]]$apcdata$Age))
  if(age_futuredata != age_APCobject) {
    stop("the number of future age categories (from futuredata) are not equal to the number of estimated
         age categories (from multiAPCMC.multifit.object$np.est apcobject). This is needed because we will assume
         that the youngest age category in the futuredata object corresponds to the youngest in the training data
         and the oldest corresponds to the oldest, etc.")
  }
  # check also if age and period are sorted (ascending)
  if(is.unsorted(age_futuredata) | is.unsorted(per_futuredata)) {
    stop("sort the futuredata data by period and age (in that order), this is needed because we will assume
         that the youngest age category in the futuredata object corresponds to the youngest in the training data
         and the oldest corresponds to the oldest, etc. The same logic applies with period")
  }

  startestage <- multiAPCMC.multifit.object$startestage
  vec.link <- multiAPCMC.multifit.object$vec.link
  vec.noperiod <- multiAPCMC.multifit.object$vec.noperiod
  vec.refper <- multiAPCMC.multifit.object$vec.refper
  vec.refcoh <- multiAPCMC.multifit.object$vec.refcoh

  depth = c(length(vec.link),
            length(vec.noperiod),
            length(vec.refper),
            length(vec.refcoh)
  )
  np.pred <- multiAPCMC.reclist(depth)

  for(l in 1:length(vec.link)) {

    for(n in 1:length(vec.noperiod)) {

      noperiod <- multiAPCMC.multifit.object$np.est[[1]][[n]][[1]][[1]]$noperiod

      for(p in 1:length(vec.refper)) {

        for(c in 1:length(vec.refcoh)) {

          # perform the predictions
          np.pred[[l]][[n]][[p]][[c]] <-  multiAPCMC.singlepred(multiAPCMC.multifit.object$np.est[[l]][[n]][[p]][[c]],
                                                           futuredata,
                                                           noproj,
                                                           cuttrend)

        }

      }

    }

  }

  return(list(
    vec.link=vec.link,
    vec.noperiod=vec.noperiod,
    vec.refper=vec.refper,
    vec.refcoh=vec.refcoh,
    np.pred=np.pred))
}
