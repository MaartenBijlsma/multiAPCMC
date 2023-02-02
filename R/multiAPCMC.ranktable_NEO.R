#'
#' @title multiAPCMC.ranktable
#'
#' @description Makes a table with summary statistics model predictions compared to out of sample data. The summary statistics reported are root mean squared error (RMSE), mean absolute error (MAE), and period corrected AIC (pcAIC) which is determined only on the data to which the model was fitted, and standard error of the drift (driftSE)
#'
#' @param multiAPCMC.multifit.object The list of fitted APCMC models that is produced by \code{multiAPCMC.multifit}.
#' @param multiAPCMC.multipred.object The list of predictions of the fitted APCMC models that is produced by \code{multiAPCMC.multipred}. If there were no future predictions to validate against and the user wants to see AIC and so forth from the original model fit, then set this parameter to NULL.
#' @param rankhow which metric should be used to rank (sort) the models? Options are no sorting ('none'), pcAIC, MAE, RMSE, and driftSE. MAE or RMSE are recommended.
#' @param oos.data a dataframe containing information on future periods, cases (e.g. incidence) and person-years (exposed population) in those years. Should have the same resolution as the data on which the models were also fitted (e.g. 1 year age group by 1 year periods, or 5-year age groups by 5-year periods, and so forth) and should have the columns Period, Age, cases, and PY (person-years).
#' @param validate.what which data should be used to test the predictions against? \code{oos} stands for 'out of sample' and hence only compares fit of predictions to data from the \code{oos.data}. The option \code{not periodfitted} refers to all data that the model was not fitted on. For example, if you have 20 years of data to fit the model on, but had \code{noperiod} of 10, then the first 10 years of that data are used, in addition to the \code{oos.data}, to compare the model against (essentially, you forecast backwards). The option \code{all} compares to everything, including the period that the model was fitted on. Hence, this last option is not recommended.
#'
#' @return returns a dataframe with model information and their ranking
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
#' multipred.example <- multiAPCMC.multipred(multifit.example,
#'                                           multiAPCMC.example.futuredata,
#'                                           cuttrend=rep(0,15))
#' multiAPCMC.example.futuredata <-
#' multiAPCMC.example.futuredata[order(multiAPCMC.example.futuredata$Period,
#'                                    multiAPCMC.example.futuredata$Age),]
#' multiAPCMC.ranktable(multifit.example,
#'                      multipred.example,
#'                      rankhow='RMSE',
#'                      multiAPCMC.example.futuredata,
#'                      validate.what="oos")
#' # note that multiAPCMC.example.futuredata contains 'cases' so that
#' # we can use it both to predict
#' # and to subsequently test the fit of our model to the cases in
#' # this future data.
multiAPCMC.ranktable_NEO <- function(multiAPCMC.multifit.object,
                                 multiAPCMC.multipred.object=NULL,
                                 rankhow=c("none","pcAIC","MAE","RMSE","driftSE"),
                                 oos.data=NULL,
                                 validate.what='oos') {

  vec.link <- multiAPCMC.multifit.object$vec.link
  vec.noperiod <- multiAPCMC.multifit.object$vec.noperiod
  vec.refper <- multiAPCMC.multifit.object$vec.refper
  vec.refcoh <- multiAPCMC.multifit.object$vec.refcoh

  ranktable <- expand.grid(vec.link,vec.noperiod,vec.refper,vec.refcoh)
  names(ranktable) <- c("link","noperiod","ref_period","ref_cohort")
  ranktable <- data.frame(ranktable)

  if(!is.null(multiAPCMC.multipred.object)) {
    multiAPCMC.multipredtest.object <- multiAPCMC.multipredtest(multiAPCMC.multipred.object = multiAPCMC.multipred.object,
                                                                oos.data = oos.data,
                                                                validate.what = validate.what)
  }

  for(l in 1:length(vec.link)) {

    lname <- vec.link[l]

    for(n in 1:length(vec.noperiod)) {

      nname <- vec.noperiod[n]

      for(p in 1:length(vec.refper)) {

        pname <- vec.refper[p]

        for(c in 1:length(vec.refcoh)) {

          cname <- vec.refcoh[c]

          ranktable$refper1st[ranktable$link==lname &
                                ranktable$noperiod==nname &
                                ranktable$ref_period==pname &
                                ranktable$ref_cohort==cname] <- multiAPCMC.multifit.object$np.est[[l]][[n]][[p]][[c]]$refper1st
          ranktable$refper2nd[ranktable$link==lname &
                                ranktable$noperiod==nname &
                                ranktable$ref_period==pname &
                                ranktable$ref_cohort==cname] <- multiAPCMC.multifit.object$np.est[[l]][[n]][[p]][[c]]$refper2nd
          ranktable$refcoh1st[ranktable$link==lname &
                                ranktable$noperiod==nname &
                                ranktable$ref_period==pname &
                                ranktable$ref_cohort==cname] <- multiAPCMC.multifit.object$np.est[[l]][[n]][[p]][[c]]$refcoh1st
          ranktable$refcoh2nd[ranktable$link==lname &
                                ranktable$noperiod==nname &
                                ranktable$ref_period==pname &
                                ranktable$ref_cohort==cname] <- multiAPCMC.multifit.object$np.est[[l]][[n]][[p]][[c]]$refcoh2nd

          ranktable$is.AIC[ranktable$link==lname &
                             ranktable$noperiod==nname &
                             ranktable$ref_period==pname &
                             ranktable$ref_cohort==cname] <- unlist(multiAPCMC.multifit.object$np.est[[l]][[n]][[p]][[c]])$glm.aic

          if(!is.null(multiAPCMC.multipred.object)) {
          ranktable$oos.RMSE[ranktable$link==lname &
                               ranktable$noperiod==nname &
                               ranktable$ref_period==pname &
                               ranktable$ref_cohort==cname] <- multiAPCMC.multipredtest.object[[l]][[n]][[p]][[c]][1]

          ranktable$oos.MAE[ranktable$link==lname &
                              ranktable$noperiod==nname &
                              ranktable$ref_period==pname &
                              ranktable$ref_cohort==cname] <- multiAPCMC.multipredtest.object[[l]][[n]][[p]][[c]][2]
          }

          modelcoefs <- summary(multiAPCMC.multifit.object$np.est[[l]][[n]][[p]][[c]]$glm)$coef
          driftindex <- match(row.names(modelcoefs),"pcode",nomatch=0)
          ranktable$driftSE[ranktable$link==lname &
                              ranktable$noperiod==nname &
                              ranktable$ref_period==pname &
                              ranktable$ref_cohort==cname] <- modelcoefs[driftindex==1,2]

        }

      }

    }

  }

  if(is.null(multiAPCMC.multipred.object) | is.null(oos.data)) {
    ranktable$oos.RMSE <- 'specify multipred object'
    ranktable$oos.MAE <- 'specify multipred object'

    rankhow <- 'pcAIC'

    print('no multipred object was specified or no oos.data was specified')
    print('ranktable is filled with model fit statistics from training')
  }

  ranktable$is.pcAIC <- ranktable$is.AIC/ranktable$noperiod

  ranktable$link <- as.character(ranktable$link)
  ranktable$noperiod <- as.numeric(ranktable$noperiod)
  ranktable$ref_period <- as.character(ranktable$ref_period)
  ranktable$ref_cohort <- as.character(ranktable$ref_cohort)

  if(rankhow=='none') {
    return(ranktable)
  } else if(rankhow=='pcAIC') {
    return(ranktable[order(ranktable$is.pcAIC),])
  } else if(rankhow=='MAE') {
    return(ranktable[order(ranktable$oos.MAE),])
  } else if(rankhow=='RMSE') {
    return(ranktable[order(ranktable$oos.RMSE),])
  } else if(rankhow=='driftSE') {
    return(ranktable[order(ranktable$driftSE),])
  }
}
