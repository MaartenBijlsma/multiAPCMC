#'
#' @title multiAPCMC.multipredtest
#'
#' @description compared model predictions to out of sample data. The function \code{multiAPCMC.multipredtest} is used within the ranktable function, so users are not expected to use this function directly.
#'
#' @param multiAPCMC.multipred.object The list of predictions of the fitted APCMC models that is produced by \code{multiAPCMC.multipred}.
#' @param oos.data a dataframe containing information on future periods, cases (e.g. incidence) and person-years (exposed population) in those years. Should have the same resolution as the data on which the models were also fitted (e.g. 1 year age group by 1 year periods, or 5-year age groups by 5-year periods, and so forth) and should have the columns Period, Age, cases, and PY (person-years). In other words, this is the same as the 'futuredata' argument from multiAPCMC.multipred, except with an additional 'cases' column.
#' @param validate.what which data should be used to test the predictions against? \code{oos} stands for 'out of sample' and hence only compares fit of predictions to data from the \code{oos.data}. The option \code{not periodfitted} refers to all data that the model was not fitted on. For example, if you have 20 years of data to fit the model on, but had \code{noperiod} of 10, then the first 10 years of that data are used, in addition to the \code{oos.data}, to compare the model against (essentially, you forecast backwards). The option \code{all} compares to everything, including the period that the model was fitted on. Hence, this last option is not recommended.
#'
#' @return returns a nested list with root mean squared error (RMSE) and mean absolute error (MAE) values.
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
#' multiAPCMC.multipredtest(multipred.example,
#'                          multiAPCMC.example.futuredata,
#'                          validate.what="oos")
#' # note that multiAPCMC.example.futuredata contains 'cases' so that we can use it both to predict
#' # and to subsequently test the fit of our model to the cases in this future data.
multiAPCMC.multipredtest <- function(multiAPCMC.multipred.object,oos.data,validate.what="oos") {

  # check if number of age groups in futureAPC is equal to the number of age groups in the apcobject
  age_oos.data <- length(unique(oos.data$Age))
  age_preds <- length(unique(multiAPCMC.multipred.object$np.pred[[1]][[1]][[1]][[1]]$Age))
  per_oos.data <- length(unique(oos.data$Period))
  per_preds <- length(unique(multiAPCMC.multipred.object$np.pred[[1]][[1]][[1]][[1]]$Period))
  if(age_oos.data != age_oos.data ) {
    stop("The number of age categories in oos.data are not equal to the number of age categories in the predictions")
  }
  if(per_oos.data != per_oos.data ) {
    stop("The number of period categories in oos.data are not equal to the number of period categories in the predictions")
  }
  # check also if age and period are sorted (ascending)
  if(is.unsorted(unique(oos.data$Age)) | is.unsorted(unique(oos.data$Period))) {
    stop("sort the oos.data by period and age (in that order), this is needed because we will assume
         that the youngest age category in the futureAPC object corresponds to the youngest in the training data
         and the oldest corresponds to the oldest, etc. The same logic applies with period")
  }

  vec.link <- multiAPCMC.multipred.object$vec.link
  vec.noperiod <- multiAPCMC.multipred.object$vec.noperiod
  vec.refper <- multiAPCMC.multipred.object$vec.refper
  vec.refcoh <- multiAPCMC.multipred.object$vec.refcoh

  depth = c(length(vec.link),
            length(vec.noperiod),
            length(vec.refper),
            length(vec.refcoh)
  )
  np.valid <- multiAPCMC.reclist(depth)

  for(l in 1:length(vec.link)) {

    for(n in 1:length(vec.noperiod)) {

      for(p in 1:length(vec.refper)) {

        for(c in 1:length(vec.refcoh)) {

          if(validate.what=="oos") {

            # only select observed==0
            temp <- multiAPCMC.multipred.object$np.pred[[l]][[n]][[p]][[c]][multiAPCMC.multipred.object$np.pred[[l]][[n]][[p]][[c]]$observed==0,]
            oos.dis <- temp$pred_cases - oos.data$cases
            rmse <- sqrt(sum((oos.dis)^2)/length(oos.dis))
            mae <- sum(abs(oos.dis))/length(oos.dis)

            temp <- NULL
            np.valid[[l]][[n]][[p]][[c]] <-  c(rmse,mae)

          } else if(validate.what=="not periodfitted") {

            # only select periodfitted==0
            temp <- multiAPCMC.multipred.object$np.pred[[l]][[n]][[p]][[c]][multiAPCMC.multipred.object$np.pred[[l]][[n]][[p]][[c]]$periodfitted==0,]

            temp.dis <- temp$pred_cases[temp$observed==1 & temp$periodfitted==0] - temp$cases[temp$observed==1 & temp$periodfitted==0]
            oos.dis <- temp$pred_cases[temp$observed==0 & temp$periodfitted==0] - oos.data$cases
            dis.tot <- c(temp.dis,oos.dis)
            rmse <- sqrt(sum((dis.tot)^2)/length(dis.tot))
            mae <- sum(abs(dis.tot))/length(dis.tot)

            temp <- NULL
            np.valid[[l]][[n]][[p]][[c]] <-  c(rmse,mae)

          }

          else if(validate.what=="all") {

            # combination  of the above two
            # observed (i.e. training data)==1 + oos
            temp <- multiAPCMC.multipred.object$np.pred[[l]][[n]][[p]][[c]][multiAPCMC.multipred.object$np.pred[[l]][[n]][[p]][[c]]$observed==1,]
            temp.dis <- temp$pred_cases - temp$cases
            oos.dis <- multiAPCMC.multipred.object$np.pred[[l]][[n]][[p]][[c]][multiAPCMC.multipred.object$np.pred[[l]][[n]][[p]][[c]]$observed==0 & multiAPCMC.multipred.object$np.pred[[l]][[n]][[p]][[c]]$periodfitted==0,]$pred_cases - oos.data$cases
            dis.tot <- c(temp.dis,oos.dis)
            rmse <- sqrt(sum((dis.tot)^2)/length(dis.tot))
            mae <- sum(abs(dis.tot))/length(dis.tot)

            temp <- NULL
            np.valid[[l]][[n]][[p]][[c]] <-  c(rmse,mae)

          }

        }

      }

    }

  }

  return(np.valid)
}
