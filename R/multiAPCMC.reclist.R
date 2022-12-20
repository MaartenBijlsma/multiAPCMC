#'
#' @title multiAPCMC.reclist
#'
#' @description makes empty nested lists. This is used by various multiAPCMC functions, but is not expected to be called directly by the user.
#'
#' @param dimensions a vector containing the dimensions of the nested list.
#'
#' @return returns an empty set of nested lists.
#' @export
#'
#' @examples
#'
#' vec.link <- c('power5','poisson')
#' vec.noperiod <- c(5,10,15,21)
#' vec.refper <- c('extremes','outer','center', 'first middle','middle last',
#'                 'first second', 'penultimate last')
#' vec.refcoh <- c('extremes','outer','center', 'first middle','middle last',
#'                 'first second', 'penultimate last')
#'
#'   depth = c(length(vec.link),
#'   length(vec.noperiod),
#'   length(vec.refper),
#'   length(vec.refcoh))
#'   multiAPCMC.reclist(depth)
#'
multiAPCMC.reclist <- function(dimensions){
  if(length(dimensions) == 1){
    vector("list", dimensions)
  } else {
    lapply(1:dimensions[1], function(...) multiAPCMC.reclist(dimensions[-1]))
  }
}
