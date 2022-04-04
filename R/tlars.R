#' Executes the Terminated-LARS (tlars) algorithm
#'
#' Modifies the generic tlarsCpp object by executing the tlars algorithm and including the results in the tlarsCpp object.
#'
#' @param obj Object of the class tlarsCpp.
#' @param T_stop Number of included knockoffs after which the random experiments (i.e., forward selection processes) are stopped.
#' @param earlyStop If TRUE, then the forward selection process is stopped after T_stop knockoffs have been included. Otherwise the entire solution path is computed.
#'
#' @importFrom stats rnorm
#'
#' @export
#'
#' @examples
#' data('Gauss_data')
#' X = Gauss_data$X
#' y = c(Gauss_data$y)
#' p = ncol(X)
#' n = nrow(X)
#' knock = matrix(stats::rnorm(n * p), nrow = n, ncol = p)
#' XK = cbind(X, knock)
#' larsObj = newObj_tlars(X = XK, y = y, L_val = p)
#' tlars(obj = larsObj, T_stop = 3, earlyStop = TRUE)
#' beta = larsObj$getLastBeta()
#' beta
tlars = function(obj,
                 T_stop = 1,
                 earlyStop = TRUE) {
  obj$executeLarsStep(T_stop = T_stop, earlyStop = earlyStop)
}
