#' Creates a Terminated-LARS (tlars) object
#'
#' Creates an object of the class tlarsCpp.
#'
#' @param lars_state Variables associated with previous stopping point (necessary to restart forward selection selection process exactly where it was previously terminated). The lars_state is only required when the object (or its pointer) of the class tlarsCpp is deleted or got lost in another R session (e.g., in parallel processing).
#' @param X Real valued Predictor matrix.
#' @param y Response vector.
#' @param L_val Number of knockoffs that are appended to the predictor matrix.
#' @param verbose If TRUE progress in computations is shown.
#' @param intercept If TRUE an intercept is included.
#' @param normalize If TRUE the predictors are standardized and the response is centered.
#' @param type 'lar' for 'LARS' and 'lasso' for Lasso.
#'
#' @return Object of the class tlarsCpp.
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
#' larsObj = newObj_tlars(X = X, y = y, L_val = p)
#' larsObj
newObj_tlars = function(lars_state,
                        X,
                        y,
                        L_val,
                        verbose = FALSE,
                        intercept = FALSE,
                        normalize = TRUE,
                        type = 'lar') {
  if (missing(lars_state)) {
    tlarsCppMod = new(
      tlars::tlarsCpp,
      X = X,
      y = drop(y),
      verbose = verbose,
      intercept = intercept,
      normalize = normalize,
      L_val = L_val,
      type = type
    )
  } else{
    tlarsCppMod = new(tlars::tlarsCpp,
                      lars_state = lars_state)
  }
  return(tlarsCppMod)
}
