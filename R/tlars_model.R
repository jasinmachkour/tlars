#' Creates a Terminated-LARS (tlars) object
#'
#' Creates an object of the class tlars_cpp.
#'
#' @param lars_state Variables associated with previous stopping point (necessary to restart
#' the forward selection process exactly where it was previously terminated). The lars_state
#' is extracted from an object of class tlars_cpp via $get_all() and is only required when the
#' object (or its pointer) of class tlars_cpp is deleted or got lost in another R session (e.g.,
#' in parallel processing).
#' @param X Real valued Predictor matrix.
#' @param y Response vector.
#' @param num_knocks Number of knockoffs that are appended to the predictor matrix.
#' @param verbose Logical. If TRUE progress in computations is shown.
#' @param intercept Logical. If TRUE an intercept is included.
#' @param standardize Logical. If TRUE the predictors are standardized and the response is centered.
#' @param type 'lar' for 'LARS' and 'lasso' for Lasso.
#' @param info Logical. If TRUE information about the T-LARS step is printed.
#'
#' @return Object of the class tlars_cpp.
#'
#' @importFrom stats rnorm
#'
#' @export
#'
#' @examples
#' data('Gauss_data')
#' X = Gauss_data$X
#' y = drop(Gauss_data$y)
#' p = ncol(X)
#' n = nrow(X)
#' knocks = matrix(stats::rnorm(n * p), nrow = n, ncol = p)
#' XK = cbind(X, knocks)
#' mod_tlars = tlars_model(X = XK, y = y, num_knocks = ncol(knocks))
#' mod_tlars
tlars_model = function(lars_state,
                       X,
                       y,
                       num_knocks,
                       verbose = FALSE,
                       intercept = FALSE,
                       standardize = TRUE,
                       type = 'lar',
                       info = TRUE) {
  if (missing(lars_state)) {
    mod_tlars = new(
      tlars::tlars_cpp,
      X = X,
      y = drop(y),
      verbose = verbose,
      intercept = intercept,
      standardize = standardize,
      num_knocks = num_knocks,
      type = type
    )
  } else{
    mod_tlars = new(tlars::tlars_cpp,
                    lars_state = lars_state)
  }

  if (info) {
    # Get name of predictor matrix passed to function argument "X"
    pred_mat_name = deparse(substitute(X))

    # Print information about the generated T-LARS model
    mes = paste0(
      "Created an object of class tlars_cpp... \n",
      "\t\t The first ",
      ncol(X) - num_knocks,
      " predictors in '",
      pred_mat_name,
      "' are the original predictors and \n",
      "\t\t the last ",
      num_knocks,
      " predictors are knockoffs."
    )
    message(mes)
  }

  return(mod_tlars)
}
