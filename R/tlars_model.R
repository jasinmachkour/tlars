#' Creates a Terminating-LARS (T-LARS) object
#'
#' Creates an object of the class tlars_cpp.
#'
#' @param lars_state List of variables associated with previous T-LARS step (necessary to restart
#' the forward selection process exactly where it was previously terminated). The lars_state
#' is extracted from an object of class tlars_cpp via get_all() and is only required when the
#' object (or its pointer) of class tlars_cpp is deleted or got lost in another R session (e.g.,
#' in parallel processing).
#' @param X Real valued predictor matrix.
#' @param y Response vector.
#' @param num_dummies Number of dummies that are appended to the predictor matrix.
#' @param verbose Logical. If TRUE progress in computations is shown when performing T-LARS steps on the created model.
#' @param intercept Logical. If TRUE an intercept is included.
#' @param standardize Logical. If TRUE the predictors are standardized and the response is centered.
#' @param type 'lar' for 'LARS' and 'lasso' for Lasso.
#' @param info Logical. If TRUE and object is not recreated from previous T-LARS state, then information about the created object is printed.
#'
#' @return Object of the class tlars_cpp.
#'
#' @importFrom stats rnorm
#'
#' @export
#'
#' @examples
#' data("Gauss_data")
#' X <- Gauss_data$X
#' y <- drop(Gauss_data$y)
#' p <- ncol(X)
#' n <- nrow(X)
#' num_dummies <- p
#' dummies <- matrix(stats::rnorm(n * p), nrow = n, ncol = num_dummies)
#' XD <- cbind(X, dummies)
#' mod_tlars <- tlars_model(X = XD, y = y, num_dummies = num_dummies)
#' mod_tlars
tlars_model <- function(lars_state,
                        X,
                        y,
                        num_dummies,
                        verbose = FALSE,
                        intercept = FALSE,
                        standardize = TRUE,
                        type = "lar",
                        info = TRUE) {
  # Error control
  if (!missing(lars_state)) {
    if (length(lars_state) != 4) {
      stop(
        "'lars_state' has to be a list containing the state variables of an object of class tlars_cpp. It has to be obtained via model$get_all(), where 'model' is the object from which the state variables are extracted."
      )
    } else {
      if (length(lars_state$l1) != 20 ||
        length(lars_state$l2) != 20 ||
        length(lars_state$l3) != 20 ||
        length(lars_state$l4) != 2) {
        stop(
          "'lars_state' has to be a list containing the state variables of an object of class tlars_cpp. It has to be obtained via model$get_all(), where 'model' is the object from which the state variables are extracted."
        )
      }
    }
  } else {
    if (!is.matrix(X)) {
      stop("'X' must be a matrix.")
    }

    if (!is.numeric(X)) {
      stop("'X' only allows numerical values.")
    }

    if (anyNA(X)) {
      stop("'X' contains NAs. Please remove or impute them before proceeding.")
    }

    if (!is.vector(drop(y))) {
      stop("'y' must be a vector.")
    }

    if (!is.numeric(y)) {
      stop("'y' only allows numerical values.")
    }

    if (anyNA(y)) {
      stop("'y' contains NAs. Please remove or impute them before proceeding.")
    }

    if (nrow(X) != length(drop(y))) {
      stop("Number of rows in X does not match length of y.")
    }

    if (num_dummies %% 1 != 0 ||
      num_dummies < 1 ||
      num_dummies > ncol(X) - 1) {
      stop(
        "'num_dummies' must be an integer larger or equal to 1 and smaller than the total number of predictors in X. This integer must be the number of dummy predictors appended to the right side of the orginal predictor matrix."
      )
    }

    if (!standardize) {
      warning(
        "'standardize' should be TRUE for the T-LARS algorithm. Since you set standardize = FALSE, we hope that you have a good reason for doing that!"
      )
    }

    if (!(type %in% c("lar", "lasso"))) {
      stop("'type' must be one of 'lar', 'lasso'.")
    }
  }

  # Create C++ object of class tlars_cpp
  if (missing(lars_state)) {
    mod_tlars <- new(
      tlars::tlars_cpp,
      X = X,
      y = drop(y),
      verbose = verbose,
      intercept = intercept,
      standardize = standardize,
      num_dummies = num_dummies,
      type = type
    )
  } else {
    mod_tlars <- new(tlars::tlars_cpp,
      lars_state = lars_state
    )
  }

  if (info && missing(lars_state)) {
    # Get name of predictor matrix passed to function argument "X"
    pred_mat_name <- deparse(substitute(X))

    # Print information about the generated T-LARS model
    mes <- paste0(
      "Created an object of class tlars_cpp... \n",
      "\t\t The first p = ",
      ncol(X) - num_dummies,
      " predictors in '",
      pred_mat_name,
      "' are the original predictors and \n",
      "\t\t the last num_dummies = ",
      num_dummies,
      " predictors are dummies"
    )
    message(mes)
  }

  return(mod_tlars)
}
