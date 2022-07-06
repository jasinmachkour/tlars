#' Prints a summary of the results stored in a C++ object of class tlars_cpp.
#'
#' Prints a summary of the results stored in a C++ object of class tlars_cpp
#' (see [tlars_cpp] for details), i.e., selected variables, computation time,
#' and number of included dummies.
#'
#' @param x Object of the class tlars_cpp. See [tlars_cpp] for details.
#' @param ... Ignored. Only added to keep structure of generic [print] function.
#'
#' @return Prints a summary of the results stored in a C++ object of class tlars_cpp.
#'
#' @importFrom stats rnorm
#' @import methods
#'
#' @export
#'
#' @seealso [tlars_cpp].
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
#' tlars(model = mod_tlars, T_stop = 3, early_stop = TRUE)
#' print(mod_tlars)
print.Rcpp_tlars_cpp <- function(x,
                                 ...) {

  # Error control
  if (!methods::is(object = x, class2 = tlars::tlars_cpp)) {
    stop("'x' must be an object of class tlars_cpp.")
  }

  # Checking whether LARS or Lasso are used.
  # Plot is only generated for LARS!
  method_type <- x$type
  stopifnot(
    "Info message is only generated for LARS, not Lasso!
            Set type = 'lar' when creating an object of class tlars_cpp!" =
      method_type == "lar"
  )

  # Numerical zero
  eps <- .Machine$double.eps

  # Retrieve data to be printed from C++ object of class tlars_cpp
  T_stop <- x$get_num_active_dummies()
  num_dummies <- x$get_num_dummies()
  beta <- x$get_beta()
  var_select_path <- x$get_actions()

  # Number of original variables (without dummies)
  p <- length(beta) - num_dummies

  # Get name of C++ object passed to function argument 'model'
  mod_name <- deparse(substitute(x))

  # Generate message to be printed
  selected <- var_select_path[var_select_path <= p]
  if (length(selected) == 0) {
    selected_var <- "No variables selected"
  } else {
    selected_var <- paste(as.character(selected), collapse = ", ")
  }
  mes <- paste0(
    "'",
    mod_name,
    "' is a C++ object of class 'tlars_cpp' ... \n",
    "\t - Number of dummies: ",
    num_dummies,
    ".\n",
    "\t - Number of included dummies: ",
    T_stop,
    ".\n",
    "\t - Selected variables: ",
    selected_var,
    "."
  )

  # Print message
  cat(mes)
}
