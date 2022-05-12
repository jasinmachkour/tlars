#' Prints a summary of the results stored in a C++ object of class tlars_cpp.
#'
#' Prints a summary of the results stored in a C++ object of class tlars_cpp
#' (see [tlars_model] for details), i.e., selected variables, computation time,
#' and number of included dummies.
#'
#' @param x Object of the class tlars_cpp. See [tlars_model] for details.
#' @param ... Ignored. Only added to keep structure of generic [print] function.
#'
#' @importFrom stats rnorm
#'
#' @export
#'
#' @seealso [tlars_model].
#'
#' @examples
#' data('Gauss_data')
#' X = Gauss_data$X
#' y = drop(Gauss_data$y)
#' p = ncol(X)
#' n = nrow(X)
#' dummies = matrix(stats::rnorm(n * p), nrow = n, ncol = p)
#' X_D = cbind(X, dummies)
#' mod_tlars = tlars_model(X = X_D, y = y, num_dummies = ncol(dummies))
#' tlars(model = mod_tlars, T_stop = 3, early_stop = TRUE)
#' print(mod_tlars)
print.Rcpp_tlars_cpp = function(x,
                                ...) {
  # Checking whether LARS or Lasso are used.
  # Plot is only generated for LARS!
  method_type = x$type
  stopifnot(
    "Info message is only generated for LARS, not Lasso!
            Set type = 'lar' when creating an object of class tlars_cpp!" =
      method_type == "lar"
  )

  # Numerical zero
  eps = .Machine$double.eps

  # Retrieve data to be printed from C++ object of class tlars_cpp
  T_stop = x$get_num_active_dummies()
  num_dummies = x$get_num_dummies()
  beta = x$get_beta()

  # Number of original variables (without dummies)
  p = length(beta) - num_dummies

  # Get name of C++ object passed to function argument "model"
  mod_name = deparse(substitute(x))

  # Generate message to be printed
  selected = which(abs(beta[seq(p)]) > eps)
  if (length(selected) == 0) {
    selected_var = "No variables selected"
  } else{
    selected_var = paste(as.character(selected), collapse = ", ")
  }
  mes = paste0(
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
