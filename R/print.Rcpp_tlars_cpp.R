#' Prints a summary of the results stored in a C++ object of class tlars_cpp.
#'
#' Prints a summary of the results stored in a C++ object of class tlars_cpp
#' (see [tlars_model] for details), i.e., selected variables, computation time,
#' and number of included knockoffs.
#'
#' @param tlars_model Object of the class tlars_cpp. See [tlars_model] for details.
#' @param ... Additional parameters to print. See [print] for details.
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
#' knocks = matrix(stats::rnorm(n * p), nrow = n, ncol = p)
#' XK = cbind(X, knocks)
#' mod_tlars = tlars_model(X = XK, y = y, num_knocks = ncol(knocks))
#' tlars(model = mod_tlars, T_stop = 3, early_stop = TRUE)
#' print(mod_tlars)
print.Rcpp_tlars_cpp = function(tlars_model, ...) {
  # Checking whether LARS or Lasso are used.
  # Plot is only generated for LARS!
  method_type = tlars_model$type
  stopifnot(
    "Info message is only generated for LARS, not Lasso!
            Set type = 'lar' when creating an object of class tlars_cpp!" =
      method_type == "lar"
  )

  # Numerical zero
  eps = .Machine$double.eps

  # Retrieve data to be printed from C++ object of class tlars_cpp
  T_stop = tlars_model$get_num_active_knocks()
  num_knocks = tlars_model$get_num_knocks()
  beta = tlars_model$get_beta()

  # Number of original variables (without knockoffs)
  p = length(beta) - num_knocks

  # Get name of C++ object passed to function argument "model"
  mod_name = deparse(substitute(tlars_model))

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
    "\t - Number of knockoffs: ",
    num_knocks,
    ".\n",
    "\t - Number of included knockoffs: ",
    T_stop,
    ".\n",
    "\t - Selected variables: ",
    selected_var,
    "."
  )

  # Print message
  cat(mes)
}
