#' Executes the Terminated-LARS (T-LARS) algorithm
#'
#' Modifies the generic tlars_cpp model by executing the T-LARS algorithm and including the results in the tlars_cpp model.
#'
#' @param model Object of the class tlars_cpp.
#' @param T_stop Number of included dummies after which the random experiments (i.e., forward selection processes) are stopped.
#' @param early_stop Logical. If TRUE, then the forward selection process is stopped after T_stop dummies have been included. Otherwise
#' the entire solution path is computed.
#' @param info If TRUE information about the T-LARS step are printed.
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
#' dummies = matrix(stats::rnorm(n * p), nrow = n, ncol = p)
#' X_D = cbind(X, dummies)
#' mod_tlars = tlars_model(X = X_D, y = y, num_dummies = ncol(dummies))
#' tlars(model = mod_tlars, T_stop = 3, early_stop = TRUE)
#' beta = mod_tlars$get_beta()
#' beta
tlars = function(model,
                 T_stop = 1,
                 early_stop = TRUE,
                 info = TRUE) {
  if (info) {
    # Print information about the T-LARS-step
    mes1 = "Executing T-LARS-step by reference..."
    message(mes1)

    # Execute and time T-LARS-step
    start = Sys.time()
    model$execute_lars_step(T_stop = T_stop, early_stop = early_stop)
    end = Sys.time()
    elapsed = end - start

    # Get name of C++ object passed to function argument "model"
    mod_name = deparse(substitute(model))

    # Print information about the executed T-LARS-step
    if (early_stop) {
      mes2 = paste0(
        "\t\t Finished T-LARS-step(s)... \n",
        "\t\t\t - The results are stored in the C++ object '",
        mod_name,
        "'.\n",
        "\t\t\t - New value of T_stop: ",
        T_stop,
        ".\n",
        "\t\t\t - Time elaped: ",
        round(elapsed, digits = 3),
        " sec."
      )
    } else{
      mes2 = paste0(
        "\t\t Finished T-LARS-step(s). No early stopping! \n",
        "\t\t\t - The results are stored in the C++ object '",
        mod_name,
        "'.\n",
        "\t\t\t - Time elaped: ",
        round(elapsed, digits = 3),
        " sec."
      )
    }
    message(mes2)
  } else{
    # Execute T-LARS-step
    model$execute_lars_step(T_stop = T_stop, early_stop = early_stop)
  }
}
