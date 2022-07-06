#' @name tlars_cpp
#'
#' @title Exposes the C++ class tlars_cpp to R
#'
#' @description Type 'tlars_cpp' in the console to see the constructors, variables, and methods of the class tlars_cpp.
#'
#' @field Constructor: new - Creates a new object of the class tlars_cpp.
#' @param X Real valued predictor matrix.
#' @param y Response vector.
#' @param verbose Logical. If TRUE progress in computations is shown.
#' @param intercept Logical. If TRUE an intercept is included.
#' @param standardize Logical. If TRUE the predictors are standardized and the response is centered.
#' @param num_dummies Number of dummies that are appended to the predictor matrix.
#' @param type Type of used algorithm (currently possible choices: 'lar' or 'lasso').
#'
#' @field Constructor: new - Re-creates an object of the class tlars_cpp based on a list of class variables that is obtained via get_all().
#' @param  lars_state Input list that was extracted from a previous tlars_cpp object using get_all().
#'
#' @field Method: execute_lars_step - Executes LARS steps until a stopping-condition is satisfied.
#' @param T_stop Number of included dummies after which the random experiments (i.e., forward selection processes) are stopped.
#' @param early_stop Logical. If TRUE, then the forward selection process is stopped after T_stop dummies have been included.
#' Otherwise the entire solution path is computed.
#'
#' @field Method: get_beta - Returns the estimate of the beta vector.
#' @field Method: get_beta_path - Returns a a matrix with the estimates of the beta vectors at all steps.
#' @field Method: get_num_active - Returns the number of active predictors.
#' @field Method: get_num_active_dummies - Returns the number of dummy variables that have been included.
#' @field Method: get_num_dummies - Returns the number of dummy predictors.
#' @field Method: get_actions - Returns the indices of added/removed variables along the solution path.
#' @field Method: get_df - Returns the degrees of freedom at each step which is given by number of active variables (+1 if intercept is true).
#' @field Method: get_R2 - Returns the R^2 statistic at each step.
#' @field Method: get_RSS - Returns the residual sum of squares at each step.
#' @field Method: get_Cp - Returns the Cp-statistic at each step.
#' @field Method: get_lambda - Returns the lambda-values (penalty parameters) at each step along the solution path.
#' @field Method: get_entry - Returns the first entry/selection steps of the predictors along the solution path.
#' @field Method: get_norm_X - Returns the L2-norm of the predictors.
#' @field Method: get_mean_X - Returns the sample means of the predictors.
#' @field Method: get_mean_y - Returns the sample mean of the response y.
#' @field Method: get_all - Returns all class variables: This list can be used as an input to the constructor to re-create an object of class tlars_cpp.
#'
#' @return No return value. Exposes the C++ class tlars_cpp to R.
#'
#' @useDynLib tlars, .registration = TRUE
#' @import methods
#' @import Rcpp
#'
#' @export tlars_cpp
#'
#' @examples
#' data("Gauss_data")
#' X <- Gauss_data$X
#' y <- drop(Gauss_data$y)
#' p <- ncol(X)
#' n <- nrow(X)
#' dummies <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
#' XD <- cbind(X, dummies)
#' mod_tlars <- tlars_model(X = XD, y = y, num_dummies = ncol(dummies))
#' tlars(model = mod_tlars, T_stop = 3, early_stop = TRUE)
#'
#' mod_tlars$get_beta()
#'
#' # mod_tlars$get_beta_path()
#' # mod_tlars$get_num_active()
#' # mod_tlars$get_num_active_dummies()
#' # mod_tlars$get_num_dummies()
#' # mod_tlars$get_actions()
#' # mod_tlars$get_df()
#' # mod_tlars$get_R2()
#' # mod_tlars$get_RSS()
#' # mod_tlars$get_Cp()
#' # mod_tlars$get_lambda()
#' # mod_tlars$get_entry()
#' # mod_tlars$get_norm_X()
#' # mod_tlars$get_mean_X()
#' # mod_tlars$get_mean_y()
#' # mod_tlars$get_all()
Rcpp::loadModule(module = "tlars_cpp", TRUE)
