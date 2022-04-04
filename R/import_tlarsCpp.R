#' tlars
#'
#' tlars package
#'
#' Imports
#' @useDynLib tlars, .registration = TRUE
#' @export tlarsCpp
#' @import methods
#' @import Rcpp
#' @import RcppArmadillo
#' @aliases tlars-package
'tlarsCpp'

Rcpp::loadModule(module = "tlarsCpp", TRUE)
