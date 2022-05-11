#' Plots the T-LARS solution path
#'
#' Plots the T-LARS solution path stored in C++ objects of class tlars_cpp
#' (see [tlars_model] for details) if the object is created with type = "lar"
#' (no plot for type = "lasso").
#'
#' @param tlars_model Object of the class tlars_cpp. See [tlars_model] for details.
#' @param include_knocks Logical. If TRUE solution paths of knockoffs are added to the plot.
#' @param actions Logical. If TRUE axis above plot with indices of added variables
#' (Knockoffs represented by 'K') along the solution path is added.
#' @param col_selected Color of lines corresponding to selected variables.
#' @param col_knocks Color of lines corresponding to included knockoffs.
#' @param xlab Label of the x-axis.
#' @param ylab Label of the y-axis.
#' @param lty_selected Line type of lines corresponding to selected variables.
#' See [par] for more details.
#' @param lty_knocks Line type of lines corresponding to included knockoffs.
#' See [par] for more details.
#' @param legend_pos Legend position. See [xy.coords] for more details.
#' @param ... Other parameter that control the appearance of the plot. See [plot] and [par].
#'
#' @importFrom stats rnorm
#' @importFrom graphics matplot axis abline mtext legend
#'
#' @export
#'
#' @seealso [tlars_model], [plot], [par], and [xy.coords].
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
#' plot(mod_tlars)
plot.Rcpp_tlars_cpp = function(tlars_model,
                               include_knocks = TRUE,
                               actions = TRUE,
                               col_selected = "black",
                               col_knocks = "red",
                               xlab = '# Included knockoffs',
                               ylab = 'Coefficients',
                               lty_selected = 'solid',
                               lty_knocks = 'dashed',
                               legend_pos = 'topleft',
                               ...) {
  # Checking whether LARS or Lasso are used.
  # Plot is only generated for LARS!
  method_type = tlars_model$type
  stopifnot(
    "Plot is only generated for LARS, not Lasso!
            Set type = 'lar' when creating an object of class tlars_cpp!" =
      method_type == "lar"
  )

  # Retrieve data to be plotted from C++ object of class tlars_cpp
  T_stop = tlars_model$get_num_active_knocks()
  num_knocks = tlars_model$get_num_knocks()
  var_select_path = tlars_model$get_actions()
  beta_path = do.call(rbind, tlars_model$get_beta_path())

  # Number of original variables (without knockoffs)
  p = ncol(beta_path) - num_knocks

  # Generate solution path plot of active variables
  knock_path = which(var_select_path > p) + 1
  knock_path_labels = seq(T_stop)
  graphics::matplot(
    beta_path[, seq(1, p)],
    col = col_selected,
    type = 'l',
    xlab = xlab,
    ylab = ylab,
    lty = lty_selected,
    xaxt = 'n',
    ...
  )
  graphics::axis(side = 1,
                 at = knock_path,
                 labels = knock_path_labels)
  graphics::abline(
    v = knock_path,
    col = col_knocks,
    lty = 1,
    lwd = 1.3
  )

  # Add knockoff solution path to plot
  if (include_knocks) {
    graphics::matlines(beta_path[, seq(p + 1, p + num_knocks)],
                       col = col_knocks,
                       type = 'l',
                       lty = lty_knocks,
                       ...)
  }

  # Add axis above plot to indicate index of added or removed variables
  # (added and removed knockoffs indicated with 'K' and "-K", respectively)
  if (actions) {
    var_select_path_positions = seq(2, length(var_select_path) + 1)
    var_select_path_labels = var_select_path
    var_select_path_labels[var_select_path_labels > p] = 'K'
    graphics::axis(side = 3,
                   at = var_select_path_positions,
                   labels = var_select_path_labels)
    graphics::mtext(
      "Index of selected variables (K indicates an included knockoff)",
      side = 3,
      line = 3
    )
    graphics::abline(v = var_select_path_positions, col = "gray", lty = 6)
  }

  # Add legend to plot if active variables and knockoffs are plotted
  if (include_knocks && !is.null(legend_pos)) {
    legend(
      legend_pos,
      legend = c('Original variables' , 'Knockoffs'),
      col = c(col_selected, col_knocks),
      lty = c(lty_selected, lty_knocks),
      lwd = rep(1, times = 2)
    )
  }
}
