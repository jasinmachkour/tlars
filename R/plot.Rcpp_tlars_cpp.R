#' Plots the T-LARS solution path
#'
#' Plots the T-LARS solution path stored in C++ objects of class tlars_cpp
#' (see [tlars_cpp] for details) if the object is created with type = "lar"
#' (no plot for type = "lasso").
#'
#' @param x Object of the class tlars_cpp. See [tlars_cpp] for details.
#' @param xlab Label of the x-axis.
#' @param ylab Label of the y-axis.
#' @param include_dummies Logical. If TRUE solution paths of dummies are added to the plot.
#' @param actions Logical. If TRUE axis above plot with indices of added variables
#' (Dummies represented by 'D') along the solution path is added.
#' @param col_selected Color of lines corresponding to selected variables.
#' @param col_dummies Color of lines corresponding to included dummies.
#' @param lty_selected Line type of lines corresponding to selected variables.
#' See [par] for more details.
#' @param lty_dummies Line type of lines corresponding to included dummies.
#' See [par] for more details.
#' @param legend_pos Legend position. See [xy.coords] for more details.
#' @param ... Ignored. Only added to keep structure of generic [plot] function.
#'
#' @return Plots the T-LARS solution path stored in C++ objects of class tlars_cpp (no plot for type = "lasso").
#'
#' @importFrom stats rnorm
#' @importFrom graphics matplot axis abline mtext legend
#' @import methods
#'
#' @export
#'
#' @seealso [tlars_cpp], [plot], [par], and [xy.coords].
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
#' plot(mod_tlars)
plot.Rcpp_tlars_cpp <- function(x,
                                xlab = "# Included dummies",
                                ylab = "Coefficients",
                                include_dummies = TRUE,
                                actions = TRUE,
                                col_selected = "black",
                                col_dummies = "red",
                                lty_selected = "solid",
                                lty_dummies = "dashed",
                                legend_pos = "topleft",
                                ...) {
  # Error control
  if (!methods::is(object = x, class2 = tlars::tlars_cpp)) {
    stop("'x' must be an object of class tlars_cpp.")
  }

  # Checking whether LARS or Lasso are used.
  # Plot is only generated for LARS!
  method_type <- x$type
  stopifnot(
    "Plot is only generated for LARS, not Lasso!
            Set type = 'lar' when creating an object of class tlars_cpp!" =
      method_type == "lar"
  )

  # Retrieve data to be plotted from C++ object of class tlars_cpp
  T_stop <- x$get_num_active_dummies()
  num_dummies <- x$get_num_dummies()
  var_select_path <- x$get_actions()
  beta_path <- do.call(rbind, x$get_beta_path())

  # Number of original variables (without dummies)
  p <- ncol(beta_path) - num_dummies

  # Generate solution path plot of active variables
  dummies_path <- which(var_select_path > p) + 1
  dummies_path_labels <- seq(T_stop)
  graphics::matplot(
    beta_path[, seq(1, p)],
    col = col_selected,
    type = "l",
    xlab = xlab,
    ylab = ylab,
    lty = lty_selected,
    xaxt = "n"
  )
  graphics::axis(
    side = 1,
    at = dummies_path,
    labels = dummies_path_labels,
    ...
  )
  graphics::abline(
    v = dummies_path,
    col = col_dummies,
    lty = 1,
    lwd = 1.3
  )

  # Add dummies solution path to plot
  if (include_dummies) {
    graphics::matlines(beta_path[, seq(p + 1, p + num_dummies)],
      col = col_dummies,
      type = "l",
      lty = lty_dummies
    )
  }

  # Add axis above plot to indicate index of added or removed variables
  # (added dummies are indicated with 'D')
  if (actions) {
    var_select_path_positions <- seq(2, length(var_select_path) + 1)
    var_select_path_labels <- var_select_path
    var_select_path_labels[var_select_path_labels > p] <- "D"
    graphics::axis(
      side = 3,
      at = var_select_path_positions,
      labels = var_select_path_labels
    )
    graphics::mtext(
      "Index of selected variables (D indicates an included dummy)",
      side = 3,
      line = 3
    )
    graphics::abline(
      v = var_select_path_positions,
      col = "gray",
      lty = 6
    )
  }

  # Add legend to plot if active variables and dummies are plotted
  if (include_dummies && !is.null(legend_pos)) {
    graphics::legend(
      legend_pos,
      legend = c("Original variables", "Dummies"),
      col = c(col_selected, col_dummies),
      lty = c(lty_selected, lty_dummies),
      lwd = rep(1, times = 2)
    )
  }
}
