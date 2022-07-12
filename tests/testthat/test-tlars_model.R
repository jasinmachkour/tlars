# 1
test_that(
  "'lars_state' contains the expected number of lists and that the sub-lists contain the correct numbers of elements",
  {
    # Setup and data generation
    data("Gauss_data")
    X <- Gauss_data$X
    y <- drop(Gauss_data$y)
    p <- ncol(X)
    n <- nrow(X)
    num_dummies <- p
    dummies <-
      matrix(stats::rnorm(n * p), nrow = n, ncol = num_dummies)
    XD <- cbind(X, dummies)

    # Create T-LARS model
    mod_tlars <- tlars_model(
      X = XD,
      y = y,
      num_dummies = num_dummies
    )

    # Execute T-LARS step
    tlars(
      model = mod_tlars,
      T_stop = 3,
      early_stop = TRUE
    )

    # Extract T-LARS state
    lars_state <- mod_tlars$get_all()

    # Corrupt T-LARS state
    lars_state <- lars_state[-c(4)]

    # Tests
    expect_error(
      tlars_model(lars_state = lars_state),
      "'lars_state' has to be a list containing the state variables of an object of class tlars_cpp. It has to be obtained via model$get_all(), where 'model' is the object from which the state variables are extracted.",
      fixed = TRUE
    )
  }
)

# 2
test_that("error control for input X works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)
  p <- ncol(X)
  n <- nrow(X)
  num_dummies <- p
  dummies <-
    matrix(stats::rnorm(n * p), nrow = n, ncol = num_dummies)
  XD <- cbind(X, dummies)
  XD_w_NA <- XD
  XD_w_NA[sample(prod(dim(XD)), size = 100)] <- NA

  # Tests
  expect_error(tlars_model(
    X = drop(XD[, 1]),
    y = y,
    num_dummies = num_dummies
  ),
  "'X' must be a matrix.",
  fixed = TRUE
  )

  expect_error(
    tlars_model(
      X = matrix(as.character(XD), ncol = ncol(XD)),
      y = y,
      num_dummies = num_dummies
    ),
    "'X' only allows numerical values.",
    fixed = TRUE
  )

  expect_error(
    tlars_model(
      X = matrix(as.factor(XD), ncol = ncol(XD)),
      y = y,
      num_dummies = num_dummies
    ),
    "'X' only allows numerical values.",
    fixed = TRUE
  )

  expect_error(
    tlars_model(
      X = XD_w_NA,
      y = y,
      num_dummies = num_dummies
    ),
    "'X' contains NAs. Please remove or impute them before proceeding.",
    fixed = TRUE
  )
})

# 3
test_that("error control for input y works", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)
  p <- ncol(X)
  n <- nrow(X)
  num_dummies <- p
  dummies <-
    matrix(stats::rnorm(n * p), nrow = n, ncol = num_dummies)
  XD <- cbind(X, dummies)
  y_w_NA <- y
  y_w_NA[sample(length(y), size = 10)] <- NA

  # Tests
  expect_error(
    tlars_model(
      X = XD,
      y = cbind(y, y),
      num_dummies = num_dummies
    ),
    "'y' must be a vector.",
    fixed = TRUE
  )

  expect_error(
    tlars_model(
      X = XD,
      y = as.character(y),
      num_dummies = num_dummies
    ),
    "'y' only allows numerical values.",
    fixed = TRUE
  )

  expect_error(
    tlars_model(
      X = XD,
      y = matrix(as.factor(y), ncol = 1),
      num_dummies = num_dummies
    ),
    "'y' only allows numerical values.",
    fixed = TRUE
  )

  expect_error(
    tlars_model(
      X = XD,
      y = y_w_NA,
      num_dummies = num_dummies
    ),
    "'y' contains NAs. Please remove or impute them before proceeding.",
    fixed = TRUE
  )

  expect_error(
    tlars_model(
      X = rbind(X, dummies),
      y = y,
      num_dummies = num_dummies
    ),
    "Number of rows in X does not match length of y.",
    fixed = TRUE
  )
})

# 4
test_that(
  "input value for 'num_dummies' is an integer larger or equal to 1 and smaller than the total number of predictors in X",
  {
    # Setup and data generation
    data("Gauss_data")
    X <- Gauss_data$X
    y <- drop(Gauss_data$y)
    p <- ncol(X)
    n <- nrow(X)
    num_dummies <- p
    dummies <-
      matrix(stats::rnorm(n * p), nrow = n, ncol = num_dummies)
    XD <- cbind(X, dummies)

    # Tests
    expect_error(
      tlars_model(
        X = XD,
        y = y,
        num_dummies = num_dummies + 1e-4
      ),
      "'num_dummies' must be an integer larger or equal to 1 and smaller than the total number of predictors in X. This integer must be the number of dummy predictors appended to the right side of the orginal predictor matrix.",
      fixed = TRUE
    )

    expect_error(
      tlars_model(
        X = XD,
        y = y,
        num_dummies = 0
      ),
      "'num_dummies' must be an integer larger or equal to 1 and smaller than the total number of predictors in X. This integer must be the number of dummy predictors appended to the right side of the orginal predictor matrix.",
      fixed = TRUE
    )

    expect_error(
      tlars_model(
        X = XD,
        y = y,
        num_dummies = p + num_dummies
      ),
      "'num_dummies' must be an integer larger or equal to 1 and smaller than the total number of predictors in X. This integer must be the number of dummy predictors appended to the right side of the orginal predictor matrix.",
      fixed = TRUE
    )
  }
)

# 5
test_that("user is warned when setting standardize = FALSE", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)
  p <- ncol(X)
  n <- nrow(X)
  num_dummies <- p
  dummies <-
    matrix(stats::rnorm(n * p), nrow = n, ncol = num_dummies)
  XD <- cbind(X, dummies)

  # Tests
  expect_warning(
    tlars_model(
      X = XD,
      y = y,
      num_dummies = num_dummies,
      standardize = FALSE
    ),
    "'standardize' should be TRUE for the T-LARS algorithm. Since you set standardize = FALSE, we hope that you have a good reason for doing that!",
    fixed = TRUE
  )
})

# 6
test_that("'type' is either 'lar' or 'lasso'", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)
  p <- ncol(X)
  n <- nrow(X)
  num_dummies <- p
  dummies <-
    matrix(stats::rnorm(n * p), nrow = n, ncol = num_dummies)
  XD <- cbind(X, dummies)

  # Tests
  expect_error(
    tlars_model(
      X = XD,
      y = y,
      num_dummies = num_dummies,
      type = "method"
    ),
    "'type' must be one of 'lar', 'lasso'.",
    fixed = TRUE
  )
})

# 7
test_that("output is a C++ object of class tlars_cpp", {
  # Setup and data generation
  data("Gauss_data")
  X <- Gauss_data$X
  y <- drop(Gauss_data$y)
  p <- ncol(X)
  n <- nrow(X)
  num_dummies <- p
  dummies <-
    matrix(stats::rnorm(n * p), nrow = n, ncol = num_dummies)
  XD <- cbind(X, dummies)

  # Create T-LARS model
  mod_tlars <- tlars_model(
    X = XD,
    y = y,
    num_dummies = num_dummies
  )

  # Tests
  expect_true(methods::is(object = mod_tlars, class2 = tlars::tlars_cpp))
})

# 8
test_that("creating a T-LARS model also works for low-dimensional data (i.e., fewer variables than samples)", {
  # Setup and data generation
  n <- 300
  p <- 100
  X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  beta <- c(rep(5, times = 3), rep(0, times = p - 3))
  y <- X %*% beta + stats::rnorm(n)
  num_dummies <- p
  dummies <-
    matrix(stats::rnorm(n * p), nrow = n, ncol = num_dummies)
  XD <- cbind(X, dummies)

  # Tests
  expect_error(
    tlars_model(
      X = XD,
      y = y,
      num_dummies = num_dummies
    ),
    NA
  )
})
