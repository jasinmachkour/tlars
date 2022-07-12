# 1
test_that(
  "T-LARS model is an object of class tlars_cpp and stays an object of the same class after T-LARS steps",
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

    # Tests
    expect_true(methods::is(object = mod_tlars, class2 = tlars::tlars_cpp))

    # Execute T-LARS step
    tlars(
      model = mod_tlars,
      T_stop = 3,
      early_stop = TRUE
    )

    # Tests
    expect_true(methods::is(object = mod_tlars, class2 = tlars::tlars_cpp))
  }
)

# 2
test_that("the input value of 'T_stop' is valid", {
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
  expect_error(
    tlars(
      model = mod_tlars,
      T_stop = 0
    ),
    paste0(
      "Value of 'T_stop' not valid. 'T_stop' must be an integer from 1 to ",
      num_dummies,
      "."
    ),
    fixed = TRUE
  )

  expect_error(
    tlars(
      model = mod_tlars,
      T_stop = num_dummies + 1
    ),
    paste0(
      "Value of 'T_stop' not valid. 'T_stop' must be an integer from 1 to ",
      num_dummies,
      "."
    ),
    fixed = TRUE
  )
})

# 3
test_that("the user is informed that the entire solution path is computed if early_stop = FALSE", {
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
  expect_message(
    tlars(
      model = mod_tlars,
      T_stop = 3,
      early_stop = FALSE
    ),
    "'T_stop' is ignored. Computing the entire solution path...",
    fixed = TRUE
  )
})

# 4
test_that("running T-LARS also works for low-dimensional data (i.e., fewer variables than samples)", {
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

  # Create T-LARS model
  mod_tlars <- tlars_model(
    X = XD,
    y = y,
    num_dummies = num_dummies
  )

  # Tests
  expect_error(
    tlars(
      model = mod_tlars,
      T_stop = 3,
      early_stop = TRUE
    ),
    NA
  )
})
