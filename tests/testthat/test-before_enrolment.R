test_that("estimate_lambda_e works", {

  n_e           <- 0 # Enrolment at exactly one year old.
  yearly_lambda <- rep(0.1, 26)
  start_index   <- 1
  init_weight   <- 0
  final_weight  <- 0

  actual <- estimate_lambda_e(n_e, yearly_lambda, start_index,
                              init_weight, final_weight)

  expected <- 0

  expect_equal(actual, expected)

  n_e          <- 1 # Enrolment at exactly two years old.
  init_weight  <- 0 # This individual became susceptible on 1st Jan
  final_weight <- 0 # Enrolment on 31st Dec

  actual <- estimate_lambda_e(n_e, yearly_lambda, start_index,
                              init_weight, final_weight)


  expected    <- 0.1

  expect_equal(actual, expected)
  #-----------------------------------------------------------------------------
  n_e          <- 1 # Enrolment at exactly two years old.
  init_weight  <- 0.3 # This individual became susceptible on 20th Apr
  final_weight <- 0.3 # Enrolment on 25th Nov

  actual <- estimate_lambda_e(n_e, yearly_lambda, start_index,
                              init_weight, final_weight)


  expected    <- 0.1 * 0.4

  expect_equal(actual, expected)
  #-----------------------------------------------------------------------------
  n_e          <- 2 # Enrolment at exactly three years old.
  init_weight  <- 0 # This individual became susceptible on 1st Jan
  final_weight <- 0 # Enrolment on 31st Dec

  actual <- estimate_lambda_e(n_e, yearly_lambda, start_index,
                              init_weight, final_weight)


  expected    <- 0.2

  expect_equal(actual, expected)

  #-----------------------------------------------------------------------------
  n_e          <- 3    # Enrolment at exactly two years old.
  init_weight  <- 0.25 # This individual became susceptible on 1st Jan
  final_weight <- 0.75 # Enrolment on 31st Dec

  actual <- estimate_lambda_e(n_e, yearly_lambda, start_index,
                              init_weight, final_weight)


  expected    <- 0.2

  expect_equal(actual, expected)
})

test_that("estimate_lambda_e works with one year of enrolment", {

  actual <- estimate_lambda_e(n_e           = 1,
                              yearly_lambda = rep(0.14, 26),
                              start_index   = 4,
                              init_weight   = 0,
                              final_weight  = 0)

  expected <- 0.14

  expect_equal(actual, expected)
})
