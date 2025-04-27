test_that("sim_infection_seroneg() works", {

  follow_up_times <- c(16, 17, 18, 19, 20, 21, 22, 24, 25)
  enrolment_year  <- 15

  weights_fu       <- c(0.6684932, 0.6109589, 0.6273973, 0.6147541, 0.6821918,
                        0.6328767, 0.9013699, 0.4027397, 0.8712329)

  enrolment_weight <- 0.295082

  set.seed(123)

  actual <- sim_infection_seroneg(follow_up_times  = follow_up_times,
                                  enrolment_year     = enrolment_year,
                                  lambda             = rep(0.1, 26),
                                  total_buckets      = 4,
                                  weights_fu         = weights_fu,
                                  enrolment_weight   = enrolment_weight,
                                  rho                = 0.05,
                                  rho_v              = 0,
                                  is_vaccinated      = 0,
                                  vac_year           = 15)

  expected <- c(0, 1, 0, 1, 1, 0, 0, 0, 0)

  expect_equal(actual, expected)
})

test_that("simulate_single_year() returns the expected names", {

  output <- simulate_single_year(yr                   = 2,
                                 n_filled_buckets_nat = 0,
                                 n_filled_buckets_vac = 0,
                                 rho                  = 0,
                                 rho_v                = 0,
                                 current_idx          = 1,
                                 cml_lambda           = 0.1,
                                 follow_up_times      = 2,
                                 weights_fu           = 0,
                                 lambda               = rep(0.1, 2),
                                 total_buckets        = 4,
                                 is_vaccinated        = 0,
                                 vac_year             = 15,
                                 inf_counter          = 0)

  actual   <- names(output)
  expected <- c("n_filled_buckets_nat", "n_filled_buckets_vac", "flag_report",
                "inf_counter", "leftover_lambda", "sim_y")

  expect_equal(actual, expected)
})

test_that("simulate_single_year() works", {

  n_filled_buckets_nat <- 0.95
  n_filled_buckets_vac <- 0
  rho                  <- 0.05
  rho_v                <- 0
  current_idx          <- 1
  cml_lambda           <- 0
  inf_counter          <- 0

  years           <- 16:18
  follow_up_times <- 16:18
  vac_year        <- 15
  weights_fu      <- rep(0, 3)
  total_buckets   <- 4
  is_vaccinated   <- 0

  lambda <- rep(0, 26)

  actual <- vector(mode = "numeric", length = length(years))

  for(yr in years)
  {
    output <- simulate_single_year(yr, n_filled_buckets_nat,
                                   n_filled_buckets_vac, rho,
                                   rho_v, current_idx, cml_lambda,
                                   follow_up_times, weights_fu,
                                   lambda, total_buckets, is_vaccinated,
                                   vac_year,
                                   inf_counter = inf_counter)

    n_filled_buckets_nat <- output[[1]]

    actual[yr - min(years) + 1] <- n_filled_buckets_nat
  }

  expected <- c(0.90, 0.85, 0.8)

  expect_equal(actual, expected)

  #-----------------------------------------------------------------------------
  set.seed(1613)

  years           <- 16:23
  follow_up_times <- 16:23
  weights_fu      <- rep(0, 8)

  n_filled_buckets_nat <- 0
  lambda               <- rep(0.1, 26)
  cml_lambda           <- 0.1
  current_idx          <- 1
  inf_counter          <- 0

  actual             <- vector(mode = "numeric", length = length(years))
  actual_inf_counter <- vector(mode = "numeric", length = length(years))

  for(yr in years)
  {
    output <- simulate_single_year(yr, n_filled_buckets_nat,
                                   n_filled_buckets_vac, rho,
                                   rho_v, current_idx, cml_lambda,
                                   follow_up_times, weights_fu, lambda,
                                   total_buckets, is_vaccinated, vac_year,
                                   inf_counter = inf_counter)

    n_filled_buckets_nat <- output[[1]]
    flag_report          <- output[[3]]
    if(flag_report == 1) current_idx <- current_idx + 1
    inf_counter          <- output[[4]]

    actual[yr - min(years) + 1]             <- n_filled_buckets_nat
    actual_inf_counter[yr - min(years) + 1] <- inf_counter
  }

  expected <- c(0.95, 0.90, 0.85, 0.80, 0.75, 1.70, 1.65, 1.60)

  expect_equal(actual, expected)

  expected_inf_counter <- c(1, 1, 1, 1, 1, 2, 2, 2)

  expect_equal(actual_inf_counter, expected_inf_counter)
})

test_that("simulate_single_year() works with individuals born after 1st Jan", {

  set.seed(1613)

  years           <- 2:5
  follow_up_times <- 2:5
  weights_fu      <- rep(0.25, length(follow_up_times))

  n_filled_buckets_nat <- 0
  n_filled_buckets_vac <- 0
  total_buckets        <- 0
  lambda               <- rep(0.1, 26)
  cml_lambda           <- 0.1 * (1 - 0.25)
  current_idx          <- 1
  inf_counter          <- 0
  rho                  <- 0
  rho_v                <- 0
  is_vaccinated        <- 0
  vac_year             <- 15

  actual_lo_lambda <- vector(mode = "numeric", length = length(years))

  for(yr in years)
  {
    output <- simulate_single_year(yr, n_filled_buckets_nat,
                                   n_filled_buckets_vac, rho,
                                   rho_v, current_idx, cml_lambda,
                                   follow_up_times, weights_fu, lambda,
                                   total_buckets, is_vaccinated, vac_year,
                                   vac_buckets = 0, is_vac_prob = 0,
                                   inf_counter)

    n_filled_buckets_nat <- output[[1]]
    flag_report          <- output[[3]]
    if(flag_report == 1) current_idx <- current_idx + 1
    inf_counter          <- output[[4]]
    lo_lambda            <- output[[5]] # leftover lambda

    actual_lo_lambda[yr - min(years) + 1] <- lo_lambda
  }

  expected_lo_lambda  <- rep(0.1 * 0.75, length(years))
  expect_equal(actual_lo_lambda, expected_lo_lambda)
})

test_that("simulate_single_year() works with vaccination", {

  set.seed(1613)

  current_idx <- 1

  is_vaccinated   <- 1
  vac_year        <- 18
  vac_buckets     <- 1

  years           <- 16:23
  follow_up_times <- 16:23
  weights_fu      <- rep(0, 8)

  n_filled_buckets_nat <- 0
  n_filled_buckets_vac <- 0

  rho   <- 0.05
  rho_v <- 0.12

  cml_lambda <- 0.1
  lambda     <- rep(0.1, 26)

  total_buckets <- 4
  inf_counter   <- 0
  is_vac_prob   <- 1


  actual             <- vector(mode = "numeric", length = length(years))
  actual_vac_buckets <- vector(mode = "numeric", length = length(years))
  actual_inf_counter <- vector(mode = "numeric", length = length(years))

  for(yr in years)
  {
    output <- simulate_single_year(yr, n_filled_buckets_nat,
                                   n_filled_buckets_vac, rho,
                                   rho_v, current_idx, cml_lambda,
                                   follow_up_times, weights_fu, lambda,
                                   total_buckets, is_vaccinated, vac_year,
                                   vac_buckets, is_vac_prob, inf_counter)

    n_filled_buckets_nat <- output[[1]]
    n_filled_buckets_vac <- output[[2]]
    flag_report          <- output[[3]]
    if(flag_report == 1) current_idx <- current_idx + 1
    inf_counter          <- output[[4]]

    actual[yr - min(years) + 1]             <- n_filled_buckets_nat
    actual_vac_buckets[yr - min(years) + 1] <- n_filled_buckets_vac
    actual_inf_counter[yr - min(years) + 1] <- inf_counter
  }

  expected_nat_buckets <- c(0.95, 0.90, 0.95, 0.9, 1.85, 1.80, 1.75, 1.70)

  expect_equal(actual, expected_nat_buckets)

  expected_vac_buckets <- c(0, 0, 0.88, 0.76, 0.64, 0.52, 0.40, 0.28)

  expect_equal(actual_vac_buckets, expected_vac_buckets)
})
