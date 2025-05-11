test_that("sim_infection_seropos() works", {

  follow_up_times <- c(16, 17, 18, 19, 20, 21, 22, 24, 25)
  enrolment_year  <- 15

  weights_fu       <- c(0.6684932, 0.6109589, 0.6273973, 0.6147541, 0.6821918,
                        0.6328767, 0.9013699, 0.4027397, 0.8712329)

  enrolment_weight <- 0.295082

  set.seed(123)

  # Number of empty buckets for each scenario (# of infections) at enrolment
  init_empty_buckets <- vector(mode = "numeric", length = 4)

  lambda_e <- 0.779394
  rho      <- 0.05
  n_e      <- 8

  B_rho <- max(0, rho * n_e * (1 - 1 / (4 * lambda_e)))

  total_buckets <- 4

  # This for-loop estimates re-susceptibility from birth until enrolment
  for(n_inf_e in 1:total_buckets)
  {
    # B_rho represents additional empty buckets due to re-susceptibility
    init_empty_buckets[n_inf_e] <- min(total_buckets,
                                       total_buckets - n_inf_e + B_rho);
  }

  prob_n_inf_e <- c(0.21863194, 0.38702796, 0.30450087, 0.08983923)

  init_empty_buckets <- c(3.2716952, 2.2716952, 1.2716952, 0.2716952)

  set.seed(987)

  init_filled_buckets <- 4 - init_empty_buckets

  actual <- sim_infection_seropos(follow_up_times     = follow_up_times,
                                  enrolment_year      = enrolment_year,
                                  lambda              = rep(0.1, 26),
                                  total_buckets       = 4,
                                  weights_fu          = weights_fu,
                                  enrolment_weight    = enrolment_weight,
                                  prob_n_inf_e        = prob_n_inf_e,
                                  init_filled_buckets = init_filled_buckets,
                                  rho                 = 0.05,
                                  rho_v               = 0.07,
                                  is_vaccinated       = 0,
                                  is_vac_prob         = 0,
                                  vac_buckets         = 1)

  expected <- c(0, 1, 0, 0, 0, 0, 0, 0, 1)

  expect_equal(actual, expected)
})

#---------calculate_infection_probability()-------------------------------------

test_that("calculate_infection_probability() works", {

  #-----------------------------------------------------------------------------
  total_buckets        <- 4
  n_filled_buckets_nat <- 1:4
  lambda_period        <- 0.1
  n_filled_buckets_vac <- c(0, 0)
  prob_n_inf_e         <- c(0.2, 0.4, 0.1, 0.3)

  actual <- calculate_infection_probability(
    total_buckets        = total_buckets,
    n_filled_buckets_nat = n_filled_buckets_nat,
    lambda_period        = lambda_period,
    n_filled_buckets_vac = n_filled_buckets_vac,
    prob_n_inf_e         = prob_n_inf_e,
    is_vac_prob          = 0,
    vac_buckets          = 0)

  probability_infection <- function(empty_buckets, lambda) {
    1 - exp(-empty_buckets * lambda)
  }

  expected <- sum(probability_infection(total_buckets - n_filled_buckets_nat,
                                        lambda_period) * prob_n_inf_e)

  expect_equal(actual, expected)

  #-----------------------------------------------------------------------------

  total_buckets        <- 4
  n_filled_buckets_nat <- c(0.680126, 1.680126, 2.680126, 3.680126)
  lambda_period        <- 0.09635752
  n_filled_buckets_vac <- c(0, 1)
  prob_n_inf_e         <- c(0.21863194, 0.38702796,
                            0.30450087, 0.08983923)

  actual <- calculate_infection_probability(
    total_buckets        = total_buckets,
    n_filled_buckets_nat = n_filled_buckets_nat,
    lambda_period        = lambda_period,
    n_filled_buckets_vac = n_filled_buckets_vac,
    prob_n_inf_e         = prob_n_inf_e,
    is_vac_prob          = 0,
    vac_buckets          = 1)

  expected <- 0.21863194 * 0.2003149 +
    0.1194247 * 0.38702796 +
    0.30450087 * 0.0303521 +
    0.08983923* 0.0303521

  expect_equal(actual, expected, tolerance = 1e-6)
})

test_that("calculate_infection_probability() works with vaccination", {

  actual <- calculate_infection_probability(
    total_buckets        = 4,
    n_filled_buckets_nat = c(1:4),
    lambda_period        = 0.1,
    n_filled_buckets_vac = c(0, 1),
    prob_n_inf_e         = c(0, 0, 1, 0),
    is_vac_prob          = 1,
    vac_buckets          = 1)

  expected <- (1 - exp(-0.1)) * 0.75

  expect_equal(actual, expected)

  actual <- calculate_infection_probability(
    total_buckets        = 4,
    n_filled_buckets_nat = c(1:4),
    lambda_period        = 0.1,
    n_filled_buckets_vac = c(0, 1),
    prob_n_inf_e         = c(0.5, 0.5, 0, 0),
    is_vac_prob          = 1,
    vac_buckets          = 1)

  expected <- 0.5 * (0.25 * (1 - exp(-0.1 * 3)) + 0.75 * (1 - exp(-0.1 * 2)) ) +
    0.5 * (0.5 * (1 - exp(-0.1 * 2)) + 0.5 * (1 - exp(-0.1 * 1)) )

  expect_equal(actual, expected)

})

test_that("calculate_infection_probability() returns a lower value when
          vaccination is probabilistic", {

            total_buckets        <- 4
            is_vac_prob          <- FALSE
            n_filled_buckets_nat <- c(1:4) # Assuming no-resusceptibility
            lambda_period        <- 0.1
            n_filled_buckets_vac <- c(0, 1)
            prob_n_inf_e         <- c(0.21863194, 0.38702796,
                                      0.30450087, 0.08983923)

            prob1 <- calculate_infection_probability(total_buckets,
                                                     n_filled_buckets_nat,
                                                     lambda_period,
                                                     n_filled_buckets_vac,
                                                     prob_n_inf_e,
                                                     is_vac_prob,
                                                     vac_buckets = 1)

            is_vac_prob <- TRUE

            prob <- calculate_infection_probability(total_buckets,
                                                    n_filled_buckets_nat,
                                                    lambda_period,
                                                    n_filled_buckets_vac,
                                                    prob_n_inf_e,
                                                    is_vac_prob,
                                                    vac_buckets = 1)

            expect_true(prob1 < prob)
          })

test_that("calculate_infection_probability() works with 2 vac buckets", {

  total_buckets        <- 4
  is_vac_prob          <- TRUE
  n_filled_buckets_nat <- c(1:4)
  lambda_period        <- 0.1
  n_filled_buckets_vac <- c(0, 1, 2) # Potential
  prob_n_inf_e         <- c(0.21863194, 0.38702796,
                            0.30450087, 0.08983923)

  actual <- calculate_infection_probability(total_buckets,
                                            n_filled_buckets_nat,
                                            lambda_period,
                                            n_filled_buckets_vac,
                                            prob_n_inf_e,
                                            is_vac_prob,
                                            vac_buckets = 2)

  expected <- 0.08819766

  expect_equal(actual, expected, tolerance = 1e-6)
})


test_that("infection_probability_given_n_inf_e works", {

  lambda_period <- 0.09635752

  actual <- infection_probability_given_n_inf_e(
    n_filled_buckets_nat = 0.680126,
    n_filled_buckets_vac = c(0, 1),
    lambda_period        = lambda_period,
    total_buckets        = 4,
    prob_sce             = c(0, 1))

  expected <- 1 - exp(-lambda_period * (4 - 1.680126))

  expect_equal(actual, expected)

  actual <- infection_probability_given_n_inf_e(
    n_filled_buckets_nat = 1.680126,
    n_filled_buckets_vac = c(0, 1),
    lambda_period        = lambda_period,
    total_buckets        = 4,
    prob_sce             = c(0, 1))

  expected <- 1 - exp(-lambda_period *(4 - 2.680126))

  expect_equal(actual, expected)

  actual <- infection_probability_given_n_inf_e(
    n_filled_buckets_nat = 2.680126,
    n_filled_buckets_vac = c(0, 1),
    lambda_period        = lambda_period,
    total_buckets        = 4,
    prob_sce             = c(0, 1))

  expected <- 1 - exp(-lambda_period *(4 - 3.680126))

  expect_equal(actual, expected)

  actual <- infection_probability_given_n_inf_e(
    n_filled_buckets_nat = 3.680126,
    n_filled_buckets_vac = c(0, 1),
    lambda_period        = lambda_period ,
    total_buckets        = 4,
    prob_sce             = c(0, 1))

  expected <- 0

  expect_equal(actual, expected)

  #----------------------------------------------------------------------------
  lambda_period <- 0.1
  total_buckets <- 4

  n_filled_buckets_nat <- 1
  n_filled_buckets_vac <- c(0, 1)

  prob_sce <- c(0.25, 0.75) # calculate_prob_per_scenario(1, 1, 4)

  actual <- infection_probability_given_n_inf_e(n_filled_buckets_nat,
                                                n_filled_buckets_vac,
                                                lambda_period,
                                                total_buckets,
                                                prob_sce)

  expected <- 0.20074738

  expect_equal(actual, expected)

  #-----------------------------------------------------------------------------
  # n_inf_e = 2
  n_filled_buckets_nat <- 2

  prob_sce <- c(0.5, 0.5) # calculate_prob_per_scenario(1, 2, 4)

  actual <- infection_probability_given_n_inf_e(n_filled_buckets_nat,
                                                n_filled_buckets_vac,
                                                lambda_period,
                                                total_buckets,
                                                prob_sce)

  expected <- 0.138215914

  expect_equal(actual, expected)

  #-----------------------------------------------------------------------------
  n_filled_buckets_nat <- 3

  prob_sce <- c(0.75, 0.25) # calculate_prob_per_scenario(1, 3, 4)

  actual <- infection_probability_given_n_inf_e(n_filled_buckets_nat,
                                                n_filled_buckets_vac,
                                                lambda_period,
                                                total_buckets,
                                                prob_sce)

  expected <- 0.071371936

  expect_equal(actual, expected)

  #-----------------------------------------------------------------------------

  n_filled_buckets_nat <- 4
  # n_inf_e = 4

  prob_sce <- c(0, 1) # calculate_prob_per_scenario(1, 4, 4)

  actual <- infection_probability_given_n_inf_e(n_filled_buckets_nat,
                                                n_filled_buckets_vac,
                                                lambda_period,
                                                total_buckets,
                                                prob_sce)

  expected <- 0

  expect_equal(actual, expected)

  #-----------------------------------------------------------------------------
  # n_inf_e = 1

  lambda_period <- 0.1
  total_buckets <- 4

  n_filled_buckets_nat <- 1
  n_filled_buckets_vac <- c(0, 1, 2)


  prob_sce <- c(0, 0.5, 0.5) # calculate_prob_per_scenario(2, 1, 4)


  actual <- infection_probability_given_n_inf_e(n_filled_buckets_nat,
                                                n_filled_buckets_vac,
                                                lambda_period,
                                                total_buckets,
                                                prob_sce)

  expected <- 0.138215914

  expect_equal(actual, expected)

  #-----------------------------------------------------------------------------

  lambda_period <- 0.1
  total_buckets <- 4

  n_filled_buckets_nat <- 2
  n_filled_buckets_vac <- c(0, 1, 2)

  prob_sce <- c(1/6, 4/6, 1/6) # calculate_prob_per_scenario(2, 2, 4)

  actual <- infection_probability_given_n_inf_e(n_filled_buckets_nat,
                                                n_filled_buckets_vac,
                                                lambda_period,
                                                total_buckets,
                                                prob_sce)

  expected <- 0.093653262

  expect_equal(actual, expected)
})

test_that("calculate_prob_per_scenario() works", {

  actual   <- calculate_prob_per_scenario(vac_buckets   = 1,
                                          n_inf_e       = 1,
                                          total_buckets = 4,
                                          is_vac_prob   = TRUE)

  expected <- c(0.25, 0.75)

  expect_equal(actual, expected)

  actual <- calculate_prob_per_scenario(vac_buckets   = 4,
                                        n_inf_e       = 3,
                                        total_buckets = 4,
                                        is_vac_prob   = FALSE)

  expected <- c(0, 1, 0, 0, 0)

  expect_equal(actual, expected)
})

test_that("calculate_prob_per_scenario() works for a vaccine that fills three
          buckets deterministically",
{
  actual <- calculate_prob_per_scenario(vac_buckets   = 3,
                                        n_inf_e       = 1,
                                        total_buckets = 4,
                                        is_vac_prob   = FALSE)

  expected <- c(0, 0, 0, 1)

  expect_equal(actual, expected)

  actual <- calculate_prob_per_scenario(vac_buckets   = 3,
                                        n_inf_e       = 4,
                                        total_buckets = 4,
                                        is_vac_prob   = FALSE)

  expected <- c(1, 0, 0, 0)

  expect_equal(actual, expected)
})
