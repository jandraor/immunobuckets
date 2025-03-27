test_that("sim_infection_seroneg() works", {

  #"NMC-321-00005"

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
                                  is_vaccinated      = 0)

  expected <- c(0, 1, 0, 1, 1, 0, 0, 1, 0)

  expect_equal(actual, expected)
})

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
                                  is_vac_prob         = 0)

  expected <- c(0, 1, 0, 0, 0, 0, 0, 0, 1)

  expect_equal(actual, expected)
})

test_that("calculate_infection_probability() works", {

  total_buckets        <- 4
  n_filled_buckets_nat <- c(0.680126, 1.680126, 2.680126, 3.680126)
  lambda_period        <- 0.09635752
  n_filled_buckets_vac <- 0
  prob_n_inf_e         <- c(0.21863194, 0.38702796,
                            0.30450087, 0.08983923)

  actual   <- calculate_infection_probability(
    total_buckets        = total_buckets,
    n_filled_buckets_nat = n_filled_buckets_nat,
    lambda_period        = lambda_period,
    n_filled_buckets_vac = n_filled_buckets_vac,
    prob_n_inf_e         = prob_n_inf_e,
    is_vac_prob          = 0)

  expected <- 0.1764751

  expect_equal(actual, expected, tolerance = 1e-6)
})

test_that("calculate_infection_probability() returns a lower value when
          vaccination is probabilistic", {

  total_buckets        <- 4
  is_vac_prob          <- FALSE
  n_filled_buckets_nat <- c(1:4) # Assuming no-resusceptibility
  lambda_period        <- 0.1
  n_filled_buckets_vac <- 1
  prob_n_inf_e         <- c(0.21863194, 0.38702796,
                            0.30450087, 0.08983923)

  prob1 <- calculate_infection_probability(total_buckets,
                                           n_filled_buckets_nat,
                                           lambda_period,
                                           n_filled_buckets_vac,
                                           prob_n_inf_e,
                                           is_vac_prob)

  is_vac_prob <- TRUE

  prob <- calculate_infection_probability(total_buckets,
                                          n_filled_buckets_nat,
                                          lambda_period,
                                          n_filled_buckets_vac,
                                          prob_n_inf_e,
                                          is_vac_prob)

  expect_true(prob1 < prob)
})

test_that("simulate_infections() works", {

  total_buckets    <- 4
  n_e              <- 8
  lambda           <- rep(0.1, 26)

  start_index  <- 8
  init_weight  <- 0.3479452
  final_weight <- 0.7068493

  serostatus       <- 1
  is_vaccinated    <- 0
  rho_v            <- 0.10
  rho              <- 0.05
  vac_buckets      <- 1
  follow_up_times  <- c(16, 17, 18, 19, 20, 21, 22, 24, 25)
  enrolment_year   <- 15
  weights_fu       <- c(0.6684932, 0.6109589, 0.6273973, 0.6147541, 0.6821918,
                        0.6328767, 0.9013699, 0.4027397, 0.8712329)
  enrolment_weight <- 0.295082
  is_vac_prob      <- 1
  switch_rho       <- 1

  set.seed(1142)
  actual <- simulate_infections(total_buckets, n_e, lambda, start_index,
                                init_weight, final_weight, serostatus,
                                is_vaccinated, rho_v, rho, vac_buckets,
                                follow_up_times, enrolment_year, weights_fu,
                                enrolment_weight, is_vac_prob, switch_rho)

  expected <- c(0, 0, 0, 0, 0, 0, 1, 0, 0)

  expect_equal(actual, expected)

  serostatus    <- 0
  is_vaccinated <- 1

  set.seed(1146)
  actual <- simulate_infections(total_buckets, n_e, lambda, start_index,
                                init_weight, final_weight, serostatus,
                                is_vaccinated, rho_v, rho, vac_buckets,
                                follow_up_times, enrolment_year, weights_fu,
                                enrolment_weight, is_vac_prob, switch_rho)

  expected <- c(1, 0, 0, 0, 0, 1, 0, 0, 1)

  expect_equal(actual, expected)
})
