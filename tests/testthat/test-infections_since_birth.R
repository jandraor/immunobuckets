test_that("simulate_infections_since_birth() works",
{
  individuals_df <- data.frame(subject_id       = 1,
                               birth_year_index = 9,
                               is_vaccinated    = 0,
                               vac_year         = 15)

  rho        <- 0.05
  stop_index <- 15
  lambda     <- rep(0.14, 15)

  set.seed(1437)

  actual <- simulate_infections_since_birth(individuals_df,
                                            lambda,
                                            rho,
                                            stop_index,
                                            rho_v = 0.12)

  expect_equal(nrow(actual), 5)

  expected <- data.frame(subject_id = 1,
                         year_index = 11:15,
                         sim_y      = c(0, 1, 0, 0, 0))

  expect_equal(actual, expected)

  individuals_df <- data.frame(subject_id       = 1:2,
                               birth_year_index = c(9, 7),
                               is_vaccinated    = 0,
                               vac_year         = 15)

  set.seed(1437)

  actual <- simulate_infections_since_birth(individuals_df,
                                            lambda,
                                            rho,
                                            stop_index,
                                            rho_v = 0.12)

  expected <- data.frame(subject_id = c(rep(1, 5), rep(2, 7)),
                         year_index = c(11:15, 9:15),
                         sim_y      = c(0, 1, 0, 0, 0,
                                        1, 0, 1, 0, 0, 1, 1))

  expect_equal(actual, expected)
})

test_that("simulate_single_individual produces a different value for different rho", {

  individual_df <- data.frame(subject_id       = 1,
                              birth_year_index = 0,
                              is_vaccinated    = 0,
                              vac_year         = 1)

  stop_index  <- 10
  lambda      <- rep(0.1, 10)
  vac_buckets <- 4
  rho         <- 0
  rho_v       <- 0

  set.seed(123)

  sim_rho_0 <- simulate_single_individual(individual_df, stop_index, lambda, rho,
                                       vac_buckets, rho_v)
  set.seed(123)
  rho         <- 0.5
  sim_rho_0.5 <- simulate_single_individual(individual_df, stop_index, lambda, rho,
                                       vac_buckets, rho_v)

  expect_false(isTRUE(all.equal(sim_rho_0, sim_rho_0.5)))
})
