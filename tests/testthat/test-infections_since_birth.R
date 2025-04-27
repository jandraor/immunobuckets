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
