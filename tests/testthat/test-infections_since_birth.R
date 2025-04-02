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
                                            stop_index)

  expected <- data.frame(subject_id = 1,
                         year_index = 10:15,
                         sim_y      = c(1, 0, 0, 0, 0, 1))

  expect_equal(actual, expected)

  individuals_df <- data.frame(subject_id       = 1:2,
                               birth_year_index = c(9, 7),
                               is_vaccinated    = 0,
                               vac_year         = 15)

  set.seed(1437)

  actual <- simulate_infections_since_birth(individuals_df,
                                            lambda,
                                            rho,
                                            stop_index)

  expected <- data.frame(subject_id = c(rep(1, 6), rep(2, 8)),
                         year_index = c(10:15, 8:15),
                         sim_y      = c(1, 0, 0, 0, 0, 1,
                                        1, 1, 0, 0, 1, 1, 0, 0))

  expect_equal(actual, expected)
})
