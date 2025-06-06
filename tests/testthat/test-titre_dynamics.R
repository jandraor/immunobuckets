test_that("add_titre_to_infection_history() works",
{
  inf_history_df <-
    data.frame(subject_id = 1,
               year_index = 11:23,
               sim_y      = 0)

  inf_history_df$sim_y[c(6, 8, 11, 12)] <- 1

  actual <- add_titre_to_infection_history(inf_history_df,
                                           long_decay_coefficient = 1,
                                           peak_scaling_factor    = 0)

  expected <- inf_history_df

  titre_after_one_year <- bi_exponential_decay(360, exp(10),
                                               1/30, 1/720, par_rho = 0.98) |>
    log()

  titre_after_two_years <- bi_exponential_decay(720, exp(10),
                                                1/30, 1/720, par_rho = 0.98) |>
    log()

  expected$cumulative_infections <- cumsum(expected$sim_y)

  expected$log_titre <- c(rep(0, 5), 10,
                          titre_after_one_year, 10,  titre_after_one_year,
                          titre_after_two_years , 10, 10,
                          titre_after_one_year)

  expect_equal(actual, expected)
})

test_that("estimate_log_peak() works", {
  actual   <- estimate_log_peak(log_peak_0 = 10, n_inf = 1, a = 2, p = 0.5)
  expected <- 10

  expect_equal(actual, expected)

  actual   <- estimate_log_peak(log_peak_0 = 10, n_inf = 2, a = 2, p = 0.5)
  expected <- 12

  expect_equal(actual, expected)
})
