test_that("add_titre_to_infection_history() works",
{
  inf_history_df <-
    data.frame(subject_id = 1,
               year_index = 11:23,
               sim_y      = 0)

  inf_history_df$sim_y[c(6, 8, 11, 12)] <- 1

  actual <- add_titre_to_infection_history(inf_history_df,
                                           long_decay_coefficient = 1,
                                           initial_rise  = 10,
                                           peak_dynamics = "power_law",
                                           peak_options  = list(p = 0.5, a = 0),
                                           ltr_0                  = 1/720)

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

test_that("add_titre_to_infection_history() works with a piecewise rate",
{
  inf_history_df <-
    data.frame(subject_id = 1,
               year_index = 11:23,
               sim_y      = 0)

  inf_history_df$sim_y[c(3, 5, 7, 8, 10)] <- 1

  actual <- add_titre_to_infection_history(
    inf_history_df,
    ltr_0                   = 1/720,
    long_term_rate_dynamics = "piecewise",
    initial_rise            = 10,
    peak_dynamics           = "power_law",
    peak_options            = list(p = 0.5, a = 0))

  expected <- inf_history_df

  first_inf <- bi_exponential_decay(c(0, 360), exp(10),
                                    1/30, 1/720, par_rho = 0.98) |> log()

  second_inf <- first_inf

  third_inf <- 10

  fourth_infection <- first_inf

  fifth_infection <- bi_exponential_decay(
    t          = 360 * 0:3, exp(10),
    par_sigma  = 1/30,
    par_lambda =  0, par_rho = 0.98) |> log()

  expected$cumulative_infections <- cumsum(expected$sim_y)

  expected$log_titre <- c(rep(0, 2), first_inf, second_inf, third_inf,
                          fourth_infection, fifth_infection)

  expect_equal(actual, expected)
})

test_that("add_titre_to_infection_history() works with a linear rises",
{
  inf_history_df <-
    data.frame(subject_id = 1,
               year_index = 11:23,
               sim_y      = 0)

  inf_history_df$sim_y[c(3, 8, 12)] <- 1

  actual <- add_titre_to_infection_history(
    inf_history_df = inf_history_df,
    long_decay_coefficient = 1,
    initial_rise  = 8,
    peak_dynamics = "linear_rise",
    peak_options  = list(slope     = -0.8,
                         intercept = 9.95),
    ltr_0                  = 1/720)

  expected <- inf_history_df

  first_inf <- bi_exponential_decay(360 * 0:5, exp(8),
                                    1/30, 1/720, par_rho = 0.98) |> log()

  last_titre <- first_inf[length(first_inf)]

  log_new_peak <- last_titre + 9.95 - 0.8 * last_titre

  second_inf <- bi_exponential_decay(360 * 0:4, exp(log_new_peak),
                                    1/30, 1/720, par_rho = 0.98) |> log()

  last_titre <- second_inf [length(second_inf)]

  log_new_peak <- last_titre + 9.95 - 0.8 * last_titre

  third_inf  <- bi_exponential_decay(360 * 0:2, exp(log_new_peak),
                                     1/30, 1/720, par_rho = 0.98) |> log()

  expected$cumulative_infections <- cumsum(expected$sim_y)

  expected$log_titre <- c(rep(0, 2),
                          first_inf[-length(first_inf)],
                          second_inf[-length(second_inf)],
                          third_inf[-length(third_inf)])

  expect_equal(actual, expected)
})

test_that("estimate_log_peak() works with the power law option", {

  peak_options <- list(a = 2, p = 0.5)

  actual   <- estimate_log_peak(log_peak_0 = 10, n_inf = 1,
                                method = "power_law",
                                peak_options = peak_options)


  expected <- 10

  expect_equal(actual, expected)

  actual   <- estimate_log_peak(log_peak_0 = 10, n_inf = 2,
                                method = "power_law",
                                peak_options = peak_options)
  expected <- 12

  expect_equal(actual, expected)
})

test_that("estimate_log_peak() works with the linear rise options",
{
  peak_options <- list(
    baseline_titre = 3,
    slope          = -0.8,
    intercept      = 9.95)

  actual   <- estimate_log_peak(log_peak_0 = 8, n_inf = 1,
                                method = "linear_rise",
                                peak_options)
  expected <- 8

  expect_equal(actual, expected)

  peak_options <- list(
    baseline_titre = 3,
    slope          = -0.8,
    intercept      = 9.95)

  actual   <- estimate_log_peak(log_peak_0 = 8, n_inf = 2,
                                method = "linear_rise",
                                peak_options)
  expected <- 3  + 9.95 - 0.8 * 3

  expect_equal(actual, expected)
})

test_that("estimate_decay_rate works()",
{
  actual <- estimate_decay_rate(4, 1/720, "piecewise", 1)

  expected <- 1/720

  expect_equal(actual, expected)

  actual <- estimate_decay_rate(5, 0, "piecewise", 1)

  expected <- 1/720
})
