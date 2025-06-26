test_that("add_titre_to_infection_history() works",
{
  inf_history_df <-
    data.frame(subject_id = 1,
               year_index = 11:23,
               age        = 2:14,
               sim_y      = 0)

  inf_history_df$sim_y[c(6, 8, 11, 12)] <- 1

  actual <- add_titre_to_infection_history(
    inf_history_df,
    long_options  = list(L = 0.22, a = -0.96, b = 18.6),
    short_options = list(intercept = 0.02,
                         slope     = 0.003342),
    rise_options  = list(intercept = 8,
                         slope     = -0.8))

  expected <- inf_history_df

  first_inf <- bi_exponential_decay(
    c(0, 360, 720), exp(8),
    par_sigma = 0.02 + 0.003342*8,
    par_gamma = (0.22 / (1 + exp(0.96 * (c(7, 7:8) - 18.6)))) / 360,
    par_rho = 0.98) |> log()

  rise <- 8 + first_inf[3] * -0.8

  second_inf <- bi_exponential_decay(
    360 * 0:3,
    exp(first_inf[3] + rise),
    0.02 + 0.003342 * rise,
    par_gamma = (0.22 / (1 + exp(0.96 * (c(9, 9:11) - 18.6)))) / 360,
    par_rho = 0.98) |> log()

  rise <- 8 + second_inf[4] * -0.8

  third_inf <- bi_exponential_decay(
    360 * 0:1,
    A0 = exp(second_inf[4] + rise),
    0.02 + 0.003342 * rise,
    par_gamma = (0.22 / (1 + exp(0.96 * (c(12, 12) - 18.6)))) / 360,
    par_rho = 0.98) |>
    log()

  rise <- 8 + third_inf[2] * -0.8

  fourth_inf <- bi_exponential_decay(
    360 * 0:1,
    A0 = exp(third_inf[2] + rise),
    0.02 + 0.003342 * rise,
    par_gamma = (0.22 / (1 + exp(0.96 * (c(13, 13) - 18.6)))) / 360,
    par_rho = 0.98) |>
    log()

  expected$cumulative_infections <- cumsum(expected$sim_y)

  expected$log_titre <- c(rep(0, 5), first_inf[1:2],
                          second_inf[1:3],
                          third_inf[1],
                          fourth_inf)

  expect_equal(actual, expected)
})

test_that("estimate_long_rate() works",
{
  long_options <- list(L = 0.22, a = -0.96, b = 18.6)

  age <- 7:9

  actual <- estimate_long_rate(age = 7:9, long_options)

  expected <- 0.22 / (1 + exp(0.96 * (7:9 - 18.6)))

  expect_equal(actual, expected)
})

test_that("estimate_rise() works",
{
  rise_options <- list(intercept = 8,
                       slope     = 0.8)

  actual <- estimate_rise(baseline_titre = 0, rise_options)

  expected <- 8

  expect_equal(actual, expected)
})

test_that("estimate_short_rate() works", {

  short_options <- list(intercept = 0.02,
                        slope     = 0.003342)

  actual <- estimate_short_rate(rise = 0, short_options)


  expected <- 0.02

  expect_equal(actual, expected)
})
