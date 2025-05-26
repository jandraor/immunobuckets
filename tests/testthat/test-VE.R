test_that("estimate_cumulative_VE() works", {

  follow_up_df <- data.frame(subject_id    = rep(1:4, 5),
                             year_index    = rep(1:5, each = 4),
                             is_vaccinated = rep(c(1, 1, 0, 0), 5),
                             infection     = c(0, 0, 1, 0,
                                               0, 1, 0, 1,
                                               1, 0, 0, 0,
                                               0, 0, 0, 0,
                                               0, 0, 0, 0))

  actual <- estimate_cumulative_VE(follow_up_df)

  expected <- data.frame(year_index = 1:5,
                         VE         = 1- c((0/2) / (1/2) ,
                                           (1/4) / (2/3),
                                           (2/5) / (2/3),
                                           (2/5) / (2/3),
                                           (2/5) / (2/3)))
  expect_equal(actual, expected)
})
