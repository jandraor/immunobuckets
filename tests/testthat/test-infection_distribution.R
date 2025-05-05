test_that("calculate_prob_n_inf() works", {

  lambda_avg <- 0.1

  actual <- calculate_prob_n_inf(lambda_avg = lambda_avg,
                                 n_e = 5,
                                 total_buckets = 4)

  expected <- c(0.135335283,
                0.41056681,
                0.352641491,
                0.095361996,
                prod(1- exp(-(1:4) *lambda_avg)) * sum(exp(- (0:4) * lambda_avg)))

  expect_equal(actual, expected)
})

test_that("estimate_conditional_prob() works", {

  prob_vector <- c(0.1, 0.2, 0.3, 0.3, 0.1)
  min_inf     <- 1

  actual   <- estimate_conditional_prob(prob_vector, min_inf)
  expected <- prob_vector[-1] / sum(prob_vector[-1])

  expect_equal(actual, expected)
})
