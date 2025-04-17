test_that("random_walk_lambda works",
{
  set.seed(321)
  actual <- random_walk_lambda(n_t   = 26,
                               dt    = 1/128,
                               sigma = 1,
                               x0    = 0.1)

  expect_type(actual, "double")

  expect_equal(length(actual), 26)
})
