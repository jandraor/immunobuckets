test_that("estimate_prob_after_vaccination() returns the expected values ", {

  prob_per_n_buckets <- c(0.35, 0.3, 0.2, 0.15)

  actual   <- estimate_prob_after_vaccination(4, 1, prob_per_n_buckets)

  expected <- c(prob_per_n_buckets[1] * 0.25,
                prob_per_n_buckets[1] * 0.75 + prob_per_n_buckets[2]  * 0.5,
                prob_per_n_buckets[2] * 0.50 + prob_per_n_buckets[3]  * 0.75,
                prob_per_n_buckets[3] * 0.25 + 0.15)

  expect_equal(actual, expected)

  actual   <- estimate_prob_after_vaccination(4, 2, prob_per_n_buckets)

  expected <- c(0,
                (3 / 6) * prob_per_n_buckets[1] + (1/ 6) * prob_per_n_buckets[2] ,
                (3 / 6) * prob_per_n_buckets[1] +
                  (4 / 6) * prob_per_n_buckets[2] +
                  (3 / 6) * prob_per_n_buckets[3],
                (1 / 6) * prob_per_n_buckets[2] +
                  (3 / 6) * prob_per_n_buckets[3] +
                  prob_per_n_buckets[4])

  expect_equal(actual, expected)
})
