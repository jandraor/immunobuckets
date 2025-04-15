test_that("draw_vac_buckets() returns the expected values ", {

  actual   <- draw_vac_buckets(inf_counter = 4, vac_buckets = 1,
                               total_buckets = 4, is_vac_prob = 1)
  expected <- 0

  expect_equal(actual,expected)

  set.seed(1629)

  actual   <- draw_vac_buckets(inf_counter = 3, vac_buckets = 1,
                               total_buckets = 4, is_vac_prob = 1)

  expected <- 0

  expect_equal(actual,expected)
})

test_that("draw_vac_buckets() works for deterministic vaccination",
{
  actual   <- draw_vac_buckets(inf_counter = 4, vac_buckets = 1,
                               total_buckets = 4, is_vac_prob = 0)
  expected <- 0

  expect_equal(actual,expected)

  actual   <- draw_vac_buckets(inf_counter = 3, vac_buckets = 1,
                               total_buckets = 4, is_vac_prob = 0)
  expected <- 1

  expect_equal(actual,expected)

  actual   <- draw_vac_buckets(inf_counter = 0, vac_buckets = 1,
                               total_buckets = 4, is_vac_prob = 0)
  expected <- 1

  expect_equal(actual,expected)

  actual   <- draw_vac_buckets(inf_counter = 1, vac_buckets = 1,
                               total_buckets = 4, is_vac_prob = 0)
  expected <- 1

  expect_equal(actual,expected)

  actual   <- draw_vac_buckets(inf_counter = 0, vac_buckets = 2,
                               total_buckets = 4, is_vac_prob = 0)
  expected <- 2

  expect_equal(actual,expected)

  actual   <- draw_vac_buckets(inf_counter = 3, vac_buckets = 2,
                               total_buckets = 4, is_vac_prob = 0)
  expected <- 1

  expect_equal(actual,expected)

  actual   <- draw_vac_buckets(inf_counter = 4, vac_buckets = 2,
                               total_buckets = 4, is_vac_prob = 0)
  expected <- 0

  expect_equal(actual,expected)
})

test_that("draw_vac_buckets() works when the number of infections is higher than
          total buckets ",
{
  actual   <- draw_vac_buckets(inf_counter = 5, vac_buckets = 1,
                               total_buckets = 4, is_vac_prob = 1)
  expected <- 0

  expect_equal(actual,expected)
})

test_that("draw_vac_buckets() works when the number of vaccination buckets is zero",
{
  actual   <- draw_vac_buckets(inf_counter = 0, vac_buckets = 0,
                               total_buckets = 4)

  expected <- 0

  expect_equal(actual, expected)

})
