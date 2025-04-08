test_that("sim_infection_seroneg() works", {

  follow_up_times <- c(16, 17, 18, 19, 20, 21, 22, 24, 25)
  enrolment_year  <- 15

  weights_fu       <- c(0.6684932, 0.6109589, 0.6273973, 0.6147541, 0.6821918,
                        0.6328767, 0.9013699, 0.4027397, 0.8712329)

  enrolment_weight <- 0.295082

  set.seed(123)

  actual <- sim_infection_seroneg(follow_up_times  = follow_up_times,
                                  enrolment_year     = enrolment_year,
                                  lambda             = rep(0.1, 26),
                                  total_buckets      = 4,
                                  weights_fu         = weights_fu,
                                  enrolment_weight   = enrolment_weight,
                                  rho                = 0.05,
                                  rho_v              = 0,
                                  is_vaccinated      = 0,
                                  vac_year           = 15)

  expected <- c(0, 1, 0, 1, 1, 0, 0, 0, 0)

  expect_equal(actual, expected)
})

# test_that("sim_infection_seroneg() works", {
#
#   follow_up_times  <- c(2:4)
#   enrolment_year   <- 1 # enrolled as soon as becoming susceptible
#   lambda           <- rep(0.1, 4)
#   weights_fu       <- rep(1, 3)
#   enrolment_weight <- 0
#   rho              <- 0
#   rho_v            <- 0
#   is_vaccinated    <- 1
#   vac_year         <- 4
#   total_buckets    <- 4
#   lambda           <- rep(0.1, 26)
#   vac_buckets      <- 1
#
#   set.seed(1048)
#
#   actual <- sim_infection_seroneg(follow_up_times    = follow_up_times,
#                                   enrolment_year     = enrolment_year,
#                                   lambda             = lambda,
#                                   total_buckets      = total_buckets,
#                                   weights_fu         = weights_fu,
#                                   enrolment_weight   = enrolment_weight,
#                                   rho                = rho,
#                                   rho_v              = rho_v,
#                                   is_vaccinated      = is_vaccinated,
#                                   vac_year           = vac_year,
#                                   vac_buckets        = vac_buckets)
#
#
#
# })
