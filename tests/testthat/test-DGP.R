test_that("simulate_infections() works", {

  m_specs <- list(total_buckets = 4,
                  vac_buckets   = 0,
                  switch_rho    = 1,
                  is_vac_prob   = 0)

  rho   <- 0.05
  rho_v <- 0.17 # should have no effect

  lambda <- rep(0.14, 26)

  ind_par_list <- list(
    n_e              = 8,
    start_index      = 8,
    init_weight      = 0.6821918,
    final_weight     = 0.2684932,
    serostatus       = 0,
    is_vaccinated    = 1,
    enrolment_year   = 15,
    enrolment_weight = 0.2704918,
    follow_up_times  = c(15, 16, 17, 18, 19, 20, 21, 22, 23, 25),
    weights_fu       = c(0.99, 0.68, 0.63, 0.64, 0.64, 0.69, 0.6410959, 0.91,
                         0.94, 0.74))


  actual <- simulate_infections(
    n_e              = ind_par_list$n_e,
    start_index      = ind_par_list$start_index,
    init_weight      = ind_par_list$init_weight,
    final_weight     = ind_par_list$final_weight,
    serostatus       = ind_par_list$serostatus,
    is_vaccinated    = ind_par_list$is_vaccinated,
    follow_up_times  = ind_par_list$follow_up_times,
    enrolment_year   = ind_par_list$enrolment_year,
    weights_fu       = ind_par_list$weights_fu,
    enrolment_weight = ind_par_list$enrolment_weight,
    vac_year         = ind_par_list$enrolment_year,
    # Global parameters
    avg_lambda       = 0.14,
    lambda           = lambda,
    rho_v            = rho_v,
    rho              = rho,
    # Model specs
    total_buckets    = m_specs$total_buckets,
    vac_buckets      = m_specs$vac_buckets,
    is_vac_prob      = m_specs$is_vac_prob,
    switch_rho       = m_specs$switch_rho)


  expect_type(actual , "double")
})

