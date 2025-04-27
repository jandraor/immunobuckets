sim_infection_seroneg <- function(follow_up_times,
                                  enrolment_year,
                                  lambda,
                                  total_buckets,
                                  weights_fu,
                                  enrolment_weight,
                                  rho,
                                  rho_v,
                                  is_vaccinated,
                                  vac_buckets,
                                  vac_year,
                                  is_vac_prob) {

  n_follow_up <- length(follow_up_times)
  y           <- vector(mode = "integer", length = n_follow_up)

  start_year <- follow_up_times[1]
  end_year   <- follow_up_times[n_follow_up]

  current_idx <- 1

  n_filled_buckets_vac <- 0

  if(vac_year <= enrolment_year)
  {
    n_filled_buckets_vac <- ifelse(is_vaccinated == 1, vac_buckets, 0)
    n_filled_buckets_vac <- max(0,
                                n_filled_buckets_vac - rho_v * (1 - enrolment_weight))
  }

  n_filled_buckets_nat <- 0

  inf_counter <- 0

  cml_lambda <- (1 - enrolment_weight) * lambda[enrolment_year]

  for(yr in start_year:end_year)
  {
    output <- simulate_single_year(yr                   = yr,
                                   n_filled_buckets_nat = n_filled_buckets_nat,
                                   n_filled_buckets_vac = n_filled_buckets_vac,
                                   rho                  = rho,
                                   rho_v                = rho_v,
                                   current_idx          = current_idx,
                                   cml_lambda           = cml_lambda,
                                   follow_up_times      = follow_up_times,
                                   weights_fu           = weights_fu,
                                   lambda               = lambda,
                                   total_buckets        = total_buckets,
                                   is_vaccinated        = is_vaccinated,
                                   vac_year             = vac_year,
                                   vac_buckets          = vac_buckets,
                                   is_vac_prob          = is_vac_prob,
                                   inf_counter          = inf_counter)

    n_filled_buckets_nat <- output[[1]]
    n_filled_buckets_vac <- output[[2]]
    flag_report          <- output[[3]]
    inf_counter          <- output[[4]]
    cml_lambda           <- output[[5]]
    sim_y                <- output[[6]]

    if(flag_report == 1)
    {
      y[current_idx] <- sim_y
      current_idx    <- current_idx + 1
    }
  }

  y
}

# returns value at the end of the year
simulate_single_year <- function(yr, n_filled_buckets_nat, n_filled_buckets_vac,
                                 rho, rho_v, current_idx, cml_lambda,
                                 follow_up_times, weights_fu, lambda,
                                 total_buckets, is_vaccinated, vac_year,
                                 vac_buckets, is_vac_prob, inf_counter) {

  weight <- determine_weight(yr, follow_up_times, weights_fu, current_idx)

  n_filled_buckets_nat <- max(0, n_filled_buckets_nat - rho * weight)
  n_filled_buckets_vac <- max(0, n_filled_buckets_vac - rho_v * weight)

  cml_lambda <- cml_lambda + weight * lambda[yr]

  sim_y <- 0L

  flag_report <- 0

  if(yr == follow_up_times[current_idx])
  {
    flag_report <- 1

    total_filled_buckets <- n_filled_buckets_vac + n_filled_buckets_nat

    n_empty_buckets <- max(0, total_buckets - total_filled_buckets)

    prob_inf <- 1 - exp(- n_empty_buckets * cml_lambda)

    sim_y <- stats::rbinom(1, 1, prob_inf)

    if (sim_y == 1)
    {
      n_filled_buckets_nat <- min(total_buckets, n_filled_buckets_nat + 1)
      inf_counter          <- inf_counter + 1
    }
  }

  leftover_lambda  <- (1 - weight) * lambda[yr]

  if(is_vaccinated == 1 && vac_buckets > 0 && yr == vac_year)
  {
    n_filled_buckets_nat <- inf_counter # redresses loss of immunity
    n_filled_buckets_vac <- draw_vac_buckets(inf_counter, vac_buckets,
                                             total_buckets = total_buckets,
                                             is_vac_prob)
  }

  n_filled_buckets_vac <- max(0, n_filled_buckets_vac - rho_v * (1 - weight))
  n_filled_buckets_nat <- max(0, n_filled_buckets_nat - rho * (1 - weight))

  c("n_filled_buckets_nat" = n_filled_buckets_nat,
    "n_filled_buckets_vac" = n_filled_buckets_vac,
    "flag_report"          = flag_report,
    "inf_counter"          = inf_counter,
    "leftover_lambda"      = leftover_lambda,
    "sim_y"                = sim_y)
}
