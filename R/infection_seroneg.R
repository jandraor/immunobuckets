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

  if(vac_year <= enrolment_year) {
    n_filled_buckets_vac <- ifelse(is_vaccinated == 1, vac_buckets, 0)
  }

  n_filled_buckets_nat <- 0

  inf_counter <- 0

  cml_lambda <- enrolment_weight * lambda[enrolment_year]

  for(yr in start_year:end_year)
  {
    # cat("\n Year: ", yr)
    # cat("\n n_filled_buckets_nat: ", n_filled_buckets_nat)
    # cat("\n n_filled_buckets_vac: ", n_filled_buckets_vac)

    weight <- determine_weight(yr, follow_up_times, weights_fu, current_idx)

    n_filled_buckets_nat <- max(0, n_filled_buckets_nat - rho * weight)
    n_filled_buckets_vac <- max(0, n_filled_buckets_vac - rho_v * weight)

    cml_lambda <- cml_lambda + weight * lambda[yr]

    if(yr == follow_up_times[current_idx]) {

      total_filled_buckets <- n_filled_buckets_vac + n_filled_buckets_nat

      n_empty_buckets <- max(0, total_buckets - total_filled_buckets)

      prob_inf <- 1 - exp(- n_empty_buckets * cml_lambda)

      y[current_idx] <- stats::rbinom(1, 1, prob_inf)

      # cat("\n y: ", y[current_idx], "\n")

      if (y[current_idx] == 1)
      {
        n_filled_buckets_nat <- min(total_buckets, n_filled_buckets_nat + 1)
        inf_counter          <- inf_counter + 1
      }

      current_idx <- current_idx + 1
    }

    # Vaccination occurs at the end of the period so that it does not interfere
    #  with the infection probability
    if(is_vaccinated == 1 && yr == vac_year)
    {
      n_filled_buckets_nat <- inf_counter # redresses loss of immunity
      n_filled_buckets_vac <- draw_vac_buckets(inf_counter, vac_buckets,
                                               total_buckets = total_buckets,
                                               is_vac_prob)
    }

    cml_lambda  <- (1 - weight) * lambda[yr]

    n_filled_buckets_vac <- max(0, n_filled_buckets_vac - rho_v * (1 - weight))
    n_filled_buckets_nat <- max(0, n_filled_buckets_nat - rho * (1 - weight))
  }

  y
}
