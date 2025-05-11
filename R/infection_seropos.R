sim_infection_seropos <- function(follow_up_times,
                                  enrolment_year,
                                  lambda,
                                  total_buckets,
                                  weights_fu,
                                  enrolment_weight,
                                  prob_n_inf_e,
                                  init_filled_buckets,
                                  rho,
                                  rho_v,
                                  is_vaccinated,
                                  vac_buckets,
                                  is_vac_prob) {

  n_follow_up <- length(follow_up_times)
  y           <- vector(mode = "integer", length = n_follow_up)
  start_year  <- follow_up_times[1]
  end_year    <- follow_up_times[n_follow_up]

  # Natural infection(nat)
  n_filled_buckets_nat <- init_filled_buckets

  lambda_period <- (1 - enrolment_weight) * lambda[enrolment_year]

  current_idx <- 1 # Refers to the index of the current measurement

  n_filled_buckets_vac <- rep(0, vac_buckets + 1)

  if(is_vaccinated == 1)
  {
    n_filled_buckets_vac <- 0:vac_buckets

    for(idx in 1:vac_buckets)
    {
      n_filled_buckets_vac[idx + 1] <-
        max(0, n_filled_buckets_vac[idx + 1] - rho_v * (1 - enrolment_weight))
    }

  }

  for(n_inf_e in 1:total_buckets)
  {
    n_filled_buckets_nat[n_inf_e] <-
      max(0, n_filled_buckets_nat[n_inf_e] - rho * (1 - enrolment_weight))
  }


  for(yr in start_year:end_year) # loop through all possible years
  {
    weight <- determine_weight(yr, follow_up_times, weights_fu, current_idx)

    lambda_period <- lambda_period + weight * lambda[yr]

    if(vac_buckets > 0)
    {
      for(idx in 1:vac_buckets)
      {
        n_filled_buckets_vac[idx + 1] <-
          max(0, n_filled_buckets_vac[idx + 1] - rho_v * weight)
      }
    }

    for(n_inf_e in 1:total_buckets)
    {
      n_filled_buckets_nat[n_inf_e] <-
        max(0, n_filled_buckets_nat[n_inf_e] - rho * weight)
    }

    if(yr == follow_up_times[current_idx])
    {
      prob <- calculate_infection_probability(total_buckets,
                                              n_filled_buckets_nat,
                                              lambda_period,
                                              n_filled_buckets_vac,
                                              prob_n_inf_e,
                                              is_vac_prob,
                                              vac_buckets)

      y[current_idx] <- stats::rbinom(1, 1, prob)

      if (y[current_idx] == 1) # If infection
      {
        for(n_inf_e in 1:total_buckets)
        {
          n_filled_buckets_nat[n_inf_e] <- min(total_buckets,
                                               n_filled_buckets_nat[n_inf_e] + 1);
        }
      }
      current_idx   <- current_idx + 1
    }

    lambda_period <- (1 - weight) * lambda[yr] # Leftover

    if(vac_buckets > 0)
    {
      for(idx in 1:vac_buckets) {
        n_filled_buckets_vac[idx + 1] <-
          max(0, n_filled_buckets_vac[idx + 1] - rho_v * (1 - weight));
      }
    }

    for(n_inf_e in 1:total_buckets)
    {
      n_filled_buckets_nat[n_inf_e] <-
        max(0, n_filled_buckets_nat[n_inf_e] - rho * (1- weight));
    }
  }

  y
}

calculate_infection_probability <- function(total_buckets,
                                            n_filled_buckets_nat,
                                            lambda_period,
                                            n_filled_buckets_vac,
                                            prob_n_inf_e,
                                            is_vac_prob,
                                            vac_buckets) {
  prob <- 0

  for (n_inf_e in 1:total_buckets)
  {
    prob_sce <- calculate_prob_per_scenario(vac_buckets,
                                            n_inf_e,
                                            total_buckets,
                                            is_vac_prob)

    prob_inf_given_n_inf_e <-
      infection_probability_given_n_inf_e(
        n_filled_buckets_nat[n_inf_e],
        n_filled_buckets_vac,
        lambda_period,
        total_buckets,
        prob_sce)

    prob <- prob + prob_inf_given_n_inf_e * prob_n_inf_e[n_inf_e];
  }
  prob
}

calculate_prob_per_scenario <- function(vac_buckets, n_inf_e, total_buckets,
                                        is_vac_prob) {

  if(vac_buckets == 0) return(1)

  if(!is_vac_prob)
  {
    n_sce    <- vac_buckets + 1
    prob_sce <- rep(0, n_sce)
    index    <- total_buckets - n_inf_e + 1

    if (index > n_sce) index <- n_sce # Map to the last position

    prob_sce[index] <- 1

    return(prob_sce)
  }

  if(vac_buckets == 1) {

    prob_not_filling <- n_inf_e / total_buckets
    prob_sce         <- c(prob_not_filling, 1 - prob_not_filling)
    return(prob_sce)
  }

  if(vac_buckets == 2) {

    if(n_inf_e == 1) prob_sce <- c(0, 2/4, 2/4)

    if(n_inf_e == 2) prob_sce <- c(1/6, 4/6, 1/6)

    if(n_inf_e == 3) prob_sce <- c(3/4, 1/4, 0)

    if(n_inf_e == 4) prob_sce <- c(1, 0, 0)

    return(prob_sce)
  }

  stop("Scenario not supported", call. = FALSE)
}

infection_probability_given_n_inf_e <- function(n_filled_buckets_nat,
                                                n_filled_buckets_vac,
                                                lambda_period,
                                                total_buckets,
                                                prob_sce) {

  # unconditional probability (unc_prob)
  unc_prob <- 0

  n_scenarios <- length(prob_sce)

  for(sce in 1:n_scenarios)
  {
    total_filled_buckets <- n_filled_buckets_nat + n_filled_buckets_vac[sce]

    n_empty_buckets <- max(0,
                           total_buckets - total_filled_buckets)

    local_prob <- 1 - exp(-n_empty_buckets * lambda_period)

    unc_prob <- unc_prob + local_prob * prob_sce[sce]


  }
  unc_prob
}
