determine_weight <- function(yr, enrolment_year, enrolment_weight,
                             follow_up_times, weights_fu, current_idx)
{
  weight <- 1; # Default value for years when there are no measurements.

  # In the first year, the weight represents the period from the
  #   collection of the blood sample to the end of that year.
  if(yr == enrolment_year) weight <- enrolment_weight;

  # In other years, if the current evaluated year (yr) matches the year of
  #   the current measurement (current_idx), the weight represents the period
  #   from the start of the year to the measurement's date.
  if(yr != enrolment_year && yr == follow_up_times[current_idx]) {
    weight = weights_fu[current_idx];
  }

  weight
}

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
                                  vac_year) {

  n_follow_up <- length(follow_up_times)
  y           <- vector(mode = "integer", length = n_follow_up)

  end_year <- follow_up_times[n_follow_up]

  cml_lambda <- 0 # Cumulative lambda

  current_idx <- 1

  n_filled_buckets_vac <- 0

  if(vac_year <= enrolment_year) {
    n_filled_buckets_vac <- ifelse(is_vaccinated == 1, vac_buckets, 0)
  }

  n_filled_buckets_nat <- 0

  inf_counter <- 0

  # I loop through the enrolment year to take into account the effect of the
  # FOI between the enrolment date and the 31st December of that year.
  for(yr in enrolment_year:end_year) # loop through all possible years
  {
    weight <- determine_weight(yr, enrolment_year, enrolment_weight,
                              follow_up_times, weights_fu, current_idx);

    cml_lambda <- cml_lambda + weight * lambda[yr];

    # Vaccination occurs at the start of the period
    if(is_vaccinated == 1 && yr > enrolment_year)
    {
      n_filled_buckets_nat <- inf_counter
      n_filled_buckets_vac <- draw_vac_buckets(inf_counter, vac_buckets,
                                               total_buckets = total_buckets)
    }

    if(yr == follow_up_times[current_idx]) {

      total_filled_buckets <- n_filled_buckets_vac + n_filled_buckets_nat;

      n_empty_buckets <- max(0,
                             total_buckets - total_filled_buckets);

      prob_inf <- max(0.00001, 1 - exp(- n_empty_buckets * cml_lambda))

      y[current_idx] <- stats::rbinom(1, 1, prob_inf)

      if (y[current_idx] == 1)
      {
        n_filled_buckets_nat <- n_filled_buckets_nat + 1
        inf_counter          <- inf_counter + 1
      }

      # Loss of immunity occurs at the end of the period
      n_filled_buckets_vac <- max(0, n_filled_buckets_vac - rho_v * weight);
      n_filled_buckets_nat <- max(0, n_filled_buckets_nat - rho * weight);

      cml_lambda  = (1 - weight) * lambda[yr]
      current_idx = current_idx + 1
    }
  }

  y
}


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
  end_year    <- follow_up_times[n_follow_up]

  # Natural infection(nat)
  n_filled_buckets_nat <- init_filled_buckets;

  lambda_period <- 0

  current_idx <- 1 # Refers to the index of the current measurement

  n_filled_buckets_vac <- rep(0, vac_buckets + 1)

  if(is_vaccinated == 1) n_filled_buckets_vac <- 0:vac_buckets


  for(yr in enrolment_year:end_year) # loop through all possible years
  {
    weight <- determine_weight(yr, enrolment_year, enrolment_weight,
                               follow_up_times, weights_fu, current_idx)

    lambda_period <- lambda_period + weight * lambda[yr]

    if(vac_buckets > 0) {

      for(idx in 1:vac_buckets) {
        n_filled_buckets_vac[idx + 1] <- max(0,
                                         n_filled_buckets_vac[idx + 1] - rho_v * weight);
      }
    }

    for(n_inf_e in 1:total_buckets)
    {
      n_filled_buckets_nat[n_inf_e] <-
        max(0,
            n_filled_buckets_nat[n_inf_e] - rho * weight);
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
          n_filled_buckets_nat[n_inf_e] <- max(0, n_filled_buckets_nat[n_inf_e] + 1);
        }
      }

      lambda_period <- (1 - weight) * lambda[yr] # Leftover
      current_idx   <- current_idx + 1
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

  if(vac_buckets == 1) {

    if(!is_vac_prob) {

      prob_sce <- c(0, 1)

      if(n_inf_e == 4) prob_sce <- c(1, 0)

      return(prob_sce)
    }

    prob_not_filling <- n_inf_e / total_buckets
    prob_sce         <- c(prob_not_filling, 1 - prob_not_filling)
    return(prob_sce)
  }

  if(vac_buckets == 2) {

    if(!is_vac_prob) {

      prob_sce <- c(0, 0, 1)

      if(n_inf_e == 3) prob_sce <- c(0, 1, 0)

      if(n_inf_e == 4) prob_sce <- c(1, 0, 0)

      return(prob_sce)
    }

    if(n_inf_e == 1) prob_sce <- c(0, 2/4, 2/4)

    if(n_inf_e == 2) prob_sce <- c(1/6, 4/6, 1/6)

    if(n_inf_e == 3) prob_sce <- c(3/4, 1/4, 0)

    if(n_inf_e == 4) prob_sce <- c(1, 0, 0)

    return(prob_sce)
  }

  if(vac_buckets == 3) {

    if(!is_vac_prob) {
      prob_sce <- rep(0, vac_buckets + 1)
      prob_sce[total_buckets - n_inf_e] <- 1

      return(prob_sce)
    }
  }

  if(vac_buckets == 4) {

    if(!is_vac_prob) {

      prob_sce <- rep(0, vac_buckets + 1)

      prob_sce[total_buckets - n_inf_e + 1] <- 1

      return(prob_sce)
    }
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

#' Simulate infections
#'
#' @param total_buckets A numeric value indicating the number of serotypes.
#' @param n_e A
#' @param lambda A numeric vector for the yearly FOI.
#' @param start_index A
#' @param init_weight A
#' @param final_weight A
#' @param serostatus A
#' @param is_vaccinated A
#' @param rho_v A
#' @param rho A
#' @param vac_buckets A numeric value indicating the number of buckets that
#'  a vaccine fills.
#' @param follow_up_times A
#' @param enrolment_year A
#' @param weights_fu A
#' @param enrolment_weight A numeric value between 0 and 1 representing the
#'  time elapsed from the date of the first blood draw to the year's end
#'  (31st December).
#' @param is_vac_prob A
#' @param switch_rho A numeric value,
#'  which can be either 1 (allows re-susceptibility) or 0 (no re-susceptibility).
#'  Employing the buckets metaphor, this parameter determines whether buckets
#'  have a hole of leaking immunity.
#'
#' @returns A vector of zeros and ones, where one represents infection.
#' @export
#'
#' @examples
#' simulate_infections(total_buckets    = 4,
#'                     n_e              = 2,
#'                     lambda           = rep(0.1, 4),
#'                     start_index      = 1,
#'                     init_weight      = 0,
#'                     final_weight     = 1,
#'                     serostatus       = 1,
#'                     is_vaccinated    = 1,
#'                     rho_v            = 0.1,
#'                     rho              = 0.5,
#'                     vac_buckets      = 1,
#'                     follow_up_times  = c(3, 4),
#'                     enrolment_year   = 2,
#'                     weights_fu       = c(1, 1),
#'                     enrolment_weight = 0,
#'                     is_vac_prob      = 1,
#'                     switch_rho       = 1,
#'                     vac_year         = 15)
simulate_infections <- function(total_buckets, n_e, lambda, start_index,
                                init_weight, final_weight, serostatus,
                                is_vaccinated, rho_v, rho, vac_buckets,
                                follow_up_times, enrolment_year, weights_fu,
                                enrolment_weight, is_vac_prob, switch_rho,
                                vac_year) {

  rho   <- rho * switch_rho
  rho_v <- rho_v * switch_rho

  lambda_e <- estimate_lambda_e(n_e, lambda, start_index, init_weight,
                                final_weight)

  p_e <- (1 - exp(-total_buckets * lambda_e))

  if(serostatus == 0) # if seronegative at enrolment
  {
    y <- sim_infection_seroneg(follow_up_times,
                               enrolment_year,
                               lambda,
                               total_buckets,
                               weights_fu,
                               enrolment_weight,
                               rho,
                               rho_v,
                               is_vaccinated,
                               vac_buckets,
                               vac_year)
  }

  if(serostatus == 1) {

    # `min_inf` represents the number of infections that we are certain the
    #   individual has experienced. Since the individual is seropositive, we
    #   know that the individual has experienced at least one infection.
    min_inf <- 1

    # Probability of a given number of infections at enrolment
    prob_n_inf_e <- calculate_prob_n_inf(min_inf, total_buckets, lambda_e)

    # Number of filled buckets by natural infection for each
    #   scenario (# of previous infections) at enrolment
    n_filled_buckets <- vector(mode = "numeric", length = total_buckets)

    B_rho <- max(0, rho * n_e * (1 - 1 / (total_buckets * lambda_e)));

    if(is_vaccinated == 1 && vac_buckets > 0) B_rho <- 0

    # This for-loop estimates re-susceptibility from birth until enrolment
    for(n_inf_e in 1:total_buckets)
    {
      # B_rho represents leaking immunity due to re-susceptibility
      n_filled_buckets[n_inf_e] <- max(0, n_inf_e - B_rho);
    }

    y <- sim_infection_seropos(follow_up_times     = follow_up_times,
                               total_buckets       = total_buckets,
                               enrolment_year      = enrolment_year,
                               enrolment_weight    = enrolment_weight,
                               weights_fu          = weights_fu,
                               prob_n_inf_e        = prob_n_inf_e,
                               init_filled_buckets = n_filled_buckets,
                               lambda              = lambda,
                               rho                 = rho,
                               rho_v               = rho_v,
                               is_vaccinated       = is_vaccinated,
                               vac_buckets         = vac_buckets,
                               is_vac_prob         = is_vac_prob)
  }
  y
}
