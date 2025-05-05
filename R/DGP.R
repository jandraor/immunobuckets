determine_weight <- function(yr, follow_up_times, weights_fu, current_idx)
{
  weight <- 1; # Default value for years when there are no measurements.

  # In other years, if the current evaluated year (yr) matches the year of
  #   the current measurement (current_idx), the weight represents the period
  #   from the start of the year to the measurement's date.
  if(yr == follow_up_times[current_idx]) {
    weight = weights_fu[current_idx];
  }

  weight
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
simulate_infections <- function(total_buckets, n_e, avg_lambda_e, lambda,
                                start_index, init_weight, final_weight,
                                serostatus, is_vaccinated, rho_v, rho,
                                vac_buckets, follow_up_times, enrolment_year,
                                weights_fu, enrolment_weight, is_vac_prob,
                                switch_rho, vac_year) {

  rho   <- rho * switch_rho
  rho_v <- rho_v * switch_rho

  #lambda_e: Cumulative force of infection from birth to enrolment
  # lambda_e <- estimate_lambda_e(n_e, lambda, start_index, init_weight,
  #                               final_weight)
  #p_e: Probability of being seropositive at enrolment
  # p_e <- (1 - exp(-total_buckets * lambda_e))

  prob_n_prev_inf <- calculate_prob_n_inf(lambda_avg    = avg_lambda_e,
                                          n_e           = n_e,
                                          total_buckets = total_buckets)

  p_e <- 1 - prob_n_prev_inf[[1]]

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
                               vac_year,
                               is_vac_prob = is_vac_prob)
  }

  if(serostatus == 1) {

    # `min_inf` represents the number of infections that we are certain the
    #   individual has experienced. Since the individual is seropositive, we
    #   know that the individual has experienced at least one infection.
    min_inf <- 1

    # Probability of a given number of infections at enrolment
    # prob_n_inf_e <- calculate_prob_n_inf(min_inf, total_buckets, lambda_e)
    prob_n_inf_e <- estimate_conditional_prob(prob_n_prev_inf, min_inf)

    # Number of filled buckets by natural infection for each
    #   scenario (# of previous infections) at enrolment
    n_filled_buckets <- vector(mode = "numeric", length = total_buckets)

    # B_rho <- max(0, rho * n_e * (1 - 1 / (total_buckets * lambda_e)))

    B_rho <- max(0, rho * (n_e - 1 / (total_buckets * avg_lambda_e)))

    if(is_vaccinated == 1 && vac_buckets > 0) B_rho <- 0

    # This for-loop estimates re-susceptibility from birth until enrolment
    for(n_inf_e in 1:total_buckets)
    {
      # B_rho represents leaking immunity due to re-susceptibility
      n_filled_buckets[n_inf_e] <- max(0, n_inf_e - B_rho)
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
