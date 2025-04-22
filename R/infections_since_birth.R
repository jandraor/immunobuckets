#' Simulate infections since birth
#'
#' @param individuals_df A data frame
#' @param lambda A numeric vector
#' @param rho A numeric value
#' @param stop_index An integer
#' @param vac_buckets An integer
#' @param rho_v A numeric value
#' @param is_vac_prob A boolean.
#'
#' @returns A data frame
#' @export
#'
#' @examples
#' individuals_df <- data.frame(subject_id       = 1,
#'                              birth_year_index = 9,
#'                              is_vaccinated    = 1)
#'
#' rho        <- 0.05
#' stop_index <- 15
#' lambda     <- rep(0.14, 15)
#'
#' simulate_infections_since_birth(individuals_df, lambda, rho, stop_index)
simulate_infections_since_birth <- function(individuals_df, lambda, rho,
                                            stop_index, vac_buckets, rho_v,
                                            is_vac_prob) {

  df_list <- split(individuals_df, individuals_df$subject_id)

  sim_list <- lapply(df_list, simulate_single_individual,
                     stop_index = stop_index, lambda = lambda, rho = rho,
                     vac_buckets = vac_buckets, rho_v = rho_v,
                     is_vac_prob = is_vac_prob)

  combined_df <- do.call(rbind, sim_list)

  row.names(combined_df) <- NULL

  combined_df
}

# Assumes perfect reporting
simulate_single_individual <- function(individual_df, stop_index, lambda, rho,
                                       vac_buckets, rho_v, is_vac_prob,
                                       total_buckets = 4) {

  subject_id           <- individual_df$subject_id
  is_vaccinated        <- individual_df$is_vaccinated
  birth_year_index     <- individual_df$birth_year_index
  vac_year             <- individual_df$vac_year

  duration_maternal_ab <- 1

  enrolment_year  <- birth_year_index + duration_maternal_ab
  start_index     <- enrolment_year + 1
  follow_up_times <- start_index:stop_index

  weights_fu       <- rep(1, length(follow_up_times))
  enrolment_weight <- 0

  y <- sim_infection_seroneg(follow_up_times  = follow_up_times,
                             enrolment_year   = enrolment_year,
                             lambda           = lambda,
                             total_buckets    = total_buckets,
                             weights_fu       = weights_fu,
                             enrolment_weight = enrolment_weight,
                             rho              = rho,
                             rho_v            = rho_v,
                             is_vaccinated    = is_vaccinated,
                             vac_buckets      = vac_buckets,
                             vac_year         = vac_year,
                             is_vac_prob      = is_vac_prob)

  data.frame(subject_id = subject_id ,
             year_index = follow_up_times,
             sim_y      = y)
}
