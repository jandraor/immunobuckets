#' Simulate infections since birth
#'
#' @param individuals_df A data frame
#' @param lambda A numeric vector
#' @param rho A numeric value
#' @param stop_index An integer
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
                                            stop_index) {

  df_list <- split(individuals_df, individuals_df$subject_id)

  sim_list <- lapply(df_list, simulate_single_individual,
                     stop_index = stop_index, lambda = lambda, rho = rho)

  combined_df <- do.call(rbind, sim_list)

  row.names(combined_df) <- NULL

  combined_df
}

# Assumes perfect reporting
simulate_single_individual <- function(individual_df, stop_index, lambda, rho) {

  subject_id           <- individual_df$subject_id
  is_vaccinated        <- individual_df$is_vaccinated
  birth_year_index     <- individual_df$birth_year_index
  vac_year             <- individual_df$vac_year

  duration_maternal_ab <- 1

  start_index          <- birth_year_index + duration_maternal_ab
  year_indexes         <- start_index:stop_index

  enrolment_year   <- birth_year_index
  total_buckets    <- 4
  weights_fu       <- rep(1, length(year_indexes))
  enrolment_weight <- 1

  y <- sim_infection_seroneg(follow_up_times = year_indexes,
                             enrolment_year,
                             lambda,
                             total_buckets,
                             weights_fu,
                             enrolment_weight,
                             rho,
                             rho_v         = 0,
                             is_vaccinated = is_vaccinated,
                             vac_buckets   = 0,
                             vac_year      = vac_year)

  data.frame(subject_id = subject_id ,
             year_index = year_indexes,
             sim_y      = y)
}
