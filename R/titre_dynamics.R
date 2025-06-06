#' Add titres to infection history
#'
#' @param inf_history_df A data frame
#' @param long_decay_coefficient A scalar in the interval (0, 1)
#'   that modulates the long-term decay rate, representing the proportion of the
#'   long-term rate that is effectively applied to the system.
#'
#' @returns A data frame
#' @export
#'
#' @examples
#' inf_history_df <-
#'   data.frame(subject_id = 1,
#'              year_index = 11:23,
#'              sim_y      = 0)
#'
#' inf_history_df$sim_y[c(6, 8, 11, 12)] <- 1
#'
#' add_titre_to_infection_history(inf_history_df)
add_titre_to_infection_history <- function(inf_history_df,
                                           long_decay_coefficient,
                                           peak_scaling_factor)
{
  inf_history_df$cumulative_infections <- cumsum(inf_history_df$sim_y)

  df_list <- split(inf_history_df, inf_history_df$cumulative_infections)

  n_df <- length(df_list)

  long_decay <- 1/720

  log_peak_0 <- 10

  A0 <- exp(10) # Rise

  for(i in seq_len(n_df))
  {
    df    <- df_list[[i]]
    n_inf <- unique(df$cumulative_infection)
    n_row <- nrow(df)

    if(n_inf == 0) log_titre <- rep(0, n_row)

    if(n_inf > 0)
    {
      log_peak <- estimate_log_peak(log_peak_0, n_inf, peak_scaling_factor,
                                    0.5)

      anchor_index <- min(df$year_index)
      t_indexes    <- (df$year_index - anchor_index) * 360
      titre        <- bi_exponential_decay(t          = t_indexes,
                                           A0         = exp(log_peak),
                                           par_sigma  = 1/30,
                                           par_lambda = long_decay,
                                           par_rho    = 0.98)
      log_titre    <- log(titre)

      long_decay   <- long_decay * long_decay_coefficient

    }

    df$log_titre <- log_titre
    df_list[[i]] <- df
  }

  output_df <- do.call(rbind, df_list)

  rownames(output_df) <- NULL

  output_df
}


bi_exponential_decay <- function(t, A0, par_sigma, par_lambda, par_rho)
{
  A0 * (par_rho * exp(-par_sigma * t) +
          (1 - par_rho) * exp(-par_lambda * t))
}

#' @param p A scalar corresponding to the growth exponent for the power law
#' @param a A scalar corresponding to the scaling factor
estimate_log_peak <- function(log_peak_0, n_inf, a, p)
{
  log_peak_0 + a * (n_inf - 1)^p
}
