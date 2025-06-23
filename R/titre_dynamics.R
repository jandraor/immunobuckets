#' Add titres to infection history
#'
#' @param inf_history_df A data frame
#' @param long_decay_coefficient A scalar in the interval (0, 1)
#'   that modulates the long-term decay rate, representing the proportion of the
#'   long-term rate that is effectively applied to the system. A value of 1
#'   keeps the rate constant after subsequent infection.
#' @param ltr_0
#' @param long_term_rate_dynamics
#' @param peak_scaling_factor
#' @param initial_rise A numeric value in the log scale.
#' @param peak_dynamics
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
add_titre_to_infection_history <- function(
    inf_history_df,
    ltr_0,
    long_term_rate_dynamics = "decreasing",
    long_decay_coefficient  = 1,
    initial_rise,
    peak_dynamics,
    peak_options)
{
  inf_history_df$cumulative_infections <- cumsum(inf_history_df$sim_y)

  df_list <- split(inf_history_df, inf_history_df$cumulative_infections)

  n_df <- length(df_list)

  long_decay <- ltr_0

  log_peak_0 <- initial_rise

  A0 <- exp(log_peak_0) # Rise

  baseline_titre <- 0

  for(i in seq_len(n_df))
  {
    df    <- df_list[[i]]
    n_inf <- unique(df$cumulative_infection)
    n_row <- nrow(df)

    #Long-decay rate
    long_rate <- estimate_decay_rate(n_inf, ltr_0, long_term_rate_dynamics,
                                     long_decay_coefficient)

    if(n_inf == 0) log_titre <- rep(0, n_row)

    if(n_inf > 0)
    {
      if(peak_dynamics == "linear_rise")
      {
        peak_options$baseline_titre <- baseline_titre
      }

      log_peak <- estimate_log_peak(
        log_peak_0   = log_peak_0,
        n_inf        = n_inf,
        method       = peak_dynamics,
        peak_options = peak_options)

      anchor_index <- min(df$year_index)

      year_indexes <- df$year_index - anchor_index

      year_indexes <- c(year_indexes, max(year_indexes) + 1)

      t_indexes    <- year_indexes * 360

      titre        <- bi_exponential_decay(t          = t_indexes,
                                           A0         = exp(log_peak),
                                           par_sigma  = 1/30,
                                           par_lambda = long_rate,
                                           par_rho    = 0.98)
      baseline_titre <- log(titre[length(titre)])

      log_titre    <- log(titre[-length(titre)])
    }

    df$log_titre <- log_titre
    df_list[[i]] <- df
  }

  output_df <- do.call(rbind, df_list)

  rownames(output_df) <- NULL

  output_df
}


#' Simulate bi-exponential decay
#'
#' @param t Time
#' @param A0 Peak
#' @param par_sigma Rate of short-term decay
#' @param par_lambda Rate of long-term decay
#' @param par_rho Proportion of short-term decay
#'
#' @returns A numeric vector
#' @export
#'
#' @examples
#' bi_exponential_decay(365, exp(10), 0.98, 1/720, 1/30)
bi_exponential_decay <- function(t, A0, par_sigma, par_lambda, par_rho)
{
  A0 * (par_rho * exp(-par_sigma * t) +
          (1 - par_rho) * exp(-par_lambda * t))
}

# @param p A scalar corresponding to the growth exponent for the power law
# @param a A scalar corresponding to the scaling factor
# @param peak_options A list
estimate_log_peak <- function(log_peak_0, n_inf, method, peak_options)
{
  if(n_inf == 1) return(log_peak_0)

  if(method == "power_law")
  {
    a        <- peak_options$a
    p        <- peak_options$p
    return(log_peak_0 + a * (n_inf - 1)^p)
  }

  if(method == "linear_rise")
  {
    intercept      <- peak_options$intercept
    slope          <- peak_options$slope
    baseline_titre <- peak_options$baseline_titre

    linear_rise <-  intercept + slope * baseline_titre

    return(baseline_titre + linear_rise)
  }

  stop()
}

estimate_decay_rate <- function(n_inf, ltr_0, long_term_rate_dynamics,
                                long_decay_coefficient) {

  if(long_term_rate_dynamics == "decreasing")
  {
    return(ltr_0 * long_decay_coefficient ** (n_inf - 1))
  }

  if(long_term_rate_dynamics == "piecewise")
  {
    return(ifelse(n_inf > 4, 0, ltr_0))
  }
}
