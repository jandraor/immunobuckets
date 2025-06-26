#' Add titres to infection history
#'
#' @param short_options A list of parameters to configure the short-term decay
#'  rate.
#' @param long_options A list of parameters to configure the long-term decay
#'  rate.
#' @param rise_options A list of parameters to configure the rise.
#' @param inf_history_df A data frame
#'
#' @returns A data frame
#' @export
#'
#' @examples
#' inf_history_df <-
#'   data.frame(subject_id = 1,
#'              year_index = 11:23,
#'              age        = 2:14,
#'              sim_y      = 0)
#' inf_history_df$sim_y[c(6, 8, 11, 12)] <- 1
#'
#' add_titre_to_infection_history(
#' inf_history_df,
#' long_options  = list(L = 0.22, a = -0.96, b = 18.6),
#' short_options = list(intercept = 0.02,
#' slope     = 0.003342),
#' rise_options  = list(intercept = 8,
#'                      slope     = -0.8))
add_titre_to_infection_history <- function(
    inf_history_df,
    short_options,
    long_options,
    rise_options)
{
  inf_history_df$cumulative_infections <- cumsum(inf_history_df$sim_y)

  df_list <- split(inf_history_df, inf_history_df$cumulative_infections)

  n_df <- length(df_list)

  baseline_titre <- 0

  for(i in seq_len(n_df))
  {
    df    <- df_list[[i]]
    n_inf <- unique(df$cumulative_infection)
    n_row <- nrow(df)

    if(n_inf == 0) log_titre <- rep(0, n_row)

    if(n_inf > 0)
    {
      rise     <- estimate_rise(baseline_titre, rise_options)
      log_peak <- baseline_titre + rise

      # short-term decay rate
      short_rate <- estimate_short_rate(rise, short_options)

      # Long-term decay rate
      long_rate <- estimate_long_rate(c(df$age[[1]], df$age), long_options)

      long_rate <- long_rate / 360

      anchor_index <- min(df$year_index)

      year_indexes <- df$year_index - anchor_index

      year_indexes <- c(year_indexes, max(year_indexes) + 1)

      t_indexes    <- year_indexes * 360

      titre        <- bi_exponential_decay(t          = t_indexes,
                                           A0         = exp(log_peak),
                                           par_sigma  = short_rate,
                                           par_gamma = long_rate,
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
bi_exponential_decay <- function(t, A0, par_sigma, par_gamma, par_rho)
{
  A0 * (par_rho * exp(-par_sigma * t) +
          (1 - par_rho) * exp(-par_gamma * t))
}

estimate_short_rate <- function(rise, short_options) {

  short_options$intercept + short_options$slope * rise
}

estimate_rise <- function(baseline_titre, rise_options) {
  rise_options$intercept + rise_options$slope * baseline_titre
}

estimate_long_rate <- function(age, long_options)
{
  L <- long_options$L
  a <- long_options$a
  b <- long_options$b

  L / (1 + exp(-a * (age - b)))
}


