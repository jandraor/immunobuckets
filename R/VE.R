#' Title
#'
#' @param follow_up_df A data frame
#'
#' @returns A data frame
#' @export
#'
#' @examples
#'   follow_up_df <- data.frame(subject_id    = rep(1:4, 5),
#'                              year_index    = rep(1:5, each = 4),
#'                              is_vaccinated = rep(c(1, 1, 0, 0), 5),
#'                              infection     = c(0, 0, 1, 0,
#'                                                0, 1, 0, 1,
#'                                                1, 0, 0, 0,
#'                                                0, 0, 0, 0,
#'                                                0, 0, 0, 0))
#'  estimate_cumulative_VE(follow_up_df)
estimate_cumulative_VE <- function(follow_up_df)
{
  infections_df <- follow_up_df[follow_up_df$infection == 1, ]

  earliest_inf_df <- stats::aggregate(year_index ~ subject_id,
                                      data = infections_df,
                                      FUN = min)

  colnames(earliest_inf_df)[2] <- "earliest_infection"

  VE_input <- follow_up_df |> dplyr::left_join(earliest_inf_df,
                                               by = "subject_id")

  VE_input$earliest_infection <- ifelse(is.na(VE_input$earliest_infection),
                                        Inf,
                                        VE_input$earliest_infection)

  VE_input$marker <- ifelse(VE_input$year_index > VE_input$earliest_infection,
                            0, 1)

  incidence_df <- stats::aggregate(cbind(n_inf = infection, n_people = marker) ~
                              year_index + is_vaccinated,
                            data = VE_input, FUN = sum)

  df_list <- split(incidence_df, incidence_df$is_vaccinated) |>
    lapply(\(df) {
      df$cumulative_incidence <- cumsum(df$n_inf)
      df$persontime           <- cumsum(df$n_people)
      df$attack_rate          <- df$cumulative_incidence / df$persontime
      df
  })

  attack_rate_df <- do.call("rbind", df_list)

  attack_rate_df$vac_group <- ifelse(attack_rate_df$is_vaccinated == 1,
                                     "Vaccine", "Placebo")

  long_df <- attack_rate_df[, c("year_index", "vac_group", "attack_rate")]

  wide_df <- long_df |>
    tidyr::pivot_wider(names_from = "vac_group",
                       values_from = "attack_rate")

  wide_df$VE <- 1 - wide_df$Vaccine / wide_df$Placebo

  wide_df[, c("year_index", "VE")] |> as.data.frame()
}
