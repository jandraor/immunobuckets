#' Estimate the cumulative force of infection at enrolment
#'
#' @param n_e A non-negative integer.
#' @param yearly_lambda A numeric vector
#' @param start_index A positive integer indicating the year in which the
#'  individual became susceptible after birth.
#' @param init_weight A numeric value corresponding to the fraction of the first
#'  year in which the individual was susceptible.
#' @param final_weight A numeric value representing the fraction of the calendar
#'  year from the 1st Jan to the date of enrolment
#'
#' @returns A real number
#' @export
#'
#' @examples
#' estimate_lambda_e(3, rep(0.1, 26), 1, 0.25, 0.75)
estimate_lambda_e <- function(n_e, yearly_lambda, start_index,
                              init_weight, final_weight) {

  lambda_e <- 0

  if(n_e == 0) return (0)

  if(n_e == 1) return((final_weight - init_weight) * yearly_lambda[start_index])

  for (i in 1:n_e) # for every year before enrolment
  {
    # fraction of a calendar year
    weight = 1;  # For all years except the first & last ones.

    if(i == 1)   weight = init_weight # First year

    if(i == n_e) weight = final_weight # Year of enrolment

    lambda_e = lambda_e + weight * yearly_lambda[start_index + i - 1];
  }

  lambda_e
}
