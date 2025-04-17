#' Simulate a random walk using Geometric Brownian Motion with no drift
#'
#' @param n_t An integer indicating the number of years
#' @param dt A real number indicating the time step for the Euler-Maruyama algorithm
#' @param sigma A real number corresponding to the volatility parameter
#' @param x0  A real number for the initial value
#'
#' @returns a numeric vector
#' @export
#'
#' @examples
#' random_walk_lambda(26, 0.01, 0.15, 0.1)
random_walk_lambda <- function(n_t, dt, sigma, x0)
{
  time_vector <- seq(0, n_t, by = dt)
  GBM         <- vector(mode = "numeric", length = length(time_vector))
  GBM[1]      <- x0

  # Euler–Maruyama simulation
  for (i in 2:length(time_vector))
  {
    z      <- stats::rnorm(1)
    dW     <- sqrt(dt) * z
    GBM[i] <- GBM[i - 1] + sigma * GBM[i - 1] * dW
  }

  GBM[time_vector[time_vector %in% 1:n_t]]
}
