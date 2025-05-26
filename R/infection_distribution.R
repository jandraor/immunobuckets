calculate_prob_n_inf <- function(lambda_avg, n_e, total_buckets) {

  #  Probability of acquiring infection given number of previous infections
  p_vec <- 1 - exp(-(total_buckets - 0:total_buckets) * lambda_avg)

  # states[i+1] = P(exactly i infections at the current time)
  states    <- numeric(total_buckets + 1)
  states[1] <- 1L

  for (yr in seq_len(n_e))
  {
    # prob of staying in the current state
    states_tilde <- states * (1 - p_vec)
    # adding prob of moving to the next state
    states_tilde[-1] <- states_tilde[-1] + states[-length(states)] * p_vec[-length(p_vec)]
    states <- states_tilde
  }

  states
}

estimate_conditional_prob <- function(prob_vec, min_inf) {
  unconditional_prob <- utils::tail(prob_vec, -min_inf)
  normalising_const  <- sum(unconditional_prob)

  unconditional_prob / normalising_const
}
