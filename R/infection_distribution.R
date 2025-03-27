prob_n_infections <- function(n_i, lambda_delta_t, n_slots) {

  n_perm <- choose(n_slots, n_i);

  prob_pos         <- (1 - exp(-lambda_delta_t))^n_i;
  prob_neg         <- exp(-lambda_delta_t)^(n_slots - n_i);
  prob_n_infection <- n_perm * prob_pos * prob_neg;

  prob_n_infection
}

calculate_prob_n_inf <- function(min_inf, n_slots, lambda) {

  # Unnormalised (un) probability of a given number of infections
  un_prob_inf <- vector(length = n_slots)
  prob_slots  <- vector(length = n_slots)

  nc <- 0 # normalising constant

  for(i in 1:n_slots)   {
    un_prob_inf[i] <- ifelse(min_inf > i,
                             0,
                             prob_n_infections(i, lambda, n_slots))
  }

  nc <- sum(un_prob_inf);

  if(nc == 0)
  {
    for(i in 1:n_slots) prob_slots[i] <- 0;
  }

  if(nc > 0)
  {
    for(i in 1:n_slots) prob_slots[i] = un_prob_inf[i] / nc;
  }

  prob_slots
}
