estimate_prob_after_vaccination <- function(total_buckets, vac_buckets,
                                            prob_per_n_buckets) {

  # probability of having a given number of infections post-vaccination
  prob_slots_pos <- vector("numeric", length = total_buckets)

  if(vac_buckets == 1)
  {
    for(i in 1:total_buckets)
    {
      if(i == 1) prob_slots_pos[i] <- (i * 1.0 / total_buckets) * prob_per_n_buckets[i];

      if(i > 1)
      {
        prob_slots_pos[i] <- (i * 1.0 / total_buckets) * prob_per_n_buckets[i] +
          ((total_buckets - (i - 1)) * 1.0 / total_buckets) * prob_per_n_buckets[i - 1];
      }
    }
    return(prob_slots_pos)
  }

  if(vac_buckets == 2) {

      prob_slots_pos[1] <- 0

      prob_slots_pos[2] <- (3 / 6) * prob_per_n_buckets[1] +
        (1 / 6) * prob_per_n_buckets[2]

      prob_slots_pos[3] <- (3 / 6) * prob_per_n_buckets[1] +
        (4 / 6) * prob_per_n_buckets[2] + (3 / 6) * prob_per_n_buckets[3]

      prob_slots_pos[4] <- (1 / 6) * prob_per_n_buckets[2] +
        (3 / 6) * prob_per_n_buckets[3] +
        prob_per_n_buckets[4]


  }
  prob_slots_pos
}

