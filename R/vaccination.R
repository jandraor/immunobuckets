draw_vac_buckets <- function(inf_counter, vac_buckets, total_buckets,
                             is_vac_prob) {

  if(vac_buckets == 0) return(0)

  # Handling the case when the number of infections is higher than total buckets
  inf_counter <- min(inf_counter, total_buckets)

  if(is_vac_prob == 0)
  {
    filling_buckets <- min(total_buckets - inf_counter, vac_buckets)
    return(filling_buckets)
  }

  if(vac_buckets == 1)
  {
    prob_not_filling <- inf_counter / total_buckets
    filling_buckets  <- stats::rbinom(1, 1, 1 - prob_not_filling)
  }

  if(vac_buckets == 2) stop("Work in progress", .call = TRUE)

  filling_buckets
}
