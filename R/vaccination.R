draw_vac_buckets <- function(inf_counter, vac_buckets, total_buckets) {

  # Handling the case when the number of infections is higher than total buckets
  inf_counter <- min(inf_counter, total_buckets)

  if(vac_buckets == 1) {
    prob_not_filling <- inf_counter / total_buckets
    filling_buckets <- rbinom(1, 1, 1 - prob_not_filling)
  }

  if(vac_buckets == 2) stop("Work in progress", .call = TRUE)

  filling_buckets
}
