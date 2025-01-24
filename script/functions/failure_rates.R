treshold_failure_rate <- function(failure_rate, index_break_point, alpha) {
  
  # Check fail
  length_failure_rate <- length(failure_rate)
  failure_rate_pre_break <- failure_rate[0:index_break_point]
  failure_rate_post_break <- failure_rate[index_break_point:length_failure_rate]
  
  # Initialize Indices break points
  index_left_point <- 0
  index_right_point <- length_failure_rate
  
  # Pre break
  for (i in index_break_point:1) {
    if (failure_rate[i] == 0) {
      #cat("Left bound failure rate found at", i, '\n')
      index_left_point <- i
      break
    }
  }
  
  # Post Break
  for (j in index_break_point:length_failure_rate) {
    if (failure_rate[j] == 0) {
      #cat("Right bound failure rate found at", j, '\n')
      index_right_point <- j
      break
    }
  }
  
  width_interval <- index_right_point - index_left_point
  left_distance_from_break <- index_break_point - index_left_point
  right_distance_from_break <- index_right_point - index_break_point
  
  return(list(
    left_distance_from_break = left_distance_from_break,
    right_distance_from_break = right_distance_from_break,
    width_interval = width_interval)
         )
}
