#' get a 7-day grid over ordinal days

get_tgrid <- function() {
  out <- data.frame(day_of_year = 1:366)
  out$tgrid <- floor((out$day_of_year - 1)/7) + 1
  out
}
