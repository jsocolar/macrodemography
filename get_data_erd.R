library(dggridR)

get_data <- function(data, year,
                     tgrid_min, tgrid_max,
                     min_lat, max_lat, min_lon, sg) {
  cell_data <- list()
  
  zfd <- data
  zfd <- as.data.frame(zfd)
  zfd <- zfd[which(zfd$year == year), ]
  zfd <- merge(zfd, get_tgrid(), by = "day_of_year", all = T)
  zfd <- zfd[zfd$tgrid >= tgrid_min & zfd$tgrid <= tgrid_max, ]
  zfd <- zfd[zfd$latitude > min_lat & zfd$latitude < max_lat & zfd$longitude > min_lon, ]
  zfd$cells6 <- dgGEO_to_SEQNUM(dg6, zfd$longitude, zfd$latitude)$seqnum
  zfd$cells11 <- dgGEO_to_SEQNUM(dg11, zfd$longitude, zfd$latitude)$seqnum
  
  pixels <- get_pixels(sg=sg)
  all_cells_11 <- pixels$conus11$cell[pixels$conus11$cell6 %in% zfd$cells6]
  
  t_counter <- 0
  for (t in tgrid_min:tgrid_max) {
    t_counter <- t_counter + 1
    cell_data[[t_counter]] <- list()
    names(cell_data)[t_counter] <- paste0("tgrid_", t)
    zfd_t <- zfd[zfd$tgrid == t, ]
    for (g in seq_along(all_cells_11)) {
      cd <- list()
      obs <- zfd_t[zfd_t$cells11 == all_cells_11[g], ]
      cd$n <- nrow(obs)
      cd$n_z <- sum(obs$obs_count == "0" & obs$only_presence_reported == 0)
      cd$n_x <- sum(obs$only_presence_reported)
      not_x <- as.integer(obs$obs_count[obs$only_presence_reported == 0])
      cd$n_value <- sum(not_x > 0)
      cd$mean_positive <- mean(not_x[not_x > 0])
      cell_data[[t_counter]][[g]] <- cd
      names(cell_data[[t_counter]])[g] <- paste0("cell_", all_cells_11[g])
    }
  }
  return(cell_data)
}
