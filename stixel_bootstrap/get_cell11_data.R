# species and sampling files should already be subset to complete checklists
# in the lower 48 for the years desired.
get_data <- function(sp_subset_path, sampling_subset_path, 
                     tgrid_min, tgrid_max,
                     max_lat, min_lon) {
  cell_data <- list()

  zfd <- auk_zerofill(sp_subset_path, sampling_subset_path, collapse = T)
  zfd <- zfd[zfd$protocol_type %in% c("Stationary", "Area", "Traveling"), ]
  zfd$oday <- as.POSIXlt(zfd$observation_date)$yday
  zfd <- merge(zfd, get_tgrid(), all = T)
  zfd <- zfd[zfd$tgrid >= tgrid_min & zfd$tgrid <= tgrid_max, ]
  zfd <- zfd[zfd$latitude < max_lat & zfd$longitude > min_lon, ]
  zfd$cells6 <- dgGEO_to_SEQNUM(dg6, zfd$longitude, zfd$latitude)$seqnum
  zfd$cells11 <- dgGEO_to_SEQNUM(dg11, zfd$longitude, zfd$latitude)$seqnum
  
  pixels <- get_pixels()
  all_cells_11 <- pixels$conus11$cell[pixels$conus11$cell6 %in% zfd$cells6]
  
  t_counter <- 0
  for (t in tgrid_min:tgrid_max) {
    t_counter <- t_counter + 1
    cell_data[[t_counter]] <- list()
    names(cell_data)[t_counter] <- paste0("tgrid_", t)
    zfd_t <- zfd[zfd$tgrid == t, ]
    for (g in seq_along(all_cells_11)) {
      #print(c(t,g))
      cd <- list()
      obs <- zfd_t[zfd_t$cells11 == all_cells_11[g], ]
      cd$n <- nrow(obs)
      cd$n_z <- sum(obs$observation_count == "0")
      cd$n_x <- sum(obs$observation_count == "X")
      not_x <- as.integer(obs$observation_count[obs$observation_count != "X"])
      cd$n_value <- sum(not_x > 0)
      cd$mean_positive <- mean(not_x[not_x > 0])
      cell_data[[t_counter]][[g]] <- cd
      names(cell_data[[t_counter]])[g] <- paste0("cell_", all_cells_11[g])
    }
  }
  return(cell_data)
}

get_cell11_data <- function (species, year, season, max_lat = NULL) {
  sampling_subset_path <- paste0("/Users/jacob/Desktop/ebird_files/Apr-2021/sampling/by_year/",
                                 season, "_", year, ".txt")
  sp_subset_path <- paste0("/Users/jacob/Desktop/ebird_files/Apr-2021/", 
                           species, "/by_year/", season, "_", year, ".txt")
  sp_dat <- spp_data[spp_data$species == species, ]
  
  if(is.null(max_lat)){max_lat <- sp_dat$max_lat}
  min_lon <- min(c(sp_dat$min_lon_1, sp_dat$min_lon_2), na.rm = T)
  tgrid <- get_tgrid()
  if(season == "spring"){
    tgrid_min <- tgrid$tgrid[tgrid$oday == sp_dat$min_day_spring]
    tgrid_max <- 21
  } else if (season == "fall"){
    tgrid_min <- 34
    tgrid_max <- tgrid$tgrid[tgrid$oday == sp_dat$max_day_fall]
  }
  out <- get_data(sp_subset_path, sampling_subset_path, 
                  tgrid_min, tgrid_max,
                  max_lat, min_lon)
  return(out)
}
