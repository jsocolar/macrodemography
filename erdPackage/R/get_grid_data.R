#' Get eBird summaries on cells of a small grid, tagged with membership in a large grid.
#' @param data zero-filled data (output from import_from_erd())
#' @param .year the year
#' @param tgrid_min the minimum week of year
#' @param tgrid_max the maximum week of year
#' @param min_lat the minimum latitude (decimal degrees) where the prediction is desired
#' @param max_lat the maximum latitude
#' @param min_lon the minimum longitude
#' @param large_grid the dggridR resolution of the large grid
#' @param small_grid the dggridR resolution of the small grid
#' @export
get_grid_data <- function(data, .year,
                     tgrid_min, tgrid_max,
                     min_lat, max_lat, min_lon,
                     large_grid = 6, small_grid = 11,
                     .cores = 4) {
  cell_data <- list()
  
  zfd <- data
  zfd <- zfd[year == .year]
  zfd <- zfd[latitude > min_lat & latitude < max_lat & longitude > min_lon]
  zfd <- merge(zfd, get_tgrid(), by = "day_of_year", all = F)
  zfd <- zfd[tgrid >= tgrid_min & tgrid <= tgrid_max]
  
  dg_large <- dggridR::dgconstruct(res=large_grid)
  dg_small <- dggridR::dgconstruct(res=small_grid)
  
  
  zfd$cells_large <- dggridR::dgGEO_to_SEQNUM(dg_large, zfd$longitude, zfd$latitude)$seqnum
  zfd$cells_small <- dggridR::dgGEO_to_SEQNUM(dg_small, zfd$longitude, zfd$latitude)$seqnum
  
  pixels <- get_pixels(dg_large = dg_large, dg_small = dg_small)
  all_cells_small <- pixels$conus_small$cell[pixels$conus_small$cell_large %in% zfd$cells_large]
  
  cl <- parallel::makeCluster(.cores, "FORK")
  doParallel::registerDoParallel(cl)
  tgrid <- tgrid_min:tgrid_max
  
  
  
  cell_data <- foreach (t = 1:length(tgrid)
           ) %dopar% {
             get_cell_data(zfd, tgrid[t], all_cells_small)
           }
  
  parallel::stopCluster(cl = cl)
  
  names(cell_data) <- paste0("tgrid_", tgrid)
  
  attr(cell_data, "species") <- attributes(data)$species
  attr(cell_data, "year") <- .year
  attr(cell_data, "large_grid") <- large_grid
  attr(cell_data, "small_grid") <- small_grid
  attr(cell_data, "min_lat") <- min_lat
  attr(cell_data, "max_lat") <- max_lat
  attr(cell_data, "min_lon") <- min_lon
  
  return(cell_data)
}





#' get a 7-day grid over ordinal days
#' @return a two-column data-frame that tags each day-of-year to a tgrid cell
get_tgrid <- function() {
  out <- data.frame(day_of_year = 1:366)
  out$tgrid <- floor((out$day_of_year - 1)/7) + 1
  return(out)
}


#' get a hexagonal grid over the contiguous US
#' @param dg_large large dggridR grid
#' @param dg_small small dggridR grid
get_pixels <- function(dg_large, dg_small) {

  
  conus <- spData::us_states
  conus_raster <- fasterize::fasterize(conus,
                                       raster::raster(ncol=1000, nrow = 1000, 
                                                      xmn = -125, xmx = -66, 
                                                      ymn =24, ymx = 50))
  conus_coords <- raster::as.data.frame(conus_raster, xy = T)
  conus_coords <- conus_coords[!is.na(conus_coords$layer), ]
  conus_cells <- data.frame(cell = unique(dggridR::dgGEO_to_SEQNUM(dg_large, conus_coords$x, conus_coords$y)[[1]]))
  conus_cells$lat <- dggridR::dgSEQNUM_to_GEO(dg_large, conus_cells$cell)$lat_deg
  conus_cells$lon <- dggridR::dgSEQNUM_to_GEO(dg_large, conus_cells$cell)$lon_deg
  
  conus_cells_small <- data.frame(cell = unique(dggridR::dgGEO_to_SEQNUM(dg_small, conus_coords$x, conus_coords$y)[[1]]))
  conus_cells_small$lat <- dggridR::dgSEQNUM_to_GEO(dg_small, conus_cells_small$cell)$lat_deg
  conus_cells_small$lon <- dggridR::dgSEQNUM_to_GEO(dg_small, conus_cells_small$cell)$lon_deg
  conus_cells_small$cell_large <- dggridR::dgGEO_to_SEQNUM(dg_large, conus_cells_small$lon, conus_cells_small$lat)[[1]]
  
  out <- list(conus_large = conus_cells, conus_small = conus_cells_small)
  return(out)
}


#' Get mean value for a stixel
#' @param x data for a single cell
#' @return average abundance for the cell
stixel_mean_small <- function(x){
  if (x$n == 0) {
    return(NA)
  } else if (x$n_x + x$n_value == 0) {
    return(0)
  } else if (x$n_value == 0) {
    return(weighted.mean(x = c(0,2), w = c(x$n_z, x$n_x)))
  } else {
    n_z_rescaled <- x$n_z * x$n_value / (x$n_value + x$n_x)
    return(weighted.mean(x = c(0, x$mean_positive), w = c(n_z_rescaled, x$n_value)))
  }
}




#' Get data for a stixel
#' @param zfd zero-filled data
#' @param t the tgrid element
#' @param all_cells_small the list of all cells
#' @return row vector
get_cell_data <- function(zfd, t, all_cells_small) {
  out <- list()
  zfd_t <- zfd[zfd$tgrid == t, ]
  for (g in seq_along(all_cells_small)) {
    cd <- list()
    obs <- zfd_t[zfd_t$cells_small == all_cells_small[g], ]
    cd$n <- nrow(obs)
    cd$n_z <- sum(obs$obs_count == "0" & obs$only_presence_reported == 0)
    cd$n_x <- sum(obs$only_presence_reported)
    not_x <- as.integer(obs$obs_count[obs$only_presence_reported == 0])
    cd$n_value <- sum(not_x > 0)
    cd$mean_positive <- mean(not_x[not_x > 0])
    cd$stixel_mean_small <- stixel_mean_small(cd)
    out[[g]] <- cd
    names(out)[g] <- paste0("cell_", all_cells_small[g])
  }
  return(out)
}

