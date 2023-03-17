#' sample a species given thresholds on effort, space and time
#' @param species_code six letter species code
#' @param erd_path path to erd
#' @param checklists filtered checklists
#' @param roi region of interest
#' @param effort_thresholds effort thresholds
#' @param extent_space spatial extent
#' @param extent_time temporal extent 
#' @param time_window stratify by time grid or not
#' @param small_grid resolution of small grid
#' @param large_grid resolution of large grid 
#' @param time_grid resolution of time 
#' @param .cores=4 cores for parallel computation
#' @export
sample_grid_abun <- function(species_code, erd_path, checklists, roi, effort_thresholds, extent_space, extent_time, time_window="full", small_grid=11, large_grid=6, time_grid=7, .cores=4){
  # verify input arguments
  assert_that(is.character(species_code))
  path_erd <- paste0(erd_path, "/erd.db")
  assert_that(file.exists(path_erd))
  assert_that(is.data.frame(checklists))
  assert_that(is.data.frame(effort_thresholds))
  assert_that(all(c("dist_max","time_min","time_max","cci_min") %in% colnames(effort_thresholds)), msg="missing threshold value(s) for dist_max, time_min, time_max, cci_min")
  assert_that(effort_thresholds$time_min < effort_thresholds$time_max)
  assert_that(effort_thresholds$dist_max > 0)
  assert_that(is.data.frame(extent_space))
  assert_that(all(c("min_lon","max_lon","min_lat","max_lat") %in% colnames(extent_space)), msg="missing threshold value(s) for min_lon, max_lon, min_lat, max_lat")
  assert_that(extent_space$min_lon < extent_space$max_lon)
  assert_that(extent_space$min_lat < extent_space$max_lat)
  assert_that(is.data.frame(extent_time))
  assert_that(all(extent_time$year_min <= extent_time$year_max))
  assert_that(all(extent_time$tgrid_min < extent_time$tgrid_max))
  assert_that(all(c("period","tgrid_min","tgrid_max","year_min", "year_max") %in% colnames(extent_time)), msg="missing threshold value(s) for period, tgrid_min, tgrid_max, year_min, year_max")
  assert_that(time_window %in% c("gridded", "full"))
  
  # load species data from ERD
  # takes ~ 2-3 mins for carwre (Carolina Wren)
  sp_data <- import_from_erd(
    species_code,
    erd_path = path_erd,
    checklists = checklists
  )
  
  # loop over time periods
  data_grid <- data_abun <- list()
  for (i in 1:nrow(extent_time)) {
    # loop over years
    years=extent_time$year_min[i]:extent_time$year_max[i]
    # initialize year lists
    data_grid[[i]] <- data_abun[[i]] <- list()
    for (y in seq_along(years)) {
      print(paste("grid sampling",species_code,"data for",extent_time$period[i],years[y],"..."))
      
      data_grid[[i]][[y]] <- 
        get_grid_data(
          data = sp_data, .year = years[y], 
          tgrid_min = extent_time$tgrid_min[i], 
          tgrid_max = extent_time$tgrid_max[i],
          time_window = time_window, min_lat = extent_space$min_lat, 
          max_lat = extent_space$max_lat, min_lon = extent_space$min_lon, 
          max_lon = extent_space$max_lon, large_grid=large_grid, 
          small_grid=small_grid, time_grid=time_grid, roi = roi, 
          .cores = .cores
        )
      
      data_abun[[i]][[y]] <- get_abun(data_grid[[i]][[y]], n_rep=100, roi = roi)
    }
    names(data_grid[[i]]) <- years
    names(data_abun[[i]]) <- years
  }
  names(data_grid) <- extent_time$period
  names(data_abun) <- extent_time$period
  
  output <- list(grid=data_grid, abun=data_abun)
}


#' Get log-ratios
#' @param data abundance index data
#' @param cells_all cells to compute
#' @param period period to compute
get_ratios <- function(data, cells_all, period=c("spring", "fall")){
  # verify we have all periods (spring and fall) available
  assert_that(all(period %in% names(data)))
  assert_that(length(period)==2)
  assert_that(is.character(period))
  # verify that loaded data contains no unknown grid cells
  # QUESTION: in what situation could this evaluate to FALSE?
  cells_present <- sapply(period, function(x) assert_that(all(data[[x]][[1]]$cell %in% cells_all)))
  assert_that(all(cells_present), msg="loaded data contains unknown grid cells")
  
  message(paste("calculating ratios with period1 =",period[1],"and period2 =",period[2]))
  
  # Extract cell-specific abundance data
  spring_abun_summary <- get_abun_summary(data[[period[1]]], n_small_min)
  fall_abun_summary <- get_abun_summary(data[[period[2]]], n_small_min)
  cell_timeseries <- get_cell_timeseries(cells_all,
                                         spring_abun_summary,
                                         fall_abun_summary)
  # Take the log-ratios along the timeseries to get a ratio series
  # cell_ratio_series is a summary of the bootstrap uncertainty, and
  # cell_ratio_series_full gives all bootstrap replicates.
  cells <- unique(spring_abun_summary[[1]]$cell)
  cell_ratio_series <- cell_ratio_series_full <- list()
  ratio_series_length <- 2*length(spring_abun_summary)-1
  
  for (i in 1:length(cell_timeseries)) {
    if (cells_all[i] %in% cells) {
      # calculate logarithmic ratios (lrats)
      lrats <- apply(cell_timeseries[[i]][ ,2:ncol(cell_timeseries[[i]])], 2, function(x){log(stocks::ratios(x))})
      cell_ratio_series[[i]] <- list()
      
      # cell / year / reference period (named after the latest period on which the ratio is based)
      cell_ratio_series[[i]]$cell <- rep(cells_all[i],ratio_series_length)
      cell_ratio_series[[i]]$year <- as.numeric(sapply(names(spring_abun_summary), function(x) rep(x,times=2)))[-1]
      cell_ratio_series[[i]]$period <- rep(rev(period), length.out=ratio_series_length)
      
      cell_ratio_series[[i]]$median <- apply(lrats, 1, median)
      cell_ratio_series[[i]]$avg <- apply(lrats, 1, mean)
      cell_ratio_series[[i]]$sd <- apply(lrats, 1, sd)
      cell_ratio_series[[i]]$q10 <- apply(lrats, 1, function(x){quantile(x, .1, na.rm = T)})
      cell_ratio_series[[i]]$q90 <- apply(lrats, 1, function(x){quantile(x, .9, na.rm = T)})
      
      cell_ratio_series_full[[i]] <- lrats
    } else {
      cell_ratio_series[[i]] <- cell_ratio_series_full[[i]] <- NA
    }
  }
  
  # make output a named list:
  names(cell_ratio_series) <- names(cell_timeseries)
  names(cell_ratio_series_full) <- names(cell_timeseries)
  
  return(list(summary=cell_ratio_series, replicates=cell_ratio_series_full))
}