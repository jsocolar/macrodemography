
params$always_download_weather <- TRUE

# i have defined datamet as:
params$daymet <- data.frame( label=c("tmax_winter","tmax_summer","swe"), variable=c("tmax","tmax","swe"), date_min=c("01-01","07-01","12-01"), date_max=c("02-28","08-31","03-15"), period=c("spring","fall","spring") )


flag_file <- paste0(params$output_path, "/weather/flag.rds")
# flag_file <- "flag.rds"

if (params$always_download_weather | !file.exists(weather_file) | !file.exists(flag_file)) {
  # initialize google earth engine
  ee_Initialize()
  
  # load or create flag variable to indicate whether data extraction was successful or not
  flag <- if (file.exists(flag_file)) readRDS(flag_file) else data.frame(cell = character(), year = integer(), success = logical())
  
  # check if weather file exists and load data
  data_daymet <- if (file.exists(weather_file)) readRDS(weather_file) else data.frame()
  
  # loop over cells and years
  for (cell in cells_all) {
    for (year in params$years) {
      # check if data extraction was successful for this cell and year
      if (nrow(flag) == 0 || !any(flag$success[flag$cell == cell & flag$year == year])) {
        # attempt data extraction
        tryCatch({
          print(paste("year =", year, "cell =", cell))
          data_daymet <- rbind(data_daymet, daymet_set_extract(year, cell, grid_large, params$daymet))
          
          # update flag variable to indicate success
          flag <- rbind(flag, data.frame(cell = cell, year = year, success = TRUE))
          saveRDS(flag, flag_file)
          
          # save weather data after each iteration of the inner loop
          saveRDS(data_daymet, weather_file)
        }, error = function(e) {
          # if an error occurred during data extraction, print an error message
          message(paste("Error occurred for year =", year, "cell =", cell))
          # update flag variable to indicate failure
          flag <- rbind(flag, data.frame(cell = cell, year = year, success = FALSE))
          saveRDS(flag, flag_file)
        })
      }
    }
  }
} else {
  # load weather data if always_download_weather is false and the weather file exists
  data_daymet <- readRDS(weather_file)
}
