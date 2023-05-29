
# In fact, I only need params$hexagon_area_large from all these! 

params <- list()
params$erd_path <- "~/Documents/erd/erd.db"
params$output_path <- "~/Documents/macrodemography/data/"
params$years <- c(2006:2019)
params$extent_space <-  data.frame( min_lon=-125, max_lon=-66, min_lat=17, max_lat=65 )
params$period <- c("spring", "fall")
params$time_grid <- 7
params$tgrid_min <- c(17, 32)   # Please note, I have changed this according to pre-post breeding of E.Phoebe!
params$tgrid_max <- c(20, 35)   # i.e., 19Apr-17May = weeks17-20; 03Aug-31Aug=weeks32-35  
params$max_altitude <- 2000
params$max_altitude_above_lat42 <- 1500
params$effort_thresholds <- data.frame( dist_max=3, time_min=5/60, time_max=1, cci_min=0 )
params$hexagon_area_large <- 70000
params$hexagon_area_small <- 300
params$n_small_min <- 10
params$n_year_min <- 5
params$daymet <- data.frame( label=c("tmax_winter","tmax_summer","swe"), variable=c("tmax","tmax","swe"), date_min=c("01-01","07-01","12-01"), date_max=c("02-28","08-31","03-15"), period=c("spring","fall","spring") )
params$species_to_process <- c("easpho")
params$always_import_checklists  <- FALSE
params$always_filter_checklists  <- FALSE
params$always_resample_bootstrap <- FALSE
params$always_run_variance_test  <- FALSE
params$always_download_weather   <- TRUE
params$always_run_regressions    <- FALSE
params$quiet                     <- TRUE
params$region <- "north_america"
params$plotting_xlim <- c(-126, -65)
params$extent_time <-
  data.frame(
    period = params$period,
    tgrid_min = params$tgrid_min,
    tgrid_max = params$tgrid_max,
    year_min = min(params$years),
    year_max = max(params$years)
  )



# Install an older version of dggridR from source. See:
# https://github.com/r-barnes/dggridR/issues/63#issuecomment-1454929653
remotes::install_github("r-barnes/dggridR", ref = "ec2a040")
# library(macrodemography)
library(data.table) # don't worry about openMP warnings on Mac.
library(brms)
library(ggplot2)
library(dplyr)
library(magrittr)
library(sf)
library(dtplyr) # enable dplyr for data.table
library(assertthat)
library(ebirdst)
library(fasterize)
library(dggridR)
library(rgee)
library(tidyr)


# for an area of ~ 70000 km^2, we get a resolution of 6:
grid_large <- dggridR::dgconstruct(area = params$hexagon_area_large)
gl <- dggridR::dgearthgrid(grid_large, frame = FALSE)
gl.cell <- data.frame(cell = sp::getSpPPolygonsIDSlots(gl), 
                      row.names = sp::getSpPPolygonsIDSlots(gl))

global <- sf::st_as_sf(SpatialPolygonsDataFrame(gl, gl.cell))

cells_all <-
  c(726, 529, 622, 363, 526, 3727, 618, 610, 621, 471, 595, 3672, 
    590, 553, 591, 648, 646, 443, 390, 724, 723, 498, 701, 645, 588, 
    593, 417, 643, 527, 676, 675, 567, 561, 620, 563, 583, 649, 696, 
    554, 670, 669, 565, 702, 391, 562, 592, 470, 444, 3700, 501, 
    612, 3673, 673, 729, 594, 619, 611, 697, 642, 564, 668, 951, 
    421, 558, 568, 728, 555, 534, 730, 528, 447, 419, 479, 639, 566, 
    617, 3670, 703, 644, 725, 507, 3751, 473, 585, 647, 589, 672, 
    734, 418, 557, 364, 727, 586, 736, 616, 615, 540, 691, 392, 695, 
    582, 671, 674, 641, 733, 422, 732, 513, 525, 894, 3699, 535, 
    3778, 667, 3671, 420, 389, 3754, 506, 614, 698, 499, 446, 952, 
    509, 699, 393, 3808, 536, 738, 533, 3690, 923, 640, 541, 812, 
    867, 581, 530, 979, 560, 474, 587, 416, 3719, 3663, 505, 893, 
    366, 449, 336, 365, 362, 839, 475, 445, 508, 472, 700, 3781, 
    450, 556, 311, 559, 613, 504, 285, 636, 502, 532, 785, 477, 635, 
    731, 531, 537, 451, 3692, 309, 335, 584, 476, 813, 424, 398, 
    500, 503, 538, 510, 3720, 3748, 338, 3667, 3750, 3753, 922, 3752, 
    3691, 720, 608, 638, 840, 868, 719, 478, 841, 448, 480, 3726, 
    866, 337, 310, 339, 3749, 256, 255, 453, 312, 370, 396, 369, 
    455, 3721, 664, 637, 722, 3693, 3718, 665, 3664, 3806, 3777, 
    394, 539, 512, 483, 395, 484, 481, 663, 257, 486, 609, 514, 423, 
    316, 452, 895, 3723, 666, 367, 735, 758, 284, 511, 314, 3776, 
    375, 3666, 3665, 950, 3779, 896, 485, 924, 721, 921, 3747, 431, 
    454, 430, 717, 690, 320, 3724, 692, 282, 482, 458, 343, 403, 
    374, 373, 580, 342, 286, 693, 348, 3697, 397, 634, 425, 283, 
    402, 368, 346, 345, 340, 3696, 401, 313, 497, 694, 457, 456, 
    737, 347, 427, 308)


##########################


library(tictoc)

tic()
# Load required packages
library(rgee)
library(sf)

# Authenticate GEE account
ee_Initialize()

# Define study area polygon (hexagon cells)
study_area_sf <- global[global$cell %in% cells_all, ]
study_area_sf <- study_area_sf[1,]

# Define date range
focus_period <- seq(as.Date("2006-03-01"), as.Date("2019-05-31"), by = "day")

# Filter the dates to include only spring months
daterange <- focus_period[format(focus_period, "%m") %in% c("03", "04", "05")]
daterange <- seq(as.Date("1980-01-01"), as.Date("1980-02-01"), by = "day")

# Checkpoint file path
checkpoint_file <- "checkpoint_spring_t.rds"

# Check if checkpoint file exists
if (file.exists(checkpoint_file)) {
  # If checkpoint file exists, load the last successfully processed cell and date
  checkpoint <- readRDS(checkpoint_file)
  start_cell <- checkpoint$cell # start from the last cell
  start_date <- as.Date(checkpoint$date) + 1  # start from the next date
} else {
  # If checkpoint file does not exist, start from the beginning
  start_cell <- 1
  start_date <- as.Date("1980-01-01")
}

# Load previously saved data if it exists
if (file.exists("test_data.rds")) {
  df <- readRDS("test_data.rds")
} else {
  # If dataframe does not exist, initialize an empty one
  df <- data.frame(
    cell_id = character(),
    date = character(),
    tmin = numeric(),
    tmax = numeric(),
    stringsAsFactors = FALSE
  )
}

# Extract daily weather data for each hexagon cell and save results after each iteration
for (i in start_cell:nrow(study_area_sf)) {
  for (date in as.character(daterange[daterange >= start_date])) {
    
    # sf as ee object
    polygon <- sf_as_ee(study_area_sf[i, ])
    
    collection <- ee$ImageCollection("NASA/ORNL/DAYMET_V4")$
      filterDate(date)$
      filterBounds(polygon)$
      select(c("tmin", "tmax"))$
      mean()$
      reduceRegion(
        reducer = ee$Reducer$median(),
        geometry = polygon,
        scale = 1000,
        maxPixels = 1e13,
        tileScale = 4
      )$getInfo()
    
    tmin <- collection$tmin
    tmax <- collection$tmax
    
    # Check for NA values before appending to the dataframe
    if (is.na(tmin) || is.na(tmax)) {
      print(paste("Missing or invalid values for cell", study_area_sf$cell[i], "and date", as.character(date)))
      print(paste("tmin:", tmin))
      print(paste("tmax:", tmax))
    } else {
      df <- rbind(df, data.frame(cell_id = study_area_sf$cell[i], date = as.character(date), tmin = tmin, tmax = tmax, stringsAsFactors = FALSE))
    }
    
    # Print progress of the extraction
    cat(paste("Cell:", study_area_sf$cell[i], "- Date:", as.character(date), "\n"))
    
    # Save the last successfully processed cell and date as checkpoint
    checkpoint <- list(cell = i, date = as.character(date))
    saveRDS(checkpoint, checkpoint_file)
    
    # Save dataframe as .rds after each iteration to avoid data loss due to sudden failure
    saveRDS(df, file = "test_data.rds")
    
  }
}

toc()
