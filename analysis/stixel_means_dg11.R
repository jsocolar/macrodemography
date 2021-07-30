# Average eBird frequencies over stixels, stratified on a 7-day window and
# 20-km grid
source(paste0("/Users/jacob/Dropbox/Work/Code/macrodemography/data_format/",
              "ebird_raw_to_seasonal.R"))
years <- 2010:2019
library(dggridR)
library(sf)
library(ggplot2)
library(data.table)

spp_data <- read.csv("/Users/jacob/Dropbox/Work/macrodemography/include_by_spp.csv")
s_app <- read_sf("/Users/jacob/Dropbox/Work/macrodemography/masks/BLBW_mask.kml")
b_app <- read_sf("/Users/jacob/Dropbox/Work/macrodemography/masks/s_app_broad_mask.kml")
nowa <- read_sf("/Users/jacob/Dropbox/Work/macrodemography/masks/NOWA.kml")

dg11 <- dgconstruct(res = 11)

# A hack to determine all grid cells that overlap the continuous US.
# Surely there's a better way...
conus <- spData::us_states
conus_raster <- fasterize::fasterize(conus,
                                     raster::raster(ncol=1000, nrow = 1000, 
                                                    xmn = -125, xmx = -66, 
                                                    ymn =24, ymx = 50))
conus_coords <- raster::as.data.frame(conus_raster, xy = T)
conus_coords <- conus_coords[!is.na(conus_coords$layer), ]
conus_coords <- conus_coords[conus_coords$x > -111 & conus_coords$y > 29, ]
conus_cells <- unique(dgGEO_to_SEQNUM(dg11, conus_coords$x, conus_coords$y)[[1]])
length(conus_cells)

##### get_cell_data #####
seasons <- c("spring", "fall")
cell_data <- list()
for (i in seq_along(spp)) {
  cell_data[[i]] <- list()
  sp_data <- spp_data[i,]
  names(cell_data)[i] <- spp_code[i]
  output_path_2 <- paste0(output_path, ebd_month, "/", spp_code[i], "/by_year")
  sampling_path_2 <- paste0(output_path, ebd_month, "/sampling/by_year")
  for (j in seq_along(years)) {
    cell_data[[i]][[j]] <- list()
    names(cell_data[[i]])[j] <- paste0("year_", years[j])
    for (s in seq_along(seasons)) {
      cell_data[[i]][[j]][[s]] <- list()
      names(cell_data[[i]][[j]])[s] <- seasons[s]
      zfd <- auk_zerofill(paste0(output_path_2, "/", seasons[s], "_", years[j], 
                                 ".txt"),
                          paste0(sampling_path_2, "/", seasons[s], "_", 
                                 years[j], ".txt"), collapse = T)
      zfd <- as.data.table(zfd)
      lat_lon <- st_as_sf(zfd[, c("longitude", "latitude")], 
                          coords = c("longitude", "latitude"),
                          crs = 4326)
      if (sp_data$mask == "S") {
        zfd <- zfd[as.vector(!st_intersects(lat_lon , s_app, sparse = F)), ]
      } else if (sp_data$mask == "B") {
        zfd <- zfd[as.vector(!st_intersects(lat_lon , b_app, sparse = F)), ]
      } else if (sp_data$mask == "NOWA") {
        zfd <- zfd[as.vector(!st_intersects(lat_lon , nowa, sparse = F)), ]
      }
      
      zfd <- zfd[protocol_type %in% c("Stationary", "Area", "Traveling"), ]
      zfd <- zfd[latitude > sp_data$min_lat & latitude < sp_data$max_lat, ]
      zfd$days <- as.POSIXlt(zfd$observation_date)$yday
      zfd$cells11 <- dgGEO_to_SEQNUM(dg11, zfd$longitude, zfd$latitude)$seqnum
      
      if (seasons[s] == "spring") {
        zfd <- zfd[days >= sp_data$min_day_spring]
        zfd$tgrid <- floor(zfd$days/7) - 3
        tgrid_max <- 19
      } else if (seasons[s] == "fall") {
        zfd <- zfd[days <= sp_data$max_day_fall]
        zfd$tgrid <- floor(zfd$days/7) - 29
        tgrid_max <- 18
      }
      
      for (t in 1:tgrid_max) {
        cell_data[[i]][[j]][[s]][[t]] <- list()
        names(cell_data[[i]][[j]][[s]])[t] <- paste0("tgrid_", t)
        obs <- zfd[tgrid == t, ]
        
        obs_summary <- obs[,.(.N, 
                              n_z = sum(observation_count == "0"),
                              n_x = sum(observation_count == "X"),
                              n_value = sum((observation_count != "X") > 0),
                              n_positive = sum((!(observation_count %in% c("X", "0"))) > 0),
                              mean_positive = mean(as.integer(observation_count[!(observation_count %in% c("X", "0"))]), na.rm = T)
                              ), by=cells11]
        print(c(i,j,s,t))
        cell_data[[i]][[j]][[s]][[t]] <- obs_summary
      }
    }
    saveRDS(cell_data, "/Users/jacob/Dropbox/Work/macrodemography/cell_data/cell_data_dg11_22Jun21.RDS")
  }
}




