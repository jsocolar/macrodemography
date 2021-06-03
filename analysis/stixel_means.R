# Average eBird frequencies over stixels. Currently, a hexagonal grid with
# 285 km spacing and a 5-day window.
source(paste0("/Users/jacob/Dropbox/Work/Code/macrodemography/data_format/",
              "ebird_raw_to_seasonal.R"))

library(dggridR)
library(sf)
library(ggplot2)
dg6 <- dgconstruct(res=6)

# A hack to determine all grid cells that overlap the continuous US.
# Surely there's a better way...
conus <- spData::us_states
conus_raster <- fasterize::fasterize(conus,
                                     raster::raster(ncol=1000, nrow = 1000, 
                                                    xmn = -125, xmx = -66, 
                                                    ymn =24, ymx = 50))
conus_coords <- raster::as.data.frame(conus_raster, xy = T)
conus_coords <- conus_coords[!is.na(conus_coords$layer), ]
conus_cells <- unique(dgGEO_to_SEQNUM(dg6, conus_coords$x, conus_coords$y)[[1]])
length(conus_cells)

# Need to submit an issue against dggridR for the below:
grid <- dgcellstogrid(dg6, conus_cells)
p <- ggplot() + geom_sf(data = conus) + 
                geom_path(data = grid, aes(x = long, y = lat))
p

# This hack works for now
p <- ggplot() + geom_sf(data = conus, fill = "gray50", col = "gray50")
pb <- txtProgressBar(min = 1, max = length(conus_cells), initial = 1, style = 3)
for (i in seq_along(conus_cells)) {
  setTxtProgressBar(pb, i)
  grid <- dgcellstogrid(dg6, conus_cells[i])
  p <- p + geom_path(data = grid, aes(x = long, y = lat), 
                     col = "black", alpha = .3) +
           geom_polygon(data = grid, aes(x = long, y = lat), 
                     col = "goldenrod", alpha = .3)
    
}
p

##### get_cell_data #####
seasons <- c("spring", "fall")
cell_data <- list()
for (i in seq_along(spp)) {
  cell_data[[i]] <- list()
  names(cell_data)[i] <- spp_code[i]
  output_path_2 <- paste0(output_path, ebd_month, "/", spp_code[i], "/by_year")
  sampling_path_2 <- paste0(output_path, ebd_month, "/sampling/by_year")
  for (j in seq_along(years)) {
    cell_data[[i]][[j]] <- list()
    names(cell_data[[i]])[j] <- paste0("year_", years[j])
    for (s in seq_along(seasons)) {
      cell_data[[i]][[j]][[s]] <- list()
      names(cell_data[[i]][[j]])[s] <- seasons[s]
      zfd <- do.call(cbind, auk_zerofill(paste0(output_path_2, "/", seasons[s], 
                                                "_", years[j], ".txt"),
                                         paste0(sampling_path_2, "/", seasons[s], 
                                                "_", years[j], ".txt")))
      zero_filled_data <- zfd[zfd$sampling_events.protocol_type %in%
                                c("Stationary", "Area", "Traveling"), ]
      days <- as.POSIXlt(zero_filled_data$sampling_events.observation_date)$yday
      
      cells6 <- dgGEO_to_SEQNUM(dg6, zero_filled_data$sampling_events.longitude, 
                                zero_filled_data$sampling_events.latitude)$seqnum
      if (seasons[s] == "spring") {
        tgrid <- round(days/5) - 5
        tgrid_max <- 27
      } else if (seasons[s] == "fall") {
        tgrid <- round(days/5) - 41
        tgrid_max <- 22
      }
      for (t in 1:tgrid_max) {
        cell_data[[i]][[j]][[s]][[t]] <- list()
        names(cell_data[[i]][[j]][[s]])[t] <- paste0("tgrid_", t)
        for (g in seq_along(conus_cells)) {
          print(c(i,j,s,t,g))
          cd <- list()
          obs <- zero_filled_data[cells6 == conus_cells[g] & tgrid == t, ]
          cd$n <- nrow(obs)
          cd$n_z <- sum(obs$observations.observation_count == "0")
          cd$n_x <- sum(obs$observations.observation_count == "X")
          not_x <- as.integer(obs$observations.observation_count[obs$observations.observation_count != "X"])
          cd$n_value <- sum(not_x > 0)
          cd$mean_positive <- mean(not_x[not_x > 0])
          cell_data[[i]][[j]][[s]][[t]][[g]] <- cd
          names(cell_data[[i]][[j]][[s]][[t]])[g] <- paste0("cell_", conus_cells[g])
        }
      }
    }
  }
  saveRDS(cell_data, "/Users/jacob/Dropbox/Work/macrodemography/cell_data/cell_data_3Jun21.RDS")
}




