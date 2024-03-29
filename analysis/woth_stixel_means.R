# Average eBird frequencies over stixels. Currently, a hexagonal grid with
# 285 km spacing and a 5-day window.
# source(paste0("/Users/jacob/Dropbox/Work/Code/macrodemography/data_format/",
# "ebird_raw_to_seasonal.R"))
years <- 2010:2020
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

##### get_cell_data #####
library(auk)

ebd_path <- "/Users/jacob/Desktop/ebd/"  # Location of eBird downloads from Cornell
ebd_month <- "Jun-2021"  # Download version (must be same as sampling file version)
output_path <- "/Users/jacob/Desktop/ebird_files/"  # Location to write output (and intermediates)
sampling_path <- paste0(output_path, ebd_month, "/sampling")  # Location to write intermediate outputted sampling files (prior to subdivision by species)
# dir.create(sampling_path)

states <- ebird_states$state_code[which(ebird_states$country_code == "US")]
L48 <- states[!(states %in% c("US-AK", "US-HI"))]

years <- 2006:2020
spp <- c("woothr")
spp_code <- c("WOTH")



seasons <- c("summer")
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
      zfd <- auk_zerofill(paste0(output_path_2, "/", seasons[s], "_", years[j], 
                                 ".txt"),
                          paste0(sampling_path_2, "/", seasons[s], "_", 
                                 years[j], ".txt"), collapse = T)
      zfd <- zfd[zfd$protocol_type %in% c("Stationary", "Area", "Traveling"), ]
      days <- as.POSIXlt(zfd$observation_date)$yday
      cells6 <- dgGEO_to_SEQNUM(dg6, zfd$longitude, zfd$latitude)$seqnum
      if (seasons[s] == "summer") {
        tgrid <- floor(days/7) - 20
        tgrid_max <- 10
      } 
      for (t in 1:tgrid_max) {
        cell_data[[i]][[j]][[s]][[t]] <- list()
        names(cell_data[[i]][[j]][[s]])[t] <- paste0("tgrid_", t)
        for (g in seq_along(conus_cells)) {
          print(c(i,j,s,t,g))
          cd <- list()
          obs <- zfd[cells6 == conus_cells[g] & tgrid == t, ]
          cd$n <- nrow(obs)
          cd$n_z <- sum(obs$observation_count == "0")
          cd$n_x <- sum(obs$observation_count == "X")
          not_x <- as.integer(obs$observation_count[obs$observation_count != "X"])
          cd$n_value <- sum(not_x > 0)
          cd$mean_positive <- mean(not_x[not_x > 0])
          cell_data[[i]][[j]][[s]][[t]][[g]] <- cd
          names(cell_data[[i]][[j]][[s]][[t]])[g] <- paste0("cell_", conus_cells[g])
        }
      }
    }
    saveRDS(cell_data, "/Users/jacob/Dropbox/Work/macrodemography/cell_data/cell_data_woth_summer.RDS")
  }
}




