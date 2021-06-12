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
conus_cells <- data.frame(cell = unique(dgGEO_to_SEQNUM(dg6, conus_coords$x, conus_coords$y)[[1]]))
conus_cells$lat <- dgSEQNUM_to_GEO(dg6, conus_cells$cell)$lat_deg
conus_cells$lon <- dgSEQNUM_to_GEO(dg6, conus_cells$cell)$lon_deg


stixel_mean <- function(cell_data_st){
  if (cell_data_st$n == 0) {
    return(NA)
  } else if (cell_data_st$n_x + cell_data_st$n_value == 0) {
    return(0)
  } else if (cell_data_st$n_value == 0) {
    return(weighted.mean(x = c(0,2), w = c(cell_data_st$n_z, cell_data_st$n_x)))
  } else {
    n_z_rescaled <- cell_data_st$n_z * cell_data_st$n_value / (cell_data_st$n_value + cell_data_st$n_x)
    return(weighted.mean(x = c(0, cell_data_st$mean_positive), w = c(n_z_rescaled, cell_data_st$n_value)))
  }
}

mean_of_stixels <- function(sp_yr, season, cells) {
  vec <- vector()
  counter <- 0
  for(i in seq_along(cells)) {
    if (season == "spring") {
      tgrid_max <- 20
      sp_yr_s <- sp_yr[[1]]
    } else if (season == "fall") {
      tgrid_max <- 19
      sp_yr_s <- sp_yr[[2]]
    }
    for(j in 1:tgrid_max){
      sp_yr_s_t <- sp_yr_s[[j]]
      for(k in seq_along(cells)){
        counter <- counter + 1
        vec[counter] <- stixel_mean(sp_yr_s_t[[which(conus_cells == cells[k])]])
      }
    }
  }
  return(mean(vec, na.rm = T))
}

get_means <- function(cell_data, sp_code, cells) {
  the_data <- cell_data[[which(spp_code == sp_code)]]
  spring_raw <- fall_raw <- vector()
  for (i in seq_along(years)) {
    spring_raw[i] <- mean_of_stixels(the_data[[i]], "spring", cells)
    fall_raw[i] <- mean_of_stixels(the_data[[i]], "fall", cells)
  }
  out <- data.frame(means = c(rbind(spring_raw, fall_raw)), 
                              yr = .5*(1:(length(spring_raw) + length(fall_raw))),
                              season = rep(c("spring", "fall"), length(spring_raw)))
  out$ratios <- log(c(stocks::ratios(out$means), NA))
  out$year_ratios <- NA
  for (i in seq_len(nrow(out))) {
    if (i <= (nrow(out) - 1)) {
      out$year_ratios[i] <- out$ratios[i] + out$ratios[i+1]
    }
  }
  out
}

years <- 2010:2019
cell_data <- readRDS("/Users/Jacob/Dropbox/Work/macrodemography/cell_data/cell_data_3Jun21.RDS")


sp_summaries <- vector(mode = "list", length = length(spp_code))
names(sp_summaries) <- spp_code
for (i in seq_along(spp_code)) {
  sp_code <- spp_code[i]
  print(sp_code)
  sp_data <- spp_data[spp_data$species == sp_code, ]
  if (sp_data$lat_divide) {
    out_list <- list()
    split_point <- sp_data$min_lat + (sp_data$max_lat - sp_data$min_lat)/2
    cells_n <- conus_cells$cell[conus_cells$lat > split_point]
    cells_s <- conus_cells$cell[conus_cells$lat < split_point]
    out_list$north <- get_means(cell_data, sp_code, cells_n)
    out_list$south <- get_means(cell_data, sp_code, cells_s)
    sp_summaries[[i]] <- out_list
  } else {
    sp_summaries[[i]] <- get_means(cell_data, sp_code, conus_cells$cell)
  }
}



north_means <- data.frame(means = c(rbind(mean_north_spring, mean_north_fall)), 
                          yr = .5*(1:(length(mean_north_spring) + length(mean_north_fall))),
                          season = rep(c("purple", "orange"), length(mean_north_spring)))
plot(means ~ yr, data = north_means, col = season, pch = 16, ylim = c(0,.5))

south_means <- data.frame(means = c(rbind(mean_south_spring, mean_south_fall)), 
                          yr = .5*(1:(length(mean_south_spring) + length(mean_south_fall))),
                          season = rep(c("blue", "red"), length(mean_south_spring)))
points(means ~ yr, data = south_means, col = season, pch = 16)

north_means$ratios <- log(c(stocks::ratios(north_means$means), NA))
north_mean_spring_ratio <- mean(north_means$ratios[north_means$season == "purple"])
north_means$ratios <- north_means$ratios - north_mean_spring_ratio
south_means$ratios <- log(c(stocks::ratios(south_means$means), NA))
south_mean_spring_ratio <- mean(south_means$ratios[south_means$season == "blue"])
south_means$ratios <- south_means$ratios - south_mean_spring_ratio

plot(ratios ~ yr, data = north_means, col = season, pch = 16, ylim = c(-5,5))
points(ratios ~ yr, data = south_means, col = season, pch = 16)

summary(lm(south_means$ratios[south_means$season == "red"] ~ north_means$ratios[north_means$season == "orange"]))
summary(lm(south_means$ratios[south_means$season == "blue"] ~ north_means$ratios[north_means$season == "purple"]))

