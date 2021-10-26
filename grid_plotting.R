library(dggridR)
library(sf)
library(ggplot2)
library(data.table)
dg6 <- dgconstruct(res=6)
dg11 <- dgconstruct(res=11)

years <- 2010:2019
cell_data <- readRDS("/Users/JacobSocolar/Dropbox/Work/macrodemography/cell_data/cell_data_dg11_22Jun21.RDS")
spp_data <- read.csv("/Users/jacobSocolar/Dropbox/Work/macrodemography/include_by_spp.csv")
spp_code <- spp_data$species

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

conus_cells11 <- data.frame(cell = unique(dgGEO_to_SEQNUM(dg11, conus_coords$x, conus_coords$y)[[1]]))
conus_cells11$lat <- dgSEQNUM_to_GEO(dg11, conus_cells11$cell)$lat_deg
conus_cells11$lon <- dgSEQNUM_to_GEO(dg11, conus_cells11$cell)$lon_deg
conus_cells11$cell6 <- dgGEO_to_SEQNUM(dg6, conus_cells11$lon, conus_cells11$lat)[[1]]

cc11 <- conus_cells11[conus_cells11$cell6 == 700,]


# This hack works for now
p <- ggplot() + geom_sf(data = conus, fill = "gray50", col = "gray50") + scale_x_continuous(limits = c(-92, -67))
poly_colors <- c("goldenrod", "darkmagenta")
pb <- txtProgressBar(min = 1, max = length(conus_cells$cell), initial = 1, style = 3)
for (i in seq_along(conus_cells$cell)) {
  setTxtProgressBar(pb, i)
  grid <- dgcellstogrid(dg6, conus_cells$cell[i])
  meanlat <- mean(grid$lat)
  meanlon <- mean(grid$long)
  if(meanlat > 29 & meanlat < 44 & meanlon > -100){
    p <- p + geom_path(data = grid, aes(x = long, y = lat),
                       col = "black", alpha = .3) +
      geom_polygon(data = grid, aes(x = long, y = lat),
                   col = poly_colors[(meanlat>37)+1], alpha = .3)
  }
}
p

pb <- txtProgressBar(min = 1, max = length(cc11$cell), initial = 1, style = 3)
for (i in seq_along(cc11$cell)) {
  setTxtProgressBar(pb, i)
  grid <- dgcellstogrid(dg11, cc11$cell[i])
  p <- p + geom_path(data = grid, aes(x = long, y = lat),
                     col = "black", alpha = .3) +
    geom_polygon(data = grid, aes(x = long, y = lat),
                 col = "dodgerblue", alpha = .3)
  
}
p
