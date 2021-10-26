#' get a hexagonal grid over the contiguous US

get_pixels <- function(sg=11) {
  dg6 <- dgconstruct(res=6)
  dg11 <- dgconstruct(res=sg)
  
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
  
  out <- list(conus6 = conus_cells, conus11 = conus_cells11)
  out
}
