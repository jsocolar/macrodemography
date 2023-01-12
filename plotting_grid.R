library(sf)
library(ggplot2)
library(magrittr)

socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
if(socolar.desktop){
  dir.path <- "/Users/JacobSocolar/Dropbox/Work/macrodemography"
}else if(socolar.laptop){
  dir.path <- "/Users/jacob/Dropbox/Work/macrodemography"
}

setwd(dir.path)

grid6.1 <- dggridR::dgconstruct(res = 6) |> 
  dggridR::dgearthgrid(frame = FALSE) |> 
  st_as_sf()


NorthAm <- rnaturalearth::ne_states(c("united states of america", "canada")) |>
  st_as_sf() |>
  st_union() |> 
  st_transform(st_crs(grid6.1)) |>
  st_crop(xmin = -110, ymin = 20, xmax = -30, ymax = 55)

grid6 <- grid6.1 %>%
  dplyr::filter(as.logical(rowSums(st_intersects(., NorthAm, sparse = F))))

grid11.1 <- dggridR::dgconstruct(res = 11)
  
  
grid11 <- grid6[st_intersects(grid6, st_point(c(-84,35)), sparse = F),] |>
  st_transform(st_crs("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs")) |>
  st_sample(size = 5000, type = "hexagonal") |>
  st_transform(st_crs(grid6)) |>
  as_Spatial() %>%
  {dggridR::dgGEO_to_SEQNUM(dggs = grid11.1, in_lon_deg = .$coords.x1, in_lat_deg = .$coords.x2)$seqnum} |>
  unique() %>%
  dggridR::dgcellstogrid(grid11.1, ., frame = FALSE) |>
  st_as_sf()

ggplot(grid11) + geom_sf()

p <- ggplot(NorthAm) + geom_sf(fill = "gray 50", color = "gray50") +
  geom_sf(data=grid6, fill = NA, color = "goldenrod") + 
  xlim(c(-110, -67)) + ylim(c(25, 49)) +
  geom_sf(data = grid11, color = "#7393B3", fill = "gray50")
p

p <- ggplot(NorthAm) + geom_sf(fill = "gray 50", color = "gray50") +
  geom_sf(data=grid6, fill = NA, color = "goldenrod") + 
  xlim(c(-90, -80)) + ylim(c(30, 38)) +
  geom_sf(data = grid11, color = "dodgerblue", fill = "gray50")
p
