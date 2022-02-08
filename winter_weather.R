library(ebirdst)
library(raster)
library(stars)
library(rgee)
ee_Initialize(gcs = T)
sp_path <- ebirdst_download(species = "woothr")
abd_lr <- load_raster(sp_path, "abundance", resolution = "lr")

abd_lr_winter <- calc(abd_lr, fun = function(x){mean(x[c(1:8, 49:52)])})

winter_stars <- stars::st_as_stars(abd_lr_winter)
newgrid <- st_as_stars(st_bbox(st_transform(winter_stars, "wgs84")), dx = .2, dy = .2)

winter_reproj <- st_warp(winter_stars, newgrid)

winter_crop <- st_crop(winter_reproj, st_bbox(c(xmin = -100,
                                                ymin = 6,
                                                xmax = -78,
                                                ymax = 23), crs = "wgs84"))
plot(winter_crop)

ee_winter_crop_ID <- sprintf("%s/winter_crop", ee_get_assethome())


ee_winter_crop <- stars_as_ee(x = winter_crop, 
                              assetId = ee_winter_crop_ID, 
                              bucket = "macrodemography_gee",
                              overwrite = TRUE)

