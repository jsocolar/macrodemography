

socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
if(socolar.desktop){
  dir.path <- "/Users/JacobSocolar/Dropbox/Work"
}else if(socolar.laptop){
  dir.path <- "/Users/jacob/Dropbox/Work"
}

setwd(dir.path)

library(dggridR)
dg6 <- dgconstruct(res=6)
dg11 <- dgconstruct(res=10)

erd_path <- "macrodemography/erd/erd.db"

species <- c("PHVI", "BBWA")
sp_codes <- c("phivir", "babwar")
years <- c(2012:2019)
seasons <- c("spring", "fall")
bands <- c("north", "south")

source("Code/macrodemography/stixel_bootstrap/get_tgrid.R")
source("Code/macrodemography/stixel_bootstrap/get_pixels.R")
source("Code/macrodemography/erd_workflow/get_data_erd.R")
source("Code/macrodemography/erd_workflow/R_package/import_from_erd.R")

for(i in seq_along(species)){
  if (i == -1) { # just to turn this off
    c_arg <- checklists
    o_arg <- obs
  } else {
    c_arg <- o_arg <- NULL
  }
  zf <- import_from_erd(sp_codes[i], erd_path, checklists = c_arg, obs = o_arg)
  zf <- zf[zf$cci > 0, ]
  cell_11_data <- list()
  for(j in seq_along(years)){
    
    cell_11_data[[j]] <- list()
    for(k in seq_along(seasons)) {
      cell_11_data[[j]][[k]] <- list()
      if (k == 1) {
        tgrid_min <- 10
        tgrid_max <- 21
      } else {
        tgrid_min <- 33
        tgrid_max <- 44
      }
      for (L in seq_along(bands)) {
        print(c(i,j, k, L))
        if (L == 1) {
          min_lat <- 29
          max_lat <- 37
          min_lon <- -100
        } else {
          min_lat <- 37
          max_lat <- 44
          min_lon <- -100
        }
        cell_11_data[[j]][[k]][[L]] <- get_data(zf, years[j], 
                                                tgrid_min = tgrid_min, tgrid_max = tgrid_max,
                                                min_lat = min_lat, max_lat = max_lat, min_lon = min_lon, sg=10)
      }
    }
  }
  saveRDS(cell_11_data, paste0("macrodemography/cell_10_data/", species[i], ".RDS"))
}

