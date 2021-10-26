species <- c("BBCU", "BLBW", "BOBO", "BTBW", "CAWA", "CERW", "CSWA",
             "GWWA", "MAWA", "MOWA", "NAWA", "NOWA", "RBGR", "SCTA")
years <- c(2012:2019)
seasons <- c("spring", "fall")

for(i in seq_along(species)){
  cell_11_data <- list()
  for(j in seq_along(years)){
    print(c(i,j))
    cell_11_data[[j]] <- list()
    for(k in seq_along(seasons)){
      cell_11_data[[j]][[k]] <- get_cell11_data(species[i], years[j], seasons[k])
    }
  }
  saveRDS(cell_11_data, paste0("/Users/jacob/Dropbox/Work/macrodemography/cell_11_data/", species[i], ".RDS"))
}
