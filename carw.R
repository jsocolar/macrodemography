socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
if(socolar.desktop){
  dir.path <- "/Users/JacobSocolar/Dropbox/Work"
}else if(socolar.laptop){
  dir.path <- "/Users/jacob/Dropbox/Work"
}
setwd(dir.path)

roxygen2::roxygenise("Code/macrodemography/erdPackage")

install.packages("Code/macrodemography/erdPackage", 
                 repos = NULL, 
                 type = "source")

library(erdPackage)
library(data.table)

years <- c(2012:2019)
cci_min <- c(-100, 0, .5)
effort_lim <- data.frame(dist_max = c(1, 3, 3),
                         time_min = c(5/60, 5/60, 15/60),
                         time_max = c(.5, 1, 1))

#checklists <- import_checklists("macrodemography/erd/erd.db")

#saveRDS(checklists, "macrodemography/erd/imported_checklists.RDS")
checklists <- readRDS("macrodemography/erd/imported_checklists.RDS")

checklists <- checklists[checklists$latitude > 25 & 
                           checklists$latitude < 50 & 
                           checklists$longitude > -102 &
                           checklists$longitude < -60 &
                           checklists$year >= min(years) &
                           checklists$year <= max(years),]

sp_data <- import_from_erd("carwre", 
                           erd_path = "macrodemography/erd/erd.db",
                           checklists = checklists)



spring_data <- fall_data <- list()
  
for (i in 1:(length(cci_min))) {
  spring_data[[i]] <- fall_data[[i]] <- list()
  spd2 <- sp_data[cci > cci_min[i]]
  for (j in 1:nrow(effort_lim)) {
    spd3 <- spd2[effort_distance_km <= effort_lim$dist_max[j] &
                   effort_hrs >= effort_lim$time_min[j] &
                   effort_hrs <= effort_lim$time_max[j]]
    
    spring_data[[i]][[j]] <- fall_data[[i]][[j]] <- list()
    for (y in seq_along(years)) {
      print(paste(i,j,y))
      spring_data[[i]][[j]][[y]] <- get_grid_data(data = spd3, .year = years[y], 
                                                  tgrid_min = 13, tgrid_max = 16, 
                                                  min_lat = 25, max_lat = 50, min_lon = -102)
      fall_data[[i]][[j]][[y]] <- get_grid_data(data = spd3, .year = years[y], 
                                                  tgrid_min = 35, tgrid_max = 38, 
                                                  min_lat = 25, max_lat = 50, min_lon = -102)
    }
  }
}

saveRDS(spring_data, "macrodemography/carw/spring_data.RDS")
saveRDS(fall_data, "macrodemography/carw/fall_data.RDS")

spring_data <- readRDS("macrodemography/carw/spring_data.RDS")
fall_data <- readRDS("macrodemography/carw/fall_data.RDS")

cci_index <- 2
effort_index <- 2

spring_abun_data <- fall_abun_data <- list()
for (y in seq_along(years)) {
  print(years[y])
  spring_abun_data[[y]] <- get_abun(spring_data[[cci_index]][[effort_index]][[y]], n_rep = 100)
  fall_abun_data[[y]] <- get_abun(fall_data[[cci_index]][[effort_index]][[y]], n_rep = 100)
}

abun_data_bycell <- function(abun_data) {
  cells <- unique(abun_data$cell)
  ad2 <- abun_data[abun_data$n_small > 10, ]
  abun_cols <- as.data.frame(matrix(nrow = length(cells), ncol = ncol(abun_data) - 5))
  names(abun_cols) <- c("average", paste0("rep_", c(1:(ncol(abun_data) - 6))))
  out <- cbind(data.frame(cell = cells), abun_cols)
  for(i in seq_along(cells)) {
    ad3 <- ad2[ad2$cell == cells[i], ]
    if (nrow(ad3) == 4) {
      out[i, 2:102] <- colMeans(ad3[6:106])
    }
  }
  return(out)
}

t1 <- abun_data_bycell(spring_abun_data[[4]])

spring_abun_summary <- fall_abun_summary <- list()
for (y in seq_along(years)) {
  print(y)
  spring_abun_summary[[y]] <- abun_data_bycell(spring_abun_data[[y]])
  fall_abun_summary[[y]] <- abun_data_bycell(fall_abun_data[[y]])
}


cells <- unique(spring_abun_summary[[1]]$cell)

cell_timeseries <- list()
for (i in seq_along(cells)) {
  print(i)
  cell_data <- spring_abun_summary[[1]][spring_abun_summary[[1]]$cell == cells[i], 2:102]
  cell_data <- rbind(cell_data, fall_abun_summary[[1]][fall_abun_summary[[1]]$cell == cells[i], 2:102])
  for(j in 2:length(years)) {
    cell_data <- rbind(cell_data, spring_abun_summary[[j]][spring_abun_summary[[j]]$cell == cells[i], 2:102])
    cell_data <- rbind(cell_data, fall_abun_summary[[j]][fall_abun_summary[[j]]$cell == cells[i], 2:102])
  }
  cell_timeseries[[i]] <- cell_data
}

cell_ratio_series <- list()
for (i in 1:length(cell_timeseries)) {
  lrats <- apply(cell_timeseries[[i]][ ,2:101], 2, function(x){log(stocks::ratios(x))})
  cell_ratio_series[[i]] <- list()
  cell_ratio_series[[i]]$median <- apply(lrats, 1, median)
  cell_ratio_series[[i]]$q10 <- apply(lrats, 1, function(x){quantile(x, .1, na.rm = T)})
  cell_ratio_series[[i]]$q90 <- apply(lrats, 1, function(x){quantile(x, .9, na.rm = T)})
}

dev.off()
for (i in 1:length(cell_ratio_series)) {
  if(sum(is.na(cell_ratio_series[[i]]$median)) < 5){
    plot(cell_ratio_series[[i]]$median, ylim = c(-2,2), pch = 16, main = cells[i], col = c("blue", rep(c("red", "blue"), 7)))
    for(j in 1:length(cell_ratio_series[[i]]$median)){
      lines(x = c(j,j), y = c(cell_ratio_series[[i]]$q10[j], cell_ratio_series[[i]]$q90[j]))
    }
  }
}



sd_holder_prod <- sd_holder_surv <- matrix(nrow = length(cells), ncol = 100)
colnames(sd_holder_prod) <- paste0("prod_sd_rep_", 1:100)
colnames(sd_holder_surv) <- paste0("surv_sd_rep_", 1:100)

cell_lrat_sd <- cbind(data.frame(cell = cells, n_prod = NA, n_surv = NA),
                      as.data.frame(sd_holder_prod), as.data.frame(sd_holder_surv))
for(i in seq_along(cells)){
  cts <- cell_timeseries[[i]]
  lrats <- log(stocks::ratios(cts$average[2:15]))
  cell_lrat_sd$n_prod[i] <- sum(!is.na(lrats[2*c(1:7)]))
  cell_lrat_sd$n_surv[i] <- sum(!is.na(lrats[-1 + 2*c(1:7)]))
  for (j in 1:100) {
    lrats <- log(stocks::ratios(cts[2:15, 1+j]))
    cell_lrat_sd[i, j + 3] <- sd(lrats[2*c(1:7)], na.rm = T)
    cell_lrat_sd[i, j + 103] <- sd(lrats[-1 + 2*c(1:7)], na.rm = T)
  }
}

n_prod_sd_greater <- avg_sd_diff <- vector()
for (i in seq_along(cells)) {
  n_prod_sd_greater[i] <- sum(t(cell_lrat_sd[i, 4:103]) > t(cell_lrat_sd[i, 104:203]), na.rm = T)
  avg_sd_diff[i] <- mean(t(cell_lrat_sd[i, 4:103]) - t(cell_lrat_sd[i, 104:203]), na.rm = T)
}
n_prod_sd_greater
summary(n_prod_sd_greater)

summary(t(cell_lrat_sd[19,4:103]))
summary(t(cell_lrat_sd[19,104:203]))
summary(t(cell_lrat_sd[19,4:103]) - t(cell_lrat_sd[19,104:203]))




plotting_data <- data.frame(cell = cells, n_prod = cell_lrat_sd$n_prod, 
                            n_surv = cell_lrat_sd$n_surv, 
                            n_prod_sd_greater = n_prod_sd_greater,
                            avg_sd_diff = avg_sd_diff,
                            avg_sd_diff_pos = as.integer(avg_sd_diff > 0))

grid6 <- dggridR::dgconstruct(res = 6)

grid <- dggridR::dgcellstogrid(grid6,plotting_data$cell,frame=TRUE,wrapcells=TRUE)
grid  <- merge(grid,plotting_data,by.x="cell")

library(ggplot2)
states <- map_data("state")

n_min <- 6
grid_2 <- grid[grid$n_prod >= n_min & grid$n_surv >= n_min, ]
p <- ggplot() + 
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_2,      aes(x=long, y=lat, group=group, fill=avg_sd_diff), alpha=0.4)    +
  geom_path   (data=grid_2,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradient(low="blue", high="red")
p
