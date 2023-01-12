##### Set up the working directory #####
socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
if(socolar.desktop){
  dir.path <- "/Users/JacobSocolar/Dropbox/Work"
}else if(socolar.laptop){
  dir.path <- "/Users/jacob/Dropbox/Work"
}
setwd(dir.path)

# Change to appropriate remotes::install_github
install.packages("Code/macrodemography/erdPackage", 
                 repos = NULL, 
                 type = "source")
library(erdPackage)
library(data.table)

##### Import ERD #####
years <- c(2006:2019)
cci_min <- c(-100, 0, .5)
effort_lim <- data.frame(dist_max = c(1, 3, 3),
                         time_min = c(5/60, 5/60, 15/60),
                         time_max = c(.5, 1, 1))

#checklists <- import_checklists("macrodemography/erd/erd.db")

#saveRDS(checklists, "macrodemography/erd/imported_checklists.RDS")
checklists <- readRDS("macrodemography/erd/imported_checklists.RDS")

checklists <- checklists[checklists$latitude > 23 & 
                           checklists$latitude < 50 & 
                           checklists$longitude > -128 &
                           checklists$longitude < -60 &
                           checklists$year >= min(years) &
                           checklists$year <= max(years) &
                           checklists$ELEV_30M_MEDIAN < 2000 &
                           ((checklists$ELEV_30M_MEDIAN < 1500) | (checklists$latitude < 42)),]

grid6 <- dggridR::dgconstruct(res = 6)
cells_all <- unique(dggridR::dgGEO_to_SEQNUM(grid6, checklists$longitude, checklists$latitude)[[1]])

species <- data.frame(four = c("bcch", "bhnu", "cach", "carw", "noca", "piwo", "tuti", 
                               "canw", "cacw", "wren", "calt", "cant", "cath", "cbth", "lbth", 
                               "btgn", "verd", "pyrr", "oati", "juti", "casj", "wosj"),
                      six = c("bkcchi", "bnhnut", "carchi", "carwre", "norcar", "pilwoo", "tuftit",
                              "canwre", "cacwre", "wrenti", "caltow", "cantow", "calthr", "cubthr", "lobthr",
                              "bktgna", "verdin", "pyrrhu", "oaktit", "juntit", "cowscj", "wooscj"))

# for(sp in 1:nrow(species)){
#   sp_data <- import_from_erd(species$six[sp],
#                              erd_path = "macrodemography/erd/erd.db",
#                              checklists = checklists)
#   dir.create(paste0("macrodemography/residents/", species$four[sp]))
#   spring_data <- fall_data <- list()
#   for (i in 1:(length(cci_min))) {
#     spring_data[[i]] <- fall_data[[i]] <- list()
#     spd2 <- sp_data[cci > cci_min[i]]
#     for (j in 1:nrow(effort_lim)) {
#       spd3 <- spd2[effort_distance_km <= effort_lim$dist_max[j] &
#                      effort_hrs >= effort_lim$time_min[j] &
#                      effort_hrs <= effort_lim$time_max[j]]
# 
#       spring_data[[i]][[j]] <- fall_data[[i]][[j]] <- list()
#       for (y in seq_along(years)) {
#         print(paste(i,j,y))
#         if (i == 2 & j == 2) {
#           spring_data[[i]][[j]][[y]] <- get_grid_data(data = spd3, .year = years[y],
#                                                       tgrid_min = 13, tgrid_max = 16, time_window = "full",
#                                                       min_lat = 23, max_lat = 50, min_lon = -128)
#           fall_data[[i]][[j]][[y]] <- get_grid_data(data = spd3, .year = years[y],
#                                                     tgrid_min = 40, tgrid_max = 43, time_window = "full",
#                                                     min_lat = 23, max_lat = 50, min_lon = -128)
#         }
#       }
#     }
#   }
# 
#   saveRDS(spring_data, paste0("macrodemography/residents/", species$four[sp], "/spring_data_fullwindow.RDS"))
#   saveRDS(fall_data, paste0("macrodemography/residents/", species$four[sp], "/fall_data_fullwindow.RDS"))
# 
#   cci_index <- 2
#   effort_index <- 2
# 
#   spring_abun_data <- fall_abun_data <- list()
#   for (y in seq_along(years)) {
#     print(years[y])
#     spring_abun_data[[y]] <- get_abun(spring_data[[cci_index]][[effort_index]][[y]], n_rep = 100)
#     fall_abun_data[[y]] <- get_abun(fall_data[[cci_index]][[effort_index]][[y]], n_rep = 100)
#   }
# 
# 
#   saveRDS(spring_abun_data, paste0("macrodemography/residents/", species$four[sp], "/spring_abun_data_fullwindow.RDS"))
#   saveRDS(fall_abun_data,  paste0("macrodemography/residents/", species$four[sp], "/fall_abun_data_fullwindow.RDS"))
# }


##### Get demographic indices #####
spring_abun_data <- readRDS("macrodemography/residents/noca/spring_abun_data_fullwindow.RDS")
fall_abun_data <- readRDS("macrodemography/residents/noca/fall_abun_data_fullwindow.RDS")

# Function to extract cell-specific abundance data
abun_data_bycell <- function(abun_data, n_small_min = 10) {
  cells <- unique(abun_data$cell)
  ad2 <- abun_data[abun_data$n_small >= n_small_min, ]
  abun_cols <- as.data.frame(matrix(nrow = length(cells), ncol = ncol(abun_data) - 5))
  names(abun_cols) <- c("average", paste0("rep_", c(1:(ncol(abun_data) - 6))))
  out <- cbind(data.frame(cell = cells), abun_cols)
  for(i in seq_along(cells)) {
    ad3 <- ad2[ad2$cell == cells[i], ]
    if (nrow(ad3) == 1) {
      out[i, 2:102] <- colMeans(ad3[6:106])
    }
  }
  return(out)
}

# Extract cell-specific abundance data
spring_abun_summary <- fall_abun_summary <- list()
for (y in seq_along(years)) {
  print(y)
  spring_abun_summary[[y]] <- abun_data_bycell(spring_abun_data[[y]])
  fall_abun_summary[[y]] <- abun_data_bycell(fall_abun_data[[y]])
}

# grid6 <- dggridR::dgconstruct(res = 6)
# cells_all <- unique(dggridR::dgGEO_to_SEQNUM(grid6, checklists$longitude, checklists$latitude)[[1]])

# Mesh spring and fall data into a timeseries
get_cell_timeseries <- function (cells_all, 
                                 spring_abun_summary, fall_abun_summary) {
  cells <- unique(spring_abun_summary[[1]]$cell)
  cell_timeseries <- list()
  for (i in seq_along(cells_all)) {
    if (cells_all[i] %in% cells) {
      print(i)
      cell_data <- spring_abun_summary[[1]][spring_abun_summary[[1]]$cell == cells_all[i], 2:102]
      cell_data <- rbind(cell_data, fall_abun_summary[[1]][fall_abun_summary[[1]]$cell == cells_all[i], 2:102])
      for(j in 2:length(years)) {
        cell_data <- rbind(cell_data, spring_abun_summary[[j]][spring_abun_summary[[j]]$cell == cells_all[i], 2:102])
        cell_data <- rbind(cell_data, fall_abun_summary[[j]][fall_abun_summary[[j]]$cell == cells_all[i], 2:102])
      }
      cell_timeseries[[i]] <- cell_data
    } else {
      cell_timeseries[[i]] <- NA
    }
    
  }
  return(cell_timeseries)
}
cell_timeseries <- get_cell_timeseries(cells_all, 
                                       spring_abun_summary,
                                       fall_abun_summary)

# Take the log-ratios along the timeseries to get a ratio series
# cell_ratio_series is a summary of the bootstrap uncertainty, and 
# cell_ratio_series_full gives all bootstrap replicates.
cells <- unique(spring_abun_summary[[1]]$cell)
cell_ratio_series <- cell_ratio_series_full <- list()
for (i in 1:length(cell_timeseries)) {
  if (cells_all[i] %in% cells) {
    lrats <- apply(cell_timeseries[[i]][ ,2:101], 2, function(x){log(stocks::ratios(x))})
    cell_ratio_series[[i]] <- list()
    cell_ratio_series[[i]]$median <- apply(lrats, 1, median)
    cell_ratio_series[[i]]$avg <- apply(lrats, 1, mean)
    
    cell_ratio_series[[i]]$q10 <- apply(lrats, 1, function(x){quantile(x, .1, na.rm = T)})
    cell_ratio_series[[i]]$q90 <- apply(lrats, 1, function(x){quantile(x, .9, na.rm = T)})
    
    cell_ratio_series_full[[i]] <- lrats
  } else {
    cell_ratio_series[[i]] <- cell_ratio_series_full[[i]] <- NA
  }
}

# Plot the cell ratio series
dev.off()
for (i in 1:length(cell_ratio_series)) {
  if (cells_all[i] %in% cells) {
    if(sum(is.na(cell_ratio_series[[i]]$median)) < 20){ # only do the plot for cells with at least a couple of years
      plot(cell_ratio_series[[i]]$median, ylim = c(-2,2), pch = 16, main = cells_all[i], col = c("blue", rep(c("red", "blue"), 7)))
      for(j in 1:length(cell_ratio_series[[i]]$median)){
        lines(x = c(j,j), y = c(cell_ratio_series[[i]]$q10[j], cell_ratio_series[[i]]$q90[j]))
      }
    }
  }
}


##### Analyze the timeseries #####
# function to determine which cell-years are analyzeable
use_cell_years <- function (ratio_series, uncertainty_high_grade = 1, 
                            inf_exclude = F, element_1_exclude = T,
                            n_min_prod = 5, n_min_surv = 5,
                            n_min_full = 5) {
  if(!identical(ratio_series, NA)){
    if(any(is.infinite(ratio_series$median))) {
      contains_inf <- TRUE
      for(i in 1:length(ratio_series)){
        ratio_series[[i]][is.infinite(ratio_series[[i]])] <- NA
      }
    } else {contains_inf <- FALSE}
    
    if(any(!is.na(ratio_series$median[1 + 2*c(0:13)])) &
       any(!is.na(ratio_series$median[2*c(1:13)]))){
      prod_dif <- max(ratio_series$median[1+2*c(1:13)], na.rm = T) -
        min(ratio_series$median[1+2*c(1:13)], na.rm = T)
      surv_dif <- max(ratio_series$median[2*c(1:13)], na.rm = T) -
        min(ratio_series$median[2*c(1:13)], na.rm = T)
      use <- rep(1, length(ratio_series$median))
      if(element_1_exclude){
        use[1] <- 0
      }else if(is.na(ratio_series$median[1])){
        use[1] <- 0
      }else{
        uncertainty <- ratio_series$q90[1] - ratio_series$q10[1]
        if(uncertainty > (uncertainty_high_grade * prod_dif)) {
          use[1] <- 0
        }
      }
      
      for (i in 1:13) {
        p_i <- 1 + 2*i
        s_i <- 2*i
        
        if(is.na(ratio_series$median[p_i])){
          use[p_i] <- 0
        } else {
          uncertainty <- ratio_series$q90[p_i] - ratio_series$q10[p_i]
          if(uncertainty > (uncertainty_high_grade * prod_dif)) {
            use[p_i] <- 0
          }
        }
        
        if(is.na(ratio_series$median[s_i])){
          use[s_i] <- 0
        } else {
          uncertainty <- ratio_series$q90[s_i] - ratio_series$q10[s_i]
          if(uncertainty > (uncertainty_high_grade * surv_dif)) {
            use[s_i] <- 0
          }
        }
      }
      
      if(inf_exclude & contains_inf) {
        use <- rep(0, length(ratio_series$median))
      }
      use <- as.logical(use)
    } else {
      use <- rep(FALSE, length(ratio_series$median))
    }
    
    if(sum(use[1+2*c(0:13)]) < n_min_prod) {
      use <- rep(FALSE, length(ratio_series$median))
    } else if(sum(use[2*c(1:13)]) < n_min_surv) {
      use <- rep(FALSE, length(ratio_series$median))
    } else if ((sum(use[1 + 2*c(1:13)] * use[2*c(1:13)])) < n_min_full) {
      use <- rep(FALSE, length(ratio_series$median))
    }
    
    return(use)
  }
}

# Get the standard deviations of the survival and productivity values
library(rstanarm)
sd_holder_prod <- sd_holder_surv <- matrix(nrow = length(cells_all), ncol = 100)
colnames(sd_holder_prod) <- paste0("prod_sd_rep_", 1:100)
colnames(sd_holder_surv) <- paste0("surv_sd_rep_", 1:100)

cell_lrat_sd <- cbind(data.frame(cell = cells_all, n_prod = NA, n_surv = NA),
                      as.data.frame(sd_holder_prod), as.data.frame(sd_holder_surv))
var_p <- rep(NA, length(cells_all))

for(i in seq_along(cells_all)){
  print(i)
  if(!identical(cell_ratio_series[[i]], NA)){
    lrats_avg <- cell_ratio_series[[i]]$avg
    lrats_avg[!use_cell_years(cell_ratio_series[[i]], inf_exclude=T)] <- NA
    if (cells_all[i] %in% cells) {
      prod_rats <- lrats_avg[1 + 2*c(1:13)]
      surv_rats <- lrats_avg[2*c(1:13)]
      cell_lrat_sd$n_prod[i] <- sum(!is.na(prod_rats))
      cell_lrat_sd$n_surv[i] <- sum(!is.na(surv_rats))
      var_p_i <- vector()
      for (j in 1:100) {
        lrats <- cell_ratio_series_full[[i]][ , j]
        lrats[!use_cell_years(cell_ratio_series[[i]], inf_exclude=T)] <- NA
        prod_rats <- lrats[1 + 2*c(1:13)]
        surv_rats <- lrats[2*c(1:13)]
        cell_lrat_sd[i, j + 3] <- sd(prod_rats, na.rm = T)
        cell_lrat_sd[i, j + 103] <- sd(surv_rats, na.rm = T)
        # if((sum(!is.na(prod_rats)) > 1) & (sum(!is.na(surv_rats)) > 1)) {
        #   prod_mod <- stan_glm(prod_rats ~ 1, prior_intercept = normal(0,1,autoscale = F), refresh = 0)
        #   surv_mod <- stan_glm(surv_rats ~ 1, prior_intercept = normal(0,1,autoscale = F), refresh = 0)
        #   prod_var <- (as.data.frame(prod_mod)$sigma)^2
        #   surv_var <- (as.data.frame(surv_mod)$sigma)^2
        #   var_p_i <- c(var_p_i, surv_var > prod_var)
        # }
      }
      # var_p[i] <- mean(var_p_i)
    }
  }
}

# saveRDS(var_p, "macrodemography/var_p_CARW.RDS")
var_p <- readRDS("macrodemography/var_p_CARW.RDS")

n_prod_sd_greater <- avg_sd_diff <- vector()
for (i in seq_along(cells_all)) {
  if (cells_all[i] %in% cells) {
    n_prod_sd_greater[i] <- sum(t(cell_lrat_sd[i, 4:103]) >= t(cell_lrat_sd[i, 104:203]), na.rm = T)
    avg_sd_diff[i] <- mean(t(cell_lrat_sd[i, 4:103]) - t(cell_lrat_sd[i, 104:203]), na.rm = T)
  } else {
    n_prod_sd_greater[i] <- avg_sd_diff[i] <- NA
  }
}
summary(n_prod_sd_greater)

plotting_data <- data.frame(cell = cells_all, n_prod = cell_lrat_sd$n_prod, 
                            n_surv = cell_lrat_sd$n_surv, 
                            n_prod_sd_greater = n_prod_sd_greater,
                            avg_sd_diff = avg_sd_diff,
                            avg_sd_diff_pos = as.integer(avg_sd_diff > 0),
                            p_surv_var_larger = var_p)
plotting_data$n_prod_sd_greater[is.na(plotting_data$avg_sd_diff)] <- NA

grid6 <- dggridR::dgconstruct(res = 6)

grid <- dggridR::dgcellstogrid(grid6,plotting_data$cell,frame=TRUE,wrapcells=TRUE)
grid  <- merge(grid,plotting_data,by.x="cell")

library(ggplot2)
states <- map_data("state")
western_states <- c("arizona", "california", "colorado", "idaho", "montana",
                    "nevada", "new mexico", "oregon", "utah", "washington",
                    "wyoming")
states <- states[!states$region %in% western_states,]

n_min <- 5

grid_2 <- grid[!is.na(grid$n_prod), ]
grid_2 <- grid_2[grid_2$n_prod >= n_min & grid_2$n_surv >= n_min, ]
names(grid_2)[names(grid_2) == "p_surv_var_larger"] <- "p(survival)"
# p <- ggplot() + 
#   geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
#   geom_polygon(data=grid_2,      aes(x=long, y=lat, group=group, fill=n_prod_sd_greater), alpha=0.4)    +
#   geom_path   (data=grid_2,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
#   viridis::scale_fill_viridis()
# p

p <- ggplot() + 
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_2,      aes(x=long, y=lat, group=group, fill=`p(survival)`), alpha=0.7)    +
  geom_path   (data=grid_2,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  viridis::scale_fill_viridis(limits = c(0.05,0.95)) +
  xlim(c(-107, -67))
p


##### Extract weather data from Earth Engine #####
library(dggridR)
library(sf)
library(rgee)
source("/Users/jacob/Dropbox/Work/Code/macrodemography/analysis/functions/daymet_extract.R")
ee_Initialize()

dg <- dgconstruct(res = 6)
# 
# janfeb_temp <- julaug_temp <- junjul_precip <- novfeb_precip <-
#    decmar_swe <- list()
# for (j in 215:length(cells_all)){ #seq_along(cells_all)) {
#   print(j)
#   janfeb_temp[[j]] <- julaug_temp[[j]] <-
#     novfeb_precip[[j]] <- junjul_precip[[j]] <- vector()
#   decmar_swe[[j]] <- vector()
#   for (i in 1:14) {
#     print(i)
#     mindate = paste0(i+2006, "-01-01")
#     maxdate = paste0(i+2006, "-02-28")
#     janfeb_temp[[j]][i] <- daymet_extract(cell = cells_all[j], variable = "tmax", mindate = mindate, maxdate = maxdate)
# 
#     mindate = paste0(i+2005, "-07-01")
#     maxdate = paste0(i+2005, "-08-31")
#     julaug_temp[[j]][i] <- daymet_extract(cell = cells_all[j], variable = "tmax", mindate = mindate, maxdate = maxdate)
# 
#     mindate = paste0(i+2005, "-11-01")
#     maxdate = paste0(i+2006, "-02-28")
#     novfeb_precip[[j]][i] <- daymet_extract(cell = cells_all[j], variable = "prcp", mindate = mindate, maxdate = maxdate)
# 
#     mindate = paste0(i+2005, "-06-01")
#     maxdate = paste0(i+2005, "-07-31")
#     junjul_precip[[j]][i] <- daymet_extract(cell = cells_all[j], variable = "prcp", mindate = mindate, maxdate = maxdate)
# 
#     mindate = paste0(i+2005, "-12-01")
#     maxdate = paste0(i+2006, "-03-15")
#     decmar_swe[[j]][i] <- daymet_extract(cell = cells_all[j], variable = "swe", mindate = mindate, maxdate = maxdate)
#   }
# }
# 
# saveRDS(janfeb_temp, "macrodemography/weather/janfeb_temp.RDS")
# saveRDS(julaug_temp, "macrodemography/weather/julaug_temp.RDS")
# saveRDS(novfeb_precip, "macrodemography/weather/novfeb_precip.RDS")
# saveRDS(junjul_precip, "macrodemography/weather/junjul_precip.RDS")
# saveRDS(decmar_swe, "macrodemography/weather/decmar_swe.RDS")


##### Analyze and summarize weather relationships #####
janfeb_temp <- readRDS("macrodemography/weather/janfeb_temp.RDS")
julaug_temp <- readRDS("macrodemography/weather/julaug_temp.RDS")
novfeb_precip <- readRDS("macrodemography/weather/novfeb_precip.RDS")
junjul_precip <- readRDS("macrodemography/weather/junjul_precip.RDS")
decmar_swe <- readRDS("macrodemography/weather/decmar_swe.RDS")

plotting_data2 <- data.frame(cell = cells_all, 
                             surv_length = NA,
                             prod_length = NA,
                             full_length = NA,
                             
                             janfeb_mean_slope = 0, 
                             janfeb_mean_r = 0,
                             janfeb_neg_and_p.1 = 0,
                             janfeb_pos_and_p.1 = 0,
                             janfeb_r_0 = 0, 
                             janfeb_r_.5 = 0,
                             janfeb_r_neg.5 = 0,
                             janfeb_r_.3 = 0,
                             janfeb_r_neg.3 = 0,
                             janfeb_mean_bayes = NA,
                             janfeb_median_bayes = NA,
                             janfeb_var_bayes = NA,
                             janfeb_skew_bayes = NA,
                             janfeb_kurt_bayes = NA,
                             janfeb_p_bayes = NA,
                             
                             julaug_mean_slope = 0, 
                             julaug_mean_r = 0,
                             julaug_neg_and_p.1 = 0,
                             julaug_pos_and_p.1 = 0,
                             julaug_r_0 = 0, 
                             julaug_r_.5 = 0,
                             julaug_r_neg.5 = 0,
                             julaug_r_.3 = 0,
                             julaug_r_neg.3 = 0,
                             
                             full_janfeb_mean_slope = 0, 
                             full_janfeb_mean_r = 0,
                             full_janfeb_neg_and_p.1 = 0,
                             full_janfeb_pos_and_p.1 = 0,
                             full_janfeb_r_0 = 0, 
                             full_janfeb_r_.5 = 0,
                             full_janfeb_r_neg.5 = 0,
                             full_janfeb_r_.3 = 0,
                             full_janfeb_r_neg.3 = 0
                             )

# Function to get p value from lm object
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

for (i in seq_along(cells_all)) {
  print(i)
  if (!identical(cell_ratio_series[[i]], NA)) {
    crsi <- cell_ratio_series[[i]]
    
    plotting_data2$surv_length[i] <- sum(is.finite(crsi$median[2*c(1:13)]))
    plotting_data2$prod_length[i] <- sum(is.finite(crsi$median[1 + 2*c(0:13)]))
    plotting_data2$full_length[i] <- sum(is.finite(crsi$median[1 + 2*c(1:13)] + crsi$median[2*c(1:13)]))
    
    janfeb_temp_slopes <- vector()
    
    
    
    for (j in 1:100) {
      crsi_surv <- crsi_prod <- crsi_full <- cell_ratio_series_full[[i]][,j]
      
      crsi_surv[!use_cell_years(cell_ratio_series[[i]], 
                                inf_exclude = F,
                                n_min_prod = 0,
                                n_min_surv = 5,
                                n_min_full = 0)] <- NA
      crsi_surv <- crsi_surv[2*c(1:13)]
      
      if (sum(!is.na(crsi_surv)) > 2) {
        themod <- lm(crsi_surv ~ unlist(janfeb_temp[[i]][1:13]))
        b <- themod$coefficients[2]
        r <- ifelse(themod$coefficients[2] > 0, 
                    sqrt(summary(themod)$r.squared),
                    -sqrt(summary(themod)$r.squared))
        
        p_pos.1 <- (lmp(themod) + (themod$coefficients[2] < 0)) < .1
        p_neg.1 <- (lmp(themod) + (themod$coefficients[2] > 0)) < .1

        plotting_data2$janfeb_mean_slope[i] <- plotting_data2$janfeb_mean_slope[i] + b/100
        plotting_data2$janfeb_mean_r[i] <- plotting_data2$janfeb_mean_r[i] + r/100
        plotting_data2$janfeb_neg_and_p.1[i] <- plotting_data2$janfeb_neg_and_p.1[i] + p_neg.1/100
        plotting_data2$janfeb_pos_and_p.1[i] <- plotting_data2$janfeb_pos_and_p.1[i] + p_pos.1/100
        plotting_data2$janfeb_r_0[i] <- plotting_data2$janfeb_r_0[i] + (r > 0)/100
        plotting_data2$janfeb_r_.5[i] <- plotting_data2$janfeb_r_.5[i] + (r > .5)/100
        plotting_data2$janfeb_r_neg.5[i] <- plotting_data2$janfeb_r_neg.5[i] + (r < -.5)/100
        plotting_data2$janfeb_r_.3[i] <- plotting_data2$janfeb_r_.3[i] + (r > .3)/100
        plotting_data2$janfeb_r_neg.3[i] <- plotting_data2$janfeb_r_neg.3[i] + (r < -.3)/100
        
        
        my_df <- data.frame(surv = crsi_surv, janfeb_temp = unlist(janfeb_temp[[i]][1:13]))
        themod2 <- stan_glm(surv ~ janfeb_temp, data = my_df,
                            prior = normal(0, .2, autoscale = F), refresh = 0)
        janfeb_temp_slopes <- c(janfeb_temp_slopes, as.data.frame(themod2)$janfeb_temp)
        
      } else {
        plotting_data2$janfeb_mean_slope[i] <- NA
        plotting_data2$janfeb_neg_and_p.1[i] <- NA
        plotting_data2$janfeb_pos_and_p.1[i] <- NA
        plotting_data2$janfeb_r_0[i] <- NA
        plotting_data2$janfeb_r_.5[i] <- NA
        plotting_data2$janfeb_r_neg.5[i] <- NA
        plotting_data2$janfeb_r_.3[i] <- NA
        plotting_data2$janfeb_r_neg.3[i] <- NA
      }

      crsi_prod[!use_cell_years(cell_ratio_series[[i]], 
                                inf_exclude = F,
                                n_min_prod = 5,
                                n_min_surv = 0,
                                n_min_full = 0)] <- NA
      crsi_prod <- crsi_prod[1 + 2*c(0:13)]
      
      if (sum(!is.na(crsi_prod)) > 2) {
        themod <- lm(crsi_prod ~ unlist(julaug_temp[[i]]))
        b <- themod$coefficients[2]
        r <- ifelse(themod$coefficients[2] > 0, 
                    sqrt(summary(themod)$r.squared),
                    -sqrt(summary(themod)$r.squared))
        
        p_pos.1 <- (lmp(themod) + (themod$coefficients[2] < 0)) < .1
        p_neg.1 <- (lmp(themod) + (themod$coefficients[2] > 0)) < .1
        
        plotting_data2$julaug_mean_slope[i] <- plotting_data2$julaug_mean_slope[i] + b/100
        plotting_data2$julaug_mean_r[i] <- plotting_data2$julaug_mean_r[i] + r/100
        plotting_data2$julaug_neg_and_p.1[i] <- plotting_data2$julaug_neg_and_p.1[i] + p_neg.1/100
        plotting_data2$julaug_pos_and_p.1[i] <- plotting_data2$julaug_pos_and_p.1[i] + p_pos.1/100
        plotting_data2$julaug_r_0[i] <- plotting_data2$julaug_r_0[i] + (r > 0)/100
        plotting_data2$julaug_r_.5[i] <- plotting_data2$julaug_r_.5[i] + (r > .5)/100
        plotting_data2$julaug_r_neg.5[i] <- plotting_data2$julaug_r_neg.5[i] + (r < -.5)/100
        plotting_data2$julaug_r_.3[i] <- plotting_data2$julaug_r_.3[i] + (r > .3)/100
        plotting_data2$julaug_r_neg.3[i] <- plotting_data2$julaug_r_neg.3[i] + (r < -.3)/100
      } else {
        plotting_data2$julaug_mean_slope[i] <- NA
        plotting_data2$julaug_neg_and_p.1[i] <- NA
        plotting_data2$julaug_pos_and_p.1[i] <- NA
        plotting_data2$julaug_r_0[i] <- NA
        plotting_data2$julaug_r_.5[i] <- NA
        plotting_data2$julaug_r_neg.5[i] <- NA
        plotting_data2$julaug_r_.3[i] <- NA
        plotting_data2$julaug_r_neg.3[i] <- NA
      }
      
      
      
      crsi_full[!use_cell_years(cell_ratio_series[[i]], 
                                inf_exclude = F,
                                n_min_prod = 0,
                                n_min_surv = 0,
                                n_min_full = 5)] <- NA
      crsi_full <- crsi_full[1 + 2*c(1:13)] + crsi_full[2*c(1:13)]
      
      if (sum(!is.na(crsi_full)) > 2) {
        themod <- lm(crsi_full ~ unlist(janfeb_temp[[i]][1:13]))
        b <- themod$coefficients[2]
        r <- ifelse(themod$coefficients[2] > 0, 
                    sqrt(summary(themod)$r.squared),
                    -sqrt(summary(themod)$r.squared))
        
        p_pos.1 <- (lmp(themod) + (themod$coefficients[2] < 0)) < .1
        p_neg.1 <- (lmp(themod) + (themod$coefficients[2] > 0)) < .1
        
        plotting_data2$full_janfeb_mean_slope[i] <- plotting_data2$full_janfeb_mean_slope[i] + b/100
        plotting_data2$full_janfeb_mean_r[i] <- plotting_data2$full_janfeb_mean_r[i] + r/100
        plotting_data2$full_janfeb_neg_and_p.1[i] <- plotting_data2$full_janfeb_neg_and_p.1[i] + p_neg.1/100
        plotting_data2$full_janfeb_pos_and_p.1[i] <- plotting_data2$full_janfeb_pos_and_p.1[i] + p_pos.1/100
        plotting_data2$full_janfeb_r_0[i] <- plotting_data2$full_janfeb_r_0[i] + (r > 0)/100
        plotting_data2$full_janfeb_r_.5[i] <- plotting_data2$full_janfeb_r_.5[i] + (r > .5)/100
        plotting_data2$full_janfeb_r_neg.5[i] <- plotting_data2$full_janfeb_r_neg.5[i] + (r < -.5)/100
        plotting_data2$full_janfeb_r_.3[i] <- plotting_data2$full_janfeb_r_.3[i] + (r > .3)/100
        plotting_data2$full_janfeb_r_neg.3[i] <- plotting_data2$full_janfeb_r_neg.3[i] + (r < -.3)/100
      } else {
        plotting_data2$full_janfeb_mean_slope[i] <- NA
        plotting_data2$full_janfeb_neg_and_p.1[i] <- NA
        plotting_data2$full_janfeb_pos_and_p.1[i] <- NA
        plotting_data2$full_janfeb_r_0[i] <- NA
        plotting_data2$full_janfeb_r_.5[i] <- NA
        plotting_data2$full_janfeb_r_neg.5[i] <- NA
        plotting_data2$full_janfeb_r_.3[i] <- NA
        plotting_data2$full_janfeb_r_neg.3[i] <- NA
      }
      
    }
    plotting_data2$janfeb_mean_bayes[i] <- mean(janfeb_temp_slopes)
    plotting_data2$janfeb_median_bayes[i] <- median(janfeb_temp_slopes)
    plotting_data2$janfeb_var_bayes[i] <- var(janfeb_temp_slopes)
    plotting_data2$janfeb_skew_bayes[i] <- moments::skewness(janfeb_temp_slopes)
    plotting_data2$janfeb_kurt_bayes[i] <- moments::kurtosis(janfeb_temp_slopes)
    plotting_data2$janfeb_p_bayes[i] <- mean(janfeb_temp_slopes > 0)
  } else {
    plotting_data2$janfeb_mean_slope[i] <- NA
    plotting_data2$janfeb_neg_and_p.1[i] <- NA
    plotting_data2$janfeb_pos_and_p.1[i] <- NA
    plotting_data2$janfeb_r_0[i] <- NA
    plotting_data2$janfeb_r_.5[i] <- NA
    plotting_data2$janfeb_r_neg.5[i] <- NA
    plotting_data2$janfeb_r_.3[i] <- NA
    plotting_data2$janfeb_r_neg.3[i] <- NA
    
    plotting_data2$julaug_mean_slope[i] <- NA
    plotting_data2$julaug_neg_and_p.1[i] <- NA
    plotting_data2$julaug_pos_and_p.1[i] <- NA
    plotting_data2$julaug_r_0[i] <- NA
    plotting_data2$julaug_r_.5[i] <- NA
    plotting_data2$julaug_r_neg.5[i] <- NA
    plotting_data2$julaug_r_.3[i] <- NA
    plotting_data2$julaug_r_neg.3[i] <- NA
    
    plotting_data2$full_janfeb_mean_slope[i] <- NA
    plotting_data2$full_janfeb_neg_and_p.1[i] <- NA
    plotting_data2$full_janfeb_pos_and_p.1[i] <- NA
    plotting_data2$full_janfeb_r_0[i] <- NA
    plotting_data2$full_janfeb_r_.5[i] <- NA
    plotting_data2$full_janfeb_r_neg.5[i] <- NA
    plotting_data2$full_janfeb_r_.3[i] <- NA
    plotting_data2$full_janfeb_r_neg.3[i] <- NA
  }
}

# saveRDS(plotting_data2, "macrodemography/plotting_data2_carw.RDS")

##### Plotting #####
plotting_data2 <- readRDS("macrodemography/plotting_data2_carw.RDS")

names(plotting_data2)[names(plotting_data2) == "janfeb_mean_bayes"] <- "mean slope"
names(plotting_data2)[names(plotting_data2) == "janfeb_p_bayes"] <- "p(positive slope)"

grid <- dggridR::dgcellstogrid(grid6,plotting_data$cell,frame=TRUE,wrapcells=TRUE)
grid_3  <- merge(grid,plotting_data2,by.x="cell")

library(tidyverse)
plotting_fun <- function(grid_data, variable, fill_lim = NULL) {
  v <- variable
  p <- ggplot() + 
    geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
    geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = !!v), alpha = 0.4)   +
    geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
    viridis::scale_fill_viridis(limits = fill_lim)
  print(p)
}

fill_lim <- NULL
for (i in 2:15) {
  if(i == 7) {fill_lim <- c(0, 1.001)}
  if(i == 14) {fill_lim <- NULL}
  if(i == 7) {fill_lim <- c(0, 1.001)}
  plotting_fun(grid_3[!is.na(grid_3$janfeb_mean_slope), ], sym(names(plotting_data2)[i]), fill_lim = fill_lim)
}

grid_data <- grid_3[!is.na(grid_3$`mean slope`), ]
fill_lim <- NULL
p <- ggplot() + 
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = `mean slope`), alpha = 2*abs(grid_data$`p(positive slope)` - 0.5))   +
  geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  viridis::scale_fill_viridis(limits = fill_lim)
p

p <- ggplot() + 
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = `p(positive slope)`), alpha = .8)   +
  geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  viridis::scale_fill_viridis(limits = fill_lim)
p


fill_lim <- NULL
for (i in 14:22) {
  if(i == 16) {fill_lim <- c(0, 1.001)}
  plotting_fun(grid_3[!is.na(grid_3$julaug_mean_slope), ], sym(names(plotting_data2)[i]), fill_lim = fill_lim)
}

fill_lim <- NULL
for (i in 23:31) {
  if(i == 25){fill_lim <- c(0, 1.001)}
  plotting_fun(grid_3[!is.na(grid_3$full_janfeb_mean_slope), ], sym(names(plotting_data2)[i]), fill_lim = fill_lim)
  
}
  


##### CAR models #####
# format data
car_data <- plotting_data2[!is.na(plotting_data2$`mean slope`), c("cell", "mean slope", "janfeb_median_bayes",
                               "janfeb_var_bayes", "janfeb_skew_bayes", "janfeb_kurt_bayes")]
car_data$slope_scaled <- scale(car_data$`mean slope`)
car_data$lat <- dgSEQNUM_to_GEO(grid6, car_data$cell)$lat_deg
car_data$lat_scaled <- scale(car_data$lat)
car_data$cell_id <- paste0("cell_", car_data$cell)
car_data$known_se <- sqrt(car_data$janfeb_var_bayes)/sd(car_data$`mean slope`)

# get adjacency matrix
adjacency_mat <- matrix(data = 0L, nrow = nrow(car_data), ncol = nrow(car_data))
row.names(adjacency_mat) <- car_data$cell_id
for (i in 1:(nrow(car_data) - 1)) {
  coord_i <- dgSEQNUM_to_GEO(grid6, car_data$cell[i])
  for (j in (i+1):nrow(car_data)) {
    coord_j <- dgSEQNUM_to_GEO(grid6, car_data$cell[j])
    cell_dist <- geosphere::distm(c(coord_i$lon_deg, coord_i$lat_deg), 
                                  c(coord_j$lon_deg, coord_j$lat_deg), 
                                  fun = geosphere::distHaversine)
    if(cell_dist < 320000) {
      adjacency_mat[i,j] <- adjacency_mat[j,i] <- 1L
    }
  }
}


library(brms)
# Exact sparse CAR
escar2_fit <- brm(slope_scaled | resp_se(known_se, sigma = TRUE) ~ lat_scaled + car(M, gr = cell_id, type = "escar"),
                    data = car_data, data2 = list(M = adjacency_mat), backend = 'cmdstanr', iter = 12000, warmup = 2000, cores = 4)
# Intrinsic CAR (bym2)
bym2_fit <- brm(slope_scaled | resp_se(known_se, sigma = TRUE) ~ lat_scaled + car(M, gr = cell_id, type = "bym2"),
                 data = car_data, data2 = list(M = adjacency_mat), backend = 'cmdstanr', iter = 12000, warmup = 2000, cores = 4)


# # Legacy code to do resp_se by hand before I knew about resp_se in brms
# escar_code <- make_stancode(janfeb_median_bayes ~ lat_scaled + car(M, gr = cell_id, type = "escar"),
#                            data = car_data, data2 = list(M = adjacency_mat))
# escar2_code <- gsub("vector\\[N\\] Y;", "vector[N] Y;  // response variable\n   vector[N] Y_var;", escar_code)
# escar2_code <- gsub("\\nmodel \\{", "\n model {\n   vector[N] sigma_new = sqrt(sigma^2 + Y_var);", escar2_code)
# escar2_code <- gsub(" mu, sigma\\);", " mu, sigma_new);", escar2_code)
# cmdstanr::write_stan_file(escar2_code, dir = "code/macrodemography/stan",
#                           basename = "escar2")
# standata_escar <- make_standata(slope_scaled ~ lat_scaled + car(M, gr = cell_id, type = "escar"),
#               data = car_data, data2 = list(M = adjacency_mat))
# standata_escar$Y_var <- car_data$known_se^2
# escar2_mod <- cmdstanr::cmdstan_model("code/macrodemography/stan/escar2.stan")
# escar2_fit <- escar2_mod$sample(data = standata_escar, 
#                                 iter_sampling = 10000,
#                                 parallel_chains = 4)

summary(escar2_fit)
summary(bym2_fit)

#####
pl <-  posterior_epred(escar2_fit, re.form = NA, incl_autocor = F)
a <- apply(pl, 2, function(x){quantile(x, .05)})
for(i in 1:9){
  a <- cbind(a, apply(pl, 2, function(x){quantile(x, i/10 + .05)}))
}
q_frame <- as.data.frame(a)
names(q_frame) <- paste0("q_", .05 + .1*c(0:9))
q_frame$latitude <- car_data$lat
q_frame$med <- apply(pl, 2, median)

q_frame <- cbind(q_frame, car_data[c("slope_scaled", "known_se")])
q_frame <- rename(q_frame, slope = slope_scaled)

q_frame <- q_frame[order(q_frame$latitude),]
q_frame$point_lower <- q_frame$slope - q_frame$known_se
q_frame$point_upper <- q_frame$slope + q_frame$known_se
certainty <- 1.2 - q_frame$known_se/max(q_frame$known_se)
set.seed(1)
q_frame$lat_jit <- q_frame$latitude + rnorm(nrow(q_frame), 0, .5)

fill_color <- "salmon2"
alpha <- .5
point_color <- "gray25"

q_frame[, ! (names(q_frame) %in% c("latitude", "lat_jit"))] <-
  q_frame[, ! (names(q_frame) %in% c("latitude", "lat_jit"))] * sd(car_data$`mean slope`) +
  mean(car_data$`mean slope`)


ggplot(q_frame) + theme_classic() +
  geom_ribbon(aes(x = latitude, ymin = q_0.05, ymax = q_0.95), fill = fill_color, alpha = alpha) +
  geom_ribbon(aes(x = latitude, ymin = q_0.15, ymax = q_0.85), fill = fill_color, alpha = alpha) +
  geom_ribbon(aes(x = latitude, ymin = q_0.25, ymax = q_0.75), fill = fill_color, alpha = alpha) +
  geom_ribbon(aes(x = latitude, ymin = q_0.35, ymax = q_0.65), fill = fill_color, alpha = alpha) +
  geom_ribbon(aes(x = latitude, ymin = q_0.45, ymax = q_0.55), fill = fill_color, alpha = alpha) +
  geom_point(aes(x = lat_jit, y = slope, alpha = certainty), color = point_color) +
  geom_segment(aes(x = lat_jit, y = point_lower, xend = lat_jit, yend = point_upper, alpha = certainty), color = point_color) +
  geom_line(aes(x = latitude, y = med))



ggplot(car_data, aes(janfeb_skew_bayes)) + geom_density() + 
  xlim(c(-.5, .5)) + xlab("skewness") + ylab("") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


car_data$excess_kurtosis <- car_data$janfeb_kurt_bayes - 3
ggplot(car_data, aes(excess_kurtosis)) + geom_density() + 
  xlim(c(-2, 2)) + xlab("excess kurtosis") + ylab("")+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())




