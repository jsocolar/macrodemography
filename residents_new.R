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
library(brms)
library(ggplot2)

#cols_bd <- c(pals::brewer.brbg(401)[1:201], "white", pals::ocean.curl(401)[201:1])
cols_bd <- c(hsv(seq(0,.17,length.out = 100),1,seq(.9,.6,length.out = 100)), hsv(seq(.45,.65, length.out = 100),1,seq(.6,1,length.out = 100)))
cols_bd2 <- c(hsv(seq(0,.17,length.out = 100),seq(1, .2, length.out = 100),.9), hsv(seq(.45,.65, length.out = 100),seq(.2, 1, length.out = 100),.9))

blank_theme <- 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
        )

states <- map_data("state")
western_states <- c("arizona", "california", "colorado", "idaho", "montana",
                    "nevada", "new mexico", "oregon", "utah", "washington",
                    "wyoming")
states <- states[!states$region %in% western_states,]

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


##### Declare species #####
species <- "CARW"

##### Get demographic indices #####
spring_abun_data <- readRDS(
  paste0("macrodemography/residents/", species, "/spring_abun_data_fullwindow.RDS")
)
fall_abun_data <- readRDS(
  paste0("macrodemography/residents/", species, "/fall_abun_data_fullwindow.RDS")
)

assertthat::assert_that(all(spring_abun_data[[1]]$cell %in% cells_all))
assertthat::assert_that(all(fall_abun_data[[1]]$cell %in% cells_all))

# Extract cell-specific abundance data
spring_abun_summary <- get_abun_summary(spring_abun_data, 10)
fall_abun_summary <- get_abun_summary(fall_abun_data, 10)
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
# Get the standard deviations of the survival and productivity values
sd_holder_prod <- sd_holder_surv <- matrix(nrow = length(cells_all), ncol = 100)
colnames(sd_holder_prod) <- paste0("prod_sd_rep_", 1:100)
colnames(sd_holder_surv) <- paste0("surv_sd_rep_", 1:100)

cell_lrat_sd <- cbind(data.frame(cell = cells_all, n_prod = NA, n_surv = NA),
                      as.data.frame(sd_holder_prod), as.data.frame(sd_holder_surv))
var_p <- var_d <- rep(NA, length(cells_all))


lrat_skews <- lrat_kurts <- vector()
diagnostic_failure <- 0
for(i in seq_along(cells_all)){
  print(i)
  if(!identical(cell_ratio_series[[i]], NA)){
    lrats_avg <- cell_ratio_series[[i]]$avg
    lrats_avg[!use_cell_years(cell_ratio_series[[i]], inf_exclude=T)] <- NA
    if (cells_all[i] %in% cells & !all(is.na(lrats_avg))) {
      assertthat::assert_that(sum(!is.na(lrats_avg)) >= 5)
      prod_rats <- lrats_avg[1 + 2*c(1:13)]
      surv_rats <- lrats_avg[2*c(1:13)]
      cell_lrat_sd$n_prod[i] <- sum(!is.na(prod_rats))
      cell_lrat_sd$n_surv[i] <- sum(!is.na(surv_rats))

      lrat_means <- rowMeans(cell_ratio_series_full[[i]])
      lrat_sds <- apply(cell_ratio_series_full[[i]], 1, sd)
      lrat_skews1 <- apply(cell_ratio_series_full[[i]], 1, moments::skewness)
      lrat_kurts1 <- apply(cell_ratio_series_full[[i]], 1, moments::kurtosis)

      prod_df <- data.frame(ratio = lrat_means[1 + 2*c(1:13)],
                            sd = lrat_sds[1 + 2*c(1:13)],
                            season = "prod")
      prod_df$ratio[is.infinite(prod_df$ratio) | is.nan(prod_df$ratio)] <- NA

      surv_df <- data.frame(ratio = lrat_means[2*c(1:13)],
                            sd = lrat_sds[2*c(1:13)],
                            season = "surv")
      surv_df$ratio[is.infinite(surv_df$ratio) | is.nan(surv_df$ratio)] <- NA


      model_df <- rbind(prod_df, surv_df)

      mod_formula <- bf(ratio | resp_se(sd, sigma = TRUE) ~ season,
                        sigma ~ season)

      if(sum(!is.na(prod_df$ratio)) > 4 & sum(!is.na(surv_df$ratio)) > 4){
        lrat_skews <- c(lrat_skews, lrat_skews1)
        lrat_kurts <- c(lrat_kurts, lrat_kurts1)
        # mod <- brm(mod_formula, data = model_df, family = gaussian(),
        #            iter = 2000, warmup = 1000, chains = 3, refresh = 0,
        #            backend = "cmdstanr")
        # diagnostics <- check_brmsfit_diagnostics(mod)
        # if (!all(diagnostics)) {
        #   adapt_delta <- .8
        #   iter <- 2000
        #   if (!diagnostics[1]) {
        #     adapt_delta <- .99
        #   }
        #   if (!diagnostics[2]) {
        #     iter <- 4000
        #   }
        #   mod <- brm(mod_formula, data = model_df, family = gaussian(),
        #              iter = iter, warmup = 1000, chains = 3, refresh = 0,
        #              adapt_delta = adapt_delta,
        #              backend = "cmdstanr")
        #   diagnostics <- check_brmsfit_diagnostics(mod)
        # }
        # if (all(diagnostics)) {
        #   d <- as_draws_df(mod)
        #   var_p[i] <- mean(d$b_sigma_seasonsurv > 0)
        #   var_d[i] <- mean(d$b_sigma_seasonsurv)
        # } else {
        #   print("diagnostic failure!")
        #   diagnostic_failure <- diagnostic_failure + 1
        # }
      }
    }
  }
}

assertthat::assert_that(!diagnostic_failure)

skews_kurts <- data.frame(
  skew = lrat_skews,
  `excess kurtosis` = lrat_kurts - 3)

ggplot(skews_kurts, aes(skew)) + 
  geom_density() + theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) + xlab("skewness")

ggplot(skews_kurts, aes(excess.kurtosis)) + 
  geom_density() + theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) + xlab("excess kurtosis")

saveRDS(var_d, paste0("macrodemography/var_d_", species, ".RDS"))
saveRDS(var_p, paste0("macrodemography/var_p_", species, ".RDS"))
var_p <- readRDS(paste0("macrodemography/var_p_", species, ".RDS"))
var_d <- readRDS(paste0("macrodemography/var_d_", species, ".RDS"))

plotting_data <- data.frame(cell = cells_all, 
                            n_prod = cell_lrat_sd$n_prod, 
                            n_surv = cell_lrat_sd$n_surv, 
                            p_surv_var_larger = var_p,
                            surv_sd_diff = var_d)

grid6 <- dggridR::dgconstruct(res = 6)

grid <- dggridR::dgcellstogrid(grid6,plotting_data$cell,frame=TRUE,wrapcells=TRUE)
grid  <- merge(grid,plotting_data,by.x="cell")

n_min <- 5

grid_2 <- grid[!is.na(grid$n_prod), ]
grid_2 <- grid_2[grid_2$n_prod >= n_min & grid_2$n_surv >= n_min, ]
names(grid_2)[names(grid_2) == "p_surv_var_larger"] <- "p(survival)"


ggplot() + blank_theme + coord_fixed() +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_2,      aes(x=long, y=lat, group=group, fill=n_prod), alpha=0.7)    +
  geom_path   (data=grid_2,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  viridis::scale_fill_viridis(limits = c(5,13)) + xlim(c(-107, -65))

ggplot() + blank_theme + coord_fixed() +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_2,      aes(x=long, y=lat, group=group, fill=n_surv), alpha=0.7)    +
  geom_path   (data=grid_2,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  viridis::scale_fill_viridis(limits = c(5,13)) + xlim(c(-107, -65))

p <- ggplot() + blank_theme + coord_fixed() +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_2,      aes(x=long, y=lat, group=group, fill=`p(survival)`), alpha=0.7)    +
  geom_path   (data=grid_2,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd2, na.value=NA, limits = c(0, 1)) + xlim(c(-107, -65))
p

fl <- max(abs(grid_2$surv_sd_diff), na.rm = T) + .1
p <- ggplot() + coord_fixed() + blank_theme +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_2, aes(x=long, y=lat, group=group, fill = surv_sd_diff), alpha = 2*abs(grid_2$`p(survival)` - 0.5))   +
  geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd, na.value=NA, limits = c(-fl, fl)) + xlim(c(-107, -65))
p

##### Extract weather data from Earth Engine #####
# library(dggridR)
# library(sf)
# library(rgee)
# source("/Users/jacob/Dropbox/Work/Code/macrodemography/analysis/functions/daymet_extract.R")
# ee_Initialize()
# 
# dg <- dgconstruct(res = 6)
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
decmar_swe <- readRDS("macrodemography/weather/decmar_swe.RDS")

plotting_data2 <- data.frame(cell = cells_all, 
                             surv_length = NA,
                             prod_length = NA,
                             full_length = NA,
                             
                             janfeb_mean_bayes = NA,
                             janfeb_median_bayes = NA,
                             janfeb_sd_bayes = NA,
                             janfeb_skew_bayes = NA,
                             janfeb_kurt_bayes = NA,
                             janfeb_p_bayes = NA,
                             
                             julaug_mean_bayes = NA,
                             julaug_median_bayes = NA,
                             julaug_sd_bayes = NA,
                             julaug_skew_bayes = NA,
                             julaug_kurt_bayes = NA,
                             julaug_p_bayes = NA,
                             
                             swe_mean_bayes = NA,
                             swe_median_bayes = NA,
                             swe_sd_bayes = NA,
                             swe_skew_bayes = NA,
                             swe_kurt_bayes = NA,
                             swe_p_bayes = NA
)

janfeb_flag <- julaug_flag <- swe_flag <- vector()
janfeb_iter <- julaug_iter <- swe_iter <- rep(2000, length(cells_all))
for (i in seq_along(cells_all)) {
  print(i)
  if(i == 233){next} # This is a cell out over open ocean
  if (!identical(cell_ratio_series[[i]], NA)) {
    crsi <- cell_ratio_series[[i]]
    
    plotting_data2$surv_length[i] <- sum(is.finite(crsi$median[2*c(1:13)]))
    plotting_data2$prod_length[i] <- sum(is.finite(crsi$median[1 + 2*c(0:13)]))
    plotting_data2$full_length[i] <- sum(is.finite(crsi$median[1 + 2*c(1:13)] + crsi$median[2*c(1:13)]))
    
    janfeb_temp_slopes <- julaug_temp_slopes <- janfeb_swe_slopes <- vector()
    
    
    lrat_means <- rowMeans(cell_ratio_series_full[[i]])
    lrat_sds <- apply(cell_ratio_series_full[[i]], 1, sd)
    
    surv_means <- prod_means <- lrat_means
    surv_sds <- prod_sds <- lrat_sds
    surv_means[!use_cell_years(cell_ratio_series[[i]], 
                               inf_exclude = F,
                               n_min_prod = 0,
                               n_min_surv = 5,
                               n_min_full = 0)] <- NA
    surv_means <- surv_means[2*c(1:13)]
    surv_sds[!use_cell_years(cell_ratio_series[[i]], 
                             inf_exclude = F,
                             n_min_prod = 0,
                             n_min_surv = 5,
                             n_min_full = 0)] <- NA
    surv_sds <- surv_sds[2*c(1:13)]
    
    prod_means[!use_cell_years(cell_ratio_series[[i]], 
                               inf_exclude = F,
                               n_min_prod = 5,
                               n_min_surv = 0,
                               n_min_full = 0)] <- NA
    prod_means <- prod_means[1 + 2*c(0:13)]
    prod_sds[!use_cell_years(cell_ratio_series[[i]], 
                             inf_exclude = F,
                             n_min_prod = 5,
                             n_min_surv = 0,
                             n_min_full = 0)] <- NA
    prod_sds <- prod_sds[1 + 2*c(0:13)]
    
    if (sum(!is.na(surv_means)) > 2) {
      surv_df <- data.frame(mean = surv_means, sd = surv_sds, 
                            janfeb_temp = unlist(janfeb_temp[[i]][1:13]),
                            decmar_swe = sqrt(unlist(decmar_swe[[i]][1:13])))
      janfeb_mod <- brm(bf(mean | resp_se(sd, sigma = TRUE) ~ janfeb_temp),
                              data = surv_df, family = gaussian(), 
                              prior = prior(std_normal(), class = "b"),
                              iter = 2000, warmup = 1000, chains = 3, refresh = 0,
                              backend = "cmdstanr")
      diagnostics <- check_brmsfit_diagnostics(janfeb_mod)
      if (!all(diagnostics)) {
        adapt_delta <- .8
        iter <- 2000
        if (!diagnostics[1]) {
          adapt_delta <- .99
        }
        if (!diagnostics[2]) {
          iter <- janfeb_iter[i] <- 4000
        }
        janfeb_mod <- brm(bf(mean | resp_se(sd, sigma = TRUE) ~ janfeb_temp),
                          data = surv_df, family = gaussian(), 
                          prior = prior(std_normal(), class = "b"),
                          iter = iter, warmup = 1000, adapt_delta = adapt_delta,
                          chains = 3, refresh = 0,
                          backend = "cmdstanr")
        diagnostics <- check_brmsfit_diagnostics(janfeb_mod)
        if(!all(diagnostics)){
          janfeb_flag <- c(janfeb_flag, i)
        }
      }
      
      janfeb_temp_slopes <- as_draws_df(janfeb_mod)$b_janfeb_temp
      
      plotting_data2$janfeb_mean_bayes[i] <- mean(janfeb_temp_slopes)
      plotting_data2$janfeb_median_bayes[i] <- median(janfeb_temp_slopes)
      plotting_data2$janfeb_sd_bayes[i] <- sd(janfeb_temp_slopes)
      plotting_data2$janfeb_skew_bayes[i] <- moments::skewness(janfeb_temp_slopes)
      plotting_data2$janfeb_kurt_bayes[i] <- moments::kurtosis(janfeb_temp_slopes)
      plotting_data2$janfeb_p_bayes[i] <- sum(janfeb_temp_slopes > 0)
      
      if (sum((surv_df$decmar_swe > 0) & (!is.na(surv_df$mean))) > 2) {
        swe_mod <- brms::brm(bf(mean | resp_se(sd, sigma = TRUE) ~ decmar_swe),
                             prior = prior(std_normal(), class = "b"),
                             data = surv_df, family = gaussian(), 
                             iter = 2000, warmup = 1000, chains = 3, refresh = 0,
                             backend = "cmdstanr")
        diagnostics <- check_brmsfit_diagnostics(swe_mod)
        if (!all(diagnostics)) {
          adapt_delta <- .8
          iter <- 2000
          if (!diagnostics[1]) {
            adapt_delta <- .99
          }
          if (!diagnostics[2]) {
            iter <- swe_iter[i] <- 4000
          }
          swe_mod <- brm(bf(mean | resp_se(sd, sigma = TRUE) ~ decmar_swe),
                            data = surv_df, family = gaussian(), 
                            prior = prior(std_normal(), class = "b"),
                            iter = iter, warmup = 1000, adapt_delta = adapt_delta,
                            chains = 3, refresh = 0,
                            backend = "cmdstanr")
          diagnostics <- check_brmsfit_diagnostics(swe_mod)
          if(!all(diagnostics)){
            swe_flag <- c(swe_flag, i)
          }
        }
        
        
        swe_slopes <- as_draws_df(swe_mod)$b_decmar_swe
        
        plotting_data2$swe_mean_bayes[i] <- mean(swe_slopes)
        plotting_data2$swe_median_bayes[i] <- median(swe_slopes)
        plotting_data2$swe_sd_bayes[i] <- sd(swe_slopes)
        plotting_data2$swe_skew_bayes[i] <- moments::skewness(swe_slopes)
        plotting_data2$swe_kurt_bayes[i] <- moments::kurtosis(swe_slopes)
        plotting_data2$swe_p_bayes[i] <- sum(swe_slopes > 0)
      }
    }
    
    
    if (sum(!is.na(prod_means)) > 2) {
      prod_df <- data.frame(mean = prod_means, sd = prod_sds, 
                            julaug_temp = unlist(julaug_temp[[i]]))
      julaug_mod <- brms::brm(bf(mean | resp_se(sd, sigma = TRUE) ~ julaug_temp),
                              data = prod_df, family = gaussian(), 
                              prior = prior(std_normal(), class = "b"),
                              iter = 2000, warmup = 1000, chains = 3, refresh = 0,
                              backend = "cmdstanr")
      diagnostics <- check_brmsfit_diagnostics(julaug_mod)
      if (!all(diagnostics)) {
        adapt_delta <- .8
        iter <- 2000
        if (!diagnostics[1]) {
          adapt_delta <- .99
        }
        if (!diagnostics[2]) {
          iter <- julaug_iter[i] <- 4000
        }
        julaug_mod <- brm(bf(mean | resp_se(sd, sigma = TRUE) ~ julaug_temp),
                       data = surv_df, family = gaussian(), 
                       prior = prior(std_normal(), class = "b"),
                       iter = iter, warmup = 1000, adapt_delta = adapt_delta,
                       chains = 3, refresh = 0,
                       backend = "cmdstanr")
        diagnostics <- check_brmsfit_diagnostics(julaug_mod)
        if(!all(diagnostics)){
          julaug_flag <- c(julaug_flag, i)
        }
      }
      julaug_temp_slopes <- as_draws_df(julaug_mod)$b_julaug_temp
      
      plotting_data2$julaug_mean_bayes[i] <- mean(julaug_temp_slopes)
      plotting_data2$julaug_median_bayes[i] <- median(julaug_temp_slopes)
      plotting_data2$julaug_sd_bayes[i] <- sd(julaug_temp_slopes)
      plotting_data2$julaug_skew_bayes[i] <- moments::skewness(julaug_temp_slopes)
      plotting_data2$julaug_kurt_bayes[i] <- moments::kurtosis(julaug_temp_slopes)
      plotting_data2$julaug_p_bayes[i] <- sum(julaug_temp_slopes > 0)
    }
  }
}

saveRDS(plotting_data2, paste0("macrodemography/plotting_data2_", species, ".RDS"))
saveRDS(janfeb_iter, paste0("macrodemography/janfeb_iter_", species, ".RDS"))
saveRDS(julaug_iter, paste0("macrodemography/julaug_iter_", species, ".RDS"))
saveRDS(swe_iter, paste0("macrodemography/swe_iter_", species, ".RDS"))

saveRDS(janfeb_flag, paste0("macrodemography/janfeb_flag_", species, ".RDS"))
saveRDS(julaug_flag, paste0("macrodemography/julaug_flag_", species, ".RDS"))
saveRDS(swe_flag, paste0("macrodemography/swe_flag_", species, ".RDS"))


##### Plotting #####
plotting_data2 <- readRDS(paste0("macrodemography/plotting_data2_", species, ".RDS"))
janfeb_iter <- readRDS(paste0("macrodemography/janfeb_iter_", species, ".RDS"))
julaug_iter <- readRDS(paste0("macrodemography/julaug_iter_", species, ".RDS"))
swe_iter <- readRDS(paste0("macrodemography/swe_iter_", species, ".RDS"))

janfeb_flag <- readRDS(paste0("macrodemography/janfeb_flag_", species, ".RDS"))
julaug_flag <- readRDS(paste0("macrodemography/julaug_flag_", species, ".RDS"))
swe_flag <- readRDS(paste0("macrodemography/swe_flag_", species, ".RDS"))


names(plotting_data2)[names(plotting_data2) == "janfeb_mean_bayes"] <- "mean slope"

plotting_data2$`p(winter temp)` <- plotting_data2$janfeb_p_bayes/(3 * (janfeb_iter - 1000))
plotting_data2$`p(summer temp)` <- plotting_data2$julaug_p_bayes/(3 * (julaug_iter - 1000))
plotting_data2$`p(winter swe)` <- plotting_data2$swe_p_bayes/(3 * (swe_iter - 1000))

plotting_data2$lon <- dggridR::dgSEQNUM_to_GEO(grid6, plotting_data2$cell)$lon_deg

plotting_data2 <- plotting_data2[plotting_data2$lon > -107,]

grid <- dggridR::dgcellstogrid(grid6,plotting_data2$cell,frame=TRUE,wrapcells=TRUE)
grid_3  <- merge(grid,plotting_data2,by.x="cell")

plotting_fun <- function(grid_data, variable, fill_lim = NULL) {
  v <- variable
  p <- ggplot() + coord_fixed() + blank_theme +
    geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
    geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = !!v), alpha = 0.8)   +
    geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
    scale_fill_gradientn(colours = cols_bd2, na.value=NA) + xlim(c(-107, -65))
  print(p)
}

grid_3 <- grid_3[grid_3$long > -105, ]

fill_lim <- NULL
for (i in 2:25) {
  plotting_fun(grid_3[!is.na(grid_3$`mean slope`), ], sym(names(plotting_data2)[i]), fill_lim = fill_lim)
}


grid_data <- grid_3[!is.na(grid_3$`mean slope`), ]
fl <- max(abs(grid_3$`mean slope`), na.rm = T) + .1
fill_lim <- c(-fl, fl)
p <- ggplot() + coord_fixed() + blank_theme + 
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = `mean slope`), alpha = 2*abs(grid_data$`p(winter temp)` - 0.5))   +
  geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd, na.value=NA, limits = fill_lim)
p

grid_data <- grid_3[!is.na(grid_3$julaug_mean_bayes), ]
fl <- max(abs(grid_3$julaug_mean_bayes), na.rm = T) + .1
fill_lim <- c(-fl, fl)
p <- ggplot() + coord_fixed() + blank_theme + 
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = julaug_mean_bayes), alpha = 2*abs(grid_data$`p(summer temp)` - 0.5))   +
  geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd, na.value=NA, limits = fill_lim)
p + labs(fill="mean slope")

grid_data <- grid_3[!is.na(grid_3$swe_mean_bayes), ]
fl <- max(abs(grid_3$swe_mean_bayes), na.rm = T) + .1
fill_lim <- c(-fl, fl)
p <- ggplot() +  coord_fixed() + blank_theme + 
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = swe_mean_bayes), alpha = 2*abs(grid_data$`p(winter temp)` - 0.5))   +
  geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd, na.value=NA, limits = fill_lim)
p + labs(fill="mean slope")


grid_data <- grid_3[!is.na(grid_3$`p(winter temp)`), ]
p <- ggplot() + blank_theme + coord_fixed() +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = `p(winter temp)`), alpha = .8)   +
  geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd2, na.value=NA, limits = c(0,1))
p + labs(fill="probability")


grid_data <- grid_3[!is.na(grid_3$`p(summer temp)`), ]
p <- ggplot() + blank_theme + coord_fixed() +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = `p(summer temp)`), alpha = .8)   +
  geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd2, na.value=NA, limits = c(0,1))
p + labs(fill="probability")

grid_data <- grid_3[!is.na(grid_3$`p(winter swe)`), ]
p <- ggplot() + blank_theme + coord_fixed() +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = `p(winter swe)`), alpha = .8)   +
  geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd2, na.value=NA, limits = c(0,1))
p + labs(fill="probability")


plotting_data2$excess_kurtosis_winter_temp <- plotting_data2$janfeb_kurt_bayes - 3
plotting_data2$excess_kurtosis_summer_temp <- plotting_data2$julaug_kurt_bayes - 3
plotting_data2$excess_kurtosis_winter_snow <- plotting_data2$swe_kurt_bayes - 3

pd3 <- dplyr::rename(plotting_data2,
                                skewness_winter_temp = janfeb_skew_bayes,
                                skewness_winter_snow = swe_skew_bayes,
                                skewness_summer_temp = julaug_skew_bayes)

pd3 <- tidyr::pivot_longer(
  pd3,
  c("skewness_winter_temp", "skewness_winter_snow", "skewness_summer_temp",
    "excess_kurtosis_winter_temp", "excess_kurtosis_winter_snow", "excess_kurtosis_summer_temp")
)

pd3$name <- factor(pd3$name, levels=c("skewness_winter_temp", "skewness_winter_snow", "skewness_summer_temp",
                          "excess_kurtosis_winter_temp", "excess_kurtosis_winter_snow", "excess_kurtosis_summer_temp"))

name_labels <- c(
  "skewness (winter temp)", "skewness (winter snow)", "skewness (summer temp)",
  "excess kurtosis (winter temp)", "excess kurtosis (winter snow)", "excess kurtosis (summer temp)")

names(name_labels) <- 
  c("skewness_winter_temp", "skewness_winter_snow", "skewness_summer_temp",
    "excess_kurtosis_winter_temp", "excess_kurtosis_winter_snow", "excess_kurtosis_summer_temp")

ggplot(pd3, aes(value)) + geom_density() + 
  facet_wrap("name",
             labeller = labeller(name = name_labels)) + 
  ylab("") +
  xlab("") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank())


##### CAR models #####
plotting_data2 <- plotting_data2[plotting_data2$cell != dggridR::dgGEO_to_SEQNUM(grid6, -102, 41)$seqnum, ]
##### janfeb #####
# format data
car_data <- plotting_data2[!is.na(plotting_data2$`mean slope`), c("cell", "mean slope", "janfeb_median_bayes",
                                                                  "janfeb_sd_bayes", "janfeb_skew_bayes", "janfeb_kurt_bayes")]
car_data$slope_scaled <- scale(car_data$`mean slope`)
car_data$lat <- dggridR::dgSEQNUM_to_GEO(grid6, car_data$cell)$lat_deg
car_data$lat_scaled <- scale(car_data$lat)
car_data$cell_id <- paste0("cell_", car_data$cell)
car_data$known_se <- car_data$janfeb_sd_bayes/sd(car_data$`mean slope`)

car_data$lon <- dggridR::dgSEQNUM_to_GEO(grid6, car_data$cell)$lon_deg
car_data <- car_data[!(car_data$lat > 38.7 & car_data$lon < -99.5), ]

# get adjacency matrix
adjacency_mat <- matrix(data = 0L, nrow = nrow(car_data), ncol = nrow(car_data))
row.names(adjacency_mat) <- car_data$cell_id
for (i in 1:(nrow(car_data) - 1)) {
  coord_i <- dggridR::dgSEQNUM_to_GEO(grid6, car_data$cell[i])
  for (j in (i+1):nrow(car_data)) {
    coord_j <- dggridR::dgSEQNUM_to_GEO(grid6, car_data$cell[j])
    cell_dist <- geosphere::distm(c(coord_i$lon_deg, coord_i$lat_deg), 
                                  c(coord_j$lon_deg, coord_j$lat_deg), 
                                  fun = geosphere::distHaversine)
    if(cell_dist < 320000) {
      adjacency_mat[i,j] <- adjacency_mat[j,i] <- 1L
    }
  }
}

# CAR model
escar2_fit <- brm(slope_scaled | resp_se(known_se, sigma = TRUE) ~ lat_scaled + car(M, gr = cell_id, type = "icar"),
                  data = car_data, data2 = list(M = adjacency_mat), backend = 'cmdstanr', iter = 12000, warmup = 2000, cores = 4,
                  adapt_delta = .8, max_treedepth = 11)
summary(escar2_fit)
params <- as_draws_df(escar2_fit)

rstan::get_bfmi(escar2_fit$fit)

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
q_frame <- dplyr::rename(q_frame, slope = slope_scaled)

q_frame <- q_frame[order(q_frame$latitude),]
q_frame$point_lower <- q_frame$slope - q_frame$known_se
q_frame$point_upper <- q_frame$slope + q_frame$known_se
certainty <- min(q_frame$known_se)/q_frame$known_se
set.seed(1)
q_frame$lat_jit <- q_frame$latitude + rnorm(nrow(q_frame), 0, .5)

fill_color <- "salmon2"
alpha <- .5
point_color <- "gray25"

q_frame[, ! (names(q_frame) %in% c("latitude", "lat_jit"))] <-
  q_frame[, ! (names(q_frame) %in% c("latitude", "lat_jit"))] * sd(car_data$`mean slope`) +
  mean(car_data$`mean slope`)


q_frame$point_lower[q_frame$point_lower < -.5] <- -.5

ggplot(q_frame) + theme_classic() +
  geom_ribbon(aes(x = latitude, ymin = q_0.05, ymax = q_0.95), fill = fill_color, alpha = alpha) +
  geom_ribbon(aes(x = latitude, ymin = q_0.15, ymax = q_0.85), fill = fill_color, alpha = alpha) +
  geom_ribbon(aes(x = latitude, ymin = q_0.25, ymax = q_0.75), fill = fill_color, alpha = alpha) +
  geom_ribbon(aes(x = latitude, ymin = q_0.35, ymax = q_0.65), fill = fill_color, alpha = alpha) +
  geom_ribbon(aes(x = latitude, ymin = q_0.45, ymax = q_0.55), fill = fill_color, alpha = alpha) +
  geom_point(aes(x = lat_jit, y = slope, alpha = certainty), color = point_color) +
  geom_segment(aes(x = lat_jit, y = point_lower, xend = lat_jit, yend = point_upper, alpha = certainty), color = point_color) +
  geom_line(aes(x = latitude, y = med)) +
  ylim(c(-.5, .5))



ggplot(car_data, aes(janfeb_skew_bayes)) + geom_density() + xlab("skewness") + ylab("") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


car_data$excess_kurtosis <- car_data$janfeb_kurt_bayes - 3
ggplot(car_data, aes(excess_kurtosis)) + geom_density() + 
  xlab("excess kurtosis") + ylab("")+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

escar2_fit <- brm(slope_scaled | resp_se(known_se, sigma = TRUE) ~ car(M, gr = cell_id, type = "icar"),
                  data = car_data, data2 = list(M = adjacency_mat), backend = 'cmdstanr', iter = 12000, 
                  warmup = 2000, cores = 4)
summary(escar2_fit)
npt <- nrow(car_data)

true_values <- matrix(nrow = 40000, ncol = npt)

for(i in 1:npt) {
  print(i)
  M <- car_data$slope_scaled[i]
  se <- car_data$known_se[i]
  LP <-  posterior_linpred(escar2_fit, re.form = NULL, incl_autocor = T)[ , i]
  sigma <- as_draws_df(escar2_fit)$sigma
  mu = (M*sigma^2 + LP * se^2)/(sigma^2 + se^2)
  scale = 1/sqrt(1/sigma^2 + 1/se^2)
  
  true_values[ , i] <- rnorm(40000, mu, scale)
}

smooth_prob <- data.frame(cell = car_data$cell, smooth_prob = apply(true_values, 2, function(x){mean(x +  mean(car_data$`mean slope`)/sd(car_data$`mean slope`) > 0)}))
smooth_mean <- data.frame(cell = car_data$cell, smooth_mean = apply(true_values, 2, function(x){mean(x) + mean(car_data$`mean slope`)/sd(car_data$`mean slope`)}))

plotting_data3 <- merge(plotting_data2, smooth_prob, by = "cell")
plotting_data3 <- merge(plotting_data3, smooth_mean, by = "cell")

grid <- dggridR::dgcellstogrid(grid6,plotting_data3$cell,frame=TRUE,wrapcells=TRUE)
grid_3  <- merge(grid,plotting_data3,by.x="cell")

grid_data <- grid_3[!is.na(grid_3$smooth_prob), ]
p <- ggplot() + coord_fixed() + blank_theme +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = smooth_prob), alpha = .8)   +
  geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd2, na.value=NA, limits = c(0,1)) + xlim(c(-107, -65))
p

fl <- max(abs(grid_data$smooth_mean), na.rm = T) + .1
p <- ggplot() + coord_fixed() + blank_theme +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = smooth_mean), alpha = 2*abs(grid_data$smooth_prob - 0.5))   +
  geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd, na.value=NA, limits = c(-fl, fl)) + xlim(c(-107, -65))
p

##### julaug #####
# format data
car_data <- plotting_data2[!is.na(plotting_data2$julaug_mean_bayes), c("cell", "julaug_mean_bayes", "julaug_median_bayes",
                                                                  "julaug_sd_bayes", "julaug_skew_bayes", "julaug_kurt_bayes")]
car_data$slope_scaled <- scale(car_data$julaug_mean_bayes)
car_data$lat <- dggridR::dgSEQNUM_to_GEO(grid6, car_data$cell)$lat_deg
car_data$lat_scaled <- scale(car_data$lat)
car_data$cell_id <- paste0("cell_", car_data$cell)
car_data$known_se <- car_data$julaug_sd_bayes/sd(car_data$julaug_mean_bayes)

car_data$lon <- dggridR::dgSEQNUM_to_GEO(grid6, car_data$cell)$lon_deg
car_data <- car_data[!(car_data$lat > 38.7 & car_data$lon < -99.5), ]

# get adjacency matrix
adjacency_mat <- matrix(data = 0L, nrow = nrow(car_data), ncol = nrow(car_data))
row.names(adjacency_mat) <- car_data$cell_id
for (i in 1:(nrow(car_data) - 1)) {
  coord_i <- dggridR::dgSEQNUM_to_GEO(grid6, car_data$cell[i])
  for (j in (i+1):nrow(car_data)) {
    coord_j <- dggridR::dgSEQNUM_to_GEO(grid6, car_data$cell[j])
    cell_dist <- geosphere::distm(c(coord_i$lon_deg, coord_i$lat_deg), 
                                  c(coord_j$lon_deg, coord_j$lat_deg), 
                                  fun = geosphere::distHaversine)
    if(cell_dist < 320000) {
      adjacency_mat[i,j] <- adjacency_mat[j,i] <- 1L
    }
  }
}


ggplot(car_data, aes(julaug_skew_bayes)) + geom_density() + 
  xlim(c(-.5, .5)) + xlab("skewness") + ylab("") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


car_data$excess_kurtosis <- car_data$julaug_kurt_bayes - 3
ggplot(car_data, aes(excess_kurtosis)) + geom_density() + 
  xlim(c(-2, 2)) + xlab("excess kurtosis") + ylab("")+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

escar2_fit <- brm(slope_scaled | resp_se(known_se, sigma = TRUE) ~ car(M, gr = cell_id, type = "icar"),
                  data = car_data, data2 = list(M = adjacency_mat), backend = 'cmdstanr', iter = 12000, warmup = 2000, cores = 4,
                  refresh = 1000)
summary(escar2_fit)
npt <- nrow(car_data)

true_values <- matrix(nrow = 40000, ncol = npt)

for(i in 1:npt) {
  print(i)
  M <- car_data$slope_scaled[i]
  se <- car_data$known_se[i]
  LP <-  posterior_linpred(escar2_fit, re.form = NULL, incl_autocor = T)[ , i]
  sigma <- as_draws_df(escar2_fit)$sigma
  mu = (M*sigma^2 + LP * se^2)/(sigma^2 + se^2)
  scale = 1/sqrt(1/sigma^2 + 1/se^2)
  
  true_values[ , i] <- rnorm(40000, mu, scale)
}

smooth_prob <- data.frame(cell = car_data$cell, smooth_prob = apply(true_values, 2, function(x){mean(x +  mean(car_data$julaug_mean_bayes)/sd(car_data$julaug_mean_bayes) > 0)}))
smooth_mean <- data.frame(cell = car_data$cell, smooth_mean = apply(true_values, 2, function(x){mean(x) + mean(car_data$julaug_mean_bayes)/sd(car_data$julaug_mean_bayes)}))

plotting_data3 <- merge(plotting_data2, smooth_prob, by = "cell")
plotting_data3 <- merge(plotting_data3, smooth_mean, by = "cell")

grid <- dggridR::dgcellstogrid(grid6,plotting_data3$cell,frame=TRUE,wrapcells=TRUE)
grid_3  <- merge(grid,plotting_data3,by.x="cell")

grid_data <- grid_3[!is.na(grid_3$smooth_prob), ]
p <- ggplot() + coord_fixed() + blank_theme +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = smooth_prob), alpha = .8)   +
  geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd2, na.value=NA, limits = c(0,1)) + xlim(c(-107, -65))
p

fl <- max(abs(grid_data$smooth_mean), na.rm = T) + .1
p <- ggplot() + coord_fixed() + blank_theme +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = smooth_mean), alpha = 2*abs(grid_data$smooth_prob - 0.5))   +
  geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd, na.value=NA, limits = c(-fl, fl)) + xlim(c(-107, -65))
p




##### swe #####
# format data
car_data <- plotting_data2[!is.na(plotting_data2$swe_mean_bayes), c("cell", "swe_mean_bayes", "swe_median_bayes",
                                                                       "swe_sd_bayes", "swe_skew_bayes", "swe_kurt_bayes")]
car_data <- car_data[car_data$swe_mean_bayes < 1 & car_data$swe_mean_bayes > -1, ]
car_data$slope_scaled <- scale(car_data$swe_mean_bayes)
car_data$lat <- dggridR::dgSEQNUM_to_GEO(grid6, car_data$cell)$lat_deg
car_data$lat_scaled <- scale(car_data$lat)
car_data$cell_id <- paste0("cell_", car_data$cell)
car_data$known_se <- car_data$swe_sd_bayes/sd(car_data$swe_mean_bayes)

car_data$lon <- dggridR::dgSEQNUM_to_GEO(grid6, car_data$cell)$lon_deg
car_data <- car_data[!(car_data$lat > 38.7 & car_data$lon < -99.5), ]

# get adjacency matrix
adjacency_mat <- matrix(data = 0L, nrow = nrow(car_data), ncol = nrow(car_data))
row.names(adjacency_mat) <- car_data$cell_id
for (i in 1:(nrow(car_data) - 1)) {
  coord_i <- dggridR::dgSEQNUM_to_GEO(grid6, car_data$cell[i])
  for (j in (i+1):nrow(car_data)) {
    coord_j <- dggridR::dgSEQNUM_to_GEO(grid6, car_data$cell[j])
    cell_dist <- geosphere::distm(c(coord_i$lon_deg, coord_i$lat_deg), 
                                  c(coord_j$lon_deg, coord_j$lat_deg), 
                                  fun = geosphere::distHaversine)
    if(cell_dist < 320000) {
      adjacency_mat[i,j] <- adjacency_mat[j,i] <- 1L
    }
  }
}


ggplot(car_data, aes(swe_skew_bayes)) + geom_density() + 
  xlab("skewness") + ylab("") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


car_data$excess_kurtosis <- car_data$swe_kurt_bayes - 3
ggplot(car_data, aes(excess_kurtosis)) + geom_density() + 
   xlab("excess kurtosis") + ylab("")+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

escar2_fit <- brm(slope_scaled | resp_se(known_se, sigma = TRUE) ~ car(M, gr = cell_id, type = "icar"),
                  data = car_data, data2 = list(M = adjacency_mat), iter = 12000, 
                  warmup = 2000, backend = 'cmdstanr', cores = 4)
summary(escar2_fit)

npt <- nrow(car_data)
true_values <- matrix(nrow = 40000, ncol = npt)

for(i in 1:npt) {
  print(i)
  M <- car_data$slope_scaled[i]
  se <- car_data$known_se[i]
  LP <-  posterior_linpred(escar2_fit, re.form = NULL, incl_autocor = T)[ , i]
  sigma <- as_draws_df(escar2_fit)$sigma
  mu = (M*sigma^2 + LP * se^2)/(sigma^2 + se^2)
  scale = 1/sqrt(1/sigma^2 + 1/se^2)
  
  true_values[ , i] <- rnorm(40000, mu, scale)
}

smooth_prob <- data.frame(cell = car_data$cell, smooth_prob = apply(true_values, 2, function(x){mean(x +  mean(car_data$swe_mean_bayes)/sd(car_data$swe_mean_bayes) > 0)}))
smooth_mean <- data.frame(cell = car_data$cell, smooth_mean = apply(true_values, 2, function(x){mean(x) + mean(car_data$swe_mean_bayes)/sd(car_data$swe_mean_bayes)}))

plotting_data3 <- merge(plotting_data2, smooth_prob, by = "cell")
plotting_data3 <- merge(plotting_data3, smooth_mean, by = "cell")

grid <- dggridR::dgcellstogrid(grid6,plotting_data3$cell,frame=TRUE,wrapcells=TRUE)
grid_3  <- merge(grid,plotting_data3,by.x="cell")

grid_data <- grid_3[!is.na(grid_3$smooth_prob), ]
p <- ggplot() + coord_fixed() + blank_theme +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = smooth_prob), alpha = .8)   +
  geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd2, na.value=NA, limits = c(0,1)) + xlim(c(-107, -65))
p

fl <- max(abs(grid_data$smooth_mean), na.rm = T) + .1
p <- ggplot() + coord_fixed() + blank_theme +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = smooth_mean), alpha = 2*abs(grid_data$smooth_prob - 0.5))   +
  geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd, na.value=NA, limits = c(-fl, fl)) + xlim(c(-107, -65))
p




average_indices <- data.frame(cell = cells_all, prod_mean = NA, surv_mean = NA)
for(i in seq_along(cells_all)){
  print(i)
  if(!identical(cell_ratio_series[[i]], NA)){
    lrats_avg <- cell_ratio_series[[i]]$avg
    lrats_avg[!use_cell_years(cell_ratio_series[[i]], inf_exclude=T)] <- NA
    if (cells_all[i] %in% cells) {
      prod_rats <- lrats_avg[1 + 2*c(1:13)]
      surv_rats <- lrats_avg[2*c(1:13)]
      average_indices$prod_mean[i] <- mean(prod_rats, na.rm = T)
      average_indices$surv_mean[i] <- mean(surv_rats, na.rm = T)
      
    }
  }
}

average_indices2 <- average_indices[(!is.na(average_indices$prod_mean)) | (!is.na(average_indices$surv_mean)), ]

grid6 <- dggridR::dgconstruct(res = 6)

grid_x <- dggridR::dgcellstogrid(grid6,average_indices2$cell,frame=TRUE,wrapcells=TRUE)
grid_x  <- merge(grid_x,average_indices2,by.x="cell")


p <- ggplot() + coord_fixed() + blank_theme +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_x, aes(x=long, y=lat, group=group, fill = prod_mean), alpha = .8)   +
  geom_path(data=grid_x, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  viridis::scale_fill_viridis()
p + labs(fill="mean index")


p <- ggplot() + coord_fixed() + blank_theme +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_x, aes(x=long, y=lat, group=group, fill = surv_mean), alpha = .8)   +
  geom_path(data=grid_x, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  viridis::scale_fill_viridis()
p + labs(fill="mean index")
