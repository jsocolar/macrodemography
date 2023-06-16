
####
#### Parameters -------------------------------------------------------------
####

params <- list()
if(Sys.info()[["user"]]=="amd427"){
  params$checklists_parquet_path <- "~/Documents/ebird2022/checklists-2022.parquet"
  params$erd_path <- "~/Documents/ebird2022/observations-2022.parquet"
  params$output_path <- "~/Dropbox/macrodemography_refactor/data/residents"
} else if(Sys.info()[["user"]]=="bg423"){
  params$checklists_parquet_path <- "~/Documents/ebird2022/checklists-2022.parquet"
  params$erd_path <- "~/Documents/ebird2022/observations-2022.parquet"
  params$output_path <- "~/Documents/macrodemography/data_parquet"
}
params$years <- c(2005:2022) # We could start from 2002, when ebird has started!
params$extent_space <-  data.frame( min_lon=-125, max_lon=-66, min_lat=17, max_lat=65 )
params$period <- c("spring", "fall")
params$time_grid <- 7
params$tgrid_min <- c(20, 31)  # May15-June07         # we vould make it specific to sp, but downstream analyses 
params$tgrid_max <- c(22, 33)  # Jul21-Aug25          # with multiple sp may become problematic!
# params$max_altitude <- 2000                         # parquet file does not have altitude, and cci!
# params$max_altitude_above_lat42 <- 1500
params$effort_thresholds <- data.frame( dist_max=3, time_min=5/60, time_max=1, cci_min=0 )
params$hexagon_area_large <- 70000
params$hexagon_area_small <- 300
params$n_small_min <- 10
params$n_year_min <- 5
params$daymet <- data.frame(label=c("tmax_jan","tmax_summer","swe"), 
                           variable=c("tmax","tmax","swe"), 
                           date_min=c("01-01","06-15","01-01"), 
                           date_max=c("01-31","07-15","01-31"), 
                           period=c("spring","fall","spring"))
params$species_to_process <- c("easpho", "amerob", "treswa", "eastow", "sonspa", "buhvir")
params$always_import_checklists  <- FALSE
params$always_filter_checklists  <- FALSE
params$always_resample_bootstrap <- FALSE
params$always_run_variance_test  <- FALSE
params$always_download_weather   <- FALSE
params$always_run_regressions    <- TRUE
params$quiet                     <- TRUE
params$region <- "north_america"
params$plotting_xlim <- c(-126, -65)
params$extent_time <-
  data.frame(
    period = params$period,
    tgrid_min = params$tgrid_min,
    tgrid_max = params$tgrid_max,
    year_min = min(params$years),
    year_max = max(params$years)
  )

####
#### Load packages -------------------------------------------------------------
####

# Install an older version of dggridR from source. See:
# https://github.com/r-barnes/dggridR/issues/63#issuecomment-1454929653
remotes::install_github("r-barnes/dggridR", ref = "ec2a040")

# we make sure we have the same version of the erdPackage as this .rmd
package_path <- paste0(strsplit(rstudioapi::getSourceEditorContext()$path, "R_scripts")[[1]][1], "macrodemography") 
install.packages(package_path, repos = NULL, type = 'source')

# package load
library(macrodemography)
library(data.table) # don't worry about openMP warnings on Mac.
library(brms)
library(ggplot2)
library(dplyr)
library(magrittr)
library(sf)
library(dtplyr) # enable dplyr for data.table
library(assertthat)
library(ebirdst)
library(fasterize)
library(dggridR)
library(rgee)
library(tidyr)

# this package also needs a working version of rcmdstan:
# see https://mc-stan.org/cmdstanr/
cmdstanr::check_cmdstan_toolchain()
assert_that(cmdstanr::cmdstan_version() >= "2.31")

####
#### Color scales / themes -------------------------------------------------------------
####
# set color scales
cols_bd <- c(hsv(seq(0,.17,length.out = 100),1,seq(.9,.6,length.out = 100)), hsv(seq(.45,.65, length.out = 100),1,seq(.6,1,length.out = 100)))
cols_bd2 <- c(hsv(seq(0,.17,length.out = 100),seq(1, .2, length.out = 100),.9), hsv(seq(.45,.65, length.out = 100),seq(.2, 1, length.out = 100),.9))

# blank ggplot2 theme
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

####
#### Regions of interest -------------------------------------------------------------
####
western_states <- c("arizona", "california", "colorado", "idaho", "montana",
                    "nevada", "new mexico", "oregon", "utah", "washington",
                    "wyoming")

if(params$region == "eastern_us"){
  region_of_interest <- spData::us_states %>% filter(!(tolower(NAME) %in% western_states))
} else if(params$region == "western_us") {
  region_of_interest <- spData::us_states %>% filter(toupper(NAME) %in% western_states)
} else if(params$region == "conus") {
  region_of_interest <- spData::us_states
} else if(params$region == "north_america") {
  region_of_interest <- rnaturalearth::ne_states(country =  c("United States of America","Canada", "Mexico"), returnclass = "sf") %>%
    filter(!region %in% "Northern Canada") %>%
    filter(!name %in% c("Hawaii", "Alaska"))
} else{
  region_of_interest <- spData::us_states %>% filter(tolower(NAME) %in% params$region)
}

# NOTE: this adds additional cropping of the region of interest
# adjust if necessary
raster_of_interest <- fasterize::fasterize(region_of_interest,
                                           raster::raster(ncol=1000, nrow = 1000,
                                                          xmn = -125, xmx = -66,
                                                          ymn =17, ymx = 65))
####
#### Define hexagon grids -------------------------------------------------------------
####

# for an area of ~ 70000 km^2, we get a resolution of 6:
grid_large <- dggridR::dgconstruct(area = params$hexagon_area_large)
# for an area of ~ 300 km^2, we get a resolution of 11:
grid_small <- dggridR::dgconstruct(area = params$hexagon_area_small)

####
#### Import checklists -------------------------------------------------------------
####

# print some species info
ebirdst::ebirdst_runs %>% filter(substr(species_code,1,6) %in% species_codes$six)

checklists_path <- paste0(params$output_path, "/checklists.RDS")
filtered_checklists_path <- paste0(params$output_path, "/filtered_checklists.RDS")


if(params$always_import_checklists | !file.exists(checklists_path)){
  checklists <- import_checklists(params$checklists_parquet_path)
  saveRDS(checklists, checklists_path)
}


if(params$always_filter_checklists | !file.exists(filtered_checklists_path)){
  checklists <- readRDS(checklists_path) %>%
    filter(latitude > params$extent_space$min_lat &
             latitude < params$extent_space$max_lat &
             longitude > params$extent_space$min_lon &
             longitude < params$extent_space$max_lon &
             year >= min(params$years) &
             year <= max(params$years) &
             #             ELEV_30M_MEDIAN < params$max_altitude &
             #             ((ELEV_30M_MEDIAN < params$max_altitude_above_lat42) | (latitude < 42)) &
             #             cci > params$effort_thresholds$cci_min &
             effort_distance_km <= params$effort_thresholds$dist_max &
             effort_hrs >= params$effort_thresholds$time_min &
             effort_hrs <= params$effort_thresholds$time_max
    ) %>%
    mutate(
      seqnum_large = dggridR::dgGEO_to_SEQNUM(
        grid_large,
        .$longitude,
        .$latitude
      )[[1]],
      seqnum_small = dggridR::dgGEO_to_SEQNUM(
        grid_small,
        .$longitude,
        .$latitude
      )[[1]]
    )
  saveRDS(checklists, filtered_checklists_path)
} else {
  checklists <- readRDS(filtered_checklists_path)
}

# extract unique large hexagons
cells_all <- unique(checklists$seqnum_large)

####
#### Bootstrap abundances (summer book-ended)  -------------------------------------------------------------
####

for(species_code in params$species_to_process){
  file_out <- paste0(params$output_path, "/abun_data/", species_code , "_summer_bookend.rds")
  
  if(params$always_resample_bootstrap | !file.exists(file_out)){
    data <-
      sample_grid_abun(
        species_code, params$erd_path, checklists, params$effort_thresholds,
        params$extent_space, params$extent_time, time_window="full",
        small_grid=grid_small$res, large_grid=grid_large$res, time_grid=7,
        roi = raster_of_interest, quiet = FALSE
      )
    if(!dir.exists(dirname(file_out))) dir.create(dirname(file_out), recursive = TRUE)
    saveRDS(data, file_out)
  }
}


####
#### Calculate summer book-ended demographic indices (spring/fall log-ratios)  --------------------------
####

for(species_code in params$species_to_process){
  ##### load abundance data
  file_species <- paste0(params$output_path, "/abun_data/", species_code , "_summer_bookend.rds")
  
  print(paste("loading data from file", file_species,"..."))
  data <- readRDS(file_species)
  
  # Get demographic indices
  cell_ratios <- get_ratios(data$abun, cells_all, n_small_min = params$n_small_min, quiet=params$quiet)
  tidy_ratios <- make_ratios_tidy(cell_ratios)
  
  # rename fall/spring ratios to productivity/recruitment:
  tidy_ratios$summary %>%
    mutate(season=ifelse(period=="fall", "prod","surv")) -> tidy_ratios$summary
  
  # save the ratio data
  saveRDS(tidy_ratios, paste0(params$output_path, "/abun_data/", species_code, "_summer_bookend_ratios.rds"))
}


####
#### Bootstrap abundances ("winter" book-ended)  -------------------------------------------------------------
####

# Re-defining a few params for "winter_bookended" abundance and indices for weather regression on overwintering survival. 

params$period <- c("winter_start", "winter_end")
params$tgrid_min <- c(48, 6)  # date_to_st_week(c("2021-12-01", "2021-12-21", "2021-02-07", "2021-02-28"))
params$tgrid_max <- c(51, 9)  # 48 51  6  9          
params$always_resample_bootstrap <- FALSE
params$extent_time <-
  data.frame(
    period = params$period,
    tgrid_min = params$tgrid_min,
    tgrid_max = params$tgrid_max,
    year_min = min(params$years),
    year_max = max(params$years)
  )

for(species_code in params$species_to_process){
  file_out <- paste0(params$output_path, "/abun_data/", species_code , "_winter_bookend.rds")
  
  if(params$always_resample_bootstrap | !file.exists(file_out)){
    data <-
      sample_grid_abun(
        species_code, params$erd_path, checklists, params$effort_thresholds,
        params$extent_space, params$extent_time, time_window="full",
        small_grid=grid_small$res, large_grid=grid_large$res, time_grid=7,
        roi = raster_of_interest, quiet = FALSE
      )
    if(!dir.exists(dirname(file_out))) dir.create(dirname(file_out), recursive = TRUE)
    saveRDS(data, file_out)
  }
}

####
#### Calculate winter book-ended demographic indices (winter_start/winter_end log-ratios)  --------------------------
####

for(species_code in params$species_to_process){
  ##### load abundance data
  file_species <- paste0(params$output_path, "/abun_data/", species_code , "_winter_bookend.rds")
  
  print(paste("loading data from file", file_species,"..."))
  data <- readRDS(file_species)
  
  # Get demographic indices
  cell_ratios <- get_ratios(data$abun, cells_all, period =c("winter_start", "winter_end"), n_small_min = params$n_small_min, quiet=params$quiet)
  tidy_ratios <- make_ratios_tidy(cell_ratios)
  
  # rename fall/spring ratios to productivity/recruitment:
  tidy_ratios$summary %>%
    mutate(season=ifelse(period=="winter_end", "surv","prod")) -> tidy_ratios$summary
  
  # save the ratio data
  saveRDS(tidy_ratios, paste0(params$output_path, "/abun_data/", species_code, "_winter_bookend_ratios.rds"))
}

####
#### Plot demographic indices ------------------------------------------------------
####

species_code <- c("easpho")
tidy_ratios <- readRDS(paste0(params$output_path, "/abun_data/", species_code, "_summer_bookend_ratios.rds"))

# select cells with at least 20 valid ratios
tidy_ratios$summary %>%
  group_by(cell) %>%
  summarize(sufficient_data = sum(is.na(median))<20) %>%
  filter(sufficient_data) %>% pull(cell) -> cells_select

# plot selected cell 20 (example)
plot_ratios(cells_select[20], data=tidy_ratios$summary)
length(cells_select) 

# Plotting demographic indices for all cells 
library(gridExtra)
# Initialize a list to hold the plots
plots <- list()
# Loop through each cell
for (i in seq_along(cells_select)) {
  # Generate a plot for the current cell
  p <- plot_ratios(cells_select[i], data=tidy_ratios$summary)+theme(legend.position = "none")
  # Add the plot to the list of plots
  plots[[length(plots) + 1]] <- p
}

# Arrange the plots in a 3x6 grid and display
p <- grid.arrange(grobs = plots, ncol = 11, nrow=13)

# ggsave("~/Documents/plots/amerob/amerob_surv_prod_flact_summer_3weeksWindow.png", plot=p, width = 80, height= 40, units = "cm")  


####
#### Verify normality assumptions ------------------------------------------------------
####

# Check the higher moments of the abundance log-ratios to verify adequacy of Gaussian approximations.
for(species_code in params$species_to_process){
  tidy_ratios <- readRDS(paste0(params$output_path, "/abun_data/", species_code, "_summer_bookend_ratios.rds"))
  # excess kurtosis = kurtosis - 3
  ggplot(tidy_ratios$summary, aes(x=kurtosis-3)) +
    geom_histogram(binwidth=.1) +
    xlab("excess kurtosis") +
    ggtitle(species_code)
  # skewness
  ggplot(tidy_ratios$summary, aes(x=skewness)) +
    geom_histogram(binwidth=.1) +
    xlab("skewness") +
    ggtitle(species_code)
}

####
#### Compare recruitment vs mortality variances  ------------------------------------------------------
####

for(species_code in params$species_to_process){
  tidy_ratios <- readRDS(paste0(params$output_path, "/abun_data/", species_code, "_summer_bookend_ratios.rds"))
  
  # compare variances in productivity and survival for each cell across years:
  var_save_path <- paste0(params$output_path, "/variance_results/", species_code,"/_summer_bookend_variance_test.rds")
  dir.create(dirname(var_save_path), recursive = TRUE, showWarnings = FALSE)
  if(params$always_run_variance_test | !file.exists(var_save_path)){
    print(paste0("processing_species ", species_code))
    data_compare_ratios <- lapply(sort(unique(tidy_ratios$summary$cell)),
                                  compare_ratio_variances, data=tidy_ratios$summary,
                                  n_ratio_min=params$n_year_min,
                                  quiet=params$quiet) %>% do.call(rbind, .)
    saveRDS(data_compare_ratios,file=var_save_path)
  } else{
    data_compare_ratios <- readRDS(file=var_save_path)
  }
}



####
#### Plot variances  ------------------------------------------------------
####

# helper function to plot cells on a map
plot_cells_on_map <- function(data, param, color_scale){
  ggplot() + blank_theme + coord_fixed() +
    geom_sf(data=region_of_interest, fill=NA, color="black") +
    geom_polygon(data=data, aes(x=long, y=lat, group=group, fill=eval(parse(text=param))), alpha=0.7) +
    geom_path(data=data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
    color_scale +
    labs(fill=param) +
    xlim(params$plotting_xlim)+
    labs(subtitle = species_code)
}

# Initialize a lists to hold the plots
plot_list <- list()

for(species_code in params$species_to_process){
  
  var_save_path <- paste0(params$output_path, "/variance_results/", species_code,"/_summer_bookend_variance_test.rds")
  cell_lrat_sd <- readRDS(var_save_path)
  
  tidy_ratios <- readRDS(paste0(params$output_path, "/abun_data/", species_code, "_summer_bookend_ratios.rds"))
  
  tidy_ratios$summary %>%
    select(cell, n_prod,n_surv,n) %>%
    group_by(cell) %>%
    filter(row_number()==1) %>%
    right_join(cell_lrat_sd, by="cell") %>%
    filter(n_prod >= params$n_year_min & n_surv >= params$n_year_min) %>%
    filter(!is.na(p_survival_variance_higher)) -> plotting_data
  
  dggridR::dgcellstogrid(
    grid_large,
    plotting_data$cell,
    frame=TRUE,
    wrapcells=TRUE
  ) %>%
    mutate(cell=as.numeric(cell)) %>%
    left_join(plotting_data, by="cell") -> grid2
  
  # define color scales
  scale_viridis <- viridis::scale_fill_viridis(limits = c(params$n_year_min, length(params$years)), 
                                               breaks=seq(params$n_year_min, length(params$years), 3))
  scale_blue_red <- scale_fill_gradientn(colours = cols_bd2, na.value=NA, limits = c(0, 1))
  fl <- max(abs(grid2$effect_size_log), na.rm = T) + .1
  

  # plot number of years with a productivity index
  p_nprod <- plot_cells_on_map(grid2, "n_prod", color_scale=scale_viridis)
  
  # plot number of years with a survival index
  p_nsurv <- plot_cells_on_map(grid2, "n_surv", color_scale=scale_viridis)
  
  # plot probability that survival variance is higher than productivity
  p_pval <- plot_cells_on_map(grid2, "p_survival_variance_higher", color_scale=scale_blue_red)+labs(fill="p-value")
  
  # plot the difference in survival and recruitment variance (log-effect size)
  # scale opacity by the probability that the survival variance is higher.
  p_slope_pval <- 
    ggplot() + coord_fixed() + blank_theme +
    geom_sf(data=region_of_interest, fill=NA, color="black") +
    geom_polygon(data=grid2, aes(x=long, y=lat, group=group, fill = effect_size_log), alpha = 2*abs(grid2$p_survival_variance_higher - 0.5))   +
    geom_path(data=grid2, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
    scale_fill_gradientn(colours = cols_bd, na.value=NA, limits = c(-6, 6), oob=scales::squish, name="slope") + # I set limit to 6, as it is most common
    xlim(params$plotting_xlim)+labs(subtitle = species_code)
  
  
  # Add the plots to the list
  plot_list[[species_code]] <- list(nprod = p_nprod, nsurv = p_nsurv, pval = p_pval, slope = p_slope_pval)
  
}
 
# Arrange and save "n_year prod" plot panels
nprod_panel <- gridExtra::grid.arrange(grobs = lapply(plot_list, function(x) x$nprod), ncol=3, nrow=2)
ggsave("~/Documents/plots/sp_combined/nyears_prod_panel.png", plot=nprod_panel, width = 40, height= 20, units = "cm")  

# Arrange and save" n_year surv" plot panel
nsurv_panel <- gridExtra::grid.arrange(grobs = lapply(plot_list, function(x) x$nsurv), ncol=3, nrow=2)
ggsave("~/Documents/plots/sp_combined/nyears_surv_panel.png", plot=nsurv_panel, width = 40, height= 20, units = "cm")  
  
# Arrange and save "p_val" plot panel 
pval_panel <- gridExtra::grid.arrange(grobs = lapply(plot_list, function(x) x$pval), ncol = 3, nrow = 2)
ggsave("~/Documents/plots/sp_combined/pval_survProd_panel.png", plot=pval_panel, width = 40, height= 20, units = "cm")  
  
# Arrange and save "slope and pval" plot panel 
slope_panel <- gridExtra::grid.arrange(grobs=lapply(plot_list, function(x) x$slope), ncol=3, nrow=2)
ggsave("~/Documents/plots/sp_combined/slope_pval_survProd_panel.png", plot=slope_panel, width = 40, height= 20, units = "cm")  
  

###
###     Smoothing recruitment and mortality variances  ----------------------------------------
###

icar_regression <- function(data, adjancency_mat, warmup=2000, cores=4, iter=12000){
  assert_that(nrow(adjacency_mat) == ncol(adjacency_mat))
  assert_that(nrow(data) == nrow(adjacency_mat))
  # ICAR model without latitude
  print("Estimating icar model ...")
  icar_fit <-
    brm(slope_scaled | resp_se(known_se, sigma = TRUE) ~
          car(M, gr = cell, type = "icar"),
        data = data,
        data2 = list(M = adjacency_mat),
        backend = 'cmdstanr',
        iter = iter, warmup = warmup, cores = cores, refresh = 0)
  # number of cells:
  npt <- nrow(data)
  # number of posterior draws:
  n_samples = cores*(iter-warmup)
  # initialize matrix to contain posterior draws for each cell:
  true_values <- matrix(nrow = n_samples, ncol = npt)
  
  print("Drawing from posterior ...")
  for(i in 1:npt) {
    print(paste("cell",i,"/",npt,"..."))
    # the (unsmoothed) posterior slope:
    M <- data$slope_scaled[i]
    # the posterior uncertainty from the weather regression:
    se <- data$known_se[i]
    # get posterior draws of the linear predictor
    # re.form=NULL (default), i.e. include all group-level effects
    # incl_autor=TRUE (default), i.e. include correlation structures in prediction
    LP <-  posterior_linpred(icar_fit, re.form = NULL, incl_autocor = T)[ , i]
    # get draws of the unobserved latent uncertainty:
    sigma <- as_draws_df(icar_fit)$sigma
    # error propagation, combining sigma form CAR model and se from weather regression
    # QUESTION: why the linear predictor here
    mu = (M*sigma^2 + LP * se^2)/(sigma^2 + se^2)
    scale = 1/sqrt(1/sigma^2 + 1/se^2)
    true_values[ , i] <- rnorm(n_samples, mu, scale)
  }
  unscale <- function(x){x + mean(data$mean)/sd(data$mean)}
  smooth_prob=apply(true_values, 2,function(x){mean(unscale(x) > 0)})
  smooth_mean=apply(true_values, 2,function(x){mean(unscale(x))})
  smooth_data=data.frame(cell=data$cell, smooth_prob=smooth_prob, smooth_mean=smooth_mean)
  output_data <- merge(data,smooth_data, by = "cell")
  output_data
}

# plots the mean probability for a nonzero effect
plot_smooth_prob <- function(grid_data, region_of_interest){
  ggplot() + coord_fixed() + blank_theme +
    geom_sf(data=region_of_interest, fill=NA, color="black")   +
    geom_polygon(data=grid_data, aes(x=long, y=lat.x, group=group, fill = smooth_prob), alpha = .8)   +
    geom_path(data=grid_data, aes(x=long, y=lat.x, group=group), alpha=0.4, color="white") +
    scale_fill_gradientn(colours = cols_bd2, na.value=NA, limits = c(0,1), oob=scales::squish) +
    xlim(params$plotting_xlim)
}

# plots the mean smoothed effect size, with opacity according to probability
plot_smooth_mean <- function(grid_data, region_of_interest){
#  fl <- max(abs(grid_data$smooth_mean), na.rm = T) + .1
  fl <- 1  # BG set it as 1 to make it comparable across species. Eastern phoebe's was max 0.7! 
  ggplot() + coord_fixed() + blank_theme +
    geom_sf(data=region_of_interest, fill=NA, color="black")   +
    geom_polygon(data=grid_data, aes(x=long, y=lat.x, group=group, fill = smooth_mean), alpha = 2*abs(grid_data$smooth_prob - 0.5))   +
    geom_path(data=grid_data, aes(x=long, y=lat.x, group=group), alpha=0.4, color="white") +
    scale_fill_gradientn(colours = cols_bd, na.value=NA, limits = c(-fl, fl), oob=scales::squish) +
    xlim(params$plotting_xlim)
}

merge_grid_to_data <- function(data, grid_large){
  grid <- dggridR::dgcellstogrid(grid_large,data$cell,frame=TRUE,wrapcells=TRUE)
  grid <- merge(grid,data,by=c("cell"))
  # remove smoothed data equal to NA
  grid[!is.na(grid$smooth_prob), ]
}

# Initialize a lists to hold the plots
plot_list <- list()

for(species_code in params$species_to_process) {
  
  var_save_path <- paste0(params$output_path, "/variance_results/", species_code,"/_summer_bookend_variance_test.rds")
  cell_lrat_sd <- readRDS(var_save_path)
  
  tidy_ratios <- readRDS(paste0(params$output_path, "/abun_data/", species_code, "_summer_bookend_ratios.rds"))
  
  tidy_ratios$summary %>%
    select(cell, n_prod,n_surv,n) %>%
    group_by(cell) %>%
    filter(row_number()==1) %>%
    right_join(cell_lrat_sd, by="cell") %>%
    filter(n_prod >= params$n_year_min & n_surv >= params$n_year_min) %>%
    filter(!is.na(p_survival_variance_higher)) -> plotting_data
  
  car_var_data <- 
    plotting_data %>% ungroup() %>% 
    mutate(lat = dggridR::dgSEQNUM_to_GEO(grid_large, cell)$lat_deg,
           lon = dggridR::dgSEQNUM_to_GEO(grid_large, cell)$lon_deg) %>%
    filter(lon > params$plotting_xlim[1]) %>%
    # remove a region where we have sparse data
    #  filter(!(lat > 38.7 & lon < -99.5)) %>%
    mutate(mean=effect_size_log,
           slope_scaled=scale(effect_size_log),
           known_se = effect_size_sd/sd(effect_size_log))
  
  adjacency_mat <- macrodemography:::get_adjacency_matrix(car_var_data$cell, grid_large)
  icar_smooth_var <- icar_regression(car_var_data, adjacency_mat)
  grid_data <- merge_grid_to_data(icar_smooth_var, grid_large)
  p1 <- plot_smooth_prob(grid_data = grid_data, region_of_interest = region_of_interest)+labs(fill="smoothed\np-value", subtitle = species_code)
  p2 <- plot_smooth_mean(grid_data = grid_data, region_of_interest = region_of_interest)+labs(fill="smoothed\nslope", subtitle = species_code)
  p <- gridExtra::grid.arrange(p1, p2, ncol=2)
  
  # Add the plots to the list
  plot_list[[species_code]] <- list(plot_prob = p1, plot_mean = p2)
  
}

# Arrange and save "smoothed prob" plot panels
smooth_pval_panel <- gridExtra::grid.arrange(grobs = lapply(plot_list, function(x) x$plot_prob), ncol=3, nrow=2)
ggsave("~/Documents/plots/sp_combined/smoothed_prob_survProd_first6sp_panel.png", plot=smooth_pval_panel, width = 40, height= 20, units = "cm")  

# Arrange and save "smoothed slope" plot panels
smooth_slope_panel <- gridExtra::grid.arrange(grobs = lapply(plot_list, function(x) x$plot_mean), ncol=3, nrow=2)
ggsave("~/Documents/plots/sp_combined/smoothed_slope_survProd_first6sp_panel.png", plot=smooth_slope_panel, width = 40, height= 20, units = "cm")  


####
#### Compare Bayesian/frequentest averages  ------------------------------------------------------
####

# Quick comparisons of inverse-sd weighted averages and the Bayesian estimates.
# Using the last generated `tidy_ratios` object for now.

# average indices over cells across years (using inverse weighting by sd)
# (note that n_prod,n_surv and n are identical within the grouping, so
# the weighted.mean just repeats the group value)
tidy_ratios$summary %>%
  group_by(cell,season) %>%
  filter(is.finite(avg)) %>%
  summarise(across(c("median","avg","q10","q90","sd","n_prod","n_surv","n", "has_inf"), \(x) weighted.mean(x, 1/sd, na.rm=TRUE)), .groups="drop_last") %>%
  mutate(has_inf=as.logical(has_inf)) -> data_cell

# join variance test results to year-averaged cell data:
data_cell <- left_join(data_cell,data_compare_ratios,by="cell")


####
#### Load weather data  ------------------------------------------------------
####
# 
# dir.create(paste0(params$output_path, "/weather"), showWarnings = FALSE)
# weather_file <- paste0(params$output_path, "/weather/",paste0(params$daymet$label, collapse = "-"),".rds")
# flag_file <- paste0(params$output_path, "/weather/flag.rds")
# 
# if(!params$quiet) print(paste("loading/processing weather file",basename(weather_file)))
# 
# if (params$always_download_weather | !file.exists(weather_file) | !file.exists(flag_file)) {
#   # initialize google earth engine
#   ee_Initialize()
# 
#   # load or create flag variable to indicate whether data extraction was successful or not
#   flag <- if (file.exists(flag_file)) readRDS(flag_file) else data.frame(cell = character(), year = integer(), success = logical())
# 
#   # check if weather file exists and load data
#   data_daymet <- if (file.exists(weather_file)) readRDS(weather_file) else data.frame()
# 
#   # loop over cells and years
#   for (cell in cells_all) {
#     for (year in params$years) {
#       # check if data extraction was successful for this cell and year
#       if (nrow(flag) == 0 || !any(flag$success[flag$cell == cell & flag$year == year])) {
#         # attempt data extraction
#         tryCatch({
#           print(paste("year =", year, "cell =", cell))
#           data_daymet <- rbind(data_daymet, daymet_set_extract(year, cell, grid_large, params$daymet))
# 
#           # update flag variable to indicate success
#           flag <- rbind(flag, data.frame(cell = cell, year = year, success = TRUE))
#           saveRDS(flag, flag_file)
# 
#           # save weather data after each iteration of the inner loop
#           saveRDS(data_daymet, weather_file)
#         }, error = function(e) {
#           # if an error occurred during data extraction, print an error message
#           message(paste("Error occurred for year =", year, "cell =", cell))
#           # update flag variable to indicate failure
#           flag <- rbind(flag, data.frame(cell = cell, year = year, success = FALSE))
#           saveRDS(flag, flag_file)
#         })
#       }
#     }
#   }
# } else {
#   # load weather data if always_download_weather is false and the weather file exists
#   data_daymet <- readRDS(weather_file)
# }


####
#### Weather regressions ------------------------------------------------------
####

# Separate regressions respective to survival (winter_bookended) and recruitment (summer_bookended)

# 1. Regressions for recruitment vs. maximum summer temperature (ERA data)

# ERA5 dataset
era_data <- read.csv("/Users/bg423/Documents/macrodemography/data_parquet/weather/ERA5_data_1964_2023.csv")

# filtering weather sampling time-window 
era_data <- era_data %>% filter(month(as.Date(date))==6 & day(as.Date(date)) >=15 |
                                month(as.Date(date))==7 & day(as.Date(date)) <=15 |
                                month(as.Date(date))==1) %>% mutate(season=ifelse(month(as.Date(date))==1, "winter", "summer"))

# Summer temperature will be used in regression with recruitment  
era_summer <- era_data %>% filter(season=="summer") %>% group_by(cell, year=year(as.Date(date))) %>% 
  summarise(tmax=mean(temperature_2m_max, na.rm = TRUE))

params$period <- c("spring", "fall")
params$era <- data.frame(label=c("tmax"),
                         variable=c("tmax"),
                            date_min=c("06-15"),
                            date_max=c("07-15"),
                            period=c("fall"))

params$always_run_regressions    <- TRUE

for(species_code in params$species_to_process){
  regression_save_path <- paste0(params$output_path, "/regression_results/", species_code)
  dir.create(regression_save_path, recursive = TRUE, showWarnings = FALSE)

  if(params$always_run_regressions | !file.exists(paste0(regression_save_path, "/summer_ERA_regressions.rds"))){
    print(paste0("processing_species ", species_code))

    file_path <- paste0(params$output_path, "/abun_data/", species_code, "_summer_bookend_ratios.rds")
    tidy_ratios <- readRDS(file_path)

    data_regression <- weather_regressions(tidy_ratios, data_daymet=era_summer,
                                           params_daymet=params$era, params$n_year_min, quiet=FALSE)

    saveRDS(data_regression, paste0(regression_save_path, "/summer_ERA_regressions.rds"))
  }
}

# 2. Regressions for overwintering survival vs. maximum winter temperature, snow depth, and their anomalies (ERA data)

# Snow and winter temperature will be used in regression with survival
era_winter <- era_data %>% filter(season=="winter") %>% group_by(cell, year=year(as.Date(date))) %>% 
  summarise(tmax=mean(temperature_2m_max, na.rm = TRUE), 
            snow_depth=mean(snow_depth, na.rm = TRUE))

params$period <- c("winter_start", "winter_end")

params$era <- data.frame(label=c("tmax", "snow_depth"),
                         variable=c("tmax", "snow_depth"),
                         date_min=c("01-01", "01-01"),
                         date_max=c("01-31", "01-31"),
                         period=c("winter_end", "winter_end"))

params$always_run_regressions    <- TRUE

for(species_code in params$species_to_process){
  regression_save_path <- paste0(params$output_path, "/regression_results/", species_code)
  dir.create(regression_save_path, recursive = TRUE, showWarnings = FALSE)

  if(params$always_run_regressions | !file.exists(paste0(regression_save_path, "/winter_ERA_regressions.rds"))){
    print(paste0("processing_species ", species_code))

    file_path <- paste0(params$output_path, "/abun_data/", species_code, "_winter_bookend_ratios.rds")
    tidy_ratios <- readRDS(file_path)

    data_regression <- weather_regressions(tidy_ratios, data_daymet=era_winter,
                                           params_daymet=params$era, params$n_year_min, quiet=FALSE)

    saveRDS(data_regression, paste0(regression_save_path, "/winter_ERA_regressions.rds"))
  }
}


####
#### Plot weather regressions------------------------------------------------------
####

plot_list <- list()

for(species_code in params$species_to_process) {

# load the data for that species:
regression_save_path <- paste0(params$output_path, "/regression_results/", species_code)
winter_regression <- readRDS(paste0(regression_save_path, "/winter_ERA_regressions.rds"))
summer_regression <- readRDS(paste0(regression_save_path, "/summer_ERA_regressions.rds"))
data_regression   <- rbind(winter_regression, summer_regression) %>% 
  mutate(label=ifelse(period=="winter_end", paste("surv", label, sep = "_"), paste("prod", label, sep = "_")))

# plot function for regression data:
plot_regression <- function(data, label_daymet, moment, x_lim, fill_lim="auto", alpha="auto", labels=FALSE, colors=cols_bd) {
  assert_that(moment %in% names(data),msg=paste("column",variable,"not found in data"))
  assert_that(label_daymet %in% unique(data_regression$label),msg=paste("daymet variable",variable,"not found in data"))
  data <- filter(data, label==label_daymet)
  # get the spatial polygons for the cells
  grid <- dggridR::dgcellstogrid(grid_large,data$cell,frame=TRUE,wrapcells=TRUE)
  # merge grid and regression data
  grid_data <- merge(grid,data,by="cell")
  # calculate average lat,long for cells:
  grid %>%
    group_by(cell) %>%
    summarize(long=mean(long), lat=mean(lat)) -> data_label
  # define opacity based on value of alpha parameter:
  if(is.number(alpha)) grid_data$alpha=alpha
  # when alpha is NULL, color according to the certainty that the coefficient is either positively or
  # negatively different from zero
  if(alpha=="auto") grid_data$alpha=2*abs(grid_data$p_value-0.5)
  
  if(identical(fill_lim,"auto")){
    max_moment = max(abs(grid_data[[moment]]), na.rm = T) + 0.1
    fill_lim = c(-max_moment,max_moment)
  }
  p <- ggplot() + coord_fixed() + blank_theme +
    geom_sf(data=region_of_interest, fill=NA, color="black")   +
    geom_polygon(data=grid_data, aes(x=long, y=lat, group=group, fill = !!sym(moment)), alpha=grid_data$alpha) +
    geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
    scale_fill_gradientn(colours = colors, na.value=NA, limits=fill_lim, oob=scales::squish) + xlim(x_lim) # +labs(subtitle = species_code)
  if(labels) p = p + geom_text(aes(x=long,y=lat,label=cell), data=data_label)
  print(p)
}

# plot excess kurtosis and skewness of the regression coefficient posterior:
data_regression %>%
  mutate(excess_kurtosis=kurtosis-3) %>%
  tidyr::pivot_longer(c(skewness, excess_kurtosis)) %>%
  group_by(label,name) %>%
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(c("name","label"))

# plot number of years of data on which the regression was based
# also add the cell number labels
surv_n <- plot_regression(data_regression, "surv_tmax", "n", params$plotting_xlim, fill_lim=c(0,18), alpha=.8, labels=FALSE, colors=viridisLite::viridis(100))
prod_n <- plot_regression(data_regression, "prod_tmax", "n", params$plotting_xlim, fill_lim=c(0,18), alpha=.8, labels=FALSE, colors=viridisLite::viridis(100))

# plot the mean regression slopes
snow_sl    <- plot_regression(data_regression, "surv_snow_depth", "mean", params$plotting_xlim, alpha=.8)
winterT_sl <- plot_regression(data_regression, "surv_tmax", "mean", params$plotting_xlim, alpha=.8)
summerT_sl <- plot_regression(data_regression, "prod_tmax", "mean", params$plotting_xlim, alpha=.8)

# plot the posterior probability that the regression coefficient is larger than zero (p_value)
snow_pval    <- plot_regression(data_regression, "surv_snow_depth", "p_value", params$plotting_xlim, alpha=.8, fill_lim=c(0,1))
winterT_pval <- plot_regression(data_regression, "surv_tmax", "p_value", params$plotting_xlim, alpha=.8, fill_lim=c(0,1))
summerT_pval <- plot_regression(data_regression, "prod_tmax", "p_value", params$plotting_xlim, alpha=.8, fill_lim=c(0,1))

# plot the mean regression slopes, with transparency denoting the probability of being non-zero
# essentially a combining the previous two types of plots in one figure
snow_sl_pval    <- plot_regression(data_regression, "surv_snow_depth", "mean", params$plotting_xlim, fill_lim = c(-1,1))
winterT_sl_pval <- plot_regression(data_regression, "surv_tmax", "mean", params$plotting_xlim, fill_lim = c(-1,1))
summerT_sl_pval <- plot_regression(data_regression, "prod_tmax", "mean", params$plotting_xlim, fill_lim = c(-1,1))

# plot the kurtosis (note for Guassian distribution the kurtosis = 3)
snow_kurt    <- plot_regression(data_regression, "surv_snow_depth", "kurtosis", params$plotting_xlim, fill_lim=c(3,9))
winterT_kurt <- plot_regression(data_regression, "surv_tmax", "kurtosis", params$plotting_xlim, fill_lim=c(3,9))
summerT_kurt <- plot_regression(data_regression, "prod_tmax", "kurtosis", params$plotting_xlim, fill_lim=c(3,9))

# plot the skewness
snow_skew    <- plot_regression(data_regression, "surv_snow_depth", "skewness", params$plotting_xlim)
winterT_skew <- plot_regression(data_regression, "surv_tmax", "skewness", params$plotting_xlim)
summerT_skew <- plot_regression(data_regression, "prod_tmax", "skewness", params$plotting_xlim)

# Add the plots to the list
plot_list[[species_code]] <- list(surv_n=surv_n, prod_n=prod_n, snow_sl=snow_sl, winterT_sl=winterT_sl, 
                                  summerT_sl=summerT_sl, snow_pval=snow_pval, winterT_pval=winterT_pval, 
                                  summerT_pval=summerT_pval, snow_sl_pval=snow_sl_pval, winterT_sl_pval=winterT_sl_pval, 
                                  summerT_sl_pval=summerT_sl_pval, snow_kurt=snow_kurt, winterT_kurt=winterT_kurt, 
                                  summerT_kurt=summerT_kurt, snow_skew=snow_skew, winterT_skew=winterT_skew, 
                                  summerT_skew=summerT_skew)
}

saveRDS(plot_list, paste0(params$output_path, "/plots/", "/weather_regression_plots.rds"))


# Arrange and save "surv_n" in plot panels
surv_n <- gridExtra::grid.arrange(grobs = lapply(plot_list, function(x) x$surv_n), ncol=3, nrow=2)
ggsave("~/Documents/plots/sp_combined/nyears_surv_first6sp_panel.png", plot=surv_n, width = 40, height= 20, units = "cm")  

# Arrange and save "winterT_sl_pval" in plot panels
winterT_sl_pval <- gridExtra::grid.arrange(grobs = lapply(plot_list, function(x) x$winterT_sl_pval), ncol=3, nrow=2)
ggsave("~/Documents/plots/sp_combined/ERA_winterT_sl_pval_first6sp_panel.png", plot=winterT_sl_pval, width = 40, height= 20, units = "cm")  

# Arrange and save "snow_sl_pval" in plot panels
snow_sl_pval <- gridExtra::grid.arrange(grobs = lapply(plot_list, function(x) x$snow_sl_pval), ncol=3, nrow=2)
ggsave("~/Documents/plots/sp_combined/ERA_snow_sl_pval_first6sp_panel.png", plot=snow_sl_pval, width = 40, height= 20, units = "cm")  

# Arrange and save "summerT_sl_pval" in plot panels
summerT_sl_pval <- gridExtra::grid.arrange(grobs = lapply(plot_list, function(x) x$summerT_sl_pval), ncol=3, nrow=2)
ggsave("~/Documents/plots/sp_combined/ERA_summerT_sl_pval_first6sp_panel.png", plot=summerT_sl_pval, width = 40, height= 20, units = "cm")  


####
#### Smoothing using ICAR/CAR models: smoothed maps----------------------------------------------
####

icar_regression <- function(data, adjancency_mat, warmup=2000, cores=4, iter=12000){
  assert_that(nrow(adjacency_mat) == ncol(adjacency_mat))
  assert_that(nrow(data) == nrow(adjacency_mat))
  # ICAR model without latitude
  print("Estimating icar model ...")
  icar_fit <-
    brm(slope_scaled | resp_se(known_se, sigma = TRUE) ~
          car(M, gr = cell, type = "icar"),
        data = data,
        data2 = list(M = adjacency_mat),
        backend = 'cmdstanr',
        iter = iter, warmup = warmup, cores = cores, refresh = 0)
  # number of cells:
  npt <- nrow(data)
  # number of posterior draws:
  n_samples = cores*(iter-warmup)
  # initialize matrix to contain posterior draws for each cell:
  true_values <- matrix(nrow = n_samples, ncol = npt)
  
  print("Drawing from posterior ...")
  for(i in 1:npt) {
    print(paste("cell",i,"/",npt,"..."))
    # the (unsmoothed) posterior slope:
    M <- data$slope_scaled[i]
    # the posterior uncertainty from the weather regression:
    se <- data$known_se[i]
    # get posterior draws of the linear predictor
    # re.form=NULL (default), i.e. include all group-level effects
    # incl_autor=TRUE (default), i.e. include correlation structures in prediction
    LP <-  posterior_linpred(icar_fit, re.form = NULL, incl_autocor = T)[ , i]
    # get draws of the unobserved latent uncertainty:
    sigma <- as_draws_df(icar_fit)$sigma
    # error propagation, combining sigma form CAR model and se from weather regression
    # QUESTION: why the linear predictor here
    mu = (M*sigma^2 + LP * se^2)/(sigma^2 + se^2)
    scale = 1/sqrt(1/sigma^2 + 1/se^2)
    true_values[ , i] <- rnorm(n_samples, mu, scale)
  }
  unscale <- function(x){x + mean(data$mean)/sd(data$mean)}
  smooth_prob=apply(true_values, 2,function(x){mean(unscale(x) > 0)})
  smooth_mean=apply(true_values, 2,function(x){mean(unscale(x))})
  smooth_data=data.frame(cell=data$cell, smooth_prob=smooth_prob, smooth_mean=smooth_mean)
  output_data <- merge(data,smooth_data, by = "cell")
  output_data
}

# plots the mean probability for a nonzero effect
plot_smooth_prob <- function(grid_data, region_of_interest){
  ggplot() + coord_fixed() + blank_theme +
    geom_sf(data=region_of_interest, fill=NA, color="black")   +
    geom_polygon(data=grid_data, aes(x=long, y=lat.x, group=group, fill = smooth_prob), alpha = .8)   +
    geom_path(data=grid_data, aes(x=long, y=lat.x, group=group), alpha=0.4, color="white") +
    scale_fill_gradientn(colours = cols_bd2, na.value=NA, limits = c(0,1), oob=scales::squish) +
    xlim(params$plotting_xlim)
}

# plots the mean smoothed effect size, with opacity according to probability
plot_smooth_mean <- function(grid_data, region_of_interest){
  fl <- max(abs(grid_data$smooth_mean), na.rm = T) + .1
  ggplot() + coord_fixed() + blank_theme +
    geom_sf(data=region_of_interest, fill=NA, color="black")   +
    geom_polygon(data=grid_data, aes(x=long, y=lat.x, group=group, fill = smooth_mean), alpha = 2*abs(grid_data$smooth_prob - 0.5))   +
    geom_path(data=grid_data, aes(x=long, y=lat.x, group=group), alpha=0.4, color="white") +
    scale_fill_gradientn(colours = cols_bd, na.value=NA, limits = c(-fl, fl), oob=scales::squish) +
    xlim(params$plotting_xlim)
}

merge_grid_to_data <- function(data, grid_large){
  grid <- dggridR::dgcellstogrid(grid_large,data$cell,frame=TRUE,wrapcells=TRUE)
  grid <- merge(grid,data,by=c("cell"))
  # remove smoothed data equal to NA
  grid[!is.na(grid$smooth_prob), ]
}


plot_list <- list()

for(species_code in params$species_to_process) {

  # load the data for that species:
  regression_save_path <- paste0(params$output_path, "/regression_results/", species_code)
  winter_regression <- readRDS(paste0(regression_save_path, "/winter_ERA_regressions.rds"))
  summer_regression <- readRDS(paste0(regression_save_path, "/summer_ERA_regressions.rds"))
  data_regression   <- rbind(winter_regression, summer_regression) %>% 
    mutate(label=ifelse(period=="winter_end", paste("surv", label, sep = "_"), paste("prod", label, sep = "_")))
  
  
# prepare CAR model data from weather regression data
data_regression %>%
  mutate(lat = dggridR::dgSEQNUM_to_GEO(grid_large, cell)$lat_deg,
         lon = dggridR::dgSEQNUM_to_GEO(grid_large, cell)$lon_deg) %>%
  filter(lon > params$plotting_xlim[1]) %>%
  # remove a region where we have sparse data
#  filter(!(lat > 38.7 & lon < -99.5)) %>%
  filter(!is.na(mean)) %>%
  group_by(label) %>%
  # rescale the data (to help with model convergence)
  mutate(slope_scaled = scale(`mean`),
         lat_scaled = scale(lat),
         known_se = `sd`/sd(`mean`)) -> car_data

# tmax_winter
car_data %>% filter(label=="surv_tmax") -> car_data_select
adjacency_mat <- macrodemography:::get_adjacency_matrix(car_data_select$cell, grid_large)
icar_smooth_tmax_winter <- icar_regression(car_data_select, adjacency_mat)
grid_data <- merge_grid_to_data(icar_smooth_tmax_winter, grid_large)
wintT_pval_sm    <- plot_smooth_prob(grid_data = grid_data, region_of_interest = region_of_interest)
wintT_sl_pval_sm <- plot_smooth_mean(grid_data = grid_data, region_of_interest = region_of_interest)

# tmax_summer
car_data %>% filter(label=="prod_tmax") -> car_data_select
adjacency_mat <- macrodemography:::get_adjacency_matrix(car_data_select$cell, grid_large)
icar_smooth_tmax_summer <- icar_regression(car_data_select, adjacency_mat)
grid_data <- merge_grid_to_data(icar_smooth_tmax_summer, grid_large)
summerT_pval_sm    <- plot_smooth_prob(grid_data = grid_data, region_of_interest = region_of_interest)
summerT_sl_pval_sm <- plot_smooth_mean(grid_data = grid_data, region_of_interest = region_of_interest)

# swe
car_data %>% filter(label=="surv_snow_depth") -> car_data_select
adjacency_mat <- macrodemography:::get_adjacency_matrix(car_data_select$cell, grid_large)
icar_smooth_swe <- icar_regression(car_data_select, adjacency_mat)
grid_data <- merge_grid_to_data(icar_smooth_swe, grid_large)
snow_pval_sm    <- plot_smooth_prob(grid_data = grid_data, region_of_interest = region_of_interest)
snow_sl_pval_sm <- plot_smooth_mean(grid_data = grid_data, region_of_interest = region_of_interest)

plot_list[[species_code]] <- list(wintT_pval_sm=wintT_pval_sm, wintT_sl_pval_sm=wintT_sl_pval_sm,
                                  summerT_pval_sm=summerT_pval_sm, summerT_sl_pval_sm=summerT_sl_pval_sm,
                                  snow_pval_sm=snow_pval_sm, snow_sl_pval_sm=snow_sl_pval_sm)
}

saveRDS(plot_list, "ERA_weather_regressions_plot_list.rds")

# Arrange and save "winterT_pval smoothed" in plot panels
wintT_pval_sm <- gridExtra::grid.arrange(grobs = lapply(plot_list, function(x) x$wintT_pval_sm), ncol=3, nrow=2)
ggsave("~/Documents/plots/sp_combined/ERA_wintT_pval_first6sp_panel_smoothed.png", plot=wintT_pval_sm, width = 40, height= 20, units = "cm")  

# Arrange and save "winterT_sl_pval smoothed" in plot panels
wintT_sl_pval_sm <- gridExtra::grid.arrange(grobs = lapply(plot_list, function(x) x$wintT_sl_pval_sm), ncol=3, nrow=2)
ggsave("~/Documents/plots/sp_combined/ERA_wintT_sl_pval_first6sp_panel_smoothed.png", plot=wintT_sl_pval_sm, width = 40, height= 20, units = "cm")  

# Arrange and save "snow_pval smoothed" in plot panels
snow_pval_sm <- gridExtra::grid.arrange(grobs = lapply(plot_list, function(x) x$snow_pval_sm), ncol=3, nrow=2)
ggsave("~/Documents/plots/sp_combined/ERA_snow_pval_first6sp_panel_smoothed.png", plot=snow_pval_sm, width = 40, height= 20, units = "cm")  

# Arrange and save "snow_sl_pval smoothed" in plot panels
snow_sl_pval_sm <- gridExtra::grid.arrange(grobs = lapply(plot_list, function(x) x$snow_sl_pval_sm), ncol=3, nrow=2)
ggsave("~/Documents/plots/sp_combined/ERA_snow_sl_pval_first6sp_panel_smoothed.png", plot=snow_sl_pval_sm, width = 40, height= 20, units = "cm")  

# Arrange and save "summerT_pval smoothed" in plot panels
summerT_pval_sm <- gridExtra::grid.arrange(grobs = lapply(plot_list, function(x) x$summerT_pval_sm), ncol=3, nrow=2)
ggsave("~/Documents/plots/sp_combined/ERA_summerT_pval_first6sp_panel_smoothed.png", plot=summerT_pval_sm, width = 40, height= 20, units = "cm")  

# Arrange and save "summerT_sl_pval smoothed" in plot panels
summerT_sl_pval_sm <- gridExtra::grid.arrange(grobs = lapply(plot_list, function(x) x$summerT_sl_pval_sm), ncol=3, nrow=2)
ggsave("~/Documents/plots/sp_combined/ERA_summerT_sl_pval_first6sp_panel_smoothed.png", plot=summerT_sl_pval_sm, width = 40, height= 20, units = "cm")  



