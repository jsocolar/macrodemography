#################################
# load libraries
#################################

# authorization token for private repo
auth_token=readLines("token")

# Install erdPackage
#if(!"erdPackage" %in% rownames(installed.packages()))
  devtools::install_github("adokter/macrodemography", subdir="erdPackage", auth_token=auth_token)

# load libraries
library(erdPackage)
library(data.table) # note openmp not enabled on Mac by default. QUESTION: does this severely affect efficiency?
library(brms)
library(ggplot2)
library(dplyr)
library(magrittr)
library(sf)
library(dtplyr) # enable dplyr for data.table
library(assertthat)
library(ebirdst)
library(fasterize)
library(rnaturalearth)
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
#- Use install_cmdstan() to install CmdStan


#################################
# paths
#################################

# set work directory
repo_path = "~/git/macrodemography"

# set work directory
setwd(repo_path)

# where to write output
output_path = ""
# path to erd
path_erd = "~/Dropbox/macrodemography/erd/erd.db"
# path to checkist file (see flag import_checklists_from_db)
path_checklists = "~/Dropbox/macrodemography/erd/imported_checklists.RDS"
path_data = "~/Dropbox/macrodemography_refactor/data/residents"
path_checklists_filtered = "~/Dropbox/macrodemography_refactor/data/filtered_checklists.RDS"
path_brms_results = "~/Dropbox/macrodemography"

#################################
# flags
#################################

# re-create checklist file set by path_checklists by loading from erd.db
import_checklists_from_db = FALSE
# re-create filtered checklist file by applying spatial/temporal extents
filter_checklists = FALSE

# re-sample species from erd on small/large grids
resample_data = FALSE

#################################
# parameters
#################################

# years to consider:
years <- c(2006:2019)
# geographic extent
extent_space=data.frame(min_lon=-128, max_lon=-60, min_lat=23, max_lat=50)
# temporal extent should contain a row for spring and fall:
extent_time = data.frame(period=c("spring","fall"), tgrid_min=c(13,40), tgrid_max=c(16,43), year_min=min(years), year_max=max(years))
max_altitude = 2000
max_altitude_above_lat42 = 1500  #QUESTION: what is this for?
# proceed with second scenario
effort_thresholds <- data.frame(dist_max=3, time_min=5/60, time_max=1, cci_min=0)
# set hexagon grid area:
hexagon_area_large <- 70000
hexagon_area_small <- 300
# define hexagon grid
# for an area of ~ 70000 km^2, we get a resolution of 6:
grid_large <- dggridR::dgconstruct(area = hexagon_area_large)
# for an area of ~ 300 km^2, we get a resolution of 11:
grid_small <- dggridR::dgconstruct(area = hexagon_area_small)

# time grid in days
time_grid <- 7
# minimum number of small cells to compute abundance index for large cell
n_small_min = 10
# minimum number of seasonal ratios required to calculate a variance:
n_ratio_min = 5
# number of bootsrap replicates:
n_replicates = 100

# species (4 and 6 letter abbreviations)
# QUESTION: why two definitions? suggest moving to eBird species codes, see e.g. ebirdst::ebirdst_runs
species <- data.frame(four = c("bcch", "bhnu", "cach", "carw", "noca", "piwo", "tuti",
                               "canw", "cacw", "wren", "calt", "cant", "cath", "cbth", "lbth",
                               "btgn", "verd", "pyrr", "oati", "juti", "casj", "wosj"),
                      six = c("bkcchi", "bnhnut", "carchi", "carwre", "norcar", "pilwoo", "tuftit",
                              "canwre", "cacwre", "wrenti", "caltow", "cantow", "calthr", "cubthr", "lobthr",
                              "bktgna", "verdin", "pyrrhu", "oaktit", "juntit", "cowscj", "wooscj"))
# species to process
species_to_process <- c("carwre","norcar")

# print more species info:
ebirdst::ebirdst_runs %>% filter(substr(species_code,1,6) %in% species$six)

# region of interest:
# option A: USA + Canada, excluding northern regions
region_of_interest <- ne_states(country =  c("United States of America","Canada"), returnclass = "sf") %>%
  filter(region != "Northern Canada") %>% 
  filter(!name %in% c("Hawaii", "Alaska"))
# option B: USA
region_of_interest <- spData::us_states
# convert to a raster
# QUESTION/TODO: make consistent with extent_space defined above
raster_of_interest <- fasterize::fasterize(region_of_interest,
                                   raster::raster(ncol=1000, nrow = 1000, 
                                                  xmn = -125, xmx = -66, 
                                                  ymn =24, ymx = 50))

#################################
# load checklists
#################################
if(import_checklists_from_db){
  checklists <- import_checklists(path_erd)
  saveRDS(checklists, path_checklists)
} else{
  if(filter_checklists) checklists <- readRDS(path_checklists)
}

#################################
# filter checklists for extent, altitude and year
#################################
if(filter_checklists){
  # filter checklists for extent, year and elevation
  checklists %>%
    filter(latitude > extent_space$min_lat &
             latitude < extent_space$max_lat &
             longitude > extent_space$min_lon &
             longitude < extent_space$max_lon &
             year >= min(years) &
             year <= max(years) &
             ELEV_30M_MEDIAN < max_altitude &
             ((ELEV_30M_MEDIAN < max_altitude_above_lat42) | (latitude < 42)) &  # QUESTION: why this selection?
             cci > effort_thresholds$cci_min &
             effort_distance_km <= effort_thresholds$dist_max &
             effort_hrs >= effort_thresholds$time_min &
             effort_hrs <= effort_thresholds$time_max
    ) -> checklists
  
  # add hexagon indices to checklists (takes ~2-3 minutes ...)
  mutate(checklists, seqnum_large=dggridR::dgGEO_to_SEQNUM(grid_large, longitude, latitude)[[1]]) -> checklists
  mutate(checklists, seqnum_small=dggridR::dgGEO_to_SEQNUM(grid_small, longitude, latitude)[[1]]) -> checklists
  # writing results to file, since above statement takes annoyingly long
  # TODO: may want to add this to all checklists once, as part of previous load_checklists section
  saveRDS(checklists, path_checklists_filtered)
  } else{
  if(import_checklists_from_db) warning("import_checklists_from_db equals TRUE, consider refiltering checklists as well by setting filter_checklists=TRUE")
  checklists <- readRDS(path_checklists_filtered)
}

# extract unique hexagons
cells_all <- unique(checklists$seqnum)


#################################
# color scales, themes, map data
#################################

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

western_states <- c("arizona", "california", "colorado", "idaho", "montana",
                    "nevada", "new mexico", "oregon", "utah", "washington",
                    "wyoming")
states <- map_data("state") %>% filter(region %in% western_states)


#################################
# functions
#################################

# sample a species given thresholds on effort, space and time
sample_grid_abun <- function(species_code, path_erd, checklists, effort_thresholds, roi, extent_space, extent_time, time_window="full", small_grid=11, large_grid=6, time_grid=7, n_rep=100, .cores=4){
  # verify input arguments
  assert_that(is.character(species_code))
  assert_that(file.exists(path_erd))
  assert_that(is.data.frame(checklists))
  assert_that(is.data.frame(effort_thresholds))
  assert_that(all(c("dist_max","time_min","time_max","cci_min") %in% colnames(effort_thresholds)), msg="missing threshold value(s) for dist_max, time_min, time_max, cci_min")
  assert_that(effort_thresholds$time_min < effort_thresholds$time_max)
  assert_that(effort_thresholds$dist_max > 0)
  assert_that(inherits(roi,"RasterLayer"))
  assert_that(is.data.frame(extent_space))
  assert_that(all(c("min_lon","max_lon","min_lat","max_lat") %in% colnames(extent_space)), msg="missing threshold value(s) for min_lon, max_lon, min_lat, max_lat")
  assert_that(extent_space$min_lon < extent_space$max_lon)
  assert_that(extent_space$min_lat < extent_space$max_lat)
  assert_that(is.data.frame(extent_time))
  assert_that(all(extent_time$year_min <= extent_time$year_max))
  assert_that(all(extent_time$tgrid_min < extent_time$tgrid_max))
  assert_that(all(c("period","tgrid_min","tgrid_max","year_min", "year_max") %in% colnames(extent_time)), msg="missing threshold value(s) for period, tgrid_min, tgrid_max, year_min, year_max")
  assert_that(time_window %in% c("gridded", "full"))

  # load species data from ERD
  # takes ~ 2-3 mins for carwre (Carolina Wren)
  sp_data <- import_from_erd(species_code,
                             erd_path = path_erd,
                             checklists = checklists)

  # loop over time periods
  data_grid <- data_abun <- list()
  for (i in 1:nrow(extent_time)) {
    # loop over years
    years=extent_time$year_min[i]:extent_time$year_max[i]
    # initialize year lists
    data_grid[[i]] <- data_abun[[i]] <- list()
    for (y in seq_along(years)) {

      print(paste("grid sampling",species_code,"data for",extent_time$period[i],years[y],"..."))

      data_grid[[i]][[y]] <- get_grid_data(data = sp_data,
                                           .year = years[y],
                                           tgrid_min = extent_time$tgrid_min[i], tgrid_max = extent_time$tgrid_max[i],
                                           time_window = time_window, roi = roi,
                                           min_lat = extent_space$min_lat, max_lat = extent_space$max_lat, min_lon = extent_space$min_lon, max_lon = extent_space$max_lon,
                                           large_grid=large_grid, small_grid=small_grid, time_grid=time_grid)

      data_abun[[i]][[y]] <- get_abun(data_grid[[i]][[y]], roi=roi, n_rep=n_rep)
    }
    names(data_grid[[i]]) <- years
    names(data_abun[[i]]) <- years
  }
  names(data_grid) <- extent_time$period
  names(data_abun) <- extent_time$period

  output <- list(grid=data_grid, abun=data_abun)
}

get_ratios <- function(data, cells_all, period=c("spring", "fall"), n_small_min=10){
  # verify we have all periods (spring and fall) available
  assert_that(all(period %in% names(data)))
  assert_that(length(period)==2)
  assert_that(is.character(period))
  # verify that loaded data contains no unknown grid cells
  # QUESTION: in what situation could this evaluate to FALSE?
  cells_present <- sapply(period, function(x) assert_that(all(data[[x]][[1]]$cell %in% cells_all)))
  assert_that(all(cells_present), msg="loaded data contains unknown grid cells")
  
  message(paste("calculating ratios with period1 =",period[1],"and period2 =",period[2]))
  
  # Extract cell-specific abundance data
  spring_abun_summary <- get_abun_summary(data[[period[1]]], n_small_min)
  fall_abun_summary <- get_abun_summary(data[[period[2]]], n_small_min)
  cell_timeseries <- get_cell_timeseries(cells_all,
                                         spring_abun_summary,
                                         fall_abun_summary)
  # Take the log-ratios along the timeseries to get a ratio series
  # cell_ratio_series is a summary of the bootstrap uncertainty, and
  # cell_ratio_series_full gives all bootstrap replicates.
  cells <- unique(spring_abun_summary[[1]]$cell)
  cell_ratio_series <- cell_ratio_series_full <- list()
  ratio_series_length <- 2*length(spring_abun_summary)-1
  
  for (i in 1:length(cell_timeseries)) {
    if (cells_all[i] %in% cells) {
      # calculate logarithmic ratios (lrats)
      lrats <- apply(cell_timeseries[[i]][ ,2:ncol(cell_timeseries[[i]])], 2, function(x){log(stocks::ratios(x))})
      cell_ratio_series[[i]] <- list()
      
      # cell / year / reference period (named after the latest period on which the ratio is based)
      cell_ratio_series[[i]]$cell <- rep(cells_all[i],ratio_series_length)
      cell_ratio_series[[i]]$year <- as.numeric(sapply(names(spring_abun_summary), function(x) rep(x,times=2)))[-1]
      cell_ratio_series[[i]]$period <- rep(rev(period), length.out=ratio_series_length)
      
      cell_ratio_series[[i]]$median <- apply(lrats, 1, median)
      cell_ratio_series[[i]]$avg <- apply(lrats, 1, mean)
      cell_ratio_series[[i]]$sd <- apply(lrats, 1, sd)
      cell_ratio_series[[i]]$q10 <- apply(lrats, 1, function(x){quantile(x, .1, na.rm = T)})
      cell_ratio_series[[i]]$q90 <- apply(lrats, 1, function(x){quantile(x, .9, na.rm = T)})
      cell_ratio_series[[i]]$skewness <- apply(lrats, 1, moments::skewness)
      cell_ratio_series[[i]]$kurtosis <- apply(lrats, 1, moments::kurtosis)
  
      cell_ratio_series_full[[i]] <- lrats
      # set rownames to year_period
      # colnames already set in format "rep_i"
      rownames(cell_ratio_series_full[[i]]) <- paste0(cell_ratio_series[[i]]$year,"_",cell_ratio_series[[i]]$period)
    } else {
      cell_ratio_series[[i]] <- cell_ratio_series_full[[i]] <- NA
    }
  }
  
  # make output a named list:
  names(cell_ratio_series) <- names(cell_timeseries)
  names(cell_ratio_series_full) <- names(cell_timeseries)
  
  return(list(summary=cell_ratio_series, replicates=cell_ratio_series_full))
}

# function to reformat the output of get_ratios() into tidy format
make_ratios_tidy <- function(cell_ratios){
  ############################
  # tidy up replicate data 
  ############################
  tidy_rep <- function(data_rep){
    assert_that(is.list(data_rep))
    # function expects a list of length one (in order to keep the list element name available)
    assert_that(length(data_rep)==1)
    output <- as_tibble(data_rep[[1]])
    replicate_cols <- colnames(output)
    output$cell=names(data_rep)[1]
    output$year <- as.numeric(substr(rownames(data_rep[[1]]),1,4))
    output$period <- substr(rownames(data_rep[[1]]),6,10^6)
    output <- tidyr::pivot_longer(output,replicate_cols,names_to = "replicate")
    output$replicate <- as.numeric(gsub("rep_","",output$replicate))
    output
  }
  data_rep <- cell_ratios$replicates[!is.na(cell_ratios$replicates)]
  output_rep <- do.call(rbind,lapply(seq_along(data_rep), function(i) tidy_rep(data_rep[i])))
  ############################
  # tidy up summary data 
  ############################ 
  # convert lists to tibble dataframe
  output_sum <- do.call(rbind,lapply(cell_ratios$summary[!is.na(cell_ratios$summary)], as_tibble))
  
  # add information on whether there are Inf values in a given cell time series
  output_sum <- left_join(output_sum, output_sum %>% 
    group_by(cell) %>%
    summarise(has_inf=sum(is.infinite(avg))>0), by="cell")
  
  # add information on number of productivity, survival and total indices per year
  # replace Inf values with NA
  
  period1=unique(output_sum$period)[1]
  period2=unique(output_sum$period)[2]
  
  left_join(output_sum, output_sum %>%
              mutate(across(everything(), ~ replace(., is.infinite(.),NA))) %>% # change infinite values to NA
              group_by(cell) %>%
              summarise(n_prod=sum(is.finite(avg) & period==period1), # typically evaluates to period=="fall
                        n_surv=sum(is.finite(avg) & period==period2), # typically evaluates to period=="fall
                        n=sum(is.finite(avg))),
            by="cell") -> output_sum
  
  # return both as a list:
  list(summary=output_sum, replicates=output_rep)
}

# compare variances between recruitment and mortality
compare_ratio_variances <- function(cell_index, data, n_ratio_min=5, warmup=1000, iter=2000, chains=3){
  assert_that(is.data.frame(data))
  assert_that(all(c("cell","year","period","season","n_prod","n_surv","has_inf","avg") %in% names(data)))
  assert_that(cell_index %in% unique(data$cell), msg = paste("no records found for cell",cell_index,"in data"))
  
  # model formula:
  mod_formula <- bf(ratio | resp_se(sd, sigma = TRUE) ~ season,
                    sigma ~ season)
  # apply data filters:
  ratio_data <- data %>% 
    filter(cell==cell_index) %>%
    filter(n_prod>=n_ratio_min) %>%
    filter(n_surv>=n_ratio_min) %>%
    filter(is.finite(avg)) %>%
    filter(!has_inf) %>%
    mutate(ratio=avg)
  
  # initialize return values
  p_survival_variance_higher=NA
  effect_size_log=NA
  
  # only fit model if we have sufficient ratios for both recruitment and productivity:
  if(nrow(ratio_data) >= n_ratio_min & all(count(ratio_data,period)$n>=n_ratio_min)){
    adapt_delta <- 0.8
    converged <- FALSE
    tries <- 0
    
    while(!converged & tries < 2){
      print(paste("starting estimation for cell",cell_index))
      mod <- brm(mod_formula, data = ratio_data, family = gaussian(),
                 iter = iter, warmup = warmup, chains = chains, refresh = 0,
                 backend = "cmdstanr")
      diagnostics <- check_brmsfit_diagnostics(mod)
      if(all(diagnostics)) converged = TRUE
      
      if (!diagnostics[1]) {
        adapt_delta <- .99
      }
      if (!diagnostics[2]) {
        iter <- iter*2
      }
      
      tries = tries + 1
      
    }
    if(converged){
      # sample the posterior distribution:
      d <- as_draws_df(mod)
      
      p_survival_variance_higher <- mean(d$b_sigma_seasonsurv > 0)
      # p_survival_variance_higher is the probability that the variance in survival
      # is larger than the variance in productivity
      
      effect_size_log <- mean(d$b_sigma_seasonsurv)
      # effect_size_log calculated as the mean estimated posterior effect
      # equals the slope (b) for sigma when season equals survival, i.e. the log-scale difference between seasons
      # equals the mean log(survival) standard deviation - mean log(productivity) standard deviation
      # effect_size_log is still on the logarithmic scale.
      
    } else{
      print(paste("diagnostic failure for cell", cell_index))
    }
  } else{
    print(paste("insufficient ratios available for cell",cell_index))
  }
  
  return(tibble(cell=cell_index, p_survival_variance_higher=p_survival_variance_higher,effect_size_log=effect_size_log))
}



#############################################
# sample small/large grids for each species #
#############################################

if(resample_data){
  for(species_code in species_to_process){
    # sample the data
    data <- sample_grid_abun(species_code, path_erd, checklists, effort_thresholds, roi = raster_of_interest, extent_space, extent_time, time_window="full", small_grid=grid_small$res, large_grid=grid_large$res, time_grid=7, n_rep=n_replicates)
    # create output filename
    file_out <- paste0(path_data, "/data_", species_code , ".rds")
    # create output directory if not present
    if(!dir.exists(dirname(file_out))) dir.create(dirname(file_out), recursive = TRUE)
    # save the data to .rds file
    saveRDS(data, file_out)
  }
}

# check we have valid values:
sapply(1:1000,function(x) data$grid$fall[["2019"]][[1]][[x]]$stixel_mean_small)
data$grid$spring$`2018`
#################################
# Calculate spring/fall log-ratios 
#################################

##### Declare species #####
species <- "carwre"  # use consistent 6-letter convention throughout

##### load abundance data #####
file_species <- list.files(paste0(path_data), pattern=".rds$", full.names=T) %>% as_tibble %>% filter(grepl(species,value)) %>% pull(value)
print(paste("loading data from file", file_species,"..."))
data <- readRDS(file_species)

##### Get demographic indices #####
cell_ratios <- get_ratios(data$abun, cells_all, n_small_min = n_small_min)
tidy_ratios <- make_ratios_tidy(cell_ratios)
  
# rename fall/spring ratios to productivity/recruitment:
tidy_ratios$summary %>%
  mutate(season=ifelse(period=="fall", "prod","surv")) -> tidy_ratios_summary

# Plot the cell ratio series
cells <- unique(data$abun$spring[[1]]$cell)

plot_ratios <- function(cell_number, data, log=T){
  if(!log){
    data$median=exp(data$median)
    data$q10=exp(data$q10)
    data$q90=exp(data$q90)
    ylabel="ratio"
  } else{
    ylabel="log-ratio"
  }
  data %>% filter(cell==cell_number) %>%
    group_by(season) %>%
    ggplot(aes(x=year,y=median, col=season)) + 
    geom_errorbar(aes(ymin=q10, ymax=q90), width=.1,
                  position=position_dodge(0.2)) + 
    geom_point() + 
    geom_line() + 
    ggtitle(paste("cell index =",cell_number)) + 
    ylab(ylabel) 
}

tidy_ratios_summary %>%
  group_by(cell) %>%
  summarize(sufficient_data = sum(is.na(median))<20) %>%  # select cells with at least 20 valid ratios
  filter(sufficient_data) %>% pull(cell) -> cells_select

# plot selected cell 20
plot_ratios(cells_select[3], data=tidy_ratios_summary)
# plot selected cell 3-5
lapply(cells_select[3:5], plot_ratios, data=tidy_ratios_summary)


##### Analyze the timeseries #####

# average indices over cells (using inverse weighting by sd)
tidy_ratios_summary %>%
  group_by(cell,season) %>%
  filter(is.finite(avg)) %>%
  summarise(across(c("median","avg","q10","q90","sd"), \(x) weighted.mean(x, 1/sd, na.rm=TRUE))) -> data_cell


cell_ratio_series=cell_ratios$summary
cell_ratio_series_full=cell_ratios$replicates

# compare variances in productivity and survival for each cell across years:
data_compare_ratios <- lapply(sort(unique(tidy_ratios_summary$cell)),compare_ratio_variances, data=tidy_ratios_summary)
tidy_ratios_summary
View(left_join(data_cell,do.call(rbind,data_compare_ratios),by="cell"))
data_cell

# dim_series equals c(number_of_ratios=(2*number_of_years-1), number_of_replicates, number_of_cells)
dim_series <- c(dim(na.omit(cell_ratio_series_full)[[1]]),length(cell_ratio_series_full))

sd_holder_prod <- sd_holder_surv <- matrix(nrow = dim_series[3], ncol = dim_series[2])
colnames(sd_holder_prod) <- paste0("prod_sd_rep_", 1:dim_series[2])
colnames(sd_holder_surv) <- paste0("surv_sd_rep_", 1:dim_series[2])
cell_lrat_sd <- cbind(data.frame(cell = cells_all, n_prod = NA, n_surv = NA),
                      as.data.frame(sd_holder_prod), as.data.frame(sd_holder_surv))
var_p <- var_d <- rep(NA, length(cells_all))

# cell for which we have data (subset of cells_all)
cells <- unique(data$abun$spring[[1]]$cell)
# we can also get this from the cell_ratio_series_full or cell_ratio_series object:
cells <- as.numeric(names(cell_ratio_series_full)[!is.na(cell_ratio_series_full)])
cells <- as.numeric(names(cell_ratio_series)[!is.na(cell_ratio_series)])




lrat_skews <- lrat_kurts <- vector()
diagnostic_failure <- 0
# QUESTION/TODO: this for-loop should be its own function acting on cell_ratios object
for(i in seq_along(cells_all)){
  print(paste("processing cell",i))
  if(!identical(cell_ratio_series[[i]], NA)){
    lrats_avg <- cell_ratio_series[[i]]$avg
    # QUESTION: can we move thresholding by use_cell_years() to get_ratios()?
    # ANSWER: no, different selections further down.
    lrats_avg[!use_cell_years(cell_ratio_series[[i]], inf_exclude=T)] <- NA
    if (cells_all[i] %in% cells & !all(is.na(lrats_avg))) {
      # QUESTION: should this threshold be user-specifiable?
      assertthat::assert_that(sum(!is.na(lrats_avg)) >= n_ratio_min)
      # QUESTION/TODO: 13 likely refers to years, should be made year-independent
      prod_rats <- lrats_avg[1 + 2*c(1:13)]
      surv_rats <- lrats_avg[2*c(1:13)]
      
      # QUESTION/TODO: move these 5 stats can be moved to get_ratios()
      cell_lrat_sd$n_prod[i] <- sum(!is.na(prod_rats))
      cell_lrat_sd$n_surv[i] <- sum(!is.na(surv_rats))
      lrat_means <- rowMeans(cell_ratio_series_full[[i]])
      lrat_sds <- apply(cell_ratio_series_full[[i]], 1, sd)
      lrat_skews1 <- apply(cell_ratio_series_full[[i]], 1, moments::skewness)
      lrat_kurts1 <- apply(cell_ratio_series_full[[i]], 1, moments::kurtosis)
      
      #QUESTION/TODO: prod_df ans surv_df and their rbind is already available in tidy_ratios above
      prod_df <- data.frame(ratio = lrat_means[1 + 2*c(1:13)],
                            sd = lrat_sds[1 + 2*c(1:13)],
                            season = "prod")
      prod_df$ratio[is.infinite(prod_df$ratio) | is.nan(prod_df$ratio)] <- NA

      surv_df <- data.frame(ratio = lrat_means[2*c(1:13)],
                            sd = lrat_sds[2*c(1:13)],
                            season = "surv")
      surv_df$ratio[is.infinite(surv_df$ratio) | is.nan(surv_df$ratio)] <- NA


      model_df <- rbind(prod_df, surv_df)

      #QUESTION: please add some rationale in comments for this main formula
      mod_formula <- bf(ratio | resp_se(sd, sigma = TRUE) ~ season,
                        sigma ~ season)
      
      #QUESTION: should this threshold of 4 be user-defineable? Why at least 4?
      if(sum(!is.na(prod_df$ratio)) >= n_ratio_min & sum(!is.na(surv_df$ratio)) >= n_ratio_min){
        # appending lrat_skews1 to lrat_skews
        lrat_skews <- c(lrat_skews, lrat_skews1)
        lrat_kurts <- c(lrat_kurts, lrat_kurts1)
        mod <- brm(mod_formula, data = model_df, family = gaussian(),
                   iter = 2000, warmup = 1000, chains = 3, refresh = 0,
                   backend = "cmdstanr")
        #QUESTION: what are these diagnostics? Should parameters .8, 2000, .99, 4000, 'warmup', 'chains' occurring below be user-defineable?
        diagnostics <- check_brmsfit_diagnostics(mod)
        if (!all(diagnostics)) {
          adapt_delta <- .8   # this is the default
          iter <- 2000
          if (!diagnostics[1]) {
            adapt_delta <- .99
          }
          if (!diagnostics[2]) {
            iter <- 4000
          }
          mod <- brm(mod_formula, data = model_df, family = gaussian(),
                     iter = iter, warmup = 1000, chains = 3, refresh = 0,
                     adapt_delta = adapt_delta,
                     backend = "cmdstanr")
          diagnostics <- check_brmsfit_diagnostics(mod)
        }
        if (all(diagnostics)) {
          d <- as_draws_df(mod)
        #QUESTION: what is b_sigma_seasonsurv, where is it defined?
          # warmup goes per chain, so no gain there
          # but more chains could speed up the sampling part.
          var_p[i] <- mean(d$b_sigma_seasonsurv > 0) # slope (b) for sigma when season equals survival. Log-scale difference between seasons
          var_d[i] <- mean(d$b_sigma_seasonsurv)
        } else {
          print("diagnostic failure!")
          diagnostic_failure <- diagnostic_failure + 1
        }
      }
    }
  }
}

######################
# REVIEWED UNTIL HERE
######################

lrat_kurts

assertthat::assert_that(!diagnostic_failure)

skews_kurts <- data.frame(
  skew = lrat_skews,
  # excess kurtosis = kurtosis - 3
  # normal distribution has kurtosis 3, excess kurtosis 0
  `excess kurtosis` = lrat_kurts - 3). 

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
# QUESTION: what are var_p and var_d?
# var_p is the probability that the variance in survival is larger than productivity
# var_p for probability, var_d for distribution, but now only saving its mean
var_p <- readRDS(paste0(path_brms_results,"/var_p_", species, ".RDS"))
# var_d is the effect size, the mean estimated posterior effect
# var_d is the mean log(survival) standard deviation - mean log(productivity) standard deviation
# var_d is still on the logarithmic scale.
var_d <- readRDS(paste0(path_brms_results,"/var_d_", species, ".RDS"))

data.frame(cell=cells_all, var_p, var_d) %>%
  left_join(tidy_ratios, by="cell")


plotting_data <- data.frame(cell = cells_all,
                            n_prod = cell_lrat_sd$n_prod,
                            n_surv = cell_lrat_sd$n_surv,
                            p_surv_var_larger = var_p,
                            surv_sd_diff = var_d)

grid_large <- dggridR::dgconstruct(res = 6)
# QUESTION: frame and wrapcells not available in dggridR version 3.0.0, please update
grid <- dggridR::dgcellstogrid(grid_large,plotting_data$cell,frame=TRUE,wrapcells=TRUE)
grid  <- merge(grid,plotting_data,by.x="cell")

# require 5 years of productivity and survival estimates
# same as na(prod_df$ratio)) > 4 above, one can be removed.
n_min <- 5

grid_2 <- grid[!is.na(grid$n_prod), ]
grid_2 <- grid_2[grid_2$n_prod >= n_min & grid_2$n_surv >= n_min, ]
# rename p_surv_var_larger:
names(grid_2)[names(grid_2) == "p_surv_var_larger"] <- "p(survival)"

# sample size productivity
ggplot() + blank_theme + coord_fixed() +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_2,      aes(x=long, y=lat, group=group, fill=n_prod), alpha=0.7)    +
  geom_path   (data=grid_2,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  viridis::scale_fill_viridis(limits = c(5,13)) + xlim(c(-107, -65))

# sample size survival
ggplot() + blank_theme + coord_fixed() +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_2,      aes(x=long, y=lat, group=group, fill=n_surv), alpha=0.7)    +
  geom_path   (data=grid_2,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  viridis::scale_fill_viridis(limits = c(5,13)) + xlim(c(-107, -65))

# probablity that the survival variance is large than the productivity variance
p <- ggplot() + blank_theme + coord_fixed() +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_2,      aes(x=long, y=lat, group=group, fill=`p(survival)`), alpha=0.7)    +
  geom_path   (data=grid_2,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd2, na.value=NA, limits = c(0, 1)) + xlim(c(-107, -65))
p

# posterior mean effect size
fl <- max(abs(grid_2$surv_sd_diff), na.rm = T) + .1
p <- ggplot() + coord_fixed() + blank_theme +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid_2, aes(x=long, y=lat, group=group, fill = surv_sd_diff), alpha = 2*abs(grid_2$`p(survival)` - 0.5))   +
  geom_path(data=grid_data, aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  scale_fill_gradientn(colours = cols_bd, na.value=NA, limits = c(-fl, fl)) + xlim(c(-107, -65))
p

##### Extract weather data from Earth Engine #####
library(dggridR)
library(sf)
library(rgee)
source("/Users/jacob/Dropbox/Work/Code/macrodemography/analysis/functions/daymet_extract.R")
ee_Initialize()

dg <- dgconstruct(res = 6)

janfeb_temp <- julaug_temp <- junjul_precip <- novfeb_precip <-
   decmar_swe <- list()
for (j in 215:length(cells_all)){ #seq_along(cells_all)) {
  print(j)
  janfeb_temp[[j]] <- julaug_temp[[j]] <-
    novfeb_precip[[j]] <- junjul_precip[[j]] <- vector()
  decmar_swe[[j]] <- vector()
  for (i in 1:14) {
    print(i)
    mindate = paste0(i+2006, "-01-01")
    maxdate = paste0(i+2006, "-02-28")
    janfeb_temp[[j]][i] <- daymet_extract(cell = cells_all[j], variable = "tmax", mindate = mindate, maxdate = maxdate)

    mindate = paste0(i+2005, "-07-01")
    maxdate = paste0(i+2005, "-08-31")
    julaug_temp[[j]][i] <- daymet_extract(cell = cells_all[j], variable = "tmax", mindate = mindate, maxdate = maxdate)

    mindate = paste0(i+2005, "-11-01")
    maxdate = paste0(i+2006, "-02-28")
    novfeb_precip[[j]][i] <- daymet_extract(cell = cells_all[j], variable = "prcp", mindate = mindate, maxdate = maxdate)

    mindate = paste0(i+2005, "-06-01")
    maxdate = paste0(i+2005, "-07-31")
    junjul_precip[[j]][i] <- daymet_extract(cell = cells_all[j], variable = "prcp", mindate = mindate, maxdate = maxdate)

    mindate = paste0(i+2005, "-12-01")
    maxdate = paste0(i+2006, "-03-15")
    decmar_swe[[j]][i] <- daymet_extract(cell = cells_all[j], variable = "swe", mindate = mindate, maxdate = maxdate)
  }
}
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
                             surv_length = NA, # length of survival timeseries
                             prod_length = NA,
                             full_length = NA, # length for both available

                             janfeb_mean_bayes = NA, # posterian mean (effect size)
                             janfeb_median_bayes = NA, # posterior median (effect size)
                             janfeb_sd_bayes = NA,
                             janfeb_skew_bayes = NA,
                             janfeb_kurt_bayes = NA,
                             janfeb_p_bayes = NA, # posterior probability that the relationship is positive

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
      # QUESTION: please add verbal rationale for this model structure
      janfeb_mod <- brm(bf(mean | resp_se(sd, sigma = TRUE) ~ janfeb_temp),
                              data = surv_df, family = gaussian(),
                              prior = prior(std_normal(), class = "b"),
                              iter = 2000, warmup = 1000, chains = 3, refresh = 0,
                              backend = "cmdstanr")
      # more chains: more sensitive to rhat convergence criterium
      # each chain will be assigned to a cpu
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

plotting_data2$lon <- dggridR::dgSEQNUM_to_GEO(grid_large, plotting_data2$cell)$lon_deg

plotting_data2 <- plotting_data2[plotting_data2$lon > -107,]

grid <- dggridR::dgcellstogrid(grid_large,plotting_data2$cell,frame=TRUE,wrapcells=TRUE)
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
plotting_data2 <- plotting_data2[plotting_data2$cell != dggridR::dgGEO_to_SEQNUM(grid_large, -102, 41)$seqnum, ]
##### janfeb #####
# format data
car_data <- plotting_data2[!is.na(plotting_data2$`mean slope`), c("cell", "mean slope", "janfeb_median_bayes",
                                                                  "janfeb_sd_bayes", "janfeb_skew_bayes", "janfeb_kurt_bayes")]
car_data$slope_scaled <- scale(car_data$`mean slope`)
car_data$lat <- dggridR::dgSEQNUM_to_GEO(grid_large, car_data$cell)$lat_deg
car_data$lat_scaled <- scale(car_data$lat)
car_data$cell_id <- paste0("cell_", car_data$cell)
car_data$known_se <- car_data$janfeb_sd_bayes/sd(car_data$`mean slope`)

car_data$lon <- dggridR::dgSEQNUM_to_GEO(grid_large, car_data$cell)$lon_deg
car_data <- car_data[!(car_data$lat > 38.7 & car_data$lon < -99.5), ]

# get adjacency matrix
adjacency_mat <- matrix(data = 0L, nrow = nrow(car_data), ncol = nrow(car_data))
row.names(adjacency_mat) <- car_data$cell_id
for (i in 1:(nrow(car_data) - 1)) {
  coord_i <- dggridR::dgSEQNUM_to_GEO(grid_large, car_data$cell[i])
  for (j in (i+1):nrow(car_data)) {
    coord_j <- dggridR::dgSEQNUM_to_GEO(grid_large, car_data$cell[j])
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

grid <- dggridR::dgcellstogrid(grid_large,plotting_data3$cell,frame=TRUE,wrapcells=TRUE)
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
car_data$lat <- dggridR::dgSEQNUM_to_GEO(grid_large, car_data$cell)$lat_deg
car_data$lat_scaled <- scale(car_data$lat)
car_data$cell_id <- paste0("cell_", car_data$cell)
car_data$known_se <- car_data$julaug_sd_bayes/sd(car_data$julaug_mean_bayes)

car_data$lon <- dggridR::dgSEQNUM_to_GEO(grid_large, car_data$cell)$lon_deg
car_data <- car_data[!(car_data$lat > 38.7 & car_data$lon < -99.5), ]

# get adjacency matrix
adjacency_mat <- matrix(data = 0L, nrow = nrow(car_data), ncol = nrow(car_data))
row.names(adjacency_mat) <- car_data$cell_id
for (i in 1:(nrow(car_data) - 1)) {
  coord_i <- dggridR::dgSEQNUM_to_GEO(grid_large, car_data$cell[i])
  for (j in (i+1):nrow(car_data)) {
    coord_j <- dggridR::dgSEQNUM_to_GEO(grid_large, car_data$cell[j])
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

grid <- dggridR::dgcellstogrid(grid_large,plotting_data3$cell,frame=TRUE,wrapcells=TRUE)
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
car_data$lat <- dggridR::dgSEQNUM_to_GEO(grid_large, car_data$cell)$lat_deg
car_data$lat_scaled <- scale(car_data$lat)
car_data$cell_id <- paste0("cell_", car_data$cell)
car_data$known_se <- car_data$swe_sd_bayes/sd(car_data$swe_mean_bayes)

car_data$lon <- dggridR::dgSEQNUM_to_GEO(grid_large, car_data$cell)$lon_deg
car_data <- car_data[!(car_data$lat > 38.7 & car_data$lon < -99.5), ]

# get adjacency matrix
adjacency_mat <- matrix(data = 0L, nrow = nrow(car_data), ncol = nrow(car_data))
row.names(adjacency_mat) <- car_data$cell_id
for (i in 1:(nrow(car_data) - 1)) {
  coord_i <- dggridR::dgSEQNUM_to_GEO(grid_large, car_data$cell[i])
  for (j in (i+1):nrow(car_data)) {
    coord_j <- dggridR::dgSEQNUM_to_GEO(grid_large, car_data$cell[j])
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

grid <- dggridR::dgcellstogrid(grid_large,plotting_data3$cell,frame=TRUE,wrapcells=TRUE)
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

grid_large <- dggridR::dgconstruct(res = 6)

grid_x <- dggridR::dgcellstogrid(grid_large,average_indices2$cell,frame=TRUE,wrapcells=TRUE)
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
