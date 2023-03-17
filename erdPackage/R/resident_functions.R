#' sample a species given thresholds on effort, space and time
#' @param species_code six letter species code
#' @param erd_path path to erd
#' @param checklists filtered checklists
#' @param roi region of interest
#' @param effort_thresholds effort thresholds
#' @param extent_space spatial extent
#' @param extent_time temporal extent
#' @param time_window stratify by time grid or not
#' @param small_grid resolution of small grid
#' @param large_grid resolution of large grid
#' @param time_grid resolution of time
#' @param .cores=4 cores for parallel computation
#' @param quiet if TRUE, suppress informational messages (warnings/errors still returned)
#' @return large list object
#' @export
sample_grid_abun <- function(
    species_code, erd_path, checklists, roi,
    effort_thresholds, extent_space, extent_time,
    time_window="full", small_grid=11, large_grid=6,
    time_grid=7, .cores=4, quiet = TRUE){
  # verify input arguments
  assertthat::assert_that(is.character(species_code))
  path_erd <- paste0(erd_path, "/erd.db")
  assertthat::assert_that(file.exists(path_erd))
  assertthat::assert_that(is.data.frame(checklists))
  assertthat::assert_that(is.data.frame(effort_thresholds))
  assertthat::assert_that(all(c("dist_max","time_min","time_max","cci_min") %in% colnames(effort_thresholds)), msg="missing threshold value(s) for dist_max, time_min, time_max, cci_min")
  assertthat::assert_that(effort_thresholds$time_min < effort_thresholds$time_max)
  assertthat::assert_that(effort_thresholds$dist_max > 0)
  assertthat::assert_that(is.data.frame(extent_space))
  assertthat::assert_that(all(c("min_lon","max_lon","min_lat","max_lat") %in% colnames(extent_space)), msg="missing threshold value(s) for min_lon, max_lon, min_lat, max_lat")
  assertthat::assert_that(extent_space$min_lon < extent_space$max_lon)
  assertthat::assert_that(extent_space$min_lat < extent_space$max_lat)
  assertthat::assert_that(is.data.frame(extent_time))
  assertthat::assert_that(all(extent_time$year_min <= extent_time$year_max))
  assertthat::assert_that(all(extent_time$tgrid_min < extent_time$tgrid_max))
  assertthat::assert_that(all(c("period","tgrid_min","tgrid_max","year_min", "year_max") %in% colnames(extent_time)), msg="missing threshold value(s) for period, tgrid_min, tgrid_max, year_min, year_max")
  assertthat::assert_that(time_window %in% c("gridded", "full"))

  # load species data from ERD
  # takes ~ 2-3 mins for carwre (Carolina Wren)
  sp_data <- import_from_erd(
    species_code,
    erd_path = path_erd,
    checklists = checklists
  )

  # loop over time periods
  data_grid <- data_abun <- list()
  for (i in 1:nrow(extent_time)) {
    # loop over years
    years=extent_time$year_min[i]:extent_time$year_max[i]
    # initialize year lists
    data_grid[[i]] <- data_abun[[i]] <- list()
    for (y in seq_along(years)) {
      if(!quiet){
        print(paste("grid sampling",species_code,"data for",extent_time$period[i],years[y],"..."))
      }
      data_grid[[i]][[y]] <-
        get_grid_data(
          data = sp_data, .year = years[y],
          tgrid_min = extent_time$tgrid_min[i],
          tgrid_max = extent_time$tgrid_max[i],
          time_window = time_window, min_lat = extent_space$min_lat,
          max_lat = extent_space$max_lat, min_lon = extent_space$min_lon,
          max_lon = extent_space$max_lon, large_grid=large_grid,
          small_grid=small_grid, time_grid=time_grid, roi = roi,
          .cores = .cores
        )

      data_abun[[i]][[y]] <- get_abun(data_grid[[i]][[y]], n_rep=100, roi = roi)
    }
    names(data_grid[[i]]) <- years
    names(data_abun[[i]]) <- years
  }
  names(data_grid) <- extent_time$period
  names(data_abun) <- extent_time$period

  output <- list(grid=data_grid, abun=data_abun)
  output
}


#' Get log-ratios
#' @param data abundance index data
#' @param cells_all cells to compute
#' @param n_small_min minimum number of small cells to compute an abundance
#' @param period period to compute
#' @return large list object
#' @param quiet if TRUE, suppress informational messages (warnings/errors still returned)
#' @export
get_ratios <- function(data, cells_all, n_small_min, period=c("spring", "fall"), quiet = TRUE){
  # verify we have all periods (spring and fall) available
  assertthat::assert_that(all(period %in% names(data)))
  assertthat::assert_that(length(period)==2)
  assertthat::assert_that(is.character(period))
  # verify that loaded data contains no unknown grid cells
  cells_present <- sapply(period, function(x) assertthat::assert_that(all(data[[x]][[1]]$cell %in% cells_all)))
  assertthat::assert_that(all(cells_present), msg="loaded data contains unknown grid cells")

  if(!quiet) message(paste("calculating ratios with period1 =",period[1],"and period2 =",period[2]))

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

      cell_ratio_series_full[[i]] <- lrats
    } else {
      cell_ratio_series[[i]] <- cell_ratio_series_full[[i]] <- NA
    }
  }

  # make output a named list:
  names(cell_ratio_series) <- names(cell_timeseries)
  names(cell_ratio_series_full) <- names(cell_timeseries)

  return(list(summary=cell_ratio_series, replicates=cell_ratio_series_full))
}

#' Get higher moments of log ratios
#' @param cell_ratios cell ratios object
#' @param cells_all names of all cells underlying cell_ratios
#' @param cells all cells over which output is desired
#' @return list giving all skews and all excess kurtoses
#' @export
get_higher_moments <- function(cell_ratios, cells_all, cells) {
  cell_ratio_series <- cell_ratios$summary
  cell_ratio_series_full <- cell_ratios$replicates
  lrat_skews <- lrat_kurts <- vector()

  for(i in seq_along(cells_all)){
    if(!identical(cell_ratio_series[[i]], NA)){
      lrats_avg <- cell_ratio_series[[i]]$avg
      lrats_avg[!use_cell_years(cell_ratio_series[[i]], inf_exclude=T)] <- NA
      if (cells_all[i] %in% cells & !all(is.na(lrats_avg))) {
        assertthat::assert_that(sum(!is.na(lrats_avg)) >= 5)
        lrat_skews1 <- apply(cell_ratio_series_full[[i]], 1, moments::skewness)
        lrat_kurts1 <- apply(cell_ratio_series_full[[i]], 1, moments::kurtosis)
        lrat_skews <- c(lrat_skews, lrat_skews1)
        lrat_kurts <- c(lrat_kurts, lrat_kurts1)
      }
    }
  }
  return(list(lrat_skews = lrat_skews, lrat_ex_kurts = lrat_kurts - 3))
}


#' perform formal test of whether productivity or survival variance is larger
#' @param cell_ratios cell ratios object
#' @param cells_all all cells in cell_ratios
#' @param cells the cells over which output is desired
#' @param years years underlying the cell ratios object
#' @param min_n minimum sample size for which to estimate a variance
#' @return dataframe. prob_surv gives posterior probability that survival variance
#'  is higher than productivity variance
#'  coef is the average log-scale effect.
#' @export
variance_test <- function(cell_ratios, cells_all, cells, years, min_n = 5){
  cell_ratio_series <- cell_ratios$summary
  cell_ratio_series_full <- cell_ratios$replicates
  years_len <- length(years) - 1

  sd_holder_prod <- sd_holder_surv <- matrix(nrow = length(cells_all), ncol = 100)
  colnames(sd_holder_prod) <- paste0("prod_sd_rep_", 1:100)
  colnames(sd_holder_surv) <- paste0("surv_sd_rep_", 1:100)

  cell_lrat_sd <- cbind(data.frame(cell = cells_all, n_prod = NA, n_surv = NA),
                        as.data.frame(sd_holder_prod), as.data.frame(sd_holder_surv),
                        var_p = NA, var_d = NA, diagnostic_pass = NA)

  pb <- txtProgressBar(min = 1, max = length(cells_all), initial = 1)
  for(i in seq_along(cells_all)){
    setTxtProgressBar(pb, i)
    if(!identical(cell_ratio_series[[i]], NA)){
      lrats_avg <- cell_ratio_series[[i]]$avg
      lrats_avg[!use_cell_years(cell_ratio_series[[i]], inf_exclude=T)] <- NA
      if (cells_all[i] %in% cells & !all(is.na(lrats_avg))) {
        # This should be guaranteed by upstream data munging
        assertthat::assert_that(sum(!is.na(lrats_avg)) >= 5)
        prod_rats <- lrats_avg[1 + 2*seq_len(years_len)]
        surv_rats <- lrats_avg[2*seq_len(years_len)]
        cell_lrat_sd$n_prod[i] <- sum(!is.na(prod_rats))
        cell_lrat_sd$n_surv[i] <- sum(!is.na(surv_rats))
        lrat_means <- rowMeans(cell_ratio_series_full[[i]])
        lrat_sds <- apply(cell_ratio_series_full[[i]], 1, sd)

        prod_df <- data.frame(ratio = lrat_means[1 + 2*seq_len(years_len)],
                              sd = lrat_sds[1 + 2*seq_len(years_len)],
                              season = "prod")
        prod_df$ratio[is.infinite(prod_df$ratio) | is.nan(prod_df$ratio)] <- NA

        surv_df <- data.frame(ratio = lrat_means[2*c(1:13)],
                              sd = lrat_sds[2*c(1:13)],
                              season = "surv")
        surv_df$ratio[is.infinite(surv_df$ratio) | is.nan(surv_df$ratio)] <- NA


        model_df <- rbind(prod_df, surv_df)

        ##### model formula #####
        # modeling the ratio with known Gaussian measurement error `| resp_se...`
        # plus an additional modeled residual term `sigma = TRUE`
        # response varies by season `~ season`
        # residual sd also varies by season `sigma ~ season`
        mod_formula <- bf(ratio | resp_se(sd, sigma = TRUE) ~ season,
                          sigma ~ season)

        if(sum(!is.na(prod_df$ratio)) >= min_n & sum(!is.na(surv_df$ratio)) >= min_n){
          capture.output(suppressWarnings(suppressMessages(
            mod <- brm(mod_formula, data = model_df, family = gaussian(),
                       iter = 2000, warmup = 1000, chains = 3, refresh = 0,
                       backend = "cmdstanr")
          )))

          # check for divergences and poor r-hat (modern version of Gelman-Rubin diagnostic)
          # that checks split, folded, rank-normalized r-hat
          # if divergences occur, increase adapt_delta to adapt smaller step-sizes in leapfrog integrator (default .8)
          # if rhat is bad, run chains for longer (default 2000)
          diagnostics <- check_brmsfit_diagnostics(mod)
          if (!all(diagnostics)) {
            adapt_delta <- .8
            iter <- 2000
            if (!diagnostics[1]) {
              adapt_delta <- .99
            }
            if (!diagnostics[2]) {
              iter <- 4000
            }
            capture.output(suppressWarnings(suppressMessages(
              mod <- brm(mod_formula, data = model_df, family = gaussian(),
                         iter = iter, warmup = 1000, chains = 3, refresh = 0,
                         adapt_delta = adapt_delta,
                         backend = "cmdstanr")
            )))
            diagnostics <- check_brmsfit_diagnostics(mod)
          }
          if (all(diagnostics)) {
            cell_lrat_sd$diagnostic_pass[i] <- TRUE

            d <- as_draws_df(mod)
            # the coefficient of the sd (log link) is b_sigma_seasonsurv
            cell_lrat_sd$prob_surv[i] <- mean(d$b_sigma_seasonsurv > 0)
            cell_lrat_sd$coef[i] <- mean(d$b_sigma_seasonsurv)
          } else {
            cell_lrat_sd$diagnostic_pass[i] <- FALSE
          }
        }
      }
    }
  }
  cell_lrat_sd
}



#' perform weather regressions. currently expects to see data (and perform regressions)
#' for winter temp, winter swe, and summer temp, but this would be easily updated by
#' setting NULL defaults for these inputs plus a bit of control flow.
#' @param cell_ratios cell ratios object
#' @param cells_all all cells in cell_ratios
#' @param years years underlying the cell ratios object
#' @param min_n minimum sample size for which to estimate a slope
#' @param janfeb_temp jan/feb temps
#' @param julaug_temp july/august temps
#' @param decmar_swe dec-mar swe
#' @return 2-element list.  prob_surv gives posterior probability that survival variance
#'  is higher than productivity variance
#'  coef is the average log-scale effect.
#' @export
weather_regressions <- function(
  cell_ratios, cells_all, years, min_n,
  janfeb_temp, julaug_temp, decmar_swe
){
  plotting_data2 <- data.frame(
    cell = cells_all,
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
    swe_p_bayes = NA,
    janfeb_flag = NA, julaug_flag = NA, swe_flag = NA
  )

  cell_ratio_series <- cell_ratios$summary
  cell_ratio_series_full <- cell_ratios$replicates
  years_len <- length(years) - 1

  pb <- txtProgressBar(min = 1, max = length(cells_all), initial = 1)
  for (i in seq_along(cells_all)) {
    setTxtProgressBar(pb, i)
    if(i == 233){next} # This is a cell out over open ocean
    if (!identical(cell_ratio_series[[i]], NA)) {
      crsi <- cell_ratio_series[[i]]
      plotting_data2$surv_length[i] <- sum(is.finite(crsi$median[2*seq_len(years_len)]))
      plotting_data2$prod_length[i] <- sum(is.finite(crsi$median[2*seq_len(years_len + 1) - 1]))
      plotting_data2$full_length[i] <- sum(is.finite(crsi$median[1 + 2*seq_len(years_len)] + crsi$median[2*seq_len(years_len)]))

      janfeb_temp_slopes <- julaug_temp_slopes <- janfeb_swe_slopes <- vector()

      lrat_means <- rowMeans(cell_ratio_series_full[[i]])
      lrat_sds <- apply(cell_ratio_series_full[[i]], 1, sd)

      surv_means <- prod_means <- lrat_means
      surv_sds <- prod_sds <- lrat_sds
      surv_means[!use_cell_years(cell_ratio_series[[i]],
                                 inf_exclude = F,
                                 n_min_prod = 0,
                                 n_min_surv = min_n,
                                 n_min_full = 0)] <- NA
      surv_means <- surv_means[2*c(1:13)]
      surv_sds[!use_cell_years(cell_ratio_series[[i]],
                               inf_exclude = F,
                               n_min_prod = 0,
                               n_min_surv = min_n,
                               n_min_full = 0)] <- NA
      surv_sds <- surv_sds[2*seq_len(years_len)]

      prod_means[!use_cell_years(cell_ratio_series[[i]],
                                 inf_exclude = F,
                                 n_min_prod = min_n,
                                 n_min_surv = 0,
                                 n_min_full = 0)] <- NA
      prod_means <- prod_means[1 + 2*c(0:13)]
      prod_sds[!use_cell_years(cell_ratio_series[[i]],
                               inf_exclude = F,
                               n_min_prod = min_n,
                               n_min_surv = 0,
                               n_min_full = 0)] <- NA
      prod_sds <- prod_sds[2*seq_len(years_len + 1) - 1]

      if (sum(!is.na(surv_means)) > min_n) {
        surv_df <- data.frame(mean = surv_means, sd = surv_sds,
                              janfeb_temp = unlist(janfeb_temp[[i]][1:13]),
                              decmar_swe = sqrt(unlist(decmar_swe[[i]][1:13])))
        capture.output(suppressWarnings(suppressMessages(
          janfeb_mod <- brm(bf(mean | resp_se(sd, sigma = TRUE) ~ janfeb_temp),
                          data = surv_df, family = gaussian(),
                          prior = prior(std_normal(), class = "b"),
                          iter = 2000, warmup = 1000, chains = 3, refresh = 0,
                          backend = "cmdstanr")
        )))
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
          capture.output(suppressWarnings(suppressMessages(
            janfeb_mod <- brm(bf(mean | resp_se(sd, sigma = TRUE) ~ janfeb_temp),
                            data = surv_df, family = gaussian(),
                            prior = prior(std_normal(), class = "b"),
                            iter = iter, warmup = 1000, adapt_delta = adapt_delta,
                            chains = 3, refresh = 0,
                            backend = "cmdstanr")
          )))
          diagnostics <- check_brmsfit_diagnostics(janfeb_mod)
          if(!all(diagnostics)){
            plotting_data2$janfeb_flag[i] <- TRUE
          }
        }

        janfeb_temp_slopes <- as_draws_df(janfeb_mod)$b_janfeb_temp

        plotting_data2$janfeb_mean_bayes[i] <- mean(janfeb_temp_slopes)
        plotting_data2$janfeb_median_bayes[i] <- median(janfeb_temp_slopes)
        plotting_data2$janfeb_sd_bayes[i] <- sd(janfeb_temp_slopes)
        plotting_data2$janfeb_skew_bayes[i] <- moments::skewness(janfeb_temp_slopes)
        plotting_data2$janfeb_kurt_bayes[i] <- moments::kurtosis(janfeb_temp_slopes)
        plotting_data2$janfeb_p_bayes[i] <- mean(janfeb_temp_slopes > 0)

        # Here we hard-code that in addition to meeting the thresholds for inclusion,
        # we also need to see at least two years with nonzero snow.
        if (sum((surv_df$decmar_swe > 0) & (!is.na(surv_df$mean))) > 2) {
          capture.output(suppressWarnings(suppressMessages(
            swe_mod <- brms::brm(bf(mean | resp_se(sd, sigma = TRUE) ~ decmar_swe),
                               prior = prior(std_normal(), class = "b"),
                               data = surv_df, family = gaussian(),
                               iter = 2000, warmup = 1000, chains = 3, refresh = 0,
                               backend = "cmdstanr")
          )))
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
            capture.output(suppressWarnings(suppressMessages(
              swe_mod <- brm(bf(mean | resp_se(sd, sigma = TRUE) ~ decmar_swe),
                           data = surv_df, family = gaussian(),
                           prior = prior(std_normal(), class = "b"),
                           iter = iter, warmup = 1000, adapt_delta = adapt_delta,
                           chains = 3, refresh = 0,
                           backend = "cmdstanr")
            )))
            diagnostics <- check_brmsfit_diagnostics(swe_mod)
            if(!all(diagnostics)){
              plotting_data2$swe_flag[i] <- TRUE
            }
          }

          swe_slopes <- as_draws_df(swe_mod)$b_decmar_swe

          plotting_data2$swe_mean_bayes[i] <- mean(swe_slopes)
          plotting_data2$swe_median_bayes[i] <- median(swe_slopes)
          plotting_data2$swe_sd_bayes[i] <- sd(swe_slopes)
          plotting_data2$swe_skew_bayes[i] <- moments::skewness(swe_slopes)
          plotting_data2$swe_kurt_bayes[i] <- moments::kurtosis(swe_slopes)
          plotting_data2$swe_p_bayes[i] <- mean(swe_slopes > 0)
        }
      }

      if (sum(!is.na(prod_means)) > min_n) {
        prod_df <- data.frame(mean = prod_means, sd = prod_sds,
                              julaug_temp = unlist(julaug_temp[[i]]))
        capture.output(suppressWarnings(suppressMessages(
          julaug_mod <- brms::brm(bf(mean | resp_se(sd, sigma = TRUE) ~ julaug_temp),
                                data = prod_df, family = gaussian(),
                                prior = prior(std_normal(), class = "b"),
                                iter = 2000, warmup = 1000, chains = 3, refresh = 0,
                                backend = "cmdstanr")
        )))
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
          capture.output(suppressWarnings(suppressMessages(
            julaug_mod <- brm(bf(mean | resp_se(sd, sigma = TRUE) ~ julaug_temp),
                            data = surv_df, family = gaussian(),
                            prior = prior(std_normal(), class = "b"),
                            iter = iter, warmup = 1000, adapt_delta = adapt_delta,
                            chains = 3, refresh = 0,
                            backend = "cmdstanr")
          )))
          diagnostics <- check_brmsfit_diagnostics(julaug_mod)
          if(!all(diagnostics)){
            plotting_data2$julaug_flag[i] <- TRUE
          }
        }
        julaug_temp_slopes <- as_draws_df(julaug_mod)$b_julaug_temp

        plotting_data2$julaug_mean_bayes[i] <- mean(julaug_temp_slopes)
        plotting_data2$julaug_median_bayes[i] <- median(julaug_temp_slopes)
        plotting_data2$julaug_sd_bayes[i] <- sd(julaug_temp_slopes)
        plotting_data2$julaug_skew_bayes[i] <- moments::skewness(julaug_temp_slopes)
        plotting_data2$julaug_kurt_bayes[i] <- moments::kurtosis(julaug_temp_slopes)
        plotting_data2$julaug_p_bayes[i] <- mean(julaug_temp_slopes > 0)
      }
    }
  }
  plotting_data2
}






