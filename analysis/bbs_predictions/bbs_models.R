##### Set up the working directory #####
socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
if(socolar.desktop){
  dir.path <- "/Users/JacobSocolar/Dropbox/Work"
}else if(socolar.laptop){
  dir.path <- "/Users/jacob/Dropbox/Work"
}
setwd(dir.path)


library(bbsBayes)
library(magrittr)
library(brms)
# fetch_bbs_data(force = TRUE)
# stratified_data <- stratify(by = "bcr")
# CARW <- prepare_data(stratified_data, 
#                      species_to_run = "Carolina Wren",
#                      min_year = 2006,
#                      max_year = 2019,
#                      min_max_route_years = 5,
#                      model = "gamye",
#                      heavy_tailed = F)


CARW <- readRDS("/Users/Jacob/Dropbox/Work/macrodemography/CARW_BBS.RDS")

carw_data <- as.data.frame(do.call(cbind, CARW[c(7:12, 15:18)]))
carw_data$count <- as.numeric(carw_data$count)
carw_data$year <- as.numeric(carw_data$year)
carw_data$row_id <- seq(nrow(carw_data))

strata <- unique(carw_data$strat_name)


out <- list()
for(i in seq_along(strata)){
  print(i)
  dd <- carw_data[carw_data$strat_name == strata[[i]], ]
  if(length(unique(dd$route)) < 6) {
    out[[i]] <- "no estimation performed"
  } else {
    out[[i]] <-
      brm(
        count ~ s(year) + (1 | row_id) + (1 | year),
        data = dd, 
        family = "poisson",
        adapt_delta = .99999,
        max_treedepth = 12,
        backend = "cmdstanr",
        cores = 4
      )
  }
}


saveRDS(out, "macrodemography/bbs_out.RDS")

linpreds <- list()
for(i in 1:length(out)){
  print(i)
  linpreds[[i]] <- matrix(nrow = 4000, ncol = 14)
  for(j in 1:14){
    linpreds[[i]][, j] <- posterior_linpred(
      out[[i]], 
      newdata = data.frame(year = j, row_id = "dummy"),
      re_formula = ~ (1 | year)
    ) |>
      rowMeans()
  }
}

change_index <- change_median <- 
  change_u95 <- change_l95 <- change_q5 <- change_q95 <- change_q25 <- change_q75 <- 
  change_sd <- change_mean <- change_skew <- change_ex_kurt <- list()
for(i in 1:16) {
  change_index[[i]] <- t(apply(linpreds[[i]], 1, diff))
  change_median[[i]] <- matrixStats::colQuantiles(change_index[[i]], probs = .5)
  change_l95[[i]] <- matrixStats::colQuantiles(change_index[[i]], probs = .025)
  change_u95[[i]] <- matrixStats::colQuantiles(change_index[[i]], probs = .975)
  change_q5[[i]] <- matrixStats::colQuantiles(change_index[[i]], probs = .05)
  change_q95[[i]] <- matrixStats::colQuantiles(change_index[[i]], probs = .95)
  change_q25[[i]] <- matrixStats::colQuantiles(change_index[[i]], probs = .25)
  change_q75[[i]] <- matrixStats::colQuantiles(change_index[[i]], probs = .75)
  change_sd[[i]] <- apply(change_index[[i]], 2, sd)
  change_skew[[i]] <- apply(change_index[[i]], 2, moments::skewness)
  change_ex_kurt[[i]] <- apply(change_index[[i]], 2, moments::kurtosis) - 3
  change_mean[[i]] <- apply(change_index[[i]], 2, mean)
}

#source("bcr_methods.R")
summer_abun_data <- readRDS(paste0("/Users/Jacob/Dropbox/Work/macrodemography/residents/", "carw", "/summer_abun_data_fullwindow.RDS"))

all(as.integer(gsub("BCR", "", strata)) %in% summer_abun_data[[1]]$cell)

overlaps <- data.frame(stratum = rep(strata, each = 13),
                       year = 1:13,
                       diff_med = NA,
                       overlap = NA,
                       bbs = NA,
                       bbs_q5 = NA,
                       bbs_q95 = NA,
                       bbs_q25 = NA,
                       bbs_q75 = NA,
                       ebird = NA,
                       ebird_q5 = NA,
                       ebird_q95 = NA,
                       bbs_mean = NA,
                       bbs_sd = NA,
                       bbs_skew = NA, 
                       bbs_ex_kurt = NA,
                       ebird_mean = NA,
                       ebird_sd = NA,
                       ebird_skew = NA,
                       ebird_ex_kurt = NA,
                       unc = NA)


for(i in 1:16){#seq_along(strata)){
  ymax <- max(change_u95[[i]])
  ymin <- min(change_l95[[i]])
  plot(change_median[[i]][3:13] ~ c(2008:2018), ylim = c(-1.5, 1.5),
       xlab = "year",
       ylab = "change index",
       main = strata[i]
       )
  for(j in 1:13){
    lines(x = c(j+2005, j+2005), y = c(change_l95[[i]][j], change_u95[[i]][j]))
  }
  bcrx <- as.integer(gsub("BCR", "", strata[i]))
  for(j in 1:13){
    sad1 <- summer_abun_data[[j]]
    sad1 <- sad1[which(sad1$cell == bcrx), ]
    sad2 <- summer_abun_data[[j + 1]]
    sad2 <- sad2[which(sad2$cell == bcrx), ]
    
    if(sad1$n_small > 100 & sad2$n_small > 100){
      lrs <- log(as.numeric(as.vector(sad2[1, 7:106]))) - log(as.numeric(as.vector(sad1[1, 7:106])))
      points(j + 2005.2, median(lrs), pch = 17)
      lines(x = c(j + 2005.2, j + 2005.2), y = c(quantile(lrs, .025), quantile(lrs, .975)))
      
      
      overlaps$diff_med[overlaps$stratum == strata[i] & overlaps$year == j] <-
        median(lrs) - change_median[[i]][j]
      overlaps$overlap[overlaps$stratum == strata[i] & overlaps$year == j] <-
        100 - max(sum(lrs > change_index[[i]][40*c(1:100), j]),
            sum(lrs < change_index[[i]][40*c(1:100), j]))
      overlaps$bbs[overlaps$stratum == strata[i] & overlaps$year == j] <-
        change_median[[i]][j]
      
      overlaps$bbs_q5[overlaps$stratum == strata[i] & overlaps$year == j] <-
        change_q5[[i]][[j]]
      overlaps$bbs_q95[overlaps$stratum == strata[i] & overlaps$year == j] <-
        change_q95[[i]][[j]]
      
      overlaps$bbs_q25[overlaps$stratum == strata[i] & overlaps$year == j] <-
        change_q25[[i]][[j]]
      overlaps$bbs_q75[overlaps$stratum == strata[i] & overlaps$year == j] <-
        change_q75[[i]][[j]]
      
      overlaps$bbs_mean[overlaps$stratum == strata[i] & overlaps$year == j] <-
        change_mean[[i]][[j]]
      
      overlaps$bbs_mean[overlaps$stratum == strata[i] & overlaps$year == j] <-
        change_mean[[i]][[j]]
      
      overlaps$bbs_sd[overlaps$stratum == strata[i] & overlaps$year == j] <-
        change_sd[[i]][[j]]
      
      overlaps$bbs_skew[overlaps$stratum == strata[i] & overlaps$year == j] <-
        change_skew[[i]][[j]]
      
      overlaps$bbs_ex_kurt[overlaps$stratum == strata[i] & overlaps$year == j] <-
        change_ex_kurt[[i]][[j]]
      
      overlaps$ebird[overlaps$stratum == strata[i] & overlaps$year == j] <- 
        median(lrs)
      
      overlaps$ebird_mean[overlaps$stratum == strata[i] & overlaps$year == j] <- 
        mean(lrs)
      
      overlaps$ebird_sd[overlaps$stratum == strata[i] & overlaps$year == j] <- 
        sd(lrs)
      
      overlaps$ebird_skew[overlaps$stratum == strata[i] & overlaps$year == j] <-
        moments::skewness(lrs)
      
      overlaps$ebird_ex_kurt[overlaps$stratum == strata[i] & overlaps$year == j] <-
        moments::kurtosis(lrs) - 3
      
      overlaps$unc[overlaps$stratum == strata[i] & overlaps$year == j] <- 
        change_u95[[i]][j] - change_l95[[i]][j] + quantile(lrs, .95) - quantile(lrs, .05)
      
    }
  }
  
}

overlaps_p <- overlaps |>
  tidyr::pivot_longer(
    c("bbs_skew", "bbs_ex_kurt", "ebird_skew", "ebird_ex_kurt")
  )

overlaps_p$name <- factor(overlaps_p$name, levels=c("bbs_skew", "bbs_ex_kurt", "ebird_skew", "ebird_ex_kurt"))


ggplot(overlaps_p, aes(value)) + geom_density() +
  facet_wrap("name") + 
  ylab("") +
  xlab("") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank())

overlaps$alpha <- (1 - overlaps$unc) / max(overlaps$unc, na.rm = T)

strata <- unique(overlaps$stratum)
strata_include <- rep(0, length(strata))
for(i in seq_along(strata)){
  overlaps2 <- overlaps %>%
    filter(stratum == strata[i]) %>%
    filter(!is.na(bbs))
  if(nrow(overlaps2) > 2){
    if(max(overlaps2$bbs_q5) > min(overlaps2$bbs_q95)){
      strata_include[i] <- 1
    } 
  }
}
use_strata2 <- strata[as.logical(strata_include)]

# brms_output <- list()
# for(i in seq_along(use_strata)){
#   the_data <- overlaps[overlaps$stratum == use_strata[[i]], ]
#   brms_output[[i]] <- brm(bbs_mean | resp_se(bbs_sd, sigma = TRUE) ~ me(ebird_mean, ebird_sd),
#                                         data = the_data, backend = 'cmdstanr', cores = 4, adapt_delta = .999,
#                                         metric = "dense_e", max_treedepth = 13, iter = 10000, warmup = 9000)
# }


# out3 <- brm(bbs_mean | resp_se(bbs_sd, sigma = TRUE) ~ me(ebird_mean, ebird_sd) + (1 + ebird_mean | stratum),
#     data = overlaps[overlaps$stratum %in% use_strata, ], backend = 'cmdstanr', cores = 4, adapt_delta = .999,
#     metric = "dense_e", max_treedepth = 13, iter = 10000, warmup = 9000)

overlaps2 <- overlaps[!is.na(overlaps$bbs) & is.finite(overlaps$ebird_mean), ]
out4 <- brm(bbs_mean | resp_se(bbs_sd, sigma = TRUE) ~ me(ebird_mean, ebird_sd) + (1 + me(ebird_mean, ebird_sd) | stratum),
            data = overlaps2, backend = 'cmdstanr', cores = 4, adapt_delta = .999,
            max_treedepth = 10, iter = 12000, warmup = 2000)

summary(out4)


out4.5 <- brm(bbs_mean | resp_se(bbs_sd, sigma = TRUE) ~ me(ebird_mean, ebird_sd) + (1 | stratum),
            data = overlaps2, backend = 'cmdstanr', cores = 4, adapt_delta = .999,
            max_treedepth = 10, iter = 12000, warmup = 2000)

summary(out4.5)


plotting_data <- conditional_effects(
  out4, 
  conditions = data.frame(
    stratum = unique(overlaps2$stratum)
    ),
  re_formula = NULL
  ) %>%
  do.call(rbind, .) %>%
  group_by(stratum) %>%
  rename(
    bbs = estimate__,
    ebird = ebird_mean)

plotting_data_2.0 <- conditional_effects(
  out4.5
) %>%
  do.call(rbind, .) %>%
  rename(
    bbs = estimate__,
    ebird = ebird_mean) %>%
  mutate(
    class = "universal slope"
  )

plotting_data_2.1 <- conditional_effects(
  out4
) %>%
  do.call(rbind, .) %>%
  rename(
    bbs = estimate__,
    ebird = ebird_mean) %>%
  mutate(
    class = "random slopes"
  )

plotting_data_2 <- rbind(plotting_data_2.0, plotting_data_2.1) %>%
  group_by(class)

library(ggplot2)

ggplot(plotting_data, aes(x = ebird, y = bbs, color = stratum, ymin = lower__, ymax = upper__)) + 
  geom_ribbon(linetype = 0, alpha = .2) + geom_line()

ggplot(plotting_data_2, aes(x = ebird, y = bbs, color = class, ymin = lower__, ymax = upper__)) + 
  geom_ribbon(linetype = 0, alpha = .2) + geom_line()

ggplot(plotting_data_2.0, aes(x = ebird, y = bbs, ymin = lower__, ymax = upper__)) + 
  geom_ribbon(linetype = 0, alpha = .2) + geom_line()


ggplot(overlaps, aes(x = ebird, y = bbs, color = stratum, alpha = alpha)) + 
  geom_point() +
  geom_smooth(method = lm, se = F)

for(i in which(use_strata %in% use_strata2)){
  print(i)
  print(conditional_effects(brms_output[[i]], prob = .8, robust = T))
}



plot(diff_med ~ year, 
     data = overlaps, 
     col = viridis::viridis(13)[as.numeric(as.factor(overlaps$stratum))],
     pch = 16,
     alpha = overlaps$overlap/50)
