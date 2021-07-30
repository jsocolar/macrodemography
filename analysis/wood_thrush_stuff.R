library(dggridR)
library(sf)
library(ggplot2)
dg6 <- dgconstruct(res=6)

# A hack to determine all grid cells that overlap the continuous US.
# Surely there's a better way...
conus <- spData::us_states
conus_raster <- fasterize::fasterize(conus,
                                     raster::raster(ncol=1000, nrow = 1000, 
                                                    xmn = -125, xmx = -66, 
                                                    ymn =24, ymx = 50))
conus_coords <- raster::as.data.frame(conus_raster, xy = T)
conus_coords <- conus_coords[!is.na(conus_coords$layer), ]
conus_cells <- data.frame(cell = unique(dgGEO_to_SEQNUM(dg6, conus_coords$x, conus_coords$y)[[1]]))
conus_cells$lat <- dgSEQNUM_to_GEO(dg6, conus_cells$cell)$lat_deg
conus_cells$lon <- dgSEQNUM_to_GEO(dg6, conus_cells$cell)$lon_deg


stixel_mean <- function(cell_data_st){
  if (cell_data_st$n == 0) {
    return(NA)
  } else if (cell_data_st$n_x + cell_data_st$n_value == 0) {
    return(0)
  } else if (cell_data_st$n_value == 0) {
    return(weighted.mean(x = c(0,2), w = c(cell_data_st$n_z, cell_data_st$n_x)))
  } else {
    n_z_rescaled <- cell_data_st$n_z * cell_data_st$n_value / (cell_data_st$n_value + cell_data_st$n_x)
    return(weighted.mean(x = c(0, cell_data_st$mean_positive), w = c(n_z_rescaled, cell_data_st$n_value)))
  }
}

mean_of_stixels <- function(sp_yr, season, cells) {
  vec <- vector()
  counter <- 0
  for(i in seq_along(cells)) {
    if (season == "summer") {
      tgrid_max <- 10
      sp_yr_s <- sp_yr[[1]]
    } 
    for(j in 1:tgrid_max){
      sp_yr_s_t <- sp_yr_s[[j]]
      for(k in seq_along(cells)){
        counter <- counter + 1
        vec[counter] <- stixel_mean(sp_yr_s_t[[which(conus_cells == cells[k])]])
      }
    }
  }
  return(mean(vec, na.rm = T))
}

get_means <- function(cell_data, sp_code, cells) {
  the_data <- cell_data[[which(spp_code == sp_code)]]
  summer_raw <- vector()
  for (i in seq_along(years)) {
    summer_raw[i] <- mean_of_stixels(the_data[[i]], "summer", cells)
  }
  out <- data.frame(means = summer_raw, 
                    yr = (1:(length(summer_raw))),
                    season = rep("summer", length(summer_raw)))
  out$ratios <- log(c(stocks::ratios(out$means), NA))
  out$year_ratios <- NA
  for (i in seq_len(nrow(out))) {
    if (i <= (nrow(out) - 1)) {
      out$year_ratios[i] <- out$ratios[i] + out$ratios[i+1]
    }
  }
  out
}

years <- 2006:2020
cell_data <- readRDS("/Users/jacob/Dropbox/Work/macrodemography/cell_data/cell_data_woth_summer.RDS")
spp_code <- "WOTH"
sp_summaries <- vector(mode = "list", length = length(spp_code))
names(sp_summaries) <- spp_code
for (i in seq_along(spp_code)) {
  sp_code <- spp_code[i]
  print(sp_code)
  
    cells_a <- conus_cells$cell
    sp_summaries[[i]] <- get_means(cell_data, sp_code, cells_a)
}


##### Plotting #####
plot((2005+sp_summaries$WOTH$yr), sp_summaries$WOTH$means/sp_summaries$WOTH$means[2])


get_color <- function(strat, summaries) {
  color_n <- c("purple", "orange")
  color_s <- c("blue", "red")
  color_a <- c("green", "coral")
  if (strat == "north") {
    return(color_n[as.integer(summaries$season == "fall") + 1])
  } else if (strat == "south") {
    return(color_s[as.integer(summaries$season == "fall") + 1])
  } else if (strat == "all") {
    return(color_a[as.integer(summaries$season == "fall") + 1])
  }
}

plot_mean_both <- function(){
  sp_sum_n <- sp_summaries[[i]]$north
  sp_sum_s <- sp_summaries[[i]]$south
  plot(means ~ yr, data = sp_sum_n, 
       col = get_color("north", sp_sum_n),
       pch = 16, ylim = c(0,max(c(sp_sum_n$means,
                                  sp_sum_s$means))),
       main = spp_code[i])
  
  points(means ~ yr, data = sp_sum_s, 
         col = get_color("south", sp_sum_s), pch = 16)
}

plot_mean_sp <- function(){
  sp_sum_n <- sp_summaries[[i]]$north[seq(1,19,2), ]
  sp_sum_s <- sp_summaries[[i]]$south[seq(1,19,2), ]
  plot(means ~ yr, data = sp_sum_n, 
       col = get_color("north", sp_sum_n),
       pch = 16, ylim = c(0,max(c(sp_sum_n$means,
                                  sp_sum_s$means))),
       main = spp_code[i])
  
  points(means ~ yr, data = sp_sum_s, 
         col = get_color("south", sp_sum_s), pch = 16)
}

plot_mean_f <- function(){
  sp_sum_n <- sp_summaries[[i]]$north[seq(2,20,2), ]
  sp_sum_s <- sp_summaries[[i]]$south[seq(2,20,2), ]
  plot(means ~ yr, data = sp_sum_n, 
       col = get_color("north", sp_sum_n),
       pch = 16, ylim = c(0,max(c(sp_sum_n$means,
                                  sp_sum_s$means))),
       main = spp_code[i])
  
  points(means ~ yr, data = sp_sum_s, 
         col = get_color("south", sp_sum_s), pch = 16)
}

plot_ratio_both <- function(){
  sp_sum_n <- sp_summaries[[i]]$north
  sp_sum_s <- sp_summaries[[i]]$south
  plot(ratios ~ yr, data = sp_sum_n, 
       col = get_color("north", sp_sum_n),
       pch = 16, ylim = c(min(c(sp_sum_n$ratios,
                                sp_sum_s$ratios)),max(c(sp_sum_n$ratios,
                                                        sp_sum_s$ratios))),
       main = spp_code[i])
  
  points(ratios ~ yr, data = sp_sum_s, 
         col = get_color("south", sp_sum_s), pch = 16)
}

plot_ratio_br <- function(){
  sp_sum_n <- sp_summaries[[i]]$north[seq(1,19,2), ]
  sp_sum_s <- sp_summaries[[i]]$south[seq(1,19,2), ]
  plot(ratios ~ yr, data = sp_sum_n, 
       col = get_color("north", sp_sum_n),
       pch = 16, ylim = c(min(c(sp_sum_n$ratios,
                                sp_sum_s$ratios)),max(c(sp_sum_n$ratios,
                                                        sp_sum_s$ratios))),
       main = spp_code[i])
  
  points(ratios ~ yr, data = sp_sum_s, 
         col = get_color("south", sp_sum_s), pch = 16)
}

plot_ratio_wi <- function(){
  sp_sum_n <- sp_summaries[[i]]$north[seq(2,20,2), ]
  sp_sum_s <- sp_summaries[[i]]$south[seq(2,20,2), ]
  plot(ratios ~ yr, data = sp_sum_n, 
       col = get_color("north", sp_sum_n),
       pch = 16, ylim = c(min(c(sp_sum_n$ratios,
                                sp_sum_s$ratios), na.rm = T),max(c(sp_sum_n$ratios,
                                                                   sp_sum_s$ratios), na.rm=T)),
       main = spp_code[i])
  
  points(ratios ~ yr, data = sp_sum_s, 
         col = get_color("south", sp_sum_s), pch = 16)
}

plot_ratio_scatter_br <- function(){
  sp_sum_n <- sp_summaries[[i]]$north
  sp_sum_s <- sp_summaries[[i]]$south
  sp_sum_n <- sp_sum_n[sp_sum_n$season == "spring", ]
  sp_sum_s <- sp_sum_s[sp_sum_s$season == "spring", ]
  
  plot(sp_sum_n$ratios ~ sp_sum_s$ratios, 
       main = spp_code[i],
       asp = 1, 
       xlab = "south productivity",
       ylab = "north productivity")
}

plot_ratio_scatter_wi <- function(){
  sp_sum_n <- sp_summaries[[i]]$north
  sp_sum_s <- sp_summaries[[i]]$south
  sp_sum_n <- sp_sum_n[sp_sum_n$season == "fall", ]
  sp_sum_s <- sp_sum_s[sp_sum_s$season == "fall", ]
  
  plot(sp_sum_n$ratios ~ sp_sum_s$ratios, 
       main = spp_code[i],
       asp = 1, 
       xlab = "south survival",
       ylab = "north survival")
}

plot_ratio_ns_spring <- function(){
  sp_sum_n <- sp_summaries[[i]]$north
  sp_sum_s <- sp_summaries[[i]]$south
  sp_sum_n <- sp_sum_n[sp_sum_n$season == "spring", "means"]
  sp_sum_s <- sp_sum_s[sp_sum_s$season == "spring", "means"]
  
  sp_sum_n2 <- sp_summaries[[i2]]$north
  sp_sum_s2 <- sp_summaries[[i2]]$south
  sp_sum_n2 <- sp_sum_n2[sp_sum_n2$season == "spring", "means"]
  sp_sum_s2 <- sp_sum_s2[sp_sum_s2$season == "spring", "means"]
  
  ratio_ns <- log(sp_sum_n/sp_sum_s)
  ratio_ns2 <- log(sp_sum_n2/sp_sum_s2)
  
  plot(ratio_ns2 ~ ratio_ns, xlab = spp_code[i], ylab = spp_code[i2], asp = 1)
}

plot_ratio_ns_fall <- function(){
  sp_sum_n <- sp_summaries[[i]]$north
  sp_sum_s <- sp_summaries[[i]]$south
  sp_sum_n <- sp_sum_n[sp_sum_n$season == "fall", "means"]
  sp_sum_s <- sp_sum_s[sp_sum_s$season == "fall", "means"]
  
  sp_sum_n2 <- sp_summaries[[i2]]$north
  sp_sum_s2 <- sp_summaries[[i2]]$south
  sp_sum_n2 <- sp_sum_n2[sp_sum_n2$season == "fall", "means"]
  sp_sum_s2 <- sp_sum_s2[sp_sum_s2$season == "fall", "means"]
  
  ratio_ns <- log(sp_sum_n/sp_sum_s)
  ratio_ns2 <- log(sp_sum_n2/sp_sum_s2)
  
  plot(ratio_ns2 ~ ratio_ns, xlab = spp_code[i], ylab = spp_code[i2], asp = 1)
}

split_spp <- c("TEWA", "BBWA", "CMWA", "PHVI", "GCTH")
ns_data_spring <- ns_data_fall <- list()
for(j in seq_along(split_spp)) {
  i <- which(spp_code == split_spp[j])
  sp_sum_n <- sp_summaries[[i]]$north
  sp_sum_s <- sp_summaries[[i]]$south
  sp_sum_n <- sp_sum_n[sp_sum_n$season == "spring", "means"]
  sp_sum_s <- sp_sum_s[sp_sum_s$season == "spring", "means"]
  ns_data_spring[[j]] <- log(sp_sum_n/sp_sum_s)
  
  sp_sum_n <- sp_summaries[[i]]$north
  sp_sum_s <- sp_summaries[[i]]$south
  sp_sum_n <- sp_sum_n[sp_sum_n$season == "fall", "means"]
  sp_sum_s <- sp_sum_s[sp_sum_s$season == "fall", "means"]
  ns_data_fall[[j]] <- log(sp_sum_n/sp_sum_s)
}
names(ns_data_spring) <- names(ns_data_fall) <- split_spp

ns_data_spring <- do.call(cbind, ns_data_spring)
ns_data_fall <- do.call(cbind, ns_data_fall)
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(ns_data_spring, lower.panel = panel.cor)
pairs(ns_data_fall, lower.panel = panel.cor)


i <- which(spp_code == "BBWA")
i2 <- which(spp_code == "CMWA")

plot_ratio_ns_spring()

plot_ratio_ns_fall()

plot_ratio_scatter_br()
plot_ratio_scatter_wi()


plot_mean_both()
lines(means ~ yr, data = sp_summaries[[i]]$north)
lines(means ~ yr, data = sp_summaries[[i]]$south, lty = 2)

plot_mean_sp()
lines(means ~ yr, data = sp_summaries[[i]]$north[seq(1, 19, 2),])
lines(means ~ yr, data = sp_summaries[[i]]$south[seq(1, 19, 2),])

plot_mean_f()
lines(means ~ yr, data = sp_summaries[[i]]$north[seq(2, 20, 2),])
lines(means ~ yr, data = sp_summaries[[i]]$south[seq(2, 20, 2),])

plot_ratio_br()
lines(ratios ~ yr, data = sp_summaries[[i]]$north[seq(1, 19, 2),])
lines(ratios ~ yr, data = sp_summaries[[i]]$south[seq(1, 19, 2),], lty = 2)

plot_ratio_wi()
lines(ratios ~ yr, data = sp_summaries[[i]]$north[seq(2, 20, 2),])
lines(ratios ~ yr, data = sp_summaries[[i]]$south[seq(2, 20, 2),])


plot(ratios ~ yr, data = north_means, col = season, pch = 16, ylim = c(-5,5))
points(ratios ~ yr, data = south_means, col = season, pch = 16)

summary(lm(south_means$ratios[south_means$season == "red"] ~ north_means$ratios[north_means$season == "orange"]))
summary(lm(south_means$ratios[south_means$season == "blue"] ~ north_means$ratios[north_means$season == "purple"]))

