library(dggridR)
library(sf)
library(ggplot2)
library(data.table)
dg6 <- dgconstruct(res=6)
dg11 <- dgconstruct(res=11)

years <- 2010:2019
cell_data <- readRDS("/Users/Jacob/Dropbox/Work/macrodemography/cell_data/cell_data_dg11_22Jun21.RDS")
spp_data <- read.csv("/Users/jacobSocolar/Dropbox/Work/macrodemography/include_by_spp.csv")
spp_code <- spp_data$species

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

conus_cells11 <- data.frame(cell = unique(dgGEO_to_SEQNUM(dg11, conus_coords$x, conus_coords$y)[[1]]))
conus_cells11$lat <- dgSEQNUM_to_GEO(dg11, conus_cells11$cell)$lat_deg
conus_cells11$lon <- dgSEQNUM_to_GEO(dg11, conus_cells11$cell)$lon_deg
conus_cells11$cell6 <- dgGEO_to_SEQNUM(dg6, conus_cells11$lon, conus_cells11$lat)[[1]]

# Get mean for a given (small) stixel
stixel_mean_11 <- function(x){
  if (x[which(names(cdijst) == "N")] == 0) {
    return(NA)
  } else if (x[which(names(cdijst) == "n_x")] + x[which(names(cdijst) == "n_positive")] == 0) {
    return(0)
  } else if (x[which(names(cdijst) == "n_positive")] == 0) {
    return(weighted.mean(x = c(0,2), w = c(x[which(names(cdijst) == "n_z")], x[which(names(cdijst) == "n_x")])))
  } else {
    n_z_rescaled <- x[which(names(cdijst) == "n_z")] * x[which(names(cdijst) == "n_positive")] / (x[which(names(cdijst) == "n_positive")] + x[which(names(cdijst) == "n_x")])
    return(weighted.mean(x = c(0, x[which(names(cdijst) == "mean_positive")]), w = c(n_z_rescaled, x[which(names(cdijst) == "n_positive")])))
  }
}

tgrid_max <- c(19, 18)

# Add the small stixel means to cell_data
for (i in seq_along(spp_code)) {
  for (j in seq_along(years)) {
    for (s in 1:2) {
      for (t in 1:tgrid_max[s]) {
        print(c(i,j,s,t))
        cdijst <- cell_data[[i]][[j]][[s]][[t]]
        if(nrow(cdijst > 0)){
          cell_data[[i]][[j]][[s]][[t]]$stixel_mean <- unlist(apply(cdijst, 1, stixel_mean_11))
        }
      }
    }
  }
}


mean_of_stixels_rep <- function(sp, yr, season, cells6) {
  vec <- vector()
  counter <- 0
  for(i in seq_along(cells6)) {
    if (season == "spring") {
      tgrid_max <- 19
    } else if (season == "fall") {
      tgrid_max <- 18
    }
    sp_yr_s <- cell_data[[sp]][[paste0("year_", yr)]][[season]]
    for(j in 1:tgrid_max){
      sp_yr_s_t <- sp_yr_s[[j]]
        if(nrow(sp_yr_s_t) != 0) {
          counter <- counter + 1
          the_rows <- sp_yr_s_t$cells11 %in% conus_cells11$cell[conus_cells11$cell6 == cells6[i]]
          vec[counter] <- weighted.mean(sp_yr_s_t$stixel_mean[the_rows],
                                        w = gtools::rdirichlet(1, rep(1, sum(sp_yr_s_t$cells11 %in% conus_cells11$cell[conus_cells11$cell6 == cells6[i]]))))
      }
    }
  }
  return(mean(vec, na.rm = T))
}



seasons <- c("spring", "fall")
cells_6_all <- unique(conus_cells11$cell6)
sp_code <- "TEWA"
sp_data <- spp_data[spp_data$species == sp_code, ]

split_point <- sp_data$min_lat + (sp_data$max_lat - sp_data$min_lat)/2
cells_n <- conus_cells$cell[conus_cells$lat > split_point &
                              conus_cells$lat < sp_data$max_lat &
                              conus_cells$lon > sp_data$min_lon_1 &
                              conus_cells$lon < sp_data$max_lon]
cells_s <- conus_cells$cell[conus_cells$lat < split_point &
                              conus_cells$lat > sp_data$min_lat &
                              conus_cells$lon > sp_data$min_lon_2 &
                              conus_cells$lon < sp_data$max_lon]
TEWA_means_n <- list()
for(y in seq_along(years)) {
  yr <- years[y]
  TEWA_means_n[[y]] <- list()
  for(s in 1:2) {
    season <- seasons[s]
    TEWA_means_n[[y]][[s]] <- vector()
    for(i in 1:100){
      print(c(y, s, i))
      TEWA_means_n[[y]][[s]][i] <- mean_of_stixels_rep(sp_code, yr, season, cells_n)
    }
  }
}

TEWA_means_s <- list()
for(y in seq_along(years)) {
  yr <- years[y]
  TEWA_means_s[[y]] <- list()
  for(s in 1:2) {
    season <- seasons[s]
    TEWA_means_s[[y]][[s]] <- vector()
    for(i in 1:100){
      print(c(y, s, i))
      TEWA_means_s[[y]][[s]][i] <- mean_of_stixels_rep(sp_code, yr, season, cells_s)
    }
  }
}

TMS <- TMN <- data.frame(rep = c(1:100))
counter <- 1
for (y in seq_along(years)) {
  for (s in 1:2) {
    counter <- counter + 1
    TMS[,counter] <- TEWA_means_s[[y]][[s]]
    TMN[,counter] <- TEWA_means_n[[y]][[s]]
  }
}

north_ratios <- log(stocks::ratios(as.numeric(TMN[1,2:21])))
for(i in 2:99) { north_ratios <- rbind(north_ratios, log(stocks::ratios(as.numeric(TMN[i,2:21]))))}

south_ratios <- log(stocks::ratios(as.numeric(TMS[1,2:21])))
for(i in 2:99) { south_ratios <- rbind(south_ratios, log(stocks::ratios(as.numeric(TMS[i,2:21]))))}

plot(1,1, xlim = c(-2,2), ylim = c(-2,2), type = 'n', xlab = "ratio: south", ylab = "ratio: north", asp =1)
for (i in 1:19) {
  if (floor(i/2)==i/2) {mc <- "blue"} else {mc <- "coral"}
  lines(x = rep(quantile(south_ratios[,i], .5), 2), y = quantile(north_ratios[,i], c(.1,.9)), col = mc)
  lines(x = quantile(south_ratios[,i], c(.1,.9)), y = rep(quantile(north_ratios[,i], .5), 2), col = mc)
}




##### Plotting #####

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

