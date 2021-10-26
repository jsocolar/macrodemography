library(dggridR)
library(ggplot2)
library(data.table)

stixel_mean_11 <- function(x){
  if (x$n == 0) {
    return(NA)
  } else if (x$n_x + x$n_value == 0) {
    return(0)
  } else if (x$n_value == 0) {
    return(weighted.mean(x = c(0,2), w = c(x$n_z, x$n_x)))
  } else {
    n_z_rescaled <- x$n_z * x$n_value / (x$n_value + x$n_x)
    return(weighted.mean(x = c(0, x$mean_positive), w = c(n_z_rescaled, x$n_value)))
  }
}

abun_rep <- function(sp_data){
  cells <- get_pixels()
  cells11 <- list()
  for(i in seq_along(cells$conus6$cell)){
    cells11 <- cells$conus11$cell[cells$conus11$cell6 == cells$conus6$cell[i]]
  }
  
  cell6_rep <- list()
  for (y in 1:length(sp_data)) {
    cell6_rep[[y]] <- list()
    for (s in 1:2) {
      cell6_rep[[y]][[s]] <- list()
      for (t in 1:length(sp_data[[1]][[s]])) {
        cell6_rep[[y]][[s]][[t]] <- rep(NA, nrow(cells$conus6))
      }
    }
  }

  for(i in seq_along(cells$conus6$cell)){
    cells11 <- cells$conus11$cell[cells$conus11$cell6 == cells$conus6$cell[i]]
    for (y in 1:length(sp_data)) {
      for (s in 1:2) {
        for (t in 1:length(sp_data[[1]][[s]])) {
          pixel_names_data <- stringr::str_extract(names(sp_data[[y]][[s]][[t]]), "[[:digit:]]+")
          which_cells_data <- which(pixel_names_data %in% cells11)
          if (length(which_cells_data) > 0) {
            data <- sp_data[[y]][[s]][[t]]$stixel_mean[which_cells_data]
            data2 <- data[!is.na(data)]
            cell6_rep[[y]][[s]][[t]][i] <- weighted.mean(data2,
                                                         w = gtools::rdirichlet(1, rep(1, length(data2))))
          }
        }
      }
    }
    
  }
  return(cell6_rep)
}

get_abun_rep <- function(sp_data){
  for (y in 1:length(sp_data)) {
    for (s in 1:2) {
      for (t in 1:length(sp_data[[1]][[s]])) {
        sp_data_yst <- sp_data[[y]][[s]][[t]]
        sp_data[[y]][[s]][[t]]$stixel_mean <- unlist(lapply(sp_data_yst, stixel_mean_11))
      }
    }
  }
  abr <- abun_rep(sp_data)
  return(abr)
}

abun_index_rep <- function(sp_data){
  abr <- get_abun_rep(sp_data)
  
  out <- list()
  out$spring <- out$fall <- list()
  
  for(y in 1:length(sp_data)){
    yspring <- abr[[y]][[1]]
    out$spring[[y]] <- mean(unlist(lapply(yspring, function(x){mean(x, na.rm=T)})), na.rm = T)
    yfall <- abr[[y]][[2]]
    out$fall[[y]] <- mean(unlist(lapply(yfall, function(x){mean(x, na.rm=T)})), na.rm = T)
  }
  return(out)
}

sr <- bbsAssistant::sauer_results
sri <- sr$`inde_best_1993-2018_expanded`
sri$Index <- as.numeric(sri$Index)
sri$X97.5.CI <- as.numeric(sri$X97.5.CI)
sri_all <- sri[sri$Region == "SUR" & sri$Year > 2011, ]

species <- c("BBCU", "BLBW", "BOBO", "BTBW", "CAWA", "CERW", "CSWA",
             "GWWA", "MAWA", "MOWA", "NAWA", "NOWA", "RBGR", "SCTA")

species_airs <- list()
for(i in seq_along(species)){
  sp <- species[i]
  sp_data <- readRDS(paste0("/Users/jacobsocolar/Dropbox/Work/macrodemography/cell_11_data/", sp, ".RDS"))
  species_airs[[i]] <- list()
  for(r in 1:20){
    print(c(i, r))
    species_airs[[i]][[r]] <- abun_index_rep(sp_data)
  }
}

plot_data <- function(sp){
  sp_AOU_code <- wildlifeR::AOU_species_codes$spp.num[wildlifeR::AOU_species_codes$alpha.code == sp]
  sri_sp <- sri_all[sri_all$AOU == sp_AOU_code, ]
  
  airs <- species_airs[[which(species == sp)]]
  
  ymax <- max(as.numeric(sri_sp$X97.5.CI))
  ymin <- min(sri_sp$X2.5.CI)
  yrs <- c(2012:2018)
  par(mfrow = c(2,2))
  a <- plot(sri_sp$Index ~ yrs, ylim = c(ymin, ymax), main = "BBS Trend 95 CI")
  for(i in 1:7){
    lines(x = c(yrs[i], yrs[i]), y = c(sri_sp$X2.5.CI[i], sri_sp$X97.5.CI[i]))
  }
  print(a)
  
  springdata <- falldata <- matrix(nrow = 7, ncol = length(airs))
  for(i in 1:length(airs)){
    springdata[,i] <- unlist(airs[[i]]$spring[1:7])
    falldata[,i] <- unlist(airs[[i]]$fall[1:7])
  }
  springM <- apply(springdata, 1, mean)
  springL <- apply(springdata, 1, function(x){quantile(x, .1)})
  springU <- apply(springdata, 1, function(x){quantile(x, .9)})
  
  fallM <- apply(falldata, 1, mean)
  fallL <- apply(falldata, 1, function(x){quantile(x, .1)})
  fallU <- apply(falldata, 1, function(x){quantile(x, .9)})
  
  ymax <- max(springU)
  ymin <- min(springL)
  a <- plot(springM ~ yrs, ylim = c(ymin, ymax), main = "eBird spring 80 CI")
  for(i in 1:7){
    lines(x = c(yrs[i], yrs[i]), y = c(springL[i], springU[i]))
  }
  print(a)
  
  
  ymax <- max(fallU)
  ymin <- min(fallL)
  a <- plot(fallM ~ yrs, ylim = c(ymin, ymax), main = "eBird fall 80 CI")
  for(i in 1:7){
    lines(x = c(yrs[i], yrs[i]), y = c(fallL[i], fallU[i]))
  }
  print(a)
  
  xmax <- max(springU)
  xmin <- min(springL)
  ymax <- max(as.numeric(sri_sp$X97.5.CI))
  ymin <- min(sri_sp$X2.5.CI)
  
  a <- plot(sri_sp$Index ~ springM, xlim= c(xmin, xmax), ylim = c(ymin, ymax), main = "BBS vs eBird spring")
  for(i in 1:7){
    lines(x = c(springM[i], springM[i]), y = c(sri_sp$X2.5.CI[i], sri_sp$X97.5.CI[i]))
    lines(x = c(springL[i], springU[i]), y = c(sri_sp$Index[i], sri_sp$Index[i]))
  }
}

species <- c("BBCU", "BLBW", "BOBO", "BTBW", "CAWA", "CERW", "CSWA",
             "GWWA", "MAWA", "MOWA", "NAWA", "NOWA", "RBGR", "SCTA")

dev.off()

plot_data("BBCU")

for(i in 1:length(species)){
  plot_data(species[i])
}





















years <- 2012:2019

abr <- get_abun_rep(sp_data)



length(BBCU)

dg6 <- dgconstruct(res=6)
dg11 <- dgconstruct(res=11)


cell_data <- readRDS("/Users/JacobSocolar/Dropbox/Work/macrodemography/cell_data/cell_data_dg11_22Jun21.RDS")
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