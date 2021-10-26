

socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
if(socolar.desktop){
  dir.path <- "/Users/JacobSocolar/Dropbox/Work/macrodemography"
}else if(socolar.laptop){
  dir.path <- "/Users/jacob/Dropbox/Work/macrodemography"
}

setwd(dir.path)

library(dggridR)
dg6 <- dgconstruct(res=6)
dg11 <- dgconstruct(res=10)


species <- c("TEWA", "GCTH", "PHVI", "BBWA")
sp_codes <- c("tenwar", "gycthr", "phivir", "babwar")
yrs <- years <- c(2012:2019)
seasons <- c("spring", "fall")
bands <- c("north", "south")

source("/Users/jacob/Dropbox/Work/Code/macrodemography/stixel_bootstrap/get_tgrid.R")
source("/Users/jacob/Dropbox/Work/Code/macrodemography/stixel_bootstrap/get_pixels.R")


library(dggridR)
library(ggplot2)
library(data.table)
sg <- 10

TEWA <- readRDS("/Users/jacob/Dropbox/Work/macrodemography/cell_10_data/TEWA.RDS")

TEWA_N <- TEWA_S <- list()
for(j in seq_along(years)){
  TEWA_N[[j]] <- TEWA_S[[j]] <- list()
  for(k in seq_along(seasons)) {
    for (L in seq_along(bands)) {
      if (L == 1) {
        TEWA_N[[j]][[k]] <- TEWA[[j]][[k]][[L]]
      } else {
        TEWA_S[[j]][[k]] <- TEWA[[j]][[k]][[L]]
      }
    }
  }
}


GCTH <- readRDS("/Users/jacob/Dropbox/Work/macrodemography/cell_10_data/GCTH.RDS")

GCTH_N <- GCTH_S <- list()
for(j in seq_along(years)){
  GCTH_N[[j]] <- GCTH_S[[j]] <- list()
  for(k in seq_along(seasons)) {
    for (L in seq_along(bands)) {
      if (L == 1) {
        GCTH_N[[j]][[k]] <- GCTH[[j]][[k]][[L]]
      } else {
        GCTH_S[[j]][[k]] <- GCTH[[j]][[k]][[L]]
      }
    }
  }
}


PHVI <- readRDS("/Users/jacob/Dropbox/Work/macrodemography/cell_10_data/PHVI.RDS")

PHVI_N <- PHVI_S <- list()
for(j in seq_along(years)){
  PHVI_N[[j]] <- PHVI_S[[j]] <- list()
  for(k in seq_along(seasons)) {
    for (L in seq_along(bands)) {
      if (L == 1) {
        PHVI_N[[j]][[k]] <- PHVI[[j]][[k]][[L]]
      } else {
        PHVI_S[[j]][[k]] <- PHVI[[j]][[k]][[L]]
      }
    }
  }
}





BBWA <- readRDS("/Users/jacob/Dropbox/Work/macrodemography/cell_10_data/BBWA.RDS")

BBWA_N <- BBWA_S <- list()
for(j in seq_along(years)){
  BBWA_N[[j]] <- BBWA_S[[j]] <- list()
  for(k in seq_along(seasons)) {
    for (L in seq_along(bands)) {
      if (L == 1) {
        BBWA_N[[j]][[k]] <- BBWA[[j]][[k]][[L]]
      } else {
        BBWA_S[[j]][[k]] <- BBWA[[j]][[k]][[L]]
      }
    }
  }
}


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



stixel_n_11 <- function(x){
  return(x$n)
}



abun_rep <- function(sp_data){
  cells <- get_pixels(sg=sg)
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



abun_mean <- function(sp_data){
  cells <- get_pixels(sg=sg)
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
            n <- sp_data[[y]][[s]][[t]]$stixel_n[which_cells_data]
            n2 <- n[!is.na(data)]
            cell6_rep[[y]][[s]][[t]][i] <- weighted.mean(data2,
                                                         w = n2)
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
        sp_data[[y]][[s]][[t]]$stixel_n <- unlist(lapply(sp_data_yst, stixel_n_11))
      }
    }
  }
  abr <- abun_rep(sp_data)
  return(abr)
}




get_abun_mean <- function(sp_data){
  for (y in 1:length(sp_data)) {
    for (s in 1:2) {
      for (t in 1:length(sp_data[[1]][[s]])) {
        sp_data_yst <- sp_data[[y]][[s]][[t]]
        sp_data[[y]][[s]][[t]]$stixel_mean <- unlist(lapply(sp_data_yst, stixel_mean_11))
        sp_data[[y]][[s]][[t]]$stixel_n <- unlist(lapply(sp_data_yst, stixel_n_11))
      }
    }
  }
  abr <- abun_mean(sp_data)
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


abun_index_mean <- function(sp_data){
  abr <- get_abun_mean(sp_data)
  
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

n_rep <- 30

TEWA_N_airs <- TEWA_S_airs <- list()
for(r in 1:n_rep){
  print(r)
  TEWA_N_airs[[r]] <- abun_index_rep(TEWA_N)
  TEWA_S_airs[[r]] <- abun_index_rep(TEWA_S)
}



GCTH_N_airs <- GCTH_S_airs <- list()
for(r in 1:n_rep){
  print(r)
  GCTH_N_airs[[r]] <- abun_index_rep(GCTH_N)
  GCTH_S_airs[[r]] <- abun_index_rep(GCTH_S)
}


PHVI_N_airs <- PHVI_S_airs <- list()
for(r in 1:n_rep){
  print(r)
  PHVI_N_airs[[r]] <- abun_index_rep(PHVI_N)
  PHVI_S_airs[[r]] <- abun_index_rep(PHVI_S)
}


BBWA_N_airs <- BBWA_S_airs <- list()
for(r in 1:n_rep){
  print(r)
  BBWA_N_airs[[r]] <- abun_index_rep(BBWA_N)
  BBWA_S_airs[[r]] <- abun_index_rep(BBWA_S)
}





TEWA_N_mean <- abun_index_mean(TEWA_N)

TEWA_S_mean <- abun_index_mean(TEWA_S)

springmean_N <- unlist(TEWA_N_mean$spring)
springmean_S <- unlist(TEWA_S_mean$spring)

fallmean_N <- unlist(TEWA_N_mean$fall)
fallmean_S <- unlist(TEWA_S_mean$fall)




plot(springmean_N ~ springmean_S)
plot(fallmean_N ~ fallmean_S)

plot(springmean_S)
plot(fallmean_S)

springdata_N <- falldata_N <- matrix(nrow = 8, ncol = n_rep)
springdata_S <- falldata_S <- matrix(nrow = 8, ncol = n_rep)
for(i in 1:n_rep){
  springdata_N[,i] <- unlist(TEWA_N_airs[[i]]$spring[1:8])
  falldata_N[,i] <- unlist(TEWA_N_airs[[i]]$fall[1:8])
  springdata_S[,i] <- unlist(TEWA_S_airs[[i]]$spring[1:8])
  falldata_S[,i] <- unlist(TEWA_S_airs[[i]]$fall[1:8])
}


spring_N_M <- apply(springdata_N, 1, mean)
spring_N_L <- apply(springdata_N, 1, function(x){quantile(x, .1)})
spring_N_U <- apply(springdata_N, 1, function(x){quantile(x, .9)})

fall_N_M <- apply(falldata_N, 1, mean)
fall_N_L <- apply(falldata_N, 1, function(x){quantile(x, .1)})
fall_N_U <- apply(falldata_N, 1, function(x){quantile(x, .9)})


spring_S_M <- apply(springdata_S, 1, mean)
spring_S_L <- apply(springdata_S, 1, function(x){quantile(x, .1)})
spring_S_U <- apply(springdata_S, 1, function(x){quantile(x, .9)})

fall_S_M <- apply(falldata_S, 1, mean)
fall_S_L <- apply(falldata_S, 1, function(x){quantile(x, .1)})
fall_S_U <- apply(falldata_S, 1, function(x){quantile(x, .9)})



ymax <- max(c(spring_N_U, spring_S_U))
ymin <- min(c(spring_N_L, spring_S_L))
spring_plot <- plot(spring_N_M ~ yrs, ylim = c(ymin, ymax), 
                    main = "TEWA eBird spring 80 CI", pch = 16,
                    xlab = "year", ylab = "abundance index")
yrs2 <- years
points(spring_S_M ~ (yrs2),pch=16)
for(i in 1:8){
  lines(x = c(yrs[i], yrs[i]), y = c(spring_N_L[i], spring_N_U[i]), col = "darkmagenta",
        lwd = 3)
  lines(x = c(yrs2[i], yrs2[i]), y = c(spring_S_L[i], spring_S_U[i]), col = "goldenrod",
        lwd = 3)
}
points(spring_N_M ~ (yrs),pch=16)
points(spring_S_M ~ (yrs2),pch=16)


ymax <- max(c(fall_N_U, fall_S_U))
ymin <- min(c(fall_N_L, fall_S_L))
fall_plot <- plot(fall_N_M ~ yrs, ylim = c(ymin, ymax), main = "TEWA eBird fall 80 CI",
                  pch = 16,
                  xlab = "year", ylab = "abundance index")
yrs2 <- years
points(fall_S_M ~ (yrs2),pch=16)
for(i in 1:8){
  lines(x = c(yrs[i], yrs[i]), y = c(fall_N_L[i], fall_N_U[i]), col = "darkmagenta",
        lwd = 3)
  lines(x = c(yrs2[i], yrs2[i]), y = c(fall_S_L[i], fall_S_U[i]), col = "goldenrod",
        lwd = 3)
}
points(fall_N_M ~ (yrs),pch=16)
points(fall_S_M ~ (yrs2),pch=16)





plot(spring_N_M ~ spring_S_M)
plot(fall_N_M ~ fall_S_M)

summary(lm(spring_N_M ~ spring_S_M))
summary(lm(fall_N_M ~ fall_S_M))


length(spring_S_M)
length(fall_S_M)

south_M <- rep(NA, 16)
south_M[1+2*c(0:7)] <- spring_S_M
south_M[2*c(1:8)] <- fall_S_M
south_ratios <- stocks::ratios(south_M)


north_M <- rep(NA, 16)
north_M[1+2*c(0:7)] <- spring_N_M
north_M[2*c(1:8)] <- fall_N_M
north_ratios <- stocks::ratios(north_M)

north_productivity <- north_ratios[1+2*(c(0:7))]
south_productivity <- south_ratios[1+2*(c(0:7))]

plot(north_productivity ~ south_productivity)


north_survival <- north_ratios[2*(c(1:7))]
south_survival <- south_ratios[2*(c(1:7))]
plot(north_survival ~ south_survival)


summary(lm(north_survival ~ south_survival))
summary(lm(north_productivity ~ south_productivity))


for(i in 1:n_rep){
  springdata_N[,i] <- unlist(TEWA_N_airs[[i]]$spring[1:8])
  falldata_N[,i] <- unlist(TEWA_N_airs[[i]]$fall[1:8])
  springdata_S[,i] <- unlist(TEWA_S_airs[[i]]$spring[1:8])
  falldata_S[,i] <- unlist(TEWA_S_airs[[i]]$fall[1:8])
}


south_survival <- north_survival <- matrix(nrow=n_rep, ncol = 7)
south_productivity <- north_productivity <- matrix(nrow=n_rep, ncol = 8)
for(i in 1:n_rep){
  south_d <- north_d <- rep(NA, 16)
  south_d[1+2*c(0:7)] <- springdata_S[,i]
  south_d[2*c(1:8)] <- falldata_S[,i]
  south_ratios <- log(stocks::ratios(south_d))
  south_productivity[i, ] <- south_ratios[1+2*c(0:7)]
  south_survival[i, ] <- south_ratios[2*c(1:7)]
  
  north_d <- north_d <- rep(NA, 16)
  north_d[1+2*c(0:7)] <- springdata_N[,i]
  north_d[2*c(1:8)] <- falldata_N[,i]
  north_ratios <- log(stocks::ratios(north_d))
  north_productivity[i, ] <- north_ratios[1+2*c(0:7)]
  north_survival[i, ] <- north_ratios[2*c(1:7)]
}

north_survival_L80 <- apply(north_survival, 2, function(x){quantile(x, .1)})
north_survival_U80 <- apply(north_survival, 2, function(x){quantile(x, .9)})
north_survival_M <- apply(north_survival, 2, median)

south_survival_L80 <- apply(south_survival, 2, function(x){quantile(x, .1)})
south_survival_U80 <- apply(south_survival, 2, function(x){quantile(x, .9)})
south_survival_M <- apply(south_survival, 2, median)


north_productivity_L80 <- apply(north_productivity, 2, function(x){quantile(x, .1)})
north_productivity_U80 <- apply(north_productivity, 2, function(x){quantile(x, .9)})
north_productivity_M <- apply(north_productivity, 2, median)

south_productivity_L80 <- apply(south_productivity, 2, function(x){quantile(x, .1)})
south_productivity_U80 <- apply(south_productivity, 2, function(x){quantile(x, .9)})
south_productivity_M <- apply(south_productivity, 2, median)

surv_min <- min(c(north_survival_L80, south_survival_L80))
surv_max <- max(c(north_survival_U80, south_survival_U80))

prod_min <- min(c(north_productivity_L80, south_productivity_L80))
prod_max <- max(c(north_productivity_U80, south_productivity_U80))

plot(north_survival_M ~ south_survival_M, 
     ylim=c(surv_min, surv_max), xlim=c(surv_min, surv_max),
     pch = 16)
for(i in 1:7) {
  lines(x = rep(south_survival_M[i], 2), 
        y = c(north_survival_L80[i], north_survival_U80[i]))
  lines(x = c(south_survival_L80[i], south_survival_U80[i]), 
        y = rep(north_survival_M[i], 2))
}



plot(north_productivity_M ~ south_productivity_M, 
     ylim=c(prod_min, prod_max), xlim=c(prod_min, prod_max),
     pch = 16)
for(i in 1:7) {
  lines(x = rep(south_productivity_M[i], 2), 
        y = c(north_productivity_L80[i], north_productivity_U80[i]))
  lines(x = c(south_productivity_L80[i], south_productivity_U80[i]), 
        y = rep(north_productivity_M[i], 2))
}

