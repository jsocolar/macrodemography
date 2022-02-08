socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
if(socolar.desktop){
  dir.path <- "/Users/JacobSocolar/Dropbox/Work"
}else if(socolar.laptop){
  dir.path <- "/Users/jacob/Dropbox/Work"
}
setwd(dir.path)

roxygen2::roxygenise("Code/macrodemography/erdPackage")

install.packages("Code/macrodemography/erdPackage", 
                 repos = NULL, 
                 type = "source")

library(erdPackage)
library(data.table)

spp <- read.csv("macrodemography/species_data.csv")
spp$breaklat <- (spp$min_lat_s + spp$max_lat_n)/2
minlat <- min(spp$min_lat_s)
maxlat <- max(spp$max_lat_n)
minlon <- min(c(spp$min_lon_n, spp$min_lon_s))
# mintsp <- min(spp$min_tgrid_spring)
# maxtsp <- max(spp$min_tgrid_spring)
# mintfa <- min(spp$min_tgrid_fall)
# maxtfa <- max(spp$min_tgrid_fall)
years <- c(2012:2019)
cci_min <- c(-100, 0, .5)
effort_lim <- data.frame(dist_max = c(1, 3, 3),
                         time_min = c(5/60, 5/60, 15/60),
                         time_max = c(.5, 1, 1))

checklists <- import_checklists("macrodemography/erd/erd.db")

saveRDS(checklists, "macrodemography/erd/imported_checklists.RDS")
checklists <- readRDS("macrodemography/erd/imported_checklists.RDS")

checklists <- checklists[checklists$latitude > minlat & 
                           checklists$latitude < maxlat & 
                           checklists$longitude > minlon &
                           checklists$longitude < -60 &
                           checklists$year >= min(years) &
                           checklists$year <= max(years),]

mp_spring_south <- mp_spring_north <- 
  mp_fall_south <- mp_fall_north <-
  avg_and_rep_spring_south <- avg_and_rep_fall_south <- 
  avg_and_rep_spring_north <- avg_and_rep_fall_north <- list()
for(s in 1:nrow(spp)) {
  mp_spring_south[[s]] <- mp_spring_north[[s]] <- 
    mp_fall_south[[s]] <- mp_fall_north[[s]] <-
    avg_and_rep_spring_south[[s]] <- avg_and_rep_fall_south[[s]] <- 
    avg_and_rep_spring_north[[s]] <- avg_and_rep_fall_north[[s]] <- list()
  
  sp_data <- import_from_erd(spp$species[s], 
                             erd_path = "macrodemography/erd/erd.db",
                             checklists = checklists)
  
  for (i in 1:(length(cci_min))) {
    mp_spring_south[[s]][[i]] <- mp_spring_north[[s]][[i]] <- 
      mp_fall_south[[s]][[i]] <- mp_fall_north[[s]][[i]] <-
      avg_and_rep_spring_south[[s]][[i]] <- avg_and_rep_fall_south[[s]][[i]] <- 
      avg_and_rep_spring_north[[s]][[i]] <- avg_and_rep_fall_north[[s]][[i]] <- 
      list()
    
    
    spd2 <- sp_data[cci > cci_min[i]]
    for (j in 1:nrow(effort_lim)) {
      spd3 <- spd2[effort_distance_km <= effort_lim$dist_max[j] &
                     effort_hrs >= effort_lim$time_min[j] &
                     effort_hrs <= effort_lim$time_max[j]]
      
      
      mp_spring_south[[s]][[i]][[j]] <- 
        median_passage(sp_data = spd3, season = "spring", 
                       min_day = spp$min_tgrid_spring[s]*7 - 6,
                       max_day = spp$max_tgrid_spring[s]*7,
                       lat1 = spp$min_lat_s[s], lat2 = spp$breaklat[s],
                       bandwidth = 2)
      mp_spring_north[[s]][[i]][[j]] <- 
        median_passage(sp_data = spd3, season = "spring", 
                       min_day = spp$min_tgrid_spring[s]*7 - 6,
                       max_day = spp$max_tgrid_spring[s]*7,
                       lat1 = spp$breaklat[s], lat2 = spp$max_lat_n[s],
                       bandwidth = 2)
      
      mp_fall_south[[s]][[i]][[j]] <- 
        median_passage(sp_data = spd3, season = "fall", 
                       min_day = spp$min_tgrid_fall[s]*7 - 6,
                       max_day = spp$max_tgrid_fall[s]*7,
                       lat1 = spp$min_lat_s[s], lat2 = spp$breaklat[s],
                       bandwidth = 2)
      mp_fall_north[[s]][[i]][[j]] <- 
        median_passage(sp_data = spd3, season = "fall", 
                       min_day = spp$min_tgrid_fall[s]*7 - 6,
                       max_day = spp$max_tgrid_fall[s]*7,
                       lat1 = spp$breaklat[s], lat2 = spp$max_lat_n[s],
                       bandwidth = 2)
      
      avg_and_rep_spring_south[[s]][[i]][[j]] <- avg_and_rep_fall_south[[s]][[i]][[j]] <-
        avg_and_rep_spring_north[[s]][[i]][[j]] <- avg_and_rep_fall_north[[s]][[i]][[j]] <-
        list()
      for (y in 1:length(years)) {
        print(c(s, i, j, y))
      #   gd <- get_grid_data(spd3, years[y],
      #                       spp$min_tgrid_spring[s], spp$max_tgrid_spring[s],
      #                       spp$min_lat_s[s], spp$breaklat[s], spp$min_lon_s[s])
      #   avg_and_rep_spring_south[[s]][[i]][[j]][[y]] <- get_abun(gd, 100)
      # 
      #   gd <- get_grid_data(spd3, years[y],
      #                       spp$min_tgrid_fall[s], spp$max_tgrid_fall[s],
      #                       spp$min_lat_s[s], spp$breaklat[s], spp$min_lon_s[s])
      #   avg_and_rep_fall_south[[s]][[i]][[j]][[y]] <- get_abun(gd, 100)
      # 
      #   gd <- get_grid_data(spd3, years[y],
      #                       spp$min_tgrid_spring[s], spp$max_tgrid_spring[s],
      #                       spp$breaklat[s], spp$max_lat_n[s], spp$min_lon_n[s])
      #   avg_and_rep_spring_north[[s]][[i]][[j]][[y]] <- get_abun(gd, 100)
      # 
      #   gd <- get_grid_data(spd3, years[y],
      #                       spp$min_tgrid_fall[s], spp$max_tgrid_fall[s],
      #                       spp$breaklat[s], spp$max_lat_n[s], spp$min_lon_n[s])
      #   avg_and_rep_fall_north[[s]][[i]][[j]][[y]] <- get_abun(gd, 100)
      }
      saveRDS(mp_spring_south, paste0("macrodemography/erd_workflow/", spp$species[s], "_mp_spring_south.RDS"))
      saveRDS(mp_fall_south, paste0("macrodemography/erd_workflow/", spp$species[s], "_mp_fall_south.RDS"))
      saveRDS(mp_spring_north, paste0("macrodemography/erd_workflow/", spp$species[s], "_mp_spring_north.RDS"))
      saveRDS(mp_fall_north, paste0("macrodemography/erd_workflow/", spp$species[s], "_mp_fall_north.RDS"))
      
      # saveRDS(avg_and_rep_spring_south, paste0("macrodemography/erd_workflow/", spp$species[s], "_avg_and_rep_spring_south.RDS"))
      # saveRDS(avg_and_rep_fall_south, paste0("macrodemography/erd_workflow/", spp$species[s], "_avg_and_rep_fall_south.RDS"))
      # saveRDS(avg_and_rep_spring_north, paste0("macrodemography/erd_workflow/", spp$species[s], "_avg_and_rep_spring_north.RDS"))
      # saveRDS(avg_and_rep_fall_north, paste0("macrodemography/erd_workflow/", spp$species[s], "_avg_and_rep_fall_north.RDS"))
    }
  }
}



avg_and_rep_spring_south <- readRDS(paste0("macrodemography/erd_workflow/babwar_avg_and_rep_spring_south.RDS"))
avg_and_rep_fall_south <- readRDS(paste0("macrodemography/erd_workflow/babwar_avg_and_rep_fall_south.RDS"))
avg_and_rep_spring_north <- readRDS(paste0("macrodemography/erd_workflow/babwar_avg_and_rep_spring_north.RDS"))
avg_and_rep_fall_north <- readRDS(paste0("macrodemography/erd_workflow/babwar_avg_and_rep_fall_north.RDS"))

mp_spring_south <- readRDS(paste0("macrodemography/erd_workflow/babwar_mp_spring_south.RDS"))
mp_fall_south <- readRDS(paste0("macrodemography/erd_workflow/babwar_mp_fall_south.RDS"))
mp_spring_north <- readRDS(paste0("macrodemography/erd_workflow/babwar_mp_spring_north.RDS"))
mp_fall_north <- readRDS(paste0("macrodemography/erd_workflow/babwar_mp_fall_north.RDS"))

spp$species
cci_min
effort_lim

species <- 1
cci <- 2
effort <- 2
dev.off()

# North
ssn <- summarize_avg_and_rep(avg_and_rep_spring_north[[species]][[cci]][[effort]], avg_and_rep_fall_north[[species]][[cci]][[effort]], 
                             year_remove=4, remove = F, min_small = 5, min_list = 1, ci = .8)
plot_one_summary(ssn, mp_spring_north[[species]][[cci]][[effort]], mp_fall_north[[species]][[cci]][[effort]])

# South
sss <- summarize_avg_and_rep(avg_and_rep_spring_south[[species]][[cci]][[effort]], avg_and_rep_fall_south[[species]][[cci]][[effort]],
                             year_remove=4, remove = F, min_small = 5, min_list = 1, ci = .8)
plot_one_summary(sss, mp_spring_south[[species]][[cci]][[effort]], mp_fall_south[[species]][[cci]][[effort]])

# Both
plot_ns_summary(ssn, sss)

summer <- readRDS("macrodemography/erd_workflow/carw_avg_and_rep_summer.RDS")

cci <- 2
effort <- 2
wth <- summarize_avg_and_rep(summer[[cci]][[effort]], summer[[cci]][[effort]], year_remove = 1)


plot(wth$mean_spring ~ c(2013:2019), xlab = "year", ylab = "index",
     ylim = c(min(wth$lci_spring), max(wth$uci_spring)),
     main = "spring index")
for(i in 1:(length(wth$mean_spring))) {
  lines(c(i+2012,i+2012), c(wth$lci_spring[i], wth$uci_spring[i]))
}





##### light at night #####
cci_min <- c(0)
light_lim <- c("dim", "bright")



mp_spring_south <- mp_spring_north <- 
  mp_fall_south <- mp_fall_north <-
  avg_and_rep_spring_south <- avg_and_rep_fall_south <- 
  avg_and_rep_spring_north <- avg_and_rep_fall_north <- list()
for(s in 1:nrow(spp)) {
  mp_spring_south[[s]] <- mp_spring_north[[s]] <- 
    mp_fall_south[[s]] <- mp_fall_north[[s]] <-
    avg_and_rep_spring_south[[s]] <- avg_and_rep_fall_south[[s]] <- 
    avg_and_rep_spring_north[[s]] <- avg_and_rep_fall_north[[s]] <- list()
  
  sp_data <- import_from_erd(spp$species[s], 
                             erd_path = "macrodemography/erd/erd.db",
                             checklists = checklists)
  
  sp_data <- sp_data[effort_distance_km <= effort_lim$dist_max[2] &
                       effort_hrs >= effort_lim$time_min[2] &
                       effort_hrs <= effort_lim$time_max[2]]
  
  for (i in 1:(length(cci_min))) {
    mp_spring_south[[s]][[i]] <- mp_spring_north[[s]][[i]] <- 
      mp_fall_south[[s]][[i]] <- mp_fall_north[[s]][[i]] <-
      avg_and_rep_spring_south[[s]][[i]] <- avg_and_rep_fall_south[[s]][[i]] <- 
      avg_and_rep_spring_north[[s]][[i]] <- avg_and_rep_fall_north[[s]][[i]] <- 
      list()
    
    
    spd2 <- sp_data[cci > cci_min[i]]
    for (j in 1:length(light_lim)) {
      if(light_lim[j] == "dim") {
        spd3 <- spd2[spd2$NTL_MEAN < 3]
      } else if (light_lim[j] == "bright") {
        spd3 <- spd2[spd2$NTL_MEAN > 3]
      }
      
      mp_spring_south[[s]][[i]][[j]] <- 
        median_passage(sp_data = spd3, season = "spring", 
                       min_day = spp$min_tgrid_spring[s]*7 - 6,
                       max_day = spp$max_tgrid_spring[s]*7,
                       lat1 = spp$min_lat_s[s], lat2 = spp$breaklat[s],
                       bandwidth = 2)
      mp_spring_north[[s]][[i]][[j]] <- 
        median_passage(sp_data = spd3, season = "spring", 
                       min_day = spp$min_tgrid_spring[s]*7 - 6,
                       max_day = spp$max_tgrid_spring[s]*7,
                       lat1 = spp$breaklat[s], lat2 = spp$max_lat_n[s],
                       bandwidth = 2)
      
      mp_fall_south[[s]][[i]][[j]] <- 
        median_passage(sp_data = spd3, season = "fall", 
                       min_day = spp$min_tgrid_fall[s]*7 - 6,
                       max_day = spp$max_tgrid_fall[s]*7,
                       lat1 = spp$min_lat_s[s], lat2 = spp$breaklat[s],
                       bandwidth = 2)
      mp_fall_north[[s]][[i]][[j]] <- 
        median_passage(sp_data = spd3, season = "fall", 
                       min_day = spp$min_tgrid_fall[s]*7 - 6,
                       max_day = spp$max_tgrid_fall[s]*7,
                       lat1 = spp$breaklat[s], lat2 = spp$max_lat_n[s],
                       bandwidth = 2)
      
      avg_and_rep_spring_south[[s]][[i]][[j]] <- avg_and_rep_fall_south[[s]][[i]][[j]] <-
        avg_and_rep_spring_north[[s]][[i]][[j]] <- avg_and_rep_fall_north[[s]][[i]][[j]] <-
        list()
      for (y in 1:length(years)) {
        print(c(s, i, j, y))
          gd <- get_grid_data(spd3, years[y],
                              spp$min_tgrid_spring[s], spp$max_tgrid_spring[s],
                              spp$min_lat_s[s], spp$breaklat[s], spp$min_lon_s[s])
          avg_and_rep_spring_south[[s]][[i]][[j]][[y]] <- get_abun(gd, 100)

          gd <- get_grid_data(spd3, years[y],
                              spp$min_tgrid_fall[s], spp$max_tgrid_fall[s],
                              spp$min_lat_s[s], spp$breaklat[s], spp$min_lon_s[s])
          avg_and_rep_fall_south[[s]][[i]][[j]][[y]] <- get_abun(gd, 100)

          gd <- get_grid_data(spd3, years[y],
                              spp$min_tgrid_spring[s], spp$max_tgrid_spring[s],
                              spp$breaklat[s], spp$max_lat_n[s], spp$min_lon_n[s])
          avg_and_rep_spring_north[[s]][[i]][[j]][[y]] <- get_abun(gd, 100)

          gd <- get_grid_data(spd3, years[y],
                              spp$min_tgrid_fall[s], spp$max_tgrid_fall[s],
                              spp$breaklat[s], spp$max_lat_n[s], spp$min_lon_n[s])
          avg_and_rep_fall_north[[s]][[i]][[j]][[y]] <- get_abun(gd, 100)
      }
      saveRDS(mp_spring_south, paste0("macrodemography/erd_workflow/light_subset/", spp$species[s], "_mp_spring_south.RDS"))
      saveRDS(mp_fall_south, paste0("macrodemography/erd_workflow/light_subset/", spp$species[s], "_mp_fall_south.RDS"))
      saveRDS(mp_spring_north, paste0("macrodemography/erd_workflow/light_subset/", spp$species[s], "_mp_spring_north.RDS"))
      saveRDS(mp_fall_north, paste0("macrodemography/erd_workflow/light_subset/", spp$species[s], "_mp_fall_north.RDS"))
      
      saveRDS(avg_and_rep_spring_south, paste0("macrodemography/erd_workflow/light_subset/", spp$species[s], "_avg_and_rep_spring_south.RDS"))
      saveRDS(avg_and_rep_fall_south, paste0("macrodemography/erd_workflow/light_subset/", spp$species[s], "_avg_and_rep_fall_south.RDS"))
      saveRDS(avg_and_rep_spring_north, paste0("macrodemography/erd_workflow/light_subset/", spp$species[s], "_avg_and_rep_spring_north.RDS"))
      saveRDS(avg_and_rep_fall_north, paste0("macrodemography/erd_workflow/light_subset/", spp$species[s], "_avg_and_rep_fall_north.RDS"))
    }
  }
}



avg_and_rep_spring_south <- readRDS(paste0("macrodemography/erd_workflow/babwar_avg_and_rep_spring_south.RDS"))
avg_and_rep_fall_south <- readRDS(paste0("macrodemography/erd_workflow/babwar_avg_and_rep_fall_south.RDS"))
avg_and_rep_spring_north <- readRDS(paste0("macrodemography/erd_workflow/babwar_avg_and_rep_spring_north.RDS"))
avg_and_rep_fall_north <- readRDS(paste0("macrodemography/erd_workflow/babwar_avg_and_rep_fall_north.RDS"))

mp_spring_south <- readRDS(paste0("macrodemography/erd_workflow/babwar_mp_spring_south.RDS"))
mp_fall_south <- readRDS(paste0("macrodemography/erd_workflow/babwar_mp_fall_south.RDS"))
mp_spring_north <- readRDS(paste0("macrodemography/erd_workflow/babwar_mp_spring_north.RDS"))
mp_fall_north <- readRDS(paste0("macrodemography/erd_workflow/babwar_mp_fall_north.RDS"))

spp$species
cci_min
effort_lim

species <- 1
cci <- 1
effort <- 1
dev.off()

# North
ssn <- summarize_avg_and_rep(avg_and_rep_spring_north[[species]][[cci]][[effort]], avg_and_rep_fall_north[[species]][[cci]][[effort]], 
                             year_remove=4, remove = F, min_small = 5, min_list = 1, ci = .8)
plot_one_summary(ssn, mp_spring_north[[species]][[cci]][[effort]], mp_fall_north[[species]][[cci]][[effort]])

# South
sss <- summarize_avg_and_rep(avg_and_rep_spring_south[[species]][[cci]][[effort]], avg_and_rep_fall_south[[species]][[cci]][[effort]],
                             year_remove=4, remove = F, min_small = 5, min_list = 1, ci = .8)
plot_one_summary(sss, mp_spring_south[[species]][[cci]][[effort]], mp_fall_south[[species]][[cci]][[effort]])

# Both
plot_ns_summary(ssn, sss)

summer <- readRDS("macrodemography/erd_workflow/carw_avg_and_rep_summer.RDS")

cci <- 2
effort <- 2
wth <- summarize_avg_and_rep(summer[[cci]][[effort]], summer[[cci]][[effort]], year_remove = 1)


plot(wth$mean_spring ~ c(2013:2019), xlab = "year", ylab = "index",
     ylim = c(min(wth$lci_spring), max(wth$uci_spring)),
     main = "spring index")
for(i in 1:(length(wth$mean_spring))) {
  lines(c(i+2012,i+2012), c(wth$lci_spring[i], wth$uci_spring[i]))
}



