  sr <- bbsAssistant::sauer_results
  sri <- sr$`inde_best_1966-2018_core`
  sri$Index <- as.numeric(sri$Index)
  sri$X2.5.CI <- as.numeric(sri$X2.5.CI)
  sri$X97.5.CI <- as.numeric(sri$X97.5.CI)
  sri_all <- sri[sri$Region == "SU1" & sri$Year > 2011, ]
  
  aou_spp <- unique(sri_all$AOU)
  
  j <- 250
  dev.off()
  for(i in j+(1:50)) {
    spd <- sri_all[sri_all$AOU == aou_spp[i], ]
    ymin <- min(spd$X2.5.CI)
    if(is.na(ymin)){ymin <- 0}
    ymax <- max(spd$X97.5.CI)
    if(is.na(ymax)){ymax <- 50}
    plot(Index ~ Year, data = spd, 
         main = unique(wildlifeR::AOU_species_codes$alpha.code[as.integer(wildlifeR::AOU_species_codes$spp.num) == as.integer(unique(spd$AOU))]),
         ylim = c(ymin, ymax))
    for(i in 1:nrow(spd)){
      lines(c(spd$Year[i], spd$Year[i]), c(spd$X2.5.CI[i], spd$X97.5.CI[i]))
    }
  }
  
  
  
