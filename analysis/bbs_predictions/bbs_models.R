library(bbsBayes)
#fetch_bbs_data()

stratified_data <- stratify(by = "bbs_usgs")

BTBW_data <- prepare_data(stratified_data, 
                               species_to_run = "Black-throated Blue Warbler",
                               min_year = 2010,
                               max_year = 2019,
                               min_max_route_years = 5,
                               model = "gamye",
                               heavy_tailed = F)


BTBW_mod <- run_model(jags_data = BTBW_data, parallel = T)
saveRDS(BTBW_mod, "/Users/jacobsocolar/Dropbox/Work/macrodemography/BBS_output/BTBW.RDS")



BBCU_data <- prepare_data(stratified_data, 
                          species_to_run = "Black-billed Cuckoo",
                          min_year = 2010,
                          max_year = 2019,
                          min_max_route_years = 5,
                          model = "gamye",
                          heavy_tailed = F)


BBCU_mod <- run_model(jags_data = BBCU_data, parallel = T)
saveRDS(BBCU_mod, "/Users/jacobsocolar/Dropbox/Work/macrodemography/BBS_output/BBCU.RDS")


CSWA_data <- prepare_data(stratified_data, 
                          species_to_run = "Chestnut-sided Warbler",
                          min_year = 2010,
                          max_year = 2019,
                          min_max_route_years = 5,
                          model = "gamye",
                          heavy_tailed = F)


CSWA_mod <- run_model(jags_data = CSWA_data, parallel = T)
saveRDS(CSWA_mod, "/Users/jacobsocolar/Dropbox/Work/macrodemography/BBS_output/CSWA.RDS")


CERW_data <- prepare_data(stratified_data, 
                          species_to_run = "Cerulean Warbler",
                          min_year = 2010,
                          max_year = 2019,
                          min_max_route_years = 5,
                          model = "gamye",
                          heavy_tailed = F)


CERW_mod <- run_model(jags_data = CERW_data, parallel = T)
saveRDS(CERW_mod, "/Users/jacobsocolar/Dropbox/Work/macrodemography/BBS_output/CERW.RDS")


CAWA_data <- prepare_data(stratified_data, 
                          species_to_run = "Canada Warbler",
                          min_year = 2010,
                          max_year = 2019,
                          min_max_route_years = 5,
                          model = "gamye",
                          heavy_tailed = F)


CAWA_mod <- run_model(jags_data = CAWA_data, parallel = T)
saveRDS(CAWA_mod, "/Users/jacobsocolar/Dropbox/Work/macrodemography/BBS_output/CAWA.RDS")


GWWA_data <- prepare_data(stratified_data, 
                          species_to_run = "Golden-winged Warbler",
                          min_year = 2010,
                          max_year = 2019,
                          min_max_route_years = 5,
                          model = "gamye",
                          heavy_tailed = F)


GWWA_mod <- run_model(jags_data = GWWA_data, parallel = T)
saveRDS(GWWA_mod, "/Users/jacobsocolar/Dropbox/Work/macrodemography/BBS_output/GWWA.RDS")


MAWA_data <- prepare_data(stratified_data, 
                          species_to_run = "Magnolia Warbler",
                          min_year = 2010,
                          max_year = 2019,
                          min_max_route_years = 5,
                          model = "gamye",
                          heavy_tailed = F)


MAWA_mod <- run_model(jags_data = MAWA_data, parallel = T)
saveRDS(MAWA_mod, "/Users/jacobsocolar/Dropbox/Work/macrodemography/BBS_output/MAWA.RDS")


MOWA_data <- prepare_data(stratified_data, 
                          species_to_run = "Mourning Warbler",
                          min_year = 2010,
                          max_year = 2019,
                          min_max_route_years = 5,
                          model = "gamye",
                          heavy_tailed = F)


MOWA_mod <- run_model(jags_data = MOWA_data, parallel = T)
saveRDS(MOWA_mod, "/Users/jacobsocolar/Dropbox/Work/macrodemography/BBS_output/MOWA.RDS")



NAWA_data <- prepare_data(stratified_data, 
                          species_to_run = "Nashville Warbler",
                          min_year = 2010,
                          max_year = 2019,
                          min_max_route_years = 5,
                          model = "gamye",
                          heavy_tailed = F)


NAWA_mod <- run_model(jags_data = NAWA_data, parallel = T)
saveRDS(NAWA_mod, "/Users/jacobsocolar/Dropbox/Work/macrodemography/BBS_output/NAWA.RDS")

NOWA_data <- prepare_data(stratified_data, 
                          species_to_run = "Northern Waterthrush",
                          min_year = 2010,
                          max_year = 2019,
                          min_max_route_years = 5,
                          model = "gamye",
                          heavy_tailed = F)


NOWA_mod <- run_model(jags_data = NOWA_data, parallel = T)
saveRDS(NOWA_mod, "/Users/jacobsocolar/Dropbox/Work/macrodemography/BBS_output/NOWA.RDS")


RBGR_data <- prepare_data(stratified_data, 
                          species_to_run = "Rose-breasted Grosbeak",
                          min_year = 2010,
                          max_year = 2019,
                          min_max_route_years = 5,
                          model = "gamye",
                          heavy_tailed = F)


RBGR_mod <- run_model(jags_data = RBGR_data, parallel = T)
saveRDS(RBGR_mod, "/Users/jacobsocolar/Dropbox/Work/macrodemography/BBS_output/RBGR.RDS")


SCTA_data <- prepare_data(stratified_data, 
                          species_to_run = "Scarlet Tanager",
                          min_year = 2010,
                          max_year = 2019,
                          min_max_route_years = 5,
                          model = "gamye",
                          heavy_tailed = F)


SCTA_mod <- run_model(jags_data = SCTA_data, parallel = T)
saveRDS(SCTA_mod, "/Users/jacobsocolar/Dropbox/Work/macrodemography/BBS_output/SCTA.RDS")