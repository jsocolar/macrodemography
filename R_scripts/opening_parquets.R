

library(arrow)
library(dplyr)



checklists <- open_dataset("~/Documents/ebird2022/checklists-2022.parquet") %>% collect() 
observations <- open_dataset("~/Documents/ebird2022/observations-2022.parquet") # %>% select(year) %>% collect() 
both <- open_dataset("~/Documents/ebird2022", unify_schemas = TRUE)
both2 <- open_dataset("~/Documents/ebird2022")

nyears <- unique(nyear)
treswa <- open_dataset("~/Documents/ebird2022", unify_schemas = TRUE) %>%
  filter(species_code == "treswa") %>%
  collect()

checklists <- arrow::open_dataset(checklists_parquet_path) %>%
  select(checklist_id, latitude, longitude, year, day_of_year, hours_of_day, protocol_id, is_stationary, 
         is_traveling, effort_hrs, effort_distance_km)  %>% collect()

carwre_obs <- arrow::open_dataset("~/Documents/ebird2022/observations-2022.parquet") %>% 
  filter(species_code=="carwre") # %>% collect()


zf <- checklists %>% left_join(carwre_obs, by="checklist_id") %>%  collect()

zf$obs_count[is.na(zf$obs_count)] <- 0
zf$only_presence_reported[is.na(zf$only_presence_reported)] <- 0

attr(zf, "species") <- sp_code

return(zf)



