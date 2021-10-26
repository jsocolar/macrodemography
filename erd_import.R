library(DBI)
library(data.table)
library(tictoc)

##### For collaborative projects--figure out what machine we're on and automatically set the working directory ####

import_sp_data <- function(sp_code, erd_path, checklists = NULL, obs = NULL) {
  db <- dbConnect(RSQLite::SQLite(), erd_path)
  if (is.null(checklists)) {
    message("querying checklists")
    tic("query checklists")
    checklists_query <- dbSendQuery(db, 
                                    "SELECT sampling_event_id, latitude, longitude, year, month, day_of_year,
                                        hours_of_day, protocol_id, is_stationary, is_traveling, 
                                        effort_hrs, effort_distance_km, cci
                                     FROM checklists")
    checklists <- dbFetch(checklists_query)
    toc()
  }
  
  if (is.null(obs)) {
    message("querying observations")
    tic("query observations")
    obs_query <- dbSendQuery(db, 
                             paste0("SELECT * FROM obs
      WHERE species_code='", sp_code, "'"))
    obs <- dbFetch(obs_query)
    toc()
  }
  
  message("converting to data table and joining")
  tic("convert to dt and join")
  checklists_dt <- data.table(checklists)
  obs_dt <- data.table(obs)
  zf <- merge.data.table(checklists_dt, obs_dt, by = "sampling_event_id", all.x = T)
  toc()
  
  zf$obs_count[is.na(zf$obs_count)] <- 0
  zf$only_presence_reported[is.na(zf$only_presence_reported)] <- 0
  return(zf)
}
