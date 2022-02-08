#install.packages("geojsonio")
library(dggridR)
library(sf)
library(rgee)

dg <- dgconstruct(res = 6)
dgGEO_to_SEQNUM(dg, -80.8, 35.2)

# 621, 622, 595
# 729, 702
swe_data <- temp_data <- vector()
for (i in 1:7) {
  print(i)
  mindate = paste0(i+2012, "-07-01")
  maxdate = paste0(i+2012, "-08-01")
 # swe_data[i] <- daymet_extract(mindate = mindate, maxdate = maxdate)
  temp_data[i] <- daymet_extract(cell = 729, variable = "tmax", mindate = mindate, maxdate = maxdate)
}
