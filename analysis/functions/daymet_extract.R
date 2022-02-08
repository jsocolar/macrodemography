daymet_extract <- function(cell = 621, res = 6, variable = "swe", mindate = "2014-01-01", maxdate = "2014-03-01"){
  dg <- dgconstruct(res = res)
  global <- dgearthgrid(dg, frame=FALSE)
  global.cell <- data.frame(cell=getSpPPolygonsIDSlots(global), row.names=getSpPPolygonsIDSlots(global))
  global <- st_as_sf(SpatialPolygonsDataFrame(global, global.cell))
  polygon = sf_as_ee(global[cell,])
  daymet <- ee$
    ImageCollection("NASA/ORNL/DAYMET_V4")$
    filterBounds(polygon)$
    select(variable)$
    filterDate(mindate, maxdate)$
    mean()$
    reduceRegion(
      reducer = ee$Reducer$median(),
      geometry = polygon,
      scale = 1000,
      maxPixels = 1e10
    )
  return(daymet$getInfo())
}
