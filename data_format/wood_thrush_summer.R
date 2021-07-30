library(auk)

ebd_path <- "/Users/jacob/Desktop/ebd/"  # Location of eBird downloads from Cornell
ebd_month <- "Jun-2021"  # Download version (must be same as sampling file version)
output_path <- "/Users/jacob/Desktop/ebird_files/"  # Location to write output (and intermediates)
sampling_path <- paste0(output_path, ebd_month, "/sampling")  # Location to write intermediate outputted sampling files (prior to subdivision by species)
# dir.create(sampling_path)

states <- ebird_states$state_code[which(ebird_states$country_code == "US")]
L48 <- states[!(states %in% c("US-AK", "US-HI"))]

years <- 2006:2020
spp <- c("woothr")
spp_code <- c("WOTH")
# Subset to lower 48 and complete checklists.  Only build the sampling file once
# (i.e. for the first species). This dramatically reduces the computation
for (i in seq_along(spp)) {
  output_path2 <- paste0(output_path, ebd_month, "/", spp_code[i])
  dir.create(output_path2)
  if (i == 1) {
    auk_ebd(paste0(ebd_path, ebd_month, "/ebd_", spp[i], "_rel", ebd_month,
                   "/ebd_", spp[i], "_rel", ebd_month, ".txt"),
            paste0(ebd_path, ebd_month, "/ebd_sampling_rel", ebd_month,
                   "/ebd_sampling_rel", ebd_month, ".txt")) %>%
      auk_complete() %>%
      auk_state(state = L48) %>%
      auk_filter(file = paste0(output_path, ebd_month, "/", spp_code[i],
                               "/L48_complete.txt"),
                 file_sampling = paste0(sampling_path, "/L48_complete.txt"),
                 overwrite = T)
  } else {
    auk_ebd(paste0(ebd_path, ebd_month, "/ebd_", spp[i], "_rel", ebd_month,
                   "/ebd_", spp[i], "_rel", ebd_month, ".txt")) %>%
      auk_complete() %>%
      auk_state(state = L48) %>%
      auk_filter(file = paste0(output_path, ebd_month, "/", spp_code[i],
                               "/L48_complete.txt"),
                 filter_sampling = F, overwrite = T)
  }
}

# Restrict to summer period
for (i in seq_along(spp)) {
  print(i)
  if (i == 1) {
    auk_ebd(paste0(output_path, ebd_month, "/", spp_code[i],
                          "/L48_complete.txt"),
                   paste0(sampling_path, "/L48_complete.txt")) %>%
      auk_date(date = c("*-06-01", "*-08-01")) %>%
      auk_filter(file = paste0(output_path, ebd_month, "/", spp_code[i],
                               "/L48_complete_summer.txt"),
                 file_sampling = paste0(sampling_path,
                                        "/L48_complete_summer.txt"),
                 overwrite = T)
  } 
}

# Extract data by year.
for (i in seq_along(spp)) {
  output_path_2 <- paste0(output_path, ebd_month, "/", spp_code[i], "/by_year")
  dir.create(output_path_2)
  sampling_path_2 <- paste0(output_path, ebd_month, "/sampling/by_year")
  if (i == 1) {dir.create(sampling_path_2)}
  for (j in seq_along(years)) {
    y <- years[j]
    print(c(i, y))
    if (i == 1) {
      auk_ebd(paste0(output_path, ebd_month, "/", spp_code[i],
                            "/L48_complete_summer.txt"),
              paste0(sampling_path, "/L48_complete_summer.txt")) %>%
        auk_year(y) %>%
        auk_filter(file = paste0(output_path_2, "/summer_", y, ".txt"),
                   file_sampling = paste0(sampling_path_2, "/summer_", y, ".txt"),
                   overwrite = T)
    } 
  }
}
