# This script reads in single-species downloads from Cornell and slices them 
# up into annualized spring and fall datasets. Currently, spring and fall are
# as "one-size-fits-all" periods that should work in the passage range of any 
# long-distance migrant. Greater care will be require if working on the fringes 
# of the passage range, where contamination from breeders or winter residents is 
# a possiblility. Likewise, care will be required if, e.g., a single
# overwintering vagrant can light up multiple stixels through an overly liberal
# early- or late-season period. The current blanket definition of Spring is 01 
# To be refined, and possibly made species-specific, in the future.

library(auk)

ebd_path <- "/Users/jacob/Desktop/ebd/"  # Location of eBird downloads from Cornell
ebd_month <- "Apr-2021"  # Download version (must be same as sampling file version)
output_path <- "/Users/jacob/Desktop/ebird_files/"  # Location to write output (and intermediates)
sampling_path <- paste0(output_path, ebd_month, "/sampling")  # Location to write intermediate outputted sampling files (prior to subdivision by species)
# dir.create(sampling_path)

states <- ebird_states$state_code[which(ebird_states$country_code == "US")]
L48 <- states[!(states %in% c("US-AK", "US-HI"))]

years <- 2006:2020
spp <- c("babwar", "bkbcuc", "bkbwar", "boboli", "btbwar", "camwar", "canwar",
         "carwre", "cerwar", "chswar", "clcspa", "foxspa", "gowwar", "gycthr", 
         "magwar", "mouwar", "naswar", "norwat", "nstspa", "palwar", "phivir",
         "robgro", "scatan", "tenwar", "wlswar")
spp_code <- c("BBWA", "BBCU", "BLBW", "BOBO", "BTBW", "CMWA", "CAWA", 
              "CARW", "CERW", "CSWA", "CCSP", "FOSP", "GWWA", "GCTH", 
              "MAWA", "MOWA", "NAWA", "NOWA", "NESP", "PAWA", "PHVI", 
              "RBGR", "SCTA", "TEWA", "WIWA")
# # Subset to lower 48 and complete checklists.  Only build the sampling file once
# # (i.e. for the first species). This dramatically reduces the computation
# for (i in seq_along(spp)) {
#   output_path2 <- paste0(output_path, ebd_month, "/", spp_code[i])
#   dir.create(output_path2)
#   if (i == 1) {
#     auk_ebd(paste0(ebd_path, ebd_month, "/ebd_US_", spp[i], "_rel", ebd_month,
#                    "/ebd_US_", spp[i], "_rel", ebd_month, ".txt"),
#             paste0(ebd_path, ebd_month, "/ebd_sampling_rel", ebd_month,
#                    "/ebd_sampling_rel", ebd_month, ".txt")) %>%
#       auk_complete() %>%
#       auk_state(state = L48) %>%
#       auk_filter(file = paste0(output_path, ebd_month, "/", spp_code[i],
#                                "/L48_complete.txt"),
#                  file_sampling = paste0(sampling_path, "/L48_complete.txt"),
#                  overwrite = T)
#   } else {
#     auk_ebd(paste0(ebd_path, ebd_month, "/ebd_US_", spp[i], "_rel", ebd_month,
#                    "/ebd_US_", spp[i], "_rel", ebd_month, ".txt")) %>%
#       auk_complete() %>%
#       auk_state(state = L48) %>%
#       auk_filter(file = paste0(output_path, ebd_month, "/", spp_code[i],
#                                "/L48_complete.txt"),
#                  filter_sampling = F, overwrite = T)
#   }
# }
# 
# # Restrict to spring and fall periods
# for (i in seq_along(spp)) {
#   print(i)
#   if (i == 1) {
#     # Spring
#     auk_ebd(paste0(output_path, ebd_month, "/", spp_code[i],
#                           "/L48_complete.txt"),
#                    paste0(sampling_path, "/L48_complete.txt")) %>%
#       auk_date(date = c("*-02-01", "*-06-10")) %>%
#       auk_filter(file = paste0(output_path, ebd_month, "/", spp_code[i],
#                                "/L48_complete_spring.txt"),
#                  file_sampling = paste0(sampling_path,
#                                         "/L48_complete_spring.txt"),
#                  overwrite = T)
#     # Fall
#     auk_ebd(paste0(output_path, ebd_month, "/", spp_code[i],
#                    "/L48_complete.txt"),
#             paste0(sampling_path, "/L48_complete.txt")) %>%
#       auk_date(date = c("*-08-01", "*-11-30")) %>%
#       auk_filter(file = paste0(output_path, ebd_month, "/", spp_code[i],
#                                "/L48_complete_fall.txt"),
#                  file_sampling = paste0(sampling_path,
#                                         "/L48_complete_fall.txt"),
#                  overwrite = T)
#   } else {
#     # Spring
#     auk_ebd(paste0(output_path, ebd_month, "/", spp_code[i],
#                           "/L48_complete.txt")) %>%
#       auk_date(date = c("*-02-01", "*-06-10")) %>%
#       auk_filter(file = paste0(output_path, ebd_month, "/", spp_code[i],
#                                "/L48_complete_spring.txt"),
#                  filter_sampling = F, overwrite = T)
#     # Fall
#     auk_ebd(paste0(output_path, ebd_month, "/", spp_code[i],
#                           "/L48_complete.txt")) %>%
#   auk_date(date = c("*-08-01", "*-11-30")) %>%
#   auk_filter(file = paste0(output_path, ebd_month, "/", spp_code[i],
#                            "/L48_complete_fall.txt"),
#              filter_sampling = F, overwrite = T)
#   }
# }
# 
# # Extract data by year.
# for (i in seq_along(spp)) {
#   output_path_2 <- paste0(output_path, ebd_month, "/", spp_code[i], "/by_year")
#   dir.create(output_path_2)
#   sampling_path_2 <- paste0(output_path, ebd_month, "/sampling/by_year")
#   if (i == 1) {dir.create(sampling_path_2)}
#   for (j in seq_along(years)) {
#     y <- years[j]
#     print(c(i, y))
#     if (i == 1) {
#       # Spring
#       auk_ebd(paste0(output_path, ebd_month, "/", spp_code[i],
#                             "/L48_complete_spring.txt"),
#               paste0(sampling_path, "/L48_complete_spring.txt")) %>%
#         auk_year(y) %>%
#         auk_filter(file = paste0(output_path_2, "/spring_", y, ".txt"),
#                    file_sampling = paste0(sampling_path_2, "/spring_", y, ".txt"),
#                    overwrite = T)
#       # Fall
#       auk_ebd(paste0(output_path, ebd_month, "/", spp_code[i],
#                      "/L48_complete_fall.txt"),
#               paste0(sampling_path, "/L48_complete_fall.txt")) %>%
#         auk_year(y) %>%
#         auk_filter(file = paste0(output_path_2, "/fall_", y, ".txt"),
#                    file_sampling = paste0(sampling_path_2, "/fall_", y, ".txt"),
#                    overwrite = T)
#     } else {
#       # Spring
#       auk_ebd(paste0(output_path, ebd_month, "/", spp_code[i],
#                      "/L48_complete_spring.txt")) %>%
#         auk_year(y) %>%
#         auk_filter(file = paste0(output_path_2, "/spring_", y, ".txt"),
#                    filter_sampling = F, overwrite = T)
#       # Fall
#       auk_ebd(paste0(output_path, ebd_month, "/", spp_code[i],
#                      "/L48_complete_fall.txt")) %>%
#         auk_year(y) %>%
#         auk_filter(file = paste0(output_path_2, "/fall_", y, ".txt"),
#                    filter_sampling = F, overwrite = T)
#     }
#   }
# }
