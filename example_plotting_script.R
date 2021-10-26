socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
if(socolar.desktop){
  dir.path <- "/Users/JacobSocolar/Dropbox/Work"
}else if(socolar.laptop){
  dir.path <- "/Users/jacob/Dropbox/Work"
}
setwd(dir.path)

# macrodemography is a private git repo.  The below line can be replaced
# remotes::install_github() but requires passing the PAT
install.packages("Code/macrodemography/erdPackage", 
                 repos = NULL, 
                 type = "source")

library(erdPackage)
library(data.table)

setwd("macrodemography")

spp <- read.csv("species_data.csv")

avg_and_rep_spring_south <- readRDS(paste0("erd_workflow/babwar_avg_and_rep_spring_south.RDS"))
avg_and_rep_fall_south <- readRDS(paste0("erd_workflow/babwar_avg_and_rep_fall_south.RDS"))
avg_and_rep_spring_north <- readRDS(paste0("erd_workflow/babwar_avg_and_rep_spring_north.RDS"))
avg_and_rep_fall_north <- readRDS(paste0("erd_workflow/babwar_avg_and_rep_fall_north.RDS"))

mp_spring_south <- readRDS(paste0("erd_workflow/babwar_mp_spring_south.RDS"))
mp_fall_south <- readRDS(paste0("erd_workflow/babwar_mp_fall_south.RDS"))
mp_spring_north <- readRDS(paste0("erd_workflow/babwar_mp_spring_north.RDS"))
mp_fall_north <- readRDS(paste0("erd_workflow/babwar_mp_fall_north.RDS"))

years <- c(2012:2019)
cci_min <- c(-100, 0, .5)
effort_lim <- data.frame(dist_max = c(1, 3, 3),
                         time_min = c(5/60, 5/60, 15/60),
                         time_max = c(.5, 1, 1))
spp$species
cci_min
effort_lim

species <- 1
cci <- 2
effort <- 1
dev.off()

# North
ssn <- summarize_avg_and_rep(avg_and_rep_spring_north[[species]][[cci]][[effort]], avg_and_rep_fall_north[[species]][[cci]][[effort]], 
                             year_remove=4, remove = T, min_small = 5, min_list = 1, ci = .8)
plot_one_summary(ssn, mp_spring_north[[species]][[cci]][[effort]], mp_fall_north[[species]][[cci]][[effort]])

# South
sss <- summarize_avg_and_rep(avg_and_rep_spring_south[[species]][[cci]][[effort]], avg_and_rep_fall_south[[species]][[cci]][[effort]],
                             year_remove=4, remove = T, min_small = 5, min_list = 1, ci = .8)
plot_one_summary(sss, mp_spring_south[[species]][[cci]][[effort]], mp_fall_south[[species]][[cci]][[effort]])

# Both
plot_ns_summary(ssn, sss)
