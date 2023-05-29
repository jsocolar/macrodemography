
library(tidyverse)
library(visdat)

spring_temp <- readRDS("~/Documents/macrodemography/data/weather/spring_tmin_tmax2006_2019.rds")
head(spring_temp)
str(spring_temp)
vis_dat(spring_temp, warn_large_data = FALSE)

df <- spring_temp %>% mutate(year=year(ymd(date)), month=month(ymd(date)), day=yday(ymd(date))) %>% 
  group_by(cell_id, date) %>% mutate(tmean=(tmin+tmax)/2) %>% ungroup()
head(df)

# tmin
ggplot(df, aes(as.Date(date), tmin, color=cell_id))+
  geom_line(linewidth=.05)+theme(legend.position = "none")+
  facet_wrap(~year, nrow = 7, scales = "free_x")

# tmax
ggplot(df, aes(as.Date(date), tmax, color=cell_id))+
  geom_line(linewidth=.05)+theme(legend.position = "none")+
  facet_wrap(~year, nrow = 7, scales = "free_x")

# tmean
ggplot(df, aes(as.Date(date), tmean, color=cell_id))+
  geom_line(linewidth=.05)+theme(legend.position = "none")+
  facet_wrap(~year, nrow = 7, scales = "free_x")




