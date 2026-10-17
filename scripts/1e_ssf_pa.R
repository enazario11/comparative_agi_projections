#libraries
library(tidyverse)
library(here)
library(sf)
library(terra)
library(rnaturalearth)
library(tidyquant)
library(amt)
set.seed(0902)

#load cleaned data (pre SSM fitting for blu, mako, swo)
alb <- readRDS("data/loc_data/processed/pre_ssm/alb_dat.rds") %>% select(x = lon, y = lat, t = date, id = id) 
blu <- readRDS("data/loc_data/processed/pre_ssm/blu_dat.rds") %>% select(x = lon, y = lat, t = date, id = id) 
mako <- readRDS("data/loc_data/processed/pre_ssm/mako_dat.rds") %>% select(x = lon, y = lat, t = date, id = id) 
swo <- readRDS("data/loc_data/processed/pre_ssm/swo_dat.rds") %>% select(x = lon, y = lat, t = date, id = id) 


#make tracks
alb_track <- make_track(alb, .x = x, .y = y, .t = t, id = id, crs = 4326) %>% transform_coords(5070) #example on amt reprojects to north american datum (epsg 5070)
summarize_sampling_rate_many(alb_track, time_unit = "day", c("id")) %>% summarise(mean = mean(mean)) #1.02

blu_track <- make_track(blu, .x = x, .y = y, .t = t, id = id, crs = 4326) %>% transform_coords(5070)
summarize_sampling_rate_many(blu_track, time_unit = "day", c("id")) %>% summarise(mean = mean(mean)) #2.4

mako_track <- make_track(mako, .x = x, .y = y, .t = t, id = id, crs = 4326) %>% transform_coords(5070)
summarize_sampling_rate_many(mako_track, time_unit = "day", c("id")) %>% summarise(mean = mean(mean)) #1.90

swo_track <- make_track(swo, .x = x, .y = y, .t = t, id = id, crs = 4326) %>% transform_coords(5070)
summarize_sampling_rate_many(swo_track, time_unit = "day", c("id")) %>% summarise(mean = mean(mean)) #5.2

#resample tracks
#albacore
alb_track_resamp <- alb_track %>% nest(data = -"id") %>% mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = days(5), tolerance = days(3)) %>% filter_min_n_burst(min_n = 3) %>% steps_by_burst()))

alb_track_resamp |> select(id, steps) |> unnest(cols = steps) |> 
  ggplot(aes(sl_/1000, fill = factor(id))) + geom_density(alpha = 0.4) + xlab("step length (km)") + theme_bw()
alb_track_resamp |> select(id, steps) |> unnest(cols = steps) |> summarise(mean_step_km = mean(sl_, na.rm = TRUE)/1000) #80.1 km avg SL (per 5 days)

#blue sharks
blu_track_resamp <- blu_track %>% nest(data = -"id") %>% mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = days(5), tolerance = days(3)) %>% filter_min_n_burst(min_n = 3) %>% steps_by_burst()))

blu_track_resamp |> select(id, steps) |> unnest(cols = steps) |> 
  ggplot(aes(sl_/1000, fill = factor(id))) + geom_density(alpha = 0.4) + xlab("step length (km)") + theme_bw()
blu_track_resamp |> select(id, steps) |> unnest(cols = steps) |> summarise(mean_step_km = mean(sl_, na.rm = TRUE)/1000) #75.3 km avg SL (per 5 days)

#mako sharks
mako_track_resamp <- mako_track %>% nest(data = -"id") %>% mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = days(5), tolerance = days(3)) %>% filter_min_n_burst(min_n = 3) %>% steps_by_burst()))

mako_track_resamp |> select(id, steps) |> unnest(cols = steps) |> 
  ggplot(aes(sl_/1000, fill = factor(id))) + geom_density(alpha = 0.4) + xlab("step length (km)") + theme_bw()
mako_track_resamp |> select(id, steps) |> unnest(cols = steps) |> summarise(mean_step_km = mean(sl_, na.rm = TRUE)/1000) #84.8 km avg SL (per 5 days)

#swo sharks
swo_track_resamp <- swo_track %>% nest(data = -"id") %>% mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = days(5), tolerance = days(3)) %>% filter_min_n_burst(min_n = 3) %>% steps_by_burst()))

swo_track_resamp |> select(id, steps) |> unnest(cols = steps) |> 
  ggplot(aes(sl_/1000, fill = factor(id))) + geom_density(alpha = 0.4) + xlab("step length (km)") + theme_bw()
swo_track_resamp |> select(id, steps) |> unnest(cols = steps) |> summarise(mean_step_km = mean(sl_, na.rm = TRUE)/1000) # 101 km avg SL (per 5 days)


# add random steps 
alb_pres_abs <- alb_track_resamp %>% random_steps(n = 10)
blu_pres_abs <- blu_track_resamp %>% random_steps(n = 10)
mako_pres_abs <- mako_track_resamp %>% random_steps(n = 10)
swo_pres_abs <- swo_track_resamp %>% random_steps(n = 10)
