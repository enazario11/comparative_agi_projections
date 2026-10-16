#libraries
library(tidyverse)
library(here)
library(sf)
library(terra)
library(rnaturalearth)
library(tidyquant)
library(amt)
set.seed(0902)

#load regularized data
alb <- readRDS("data/loc_data/alb_tag/validTagsLocns_NOAAonly_bathym_corrected.rds") %>%
  mutate(lon = ifelse(lon360 > 180, lon360 - 360, lon360), 
         dateRd = as.POSIXct(dateRd, format = "%Y-%m-%d %H:%M:%S")) %>%
  select(x = "lon", y = "lat", t = "dateRd", id = "tag") 

blu <- readRDS("data/loc_data/processed/pre_ssm/blu_dat.rds") %>% filter(lc != "P" & lc != "D") %>%
  mutate(lon = ifelse(lon > 180, lon - 360, lon), 
         date = as.POSIXct(date)) %>%
  select(x = "lon", y = "lat", t = "date", id = "tag") 

mako <- readRDS("data/loc_data/processed/pre_ssm/mako_dat.rds") %>%
  mutate(lon = ifelse(lon > 180, lon - 360, lon), 
         date = as.POSIXct(date)) %>%
  select(x = "lon", y = "lat", t = "date", id = "tag") 

swo <- readRDS("data/loc_data/processed/pre_ssm/swo_dat.rds") %>%
  mutate(lon = ifelse(lon > 180, lon - 360, lon), 
         date = as.POSIXct(date)) %>%
  select(x = "lon", y = "lat", t = "date", id = "tag") 

swo_sf <- swo %>% st_as_sf(coords = c("x", "y"), crs = 4326, remove = FALSE)

bbox <- st_polygon(list(matrix(c(
  -170, 0,
  -100, 0,
  -100, 55,
  -170, 55,
  -170, 0
), ncol = 2, byrow = TRUE))) %>% 
  st_sfc(crs = 4326) 

swo_filt <- st_intersection(swo_sf, bbox)

swo_coords <- st_coordinates(swo_filt)


#make tracks
alb_track <- make_track(alb, .x = x, .y = y, .t = t, id = id, crs = 4326) %>% transform_coords(5070) #example on amt reprojects to north american datum (epsg 5070)
summarize_sampling_rate_many(alb_track, time_unit = "day", c("id")) %>% summarise(mean = mean(mean))

blu_track <- make_track(blu, .x = x, .y = y, .t = t, id = id, crs = 4326) %>% transform_coords(5070)
summarize_sampling_rate_many(blu_track, time_unit = "day", c("id")) %>% summarise(mean = mean(mean))

mako_track <- make_track(mako, .x = x, .y = y, .t = t, id = id, crs = 4326) %>% transform_coords(5070)
summarize_sampling_rate_many(mako_track, time_unit = "day", c("id")) %>% summarise(mean = mean(mean))

#resample tracks
#albacore
alb_track_resamp <- alb_track %>% nest(data = -"id") %>% mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = days(1), tolerance = days(5)) %>% filter_min_n_burst(min_n = 3) %>% steps_by_burst()))

alb_track_resamp |> select(id, steps) |> unnest(cols = steps) |> 
  ggplot(aes(sl_/1000, fill = factor(id))) + geom_density(alpha = 0.4) + xlab("step length (km)")
alb_track_resamp |> select(id, steps) |> unnest(cols = steps) |> summarise(mean_step_km = mean(sl_, na.rm = TRUE)/1000)

#blue sharks
blu_track_resamp <- blu_track %>% nest(data = -"id") %>% mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = days(2), tolerance = days(5)) %>% filter_min_n_burst(min_n = 3) %>% steps_by_burst()))

blu_track_resamp |> select(id, steps) |> unnest(cols = steps) |> 
  ggplot(aes(sl_/1000, fill = factor(id))) + geom_density(alpha = 0.4) + xlab("step length (km)")
blu_track_resamp |> select(id, steps) |> unnest(cols = steps) |> summarise(mean_step_km = mean(sl_, na.rm = TRUE)/1000)

#mako sharks
mako_track_resamp <- mako_track %>% nest(data = -"id") %>% mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = days(2), tolerance = days(5)) %>% filter_min_n_burst(min_n = 3) %>% steps_by_burst()))

mako_track_resamp |> select(id, steps) |> unnest(cols = steps) |> 
  ggplot(aes(sl_/1000, fill = factor(id))) + geom_density(alpha = 0.4) + xlab("step length (km)")
mako_track_resamp |> select(id, steps) |> unnest(cols = steps) |> summarise(mean_step_km = mean(sl_, na.rm = TRUE)/1000)


# add random steps 
alb_pres_abs <- alb_track_resamp %>% random_steps(n = 10)
blu_pres_abs <- blu_track_resamp %>% random_steps(n = 10)
mako_pres_abs <- mako_track_resamp %>% random_steps(n = 10)
