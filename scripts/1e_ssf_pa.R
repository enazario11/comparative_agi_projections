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
blu <- readRDS("data/loc_data/processed/pre_ssm/blu_dat.rds") %>% select(x = lon, y = lat, t = date, id = id) 
swo <- readRDS("data/loc_data/processed/pre_ssm/swo_dat.rds") %>% select(x = lon, y = lat, t = date, id = id) 

#make tracks
blu_track <- make_track(blu, .x = x, .y = y, .t = t, id = id, crs = 4326) %>% transform_coords(5070)
summarize_sampling_rate_many(blu_track, time_unit = "day", c("id")) %>% summarise(mean = mean(mean)) #2.4

swo_track <- make_track(swo, .x = x, .y = y, .t = t, id = id, crs = 4326) %>% transform_coords(5070)
summarize_sampling_rate_many(swo_track, time_unit = "day", c("id")) %>% summarise(mean = mean(mean)) #5.2

#resample tracks
#blue sharks
blu_track_resamp <- blu_track %>% nest(data = -"id") %>% mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = days(7), tolerance = days(1)) %>% filter_min_n_burst(min_n = 3) %>% steps_by_burst()))

blu_track_resamp |> select(id, steps) |> unnest(cols = steps) |> 
  ggplot(aes(sl_/1000, fill = factor(id))) + geom_density(alpha = 0.4) + xlab("step length (km)") + theme_bw()
blu_track_resamp |> select(id, steps) |> unnest(cols = steps) |> summarise(mean_step_km = mean(sl_, na.rm = TRUE)/1000) # 164 km avg SL (per 7 days)

#swo sharks
swo_track_resamp <- swo_track %>% nest(data = -"id") %>% mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = days(7), tolerance = days(1)) %>% filter_min_n_burst(min_n = 3) %>% steps_by_burst()))

swo_track_resamp |> select(id, steps) |> unnest(cols = steps) |> 
  ggplot(aes(sl_/1000, fill = factor(id))) + geom_density(alpha = 0.4) + xlab("step length (km)") + theme_bw()
swo_track_resamp |> select(id, steps) |> unnest(cols = steps) |> summarise(mean_step_km = mean(sl_, na.rm = TRUE)/1000) # 196 km avg SL (per 7 days)


# add random steps 
#blue sharks
  #make tracks
blu_list <- blu %>% nest(-id) 
blu_list <- blu_list %>%
  mutate(trk = lapply(data, function(d){
    make_track(d, x, y, t, 
               crs = 4326) %>%
              transform_coords(5070)
  }))

  #view sampling rate summary stats by tag id
blu_list %>% mutate(sr = lapply(trk, summarize_sampling_rate)) %>%
  select(id, sr) %>% unnest

  #resample and generate random steps
blu_pres_abs <- blu_list %>% 
                  mutate(steps = lapply(trk, function(x) { 
                          x %>% 
                            track_resample(rate = days(7), tolerance = days(1)) %>% 
                            filter_min_n_burst(min_n = 3) %>%
                            steps_by_burst()
                          }), 
                        pseudo = lapply(steps, function(x) {
                          x %>% random_steps(n = 10)
                        })                       
                  ) %>%
                  select(id, pseudo) %>% 
                  unnest(pseudo)

saveRDS(blu_pres_abs, here("data/loc_data/processed/ssf_pa/blu_pres_abs.rds"))

#swordfish
 #make tracks
swo_list <- swo %>% nest(-id) 
swo_list <- swo_list %>%
  mutate(trk = lapply(data, function(d){
    make_track(d, x, y, t, 
               crs = 4326) %>%
              transform_coords(5070)
  }))

  #view sampling rate summary stats by tag id
swo_list %>% mutate(sr = lapply(trk, summarize_sampling_rate)) %>%
  select(id, sr) %>% unnest() %>% print(n = 30)

  #resample and generate random steps -- add safeguards for step and random step generation due to burts with low n causing previous code to fail. Removes 3 IDs with low sample sizes.
make_steps <- possibly(function(x) {
  x %>%
    track_resample(rate = days(7), tolerance = days(1)) %>%
    filter_min_n_burst(min_n = 3) %>%
    steps_by_burst()
}, otherwise = NULL)

make_pseudo <- possibly(function(x) random_steps(x, n = 10), otherwise = NULL)

swo_pres_abs <- swo_list %>%
  mutate(steps = map(trk, make_steps)) %>%
  filter(!map_lgl(steps, is.null)) %>%
  mutate(pseudo = map(steps, make_pseudo)) %>%
  filter(!map_lgl(pseudo, is.null)) %>%
  select(id, pseudo) %>%
  unnest(pseudo)

swo_check <- swo_list %>%
  mutate(steps  = map(trk, make_steps),
         pseudo = map(steps, ~ if (is.null(.x)) NULL else make_pseudo(.x)))

swo_check %>% filter(map_lgl(pseudo, is.null)) %>% pull(id)

saveRDS(swo_pres_abs, here("data/loc_data/processed/ssf_pa/swo_pres_abs.rds"))
