#libraries
library(tidyverse)
library(here)
library(sf)
library(terra)
library(rnaturalearth)
library(tidyquant)
library(aniMotum)
source(here("functions/keep_windows.R"))
set.seed(0902)

#TO DO 
#4) get pseudo-absences to a 1:1 ratio.
#5) Run albacore SSM and then generate PAs

# land file
land <- ne_countries(scale = "large", returnclass = "sf") %>% st_make_valid() 
land <- st_transform(land, crs = 4326)

#load re-routed SSMs
blu <- readRDS("data/loc_data/processed/ssm/blu_ssm.rds")
mako <- readRDS("data/loc_data/processed/ssm/mako_ssm.rds")
swo <- readRDS("data/loc_data/processed/ssm/swo_ssm.rds")

#get min/max lats across species to filter domain down to 
#albacore
alb <- readRDS("data/loc_data/processed/pre_ssm/alb_dat.rds") %>%
  mutate(lon = ifelse(lon > 180, lon - 360, lon))
min(alb$lat) - 2 #23.1
max(alb$lat) + 2 #54.6
min(alb$lon) - 2 #-182.0

#blue sharks
blu_ssm <- grab(blu, what = "predicted")
min(blu_ssm$lat) - 2 #2.2
max(blu_ssm$lat) + 2 #52.7
min(blu_ssm$lon) -2 #-163.1

#mako sharks
mako_ssm <- grab(mako, what = "predicted")
min(mako_ssm$lat) - 2 #0.8
max(mako_ssm$lat) + 2 #49.2
min(mako_ssm$lon) - 2 #-157.3

#swordfish
swo_ssm <- grab(swo, what = "predicted")
min(swo_ssm$lat) - 2 #-5.3
max(swo_ssm$lat) + 2 #51.1
min(swo_ssm$lon) - 2 #-168.4

### PA generation ####
#### blue sharks ######
blu_pa <- sim_fit(blu, what = "predicted", reps = 100)
blu_pa_filt <- sim_filter(blu_pa, keep = 0.25, flag = 1) #flag based on hazen et al., 2017 journal of applied ecology
blu_pa_r <- route_path(blu_pa_filt, centroids = TRUE)

plot(blu_pa_r[1,])
saveRDS(blu_pa_r, here("data/loc_data/processed/pa/blu_pa_routed.rds"))

#### mako sharks #####
mako_pa <- sim_fit(mako, what = "predicted", reps = 100)
mako_pa_filt <- sim_filter(mako_pa, keep = 0.25, flag = 1) #flag based on hazen et al., 2017 journal of applied ecology
mako_pa_r <- route_path(mako_pa_filt, centroids = TRUE)

plot(mako_pa_r[4,])
saveRDS(mako_pa_r, here("data/loc_data/processed/pa/mako_pa_routed.rds"))

#### swordfish#####
swo_pa <- sim_fit(swo, what = "predicted", reps = 100)
swo_pa_filt <- sim_filter(swo_pa, keep = 0.25, flag = 1) #flag based on hazen et al., 2017 journal of applied ecology
swo_pa_r <- route_path(swo_pa_filt, centroids = TRUE)

plot(swo_pa_r[1,])
saveRDS(swo_pa_r, here("data/loc_data/processed/pa/swo_pa_routed.rds"))

#### remove locations in gap windows, on land, and outside of MOM6 domain ####
#combine rerouted presences and PAs into one df, then make locs a spatial feature
#albacore 


#blue sharks
blu_locs <- readRDS(here("data/loc_data/processed/ssm/blu_ssm.rds")) %>% 
  grab(what = "rerouted") %>%
  mutate(loc_type = "presence", 
         rep = NA) %>%
  select(-c("x", "y", "x.se", "y.se"))
blu_pa <- readRDS(here("data/loc_data/processed/pa/blu_pa_routed.rds")) %>% 
  unnest(sims) %>%
  mutate(loc_type = "absence") %>%
  select(-c("x", "y", "model"))

blu <- rbind(blu_locs, blu_pa) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)

#mako sharks
mako_locs <- readRDS(here("data/loc_data/processed/ssm/mako_ssm.rds")) %>% 
  grab(what = "rerouted") %>%
  mutate(loc_type = "presence", 
         rep = NA) %>%
  select(-c("x", "y", "x.se", "y.se"))
mako_pa <- readRDS(here("data/loc_data/processed/pa/mako_pa_routed.rds")) %>% 
  unnest(sims) %>%
  mutate(loc_type = "absence") %>%
  select(-c("x", "y", "model"))

mako <- rbind(mako_locs, mako_pa) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)

#swordfish
swo_locs <- readRDS(here("data/loc_data/processed/ssm/swo_ssm.rds")) %>% 
  grab(what = "rerouted") %>%
  mutate(loc_type = "presence", 
         rep = NA) %>%
  select(-c("x", "y", "x.se", "y.se"))
swo_pa <- readRDS(here("data/loc_data/processed/pa/swo_pa_routed.rds")) %>% 
  unnest(sims) %>%
  mutate(loc_type = "absence") %>%
  select(-c("x", "y", "model"))

swo <- rbind(swo_locs, swo_pa) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)

#### Gap windows #####
#albacore


#blue sharks
blu_raw <- readRDS("data/loc_data/processed/pre_ssm/blu_dat.rds") %>% filter(lc != "P" & lc != "D")
colnames(blu_raw) <- c("id", "date", "lc", "sp", "lon", "lat")
blu_windows <- keep_windows(blu_raw)

blu_no_gaps <- blu %>%
  inner_join(bind_rows(blu_windows), 
             by = join_by(id, between(date, start_time, end_time)))

#mako sharks
mako_raw <- readRDS("data/loc_data/processed/pre_ssm/mako_dat.rds") %>% filter(lc != "P" & lc != "D")
colnames(mako_raw) <- c("id", "date", "lc", "sp", "lon", "lat")
mako_windows <- keep_windows(mako_raw)

mako_no_gaps <- mako %>%
  inner_join(bind_rows(mako_windows), 
             by = join_by(id, between(date, start_time, end_time)))


#swordfish
swo_raw <- readRDS("data/loc_data/processed/pre_ssm/swo_dat.rds") %>% filter(lc != "P" & lc != "D")
colnames(swo_raw) <- c("id", "date", "lc", "sp", "lon", "lat")
swo_windows <- keep_windows(swo_raw)

swo_no_gaps <- swo %>%
  inner_join(bind_rows(swo_windows), 
             by = join_by(id, between(date, start_time, end_time)))


#### Land #####
#albacore 


#blue sharks
blu_no_land <- st_difference(blu_no_gaps, land)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = blu_no_land, aes(color = loc_type)) + 
    coord_sf(xlim = c(-170, -100),
      ylim = c(0, 55),
      expand = FALSE) +
    theme_bw() 

#mako sharks


#swordfish


#### MOM6 #####
#get domain and turn into sf object
mom_domain <- rast(here("../MI_AGI/data/enviro/nep/temp/raw/tob.nep.full.hcast.monthly.regrid.r20250912.199301-202506"))
mom_domain <- mom_domain[1]

mom_poly <- st_as_sf(as.polygons(mom_domain))

#albacore


#blue sharks
blu_mom <- st_intersection(blu_no_land, mom_poly)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = blu_mom, aes(color = loc_type)) + 
    coord_sf(xlim = c(-170, -100),
      ylim = c(0, 55),
      expand = FALSE) +
    theme_bw() 

#mako sharks



#swordfish


#### randomly sample PAs to get 1:1 with presences #####




