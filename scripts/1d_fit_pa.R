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
# 1) finalize study domain boundary to filter locs and PAs 
#2) get pseudo-absences to a 1:1 ratio.
#3) Run albacore SSM and then generate PAs

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

#### remove locations in gap windows, on land, and outside of study domain ####
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
land_union <- st_union(land) #speeds up st_filter

#albacore 


#blue sharks
blu_no_land <- st_filter(blu_no_gaps, land_union, .predicate = st_disjoint)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = blu_no_land, aes(color = loc_type)) + 
    coord_sf(xlim = c(-170, -100),
      ylim = c(0, 55),
      expand = FALSE) +
    theme_bw() 

#mako sharks
mako_no_land <- st_filter(mako_no_gaps, land_union, .predicate = st_disjoint)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = mako_no_land, aes(color = loc_type)) + 
    coord_sf(xlim = c(-170, -100),
      ylim = c(0, 55),
      expand = FALSE) +
    theme_bw() 


#swordfish
swo_no_land <- st_filter(swo_no_gaps, land_union, .predicate = st_disjoint)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = swo_no_land, aes(color = loc_type)) + 
    coord_sf(xlim = c(-170, -100),
      ylim = c(0, 55),
      expand = FALSE) +
    theme_bw() 

#### study domain #####
#get domain and turn into sf object -- Update once have domain boundary
# mom_domain <- rast(here("../MI_AGI/data/enviro/nep/temp/raw/tob.nep.full.hcast.monthly.regrid.r20250912.199301-202506"))
# mom_domain <- mom_domain[1]

# mom_poly <- st_as_sf(as.polygons(mom_domain))

bbox <- st_bbox(c(xmin = -170, ymin = 0, xmax = -100, ymax = 55), crs = 4326)

#albacore


#blue sharks
blu_domain <- st_filter(blu_no_land, st_as_sfc(bbox), .predicate = st_within)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = blu_domain, aes(color = loc_type)) + 
    coord_sf(xlim = c(-175, -98),
      ylim = c(-5, 60),
      expand = FALSE) +
    theme_bw() 

#mako sharks
mako_domain <- st_filter(mako_no_land, st_as_sfc(bbox), .predicate = st_within)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = mako_domain, aes(color = loc_type)) + 
    coord_sf(xlim = c(-175, -98),
      ylim = c(-5, 60),
      expand = FALSE) +
    theme_bw() 

#swordfish
swo_domain <- st_filter(swo_no_land, st_as_sfc(bbox), .predicate = st_within)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = swo_domain, aes(color = loc_type)) + 
    coord_sf(xlim = c(-175, -98),
      ylim = c(-5, 60),
      expand = FALSE) +
    theme_bw() 

#### randomly sample PAs to get 1:1 with presences #####
pa_ratio <- function(sp_dat, ratio = 1){

  pres_abs_df <- data.frame()

for(i in 1:length(unique(sp_dat$id))){

  curr_id = unique(sp_dat$id)[i]
  temp_dat = sp_dat %>% filter(id == curr_id)
  num_presence = temp_dat %>% filter(loc_type == "presence") %>% summarise(n = n())
  num_presence = num_presence$n

  pres_df = temp_dat %>% filter(loc_type == "presence")
  abs_df = temp_dat %>% filter(loc_type == "absence")
  abs_df_sub = abs_df[sample(nrow(abs_df), size = num_presence*ratio, replace = FALSE), ]

  pres_abs_temp <- rbind(pres_df, abs_df_sub)

  pres_abs_df <- rbind(pres_abs_df, pres_abs_temp)

} #end for loop
  return(pres_abs_df)
} #end function

#albacore 


#blue sharks
blu_pres_abs <- pa_ratio(blu_domain)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = blu_pres_abs, aes(color = loc_type), size = 2, alpha = 0.8) + 
    coord_sf(xlim = c(-175, -98),
      ylim = c(-5, 60),
      expand = FALSE) +
    theme_bw() +
  scale_color_manual(values = c("#77ABD9", "dodgerblue4"))

saveRDS(blu_pres_abs, here("data/loc_data/processed/pres_abs/blu_pres_abs.rds"))

#mako sharks
mako_pres_abs <- pa_ratio(mako_domain)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = mako_pres_abs, aes(color = loc_type), size = 2, alpha = 0.8) + 
    coord_sf(xlim = c(-175, -98),
      ylim = c(-5, 60),
      expand = FALSE) +
    theme_bw() +
  scale_color_manual(values = c("#77ABD9", "dodgerblue4"))

saveRDS(mako_pres_abs, here("data/loc_data/processed/pres_abs/mako_pres_abs.rds"))

#swordfish
swo_pres_abs <- pa_ratio(swo_domain)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = swo_pres_abs, aes(color = loc_type), size = 2, alpha = 0.8) + 
    coord_sf(xlim = c(-175, -98),
      ylim = c(-5, 60),
      expand = FALSE) +
    theme_bw() +
  scale_color_manual(values = c("#77ABD9", "dodgerblue4"))

saveRDS(swo_pres_abs, here("data/loc_data/processed/pres_abs/swo_pres_abs.rds"))


