#libraries
library(tidyverse)
library(here)
library(sf)
library(terra)
library(rnaturalearth)
library(tidyquant)
library(aniMotum)
set.seed(0902)

#TO DO 
#1) apply keep window filter to PA locs
#2) remove presences on land and the corresponding absences by date
#3) remove presences and pseudo-absences outside of the MOM6 domain
#4) get pseudo-absences to a 1:1 ratio.
#5) Check process for albacore? Run SSM and then generate PAs?

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
blu_pa <- sim_fit(blu[1,], what = "predicted", reps = 100)
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
swo_pa <- sim_fit(swo[1,], what = "predicted", reps = 100)
swo_pa_filt <- sim_filter(swo_pa, keep = 0.25, flag = 1) #flag based on hazen et al., 2017 journal of applied ecology
swo_pa_r <- route_path(swo_pa_filt, centroids = TRUE)

plot(swo_pa_r[1,])
saveRDS(swo_pa_r, here("data/loc_data/processed/pa/swo_pa_routed.rds"))




