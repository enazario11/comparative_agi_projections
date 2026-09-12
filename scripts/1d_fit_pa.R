#libraries
library(tidyverse)
library(here)
library(sf)
library(terra)
library(rnaturalearth)
library(tidyquant)
library(aniMotum)

#TO DO 
#1) apply keep window filter to PA locs

#load re-routed SSMs
blu <- readRDS("data/loc_data/processed/ssm/blu_ssm.rds")
mako <- readRDS("data/loc_data/processed/ssm/mako_ssm.rds")
swo <- readRDS("data/loc_data/processed/ssm/swo_ssm.rds")

#get bounding polygon for MOM6 domain 
domain <- rast("../MI_AGI/data/enviro/nep/temp/processed/btemp_nep_median_rot.nc")
rast_mask <- !is.na(domain)
rast_bound <- as.polygons(rast_mask, values = FALSE, dissolve = TRUE)

  #get min/max lats across species to filter domain down to 
#albacore
alb <- readRDS("data/loc_data/alb_tag/validTagsLocns_NOAAonly_bathym_corrected.rds")
min(alb$lat) - 2 #23.1
max(alb$lat) + 2 #54.6

#blue sharks
blu_ssm <- grab(blu, what = "predicted")
min(blu_ssm$lat) - 2 #2.0
max(blu_ssm$lat) + 2 #52.7

#mako sharks
mako_ssm <- grab(mako, what = "predicted")
min(mako_ssm$lat) - 2 #0.9
max(mako_ssm$lat) + 2 #49.3

#swordfish
swo_ssm <- grab(swo, what = "predicted")
min(swo_ssm$lat) - 2 #5.1
max(swo_ssm$lat) + 2 #51.2

  #filter bounding shape
crop_lat <- ext(-180, 180, 0, 55)
bound_crop <- crop(rast_bound, crop_lat)
plot(bound_crop)
