#libraries
library(tidyverse)
library(here)
library(terra)
library(sf)

# load location data (presence and absence, 10:1 ratio of PA to presence). Filter any last locations outside of study domain.
alb <- readRDS("data/loc_data/processed/pres_abs/alb_pres_abs.rds")
blu <- readRDS("data/loc_data/processed/pres_abs/blu_pres_abs.rds") 
mako <- readRDS("data/loc_data/processed/pres_abs/mako_pres_abs.rds") 
swo <- readRDS("data/loc_data/processed/pres_abs/swo_pres_abs.rds") 

swo_test <- getvarGLORYS(loc_dat = swo, enviro_folder = "data/enviro/merc/regrid/chl", depth = "_0m", col_name = "chl")

