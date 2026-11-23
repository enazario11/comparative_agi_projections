#libraries
library(tidyverse)
library(here)
library(terra)
library(sf)
library(foreach)
library(doParallel)
source("functions/enviro_extract.R")

#load location data (presence and absence, 10:1 ratio of PA to presence). Filter any last locations outside of study domain.
#filter out locations after October 22, 2024 because the physics data are not available after that date...
alb <- readRDS("data/loc_data/processed/pres_abs/alb_pres_abs.rds") %>% filter(date <= as.Date("2024-10-22"))
blu <- readRDS("data/loc_data/processed/pres_abs/blu_pres_abs.rds") %>% filter(date <= as.Date("2024-10-22")) 
mako <- readRDS("data/loc_data/processed/pres_abs/mako_pres_abs.rds") %>% filter(date <= as.Date("2024-10-22"))
swo <- readRDS("data/loc_data/processed/pres_abs/swo_pres_abs.rds") %>% filter(date <= as.Date("2024-10-22"))

#create vector of all varaible depth combinations
depth_vect <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 150, #remove 140 b/c would duplicate the 133.08 layer
             200, 300, 400, 500, 600, 700, 800, 900, 1000)
depth_vect_2d = 0

var_vect = c("AGI", "o2")
var_vect_2d = c("T", "S", "MLD", "chl", "su", "sustr", "sv", "svstr")

depth_var_vect = expand.grid(depth = depth_vect, var = var_vect)
depth_var_vect_2d = expand.grid(depth = depth_vect_2d, var = var_vect_2d)

depth_var_vect_all = rbind(depth_var_vect, depth_var_vect_2d)

#parallelize
cl <- makeCluster(5)
clusterEvalQ(cl, {
  library(tidyverse)
  library(terra)
  library(here)
  library(sf)
})

#extract all depth x variable combinations for each species
#albacore
for (i in 1:nrow(depth_var_vect_all)){

    depth = depth_var_vect_all$depth[i]
    var = depth_var_vect_all$var[i]

    if(var != "AGI"){ 
      enviro_folder = "data/enviro/merc/regrid/"
    } else if (var == "AGI") {
      enviro_folder = "data/enviro/agi/alb"
    }

    enviro_col_out <- getvarGLORYS_par(loc_dat = alb, enviro_folder = enviro_folder, depth = depth, col_name = var)
    alb[[paste0(var,"_", depth,"m")]] <- enviro_col_out[[1]] 
}

saveRDS(alb, here("data/loc_enviro/alb_enviro.rds"))

#blue
for (i in 1:nrow(depth_var_vect_all)){

    depth = depth_var_vect_all$depth[i]
    var = depth_var_vect_all$var[i]

    if(var != "AGI"){ 
      enviro_folder = "data/enviro/merc/regrid/"
    } else if (var == "AGI") {
      enviro_folder = "data/enviro/agi/blu"
    }

    enviro_col_out <- getvarGLORYS_par(loc_dat = blu, enviro_folder = enviro_folder, depth = depth, col_name = var)
    blu[[paste0(var,"_", depth,"m")]] <- enviro_col_out[[1]] 
}

saveRDS(blu, here("data/loc_enviro/blu_enviro.rds"))

#mako 
for (i in 1:nrow(depth_var_vect_all)){

    depth = depth_var_vect_all$depth[i]
    var = depth_var_vect_all$var[i]

    if(var != "AGI"){ 
      enviro_folder = "data/enviro/merc/regrid/"
    } else if (var == "AGI") {
      enviro_folder = "data/enviro/agi/mako"
    }

    enviro_col_out <- getvarGLORYS_par(loc_dat = mako, enviro_folder = enviro_folder, depth = depth, col_name = var)
    mako[[paste0(var,"_", depth,"m")]] <- enviro_col_out[[1]] 
}

saveRDS(mako, here("data/loc_enviro/mako_enviro.rds"))

#swordfish
  #day
swo_day <- swo
for (i in 1:nrow(depth_var_vect_all)){

    depth = depth_var_vect_all$depth[i]
    var = depth_var_vect_all$var[i]

    if(var != "AGI"){ 
      enviro_folder = "data/enviro/merc/regrid/"
    } else if (var == "AGI") {
      enviro_folder = "data/enviro/agi/swo_day"
    }

    enviro_col_out <- getvarGLORYS_par(loc_dat = swo_day, enviro_folder = enviro_folder, depth = depth, col_name = var)
    swo[[paste0(var,"_", depth,"m")]] <- enviro_col_out[[1]] 
}

saveRDS(swo_day, here("data/loc_enviro/swo_day_enviro.rds"))

  #night
swo_night <- swo
for (i in 1:nrow(depth_var_vect_all)){

    depth = depth_var_vect_all$depth[i]
    var = depth_var_vect_all$var[i]

    if(var != "AGI"){ 
      enviro_folder = "data/enviro/merc/regrid/"
    } else if (var == "AGI") {
      enviro_folder = "data/enviro/agi/swo_night"
    }

    enviro_col_out <- getvarGLORYS_par(loc_dat = swo_night, enviro_folder = enviro_folder, depth = depth, col_name = var)
    swo_night[[paste0(var,"_", depth,"m")]] <- enviro_col_out[[1]] 
}

saveRDS(swo_night, here("data/loc_enviro/swo_night_enviro.rds"))

# add bathymetry and rugosity 
alb <- readRDS(here("data/loc_enviro/alb_enviro.rds"))
blu <- readRDS(here("data/loc_enviro/blu_enviro.rds"))
mako <- readRDS(here("data/loc_enviro/mako_enviro.rds"))
swo_day <- readRDS(here("data/loc_enviro/swo_day_enviro.rds"))
swo_night <- readRDS(here("data/loc_enviro/swo_night_enviro.rds"))

bathy <- rast(here("data/enviro/bathy/gebco_2026_n55.0_s0.0_w-170.0_e-100.0.nc")) %>% clamp(upper = 0, values = FALSE)
# writeCDF(bathy, here("data/enviro/bathy/bathy_sub_0m.nc"))

agg_factor <- round(0.25/res(bathy)[1])

#albacore
alb_bathy_out <- getbathy_par(loc_dat = alb, agg_factor = agg_factor)
alb[["bathy"]] <- alb_bathy_out[[1]] 
alb[["rugosity"]] <- alb_bathy_out[[2]] 

saveRDS(alb, here("data/loc_enviro/alb_enviro_bathy.rds"))

#blue
blu_bathy_out <- getbathy_par(loc_dat = blu, agg_factor = agg_factor)
blu[["bathy"]] <- blu_bathy_out[[1]] 
blu[["rugosity"]] <- blu_bathy_out[[2]] 

saveRDS(blu, here("data/loc_enviro/blu_enviro_bathy.rds"))

#mako
mako_bathy_out <- getbathy_par(loc_dat = mako, agg_factor = agg_factor)
mako[["bathy"]] <- mako_bathy_out[[1]] 
mako[["rugosity"]] <- mako_bathy_out[[2]] 

saveRDS(mako, here("data/loc_enviro/mako_enviro_bathy.rds"))

#swordfish day
swo_day_bathy_out <- getbathy_par(loc_dat = swo_day, agg_factor = agg_factor)
swo_day[["bathy"]] <- swo_day_bathy_out[[1]] 
swo_day[["rugosity"]] <- swo_day_bathy_out[[2]] 

saveRDS(swo_day, here("data/loc_enviro/swo_day_enviro_bathy.rds"))

#swordfish night
swo_night_bathy_out <- getbathy_par(loc_dat = swo_night, agg_factor = agg_factor)
swo_night[["bathy"]] <- swo_night_bathy_out[[1]] 
swo_night[["rugosity"]] <- swo_night_bathy_out[[2]] 

saveRDS(swo_night, here("data/loc_enviro/swo_night_enviro_bathy.rds"))
