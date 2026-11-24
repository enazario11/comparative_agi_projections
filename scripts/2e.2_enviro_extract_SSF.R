#libraries
library(tidyverse)
library(here)
library(terra)
library(sf)
library(foreach)
library(doParallel)
source("functions/enviro_extract_ssf.R")

#load location data (presence and absence, 10:1 ratio of PA to presence). Filter any last locations outside of study domain.
#filter out locations after October 22, 2024 because the physics data are not available after that date...
#removed the start coordinates (x1 and y1) and start time (t1). Primarily cared about where they ended up in each 7 day block.
blu <- readRDS("data/loc_data/processed/ssf_pa/blu_pres_abs.rds") %>% 
  select(-c(x1_, y1_, t1_)) %>% 
  filter(t2_ <= as.Date("2024-10-22")) %>% 
  st_as_sf(coords = c("x2_", "y2_"), crs = 5070) %>%
  st_transform(crs = 4326) %>%
  rename(date = t2_) %>%
  mutate(loc_type = ifelse(case_ == TRUE, "presence", "absence"), 
         nearest_yw_start = round_date(date, unit = "week", week_start = 4), #week anchored on Thursday to minimize week distance for nearest match when matching raster layers
         nearest_yw = format(nearest_yw_start, "%G-%V")) %>%
  select(-c(case_, nearest_yw_start))

swo <- readRDS("data/loc_data/processed/ssf_pa/swo_pres_abs.rds") %>% 
  select(-c(x1_, y1_, t1_)) %>% 
  filter(t2_ <= as.Date("2024-10-22"))  %>% 
  st_as_sf(coords = c("x2_", "y2_"), crs = 5070) %>%
  st_transform(crs = 4326) %>%
  rename(date = t2_) %>%
  mutate(loc_type = ifelse(case_ == TRUE, "presence", "absence"), 
         nearest_yw_start = round_date(date, unit = "week", week_start = 4), #week anchored on Thursday to minimize week distance for nearest match when matching raster layers
         nearest_yw = format(nearest_yw_start, "%G-%V")) %>%
  select(-c(case_, nearest_yw_start))

#confirm there are no grid cell repeats for each year-week location
#blue sharks
check_blu <- data.frame(x = blu$geometry[, 1], y = blu$geometry[, 2], iso_yw = blu$nearest_yw)
any(duplicated(check_blu))

#swordfish
check_swo <- data.frame(x = swo$geometry[, 1], y = swo$geometry[, 2], iso_yw = swo$nearest_yw)
any(duplicated(check_swo))

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
  library(ISOweek)
})

#extract all depth x variable combinations for each species
#blue
for (i in 1:nrow(depth_var_vect_all)){

    depth = depth_var_vect_all$depth[i]
    var = depth_var_vect_all$var[i]

    if(var != "AGI"){ 
      enviro_folder = "data/enviro/merc/regrid/"
    } else if (var == "AGI") {
      enviro_folder = "data/enviro/agi/blu"
    }

    enviro_col_out <- getvarGLORYS_par_ssf(loc_dat = blu, enviro_folder = enviro_folder, depth = depth, col_name = var)
    blu[[paste0(var,"_", depth,"m")]] <- enviro_col_out[[1]] 
}

saveRDS(blu, here("data/loc_enviro/SSF/blu_enviro.rds"))

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

    enviro_col_out <- getvarGLORYS_par_ssf(loc_dat = swo_day, enviro_folder = enviro_folder, depth = depth, col_name = var)
    swo[[paste0(var,"_", depth,"m")]] <- enviro_col_out[[1]] 
}

saveRDS(swo_day, here("data/loc_enviro/SSF/swo_day_enviro.rds"))

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

    enviro_col_out <- getvarGLORYS_par_ssf(loc_dat = swo_night, enviro_folder = enviro_folder, depth = depth, col_name = var)
    swo_night[[paste0(var,"_", depth,"m")]] <- enviro_col_out[[1]] 
}

saveRDS(swo_night, here("data/loc_enviro/SSF/swo_night_enviro.rds"))

# add bathymetry and rugosity 
blu <- readRDS(here("data/loc_enviro/blu_enviro.rds"))
swo_day <- readRDS(here("data/loc_enviro/swo_day_enviro.rds"))
swo_night <- readRDS(here("data/loc_enviro/swo_night_enviro.rds"))

bathy <- rast(here("data/enviro/bathy/gebco_2026_n55.0_s0.0_w-170.0_e-100.0.nc")) %>% clamp(upper = 0, values = FALSE)
# writeCDF(bathy, here("data/enviro/bathy/bathy_sub_0m.nc"))

agg_factor <- round(0.25/res(bathy)[1])

#blue
blu_bathy_out <- getbathy_par_ssf(loc_dat = blu, agg_factor = agg_factor)
blu[["bathy"]] <- blu_bathy_out[[1]] 
blu[["rugosity"]] <- blu_bathy_out[[2]] 

saveRDS(blu, here("data/loc_enviro/SSF/blu_enviro_bathy.rds"))

#swordfish day
swo_day_bathy_out <- getbathy_par_ssf(loc_dat = swo_day, agg_factor = agg_factor)
swo_day[["bathy"]] <- swo_day_bathy_out[[1]] 
swo_day[["rugosity"]] <- swo_day_bathy_out[[2]] 

saveRDS(swo_day, here("data/loc_enviro/swo_day_enviro_bathy.rds"))

#swordfish night
swo_night_bathy_out <- getbathy_par_ssf(loc_dat = swo_night, agg_factor = agg_factor)
swo_night[["bathy"]] <- swo_night_bathy_out[[1]] 
swo_night[["rugosity"]] <- swo_night_bathy_out[[2]] 

saveRDS(swo_night, here("data/loc_enviro/SSF/swo_night_enviro_bathy.rds"))
