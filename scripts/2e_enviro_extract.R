#libraries
library(tidyverse)
library(here)
library(terra)
library(sf)
library(foreach)
library(doParallel)
source("functions/enviro_extract.R")

#load location data (presence and absence, 10:1 ratio of PA to presence). Filter any last locations outside of study domain.
#filter out locations after October 22, 2024 because the physics data is not available after that date...
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

clusterExport(cl, varlist = c("depth_var_vect_all", "alb", "blu", "mako", "swo", "getvarGLORYS"))
registerDoParallel(cl)

#extract all depth x variable combinations for each species
#albacore
alb_enviro_temp <- data.frame(x = 1:nrow(alb))
clusterExport(cl, varlist = c("depth_var_vect_all", "alb", "alb_enviro_temp", "getvarGLORYS"))
registerDoParallel(cl)

alb_enviro <- foreach (i=1:nrow(depth_var_vect_all), .combine = cbind) %dopar% {

    depth = depth_var_vect_all$depth[i]
    var = depth_var_vect_all$var[i]

    cat(sprintf("Completed %d\n", i), file = "progress.txt", append = TRUE)

    if(var != "AGI"){ 
      enviro_folder = "data/enviro/merc/regrid/"
    } else if (var == "AGI") {
      enviro_folder = "data/enviro/agi/alb"
    }

    alb_enviro_col <- getvarGLORYS(loc_dat = alb, enviro_folder = enviro_folder, depth = depth, col_name = var)
    alb_enviro_temp[paste0(var, "_", depth, "m")] <- alb_enviro_col[,1]
  }

stopCluster(cl)
alb_enviro_df <- cbind(alb, alb_enviro)
saveRDS(alb_enviro_df, here("data/loc_enviro/alb_enviro.rds"))

#blue
for (i in 1:nrow(depth_var_vect_all)){

  print(depth_var_vect_all[i,])

  depth = depth_var_vect_all$depth[i]
  var = depth_var_vect_all$var[i]

  if(var != "AGI"){ 
    enviro_folder = "data/enviro/merc/regrid/"
  } else if (var == "AGI") {
    enviro_folder = "data/enviro/agi/blu"
  }

  if (i == 1){
    blu_enviro <- getvarGLORYS(loc_dat = blu, enviro_folder = enviro_folder, depth = depth, col_name = var, depth_var_num = i)
  } else {
    blu_enviro_temp <- getvarGLORYS(loc_dat = blu_enviro, enviro_folder = enviro_folder, depth = depth, col_name = var, depth_var_num = i)
    blu_enviro[paste0(var, "_", depth, "m")] <- blu_enviro_temp[,1]
    blu_enviro <- blu_enviro %>% relocate(geometry, .before = 6)
  }
}

writeRDS(blu_enviro, here("data/loc_enviro/blu_enviro.rds"))

#mako 
for (i in 1:nrow(depth_var_vect_all)){

  print(depth_var_vect_all[i,])

  depth = depth_var_vect_all$depth[i]
  var = depth_var_vect_all$var[i]

  if(var != "AGI"){ 
    enviro_folder = "data/enviro/merc/regrid/"
  } else if (var == "AGI") {
    enviro_folder = "data/enviro/agi/mako"
  }

  if (i == 1){
    mako_enviro <- getvarGLORYS(loc_dat = mako, enviro_folder = enviro_folder, depth = depth, col_name = var, depth_var_num = i)
  } else {
    mako_enviro_temp <- getvarGLORYS(loc_dat = mako_enviro, enviro_folder = enviro_folder, depth = depth, col_name = var, depth_var_num = i)
    mako_enviro[paste0(var, "_", depth, "m")] <- mako_enviro_temp[,1]
    mako_enviro <- mako_enviro %>% relocate(geometry, .before = 6)
  }
}

writeRDS(mako_enviro, here("data/loc_enviro/mako_enviro.rds"))

#swordfish
  #day
for (i in 1:nrow(depth_var_vect_all)){

  print(depth_var_vect_all[i,])

  depth = depth_var_vect_all$depth[i]
  var = depth_var_vect_all$var[i]

  if(var != "AGI"){ 
    enviro_folder = "data/enviro/merc/regrid/"
  } else if (var == "AGI") {
    enviro_folder = "data/enviro/agi/swo_day"
  }

  if (i == 1){
    swo_day_enviro <- getvarGLORYS(loc_dat = swo_day, enviro_folder = enviro_folder, depth = depth, col_name = var, depth_var_num = i)
  } else {
    swo_day_enviro_temp <- getvarGLORYS(loc_dat = swo_day_enviro, enviro_folder = enviro_folder, depth = depth, col_name = var, depth_var_num = i)
    swo_day_enviro[paste0(var, "_", depth, "m")] <- swo_day_enviro_temp[,1]
    swo_day_enviro <- swo_day_enviro %>% relocate(geometry, .before = 6)
  }
}

writeRDS(swo_day_enviro, here("data/loc_enviro/swo_day_enviro.rds"))

  #night
for (i in 1:nrow(depth_var_vect_all)){

  print(depth_var_vect_all[i,])

  depth = depth_var_vect_all$depth[i]
  var = depth_var_vect_all$var[i]

  if(var != "AGI"){ 
    enviro_folder = "data/enviro/merc/regrid/"
  } else if (var == "AGI") {
    enviro_folder = "data/enviro/agi/swo_night"
  }

  if (i == 1){
    swo_night_enviro <- getvarGLORYS(loc_dat = swo_night, enviro_folder = enviro_folder, depth = depth, col_name = var, depth_var_num = i)
  } else {
    swo_night_enviro_temp <- getvarGLORYS(loc_dat = swo_night_enviro, enviro_folder = enviro_folder, depth = depth, col_name = var, depth_var_num = i)
    swo_night_enviro[paste0(var, "_", depth, "m")] <- swo_night_enviro_temp[,1]
    swo_night_enviro <- swo_night_enviro %>% relocate(geometry, .before = 6)
  }
}

writeRDS(swo_night_enviro, here("data/loc_enviro/swo_night_enviro.rds"))
