#libraries
library(tidyverse)
library(here)
library(terra)
library(sf)
library(doParallel)
library(foreach)
source("functions/oxy_demand_functions.R")

#set AGI parameters for each species
#albacore



#blue



#mako
mako_agi <- list(
  Linf = 321, #FishBase CA
  LwA = 0.01670, #FishBase New Zealand (closest available)
  LwB = 2.847, #FishBase New Zealand (closest available)
  Tpref = 18, #Nasby-Lucas et al., 2019
  PO2_thresh = 2 #ml/L vetter et al., 2008
)

#swordfish



#calculate AGI layers per species and depth layer
#parallel
cl <- makeCluster(8)

clusterEvalQ(cl, {
  library(tidyverse)
  library(here)
  library(terra)
  library(doParallel)
  library(foreach)
})

clusterExport(cl, "OxyDemand")

registerDoParallel(cl)

#locate DO mmol/L, temp, and salinity files 
do_files <- grep(list.files(here("data/enviro/merc/regrid/o2_atm"), full.names = TRUE), pattern = "depth_140m", invert = TRUE, value = TRUE) 
temp_files <- grep(list.files(here("data/enviro/merc/regrid/T"), full.names = TRUE), pattern = "depth_140m", invert = TRUE, value = TRUE)
sal_files <- grep(list.files(here("data/enviro/merc/regrid/S"), full.names = TRUE), pattern = "depth_140m", invert = TRUE, value = TRUE)

foreach (i=1:length(do_files)) %dopar% {
  #load all files needed for atm conversion
  do_file = rast(do_files[i])
  temp_file = rast(temp_files[i])
  sal_file = rast(sal_files[i])

  #get filename to save new atm file
  curr_do_filename <- unlist(strsplit(do_files[i], "/"))[12]
  new_do_filename <- paste0("regrid_AGI", unlist(strsplit(curr_do_filename, "regrid_o2"))[2])

  #confirm matching do, temp, salinity files
  do_range <- unlist(strsplit(curr_do_filename, "regrid_o2_"))[2]
  
  curr_temp_filename <- unlist(strsplit(temp_files[i], "/"))[12]
  temp_range <- unlist(strsplit(curr_temp_filename, "regrid_T_"))[2]

  curr_sal_filename <- unlist(strsplit(sal_files[i], "/"))[12]
  sal_range <- unlist(strsplit(curr_sal_filename, "regrid_S_"))[2]

  if(do_range == temp_range & do_range == sal_range & temp_range == sal_range){
    #convert PO2_thresh to atm
    do_thresh_atm <- do_to_atm(mako_agi$PO2_thresh, mean(temp_file, na.rm = TRUE), mean(sal_file, na.rm = TRUE), thresh = TRUE)
    met_demand <- OxyDemand(Tpref = mako_agi$Tpref, PO2_thresh = do_thresh_atm, T_C = temp_file, Linf = mako_agi$Linf, LwA = mako_ag$LwA, LwB = mako_agi$LwB)
    agi <- do_file/met_demand

    #save
    writeCDF(agi, paste0("data/enviro/agi/mako", new_do_filename))
} else {paste0("ERROR DO: ", do_range, " T: ", temp_range, " S: ", sal_range)}

}

stopCluster(cl)
