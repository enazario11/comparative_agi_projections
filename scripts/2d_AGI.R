#libraries
library(tidyverse)
library(here)
library(terra)
library(sf)
source("functions/oxy_demand_functions.R")

#set AGI parameters for each species
#albacore
alb_agi <- list(
  Linf = 134, #FishBase overall (CA values doubtful so used overall Loo)
  LwA = 0.02720, #FishBase New Zealand (closest available)
  LwB = 2.8, #FishBase New Zealand (closest available)
  Tpref = 17.6, #Childers et al., 2011, fisheries oceanography
  PO2_thresh = 1.5 #ml/L, threshold of dissolved oxygen concentrations ~1.5 ml/L below which only <2.5% of albacore was caught. Bertrand et al., 2002, Canadian journal of fisheries and aquatic sciences
)

#blue
blu_agi <- list(
  Linf = 295, #FishBase CA
  LwA = 0.00248, #FishBase New Zealand (closest available, averaged male and female values)
  LwB = 3.384, #FishBase New Zealand (closest available, averaged male and female values)
  Tpref = 20, #Musyl et al., 2011, Fishery Bulletin. Daytime mean temp was 19 and night was 21. Selected 20 as middle.
  PO2_thresh = 3.5 #ml/L vedor et al., 2021, eLife
)

#mako
mako_agi <- list(
  Linf = 321, #FishBase CA
  LwA = 0.01670, #FishBase New Zealand (closest available)
  LwB = 2.847, #FishBase New Zealand (closest available)
  Tpref = 18, #Nasby-Lucas et al., 2019
  PO2_thresh = 2 #ml/L vetter et al., 2008
)

#swordfish
swo_day_agi <- list(
  Linf = 252, #FishBase overall (limited Pacific values so used overall Loo)
  LwA = 0.00776, #FishBase New Zealand (values doubtful but comparable to US Atl. values)
  LwB = 3.210, #FishBase New Zealand (values doubtful but comparable to US Atl. values)
  Tpref = 9, # 9C day temperature preferences. Abecassis et al., 2012, Marine Ecology Progress Series; Dewar et al., 2011 Fisheries Oceanography
  PO2_thresh = 1.47 #ml/L median concentration across all maximum depths. Abecassis et al., 2012, Marine Ecology Progress Series
)

swo_night_agi <- list(
  Linf = 252, #FishBase overall (limited Pacific values so used overall Loo)
  LwA = 0.00776, #FishBase New Zealand (values doubtful but comparable to US Atl. values)
  LwB = 3.210, #FishBase New Zealand (values doubtful but comparable to US Atl. values)
  Tpref = 18, # 18C night temperature preferences. Abecassis et al., 2012, Marine Ecology Progress Series; Dewar et al., 2011 Fisheries Oceanography
  PO2_thresh = 1.47 #ml/L median concentration across all maximum depths. Abecassis et al., 2012, Marine Ecology Progress Series
)

#calculate AGI layers per species and depth layer
#locate DO mmol/L, temp, and salinity files 
do_files <- list.files(here("data/enviro/merc/regrid/o2_atm"), full.names = TRUE)

do_files <- do_files[
  str_detect(do_files, "2022") | str_detect(do_files, "2023") | str_detect(do_files, "2024") | str_detect(do_files, "2025")
]

temp_files <- grep(list.files(here("data/enviro/merc/regrid/T"), full.names = TRUE), pattern = "depth_140m", invert = TRUE, value = TRUE)

temp_files <- temp_files[
  str_detect(temp_files, "2022") | str_detect(temp_files, "2023") | str_detect(temp_files, "2024") | str_detect(temp_files, "2025")
]

sal_files <- grep(list.files(here("data/enviro/merc/regrid/S"), full.names = TRUE), pattern = "depth_140m", invert = TRUE, value = TRUE)

sal_files <- sal_files[
  str_detect(sal_files, "2022") | str_detect(sal_files, "2023") | str_detect(sal_files, "2024") | str_detect(sal_files, "2025")
]

agi_rast <- function(sp_name, sp_list, do_file, temp_file, sal_file) {
    #load all files needed for atm conversion
    do_file_rast = rast(do_file)
    temp_file_rast = rast(temp_file)
    sal_file_rast = rast(sal_file)

    #get filename to save new atm file
    curr_do_filename <- unlist(strsplit(do_file, "/"))[12]
    new_do_filename <- paste0("regrid_AGI_", sp_name, unlist(strsplit(curr_do_filename, "regrid_o2"))[2])

    #confirm matching do, temp, salinity files
    do_range <- unlist(strsplit(curr_do_filename, "regrid_o2_atm_"))[2]
    
    curr_temp_filename <- unlist(strsplit(temp_file, "/"))[12]
    temp_range <- unlist(strsplit(curr_temp_filename, "regrid_T_"))[2]

    curr_sal_filename <- unlist(strsplit(sal_file, "/"))[12]
    sal_range <- unlist(strsplit(curr_sal_filename, "regrid_S_"))[2]

    #get median temp and salinity values to convert DO Thresh value
    med_temp_rast <- median(temp_file_rast)
    med_temp <-  global(med_temp_rast, FUN = "median", na.rm = TRUE)[1,1]

    med_sal_rast <- median(sal_file_rast)
    med_sal <- global(med_sal_rast, FUN = "median", na.rm = TRUE)[1,1]

    if(do_range == temp_range & do_range == sal_range & temp_range == sal_range){
      #convert PO2_thresh to atm
      do_thresh_atm <- do_to_atm(sp_list$PO2_thresh, med_temp, med_sal, thresh = TRUE)
      met_demand <- OxyDemand(Tpref = sp_list$Tpref, PO2_thresh = do_thresh_atm, T_C = temp_file_rast, Linf = sp_list$Linf, LwA = sp_list$LwA, LwB = sp_list$LwB)
      agi <- do_file_rast/met_demand

      #save
      writeCDF(agi, paste0("data/enviro/agi/", sp_name, "/", new_do_filename), overwrite = TRUE)
  } else {paste0("ERROR DO: ", do_range, " T: ", temp_range, " S: ", sal_range)}

}

#albacore
for (i in 1:length(do_files)) {
  do_file = do_files[i]
  temp_file = temp_files[i]
  sal_file = sal_files[i]

  agi_rast(sp_name = "alb", sp_list = alb_agi, do_file = do_file, temp_file = temp_file, sal_file = sal_file)
}

#blue 
for (i in 1:length(do_files)) {
  do_file = do_files[i]
  temp_file = temp_files[i]
  sal_file = sal_files[i]

  agi_rast(sp_name = "blu", sp_list = blu_agi, do_file = do_file, temp_file = temp_file, sal_file = sal_file)
}

#mako 
for (i in 1:length(do_files)) {
  do_file = do_files[i]
  temp_file = temp_files[i]
  sal_file = sal_files[i]

  agi_rast(sp_name = "mako", sp_list = mako_agi, do_file = do_file, temp_file = temp_file, sal_file = sal_file)
}

#swordfish (day and night)
for (i in 1:length(do_files)) {
  do_file = do_files[i]
  temp_file = temp_files[i]
  sal_file = sal_files[i]

  agi_rast(sp_name = "swo_day", sp_list = swo_day_agi, do_file = do_file, temp_file = temp_file, sal_file = sal_file)
}

for (i in 1:length(do_files)) {
  do_file = do_files[i]
  temp_file = temp_files[i]
  sal_file = sal_files[i]

  agi_rast(sp_name = "swo_night", sp_list = swo_night_agi, do_file = do_file, temp_file = temp_file, sal_file = sal_file)
}
