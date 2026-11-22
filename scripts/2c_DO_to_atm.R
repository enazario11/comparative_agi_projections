#libraries
library(tidyverse)
library(here)
library(terra)
source("functions/oxy_demand_functions.R")


#locate DO mmol/L, temp, and salinity files 
do_mmol_m3_files <- grep(list.files(here("data/enviro/merc/regrid/o2"), full.names = TRUE), pattern = "depth_140m", invert = TRUE, value = TRUE) 

do_mmol_m3_files <- do_mmol_m3_files[str_detect(do_mmol_m3_files, "2002")]

# do_mmol_m3_files <- do_mmol_m3_files[
#   str_detect(do_mmol_m3_files, "2022") | str_detect(do_mmol_m3_files, "2023") | str_detect(do_mmol_m3_files, "2024") | str_detect(do_mmol_m3_files, "2025")
# ]

temp_files <- grep(list.files(here("data/enviro/merc/regrid/T"), full.names = TRUE), pattern = "depth_140m", invert = TRUE, value = TRUE)

temp_files <- temp_files[str_detect(temp_files, "2002")] 

# temp_files <- temp_files[
#   str_detect(temp_files, "2022") | str_detect(temp_files, "2023") | str_detect(temp_files, "2024") | str_detect(temp_files, "2025")
# ]

sal_files <- grep(list.files(here("data/enviro/merc/regrid/S"), full.names = TRUE), pattern = "depth_140m", invert = TRUE, value = TRUE)

sal_files <- sal_files[str_detect(sal_files, "2002")]

# sal_files <- sal_files[
#   str_detect(sal_files, "2022") | str_detect(sal_files, "2023") | str_detect(sal_files, "2024") | str_detect(sal_files, "2025")
# ]

#loop through convert to atm 
for (i in 1:length(do_mmol_m3_files)) {
  #load all files needed for atm conversion
  do_file = rast(do_mmol_m3_files[i])
  temp_file = rast(temp_files[i])
  sal_file = rast(sal_files[i])

  #get filename to save new atm file
  curr_do_filename <- unlist(strsplit(do_mmol_m3_files[i], "/"))[12]
  new_do_filename <- paste0("regrid_o2_atm", unlist(strsplit(curr_do_filename, "regrid_o2"))[2])

  #confirm matching do, temp, salinity files
  do_range <- unlist(strsplit(curr_do_filename, "regrid_o2_"))[2]
  
  curr_temp_filename <- unlist(strsplit(temp_files[i], "/"))[12]
  temp_range <- unlist(strsplit(curr_temp_filename, "regrid_T_"))[2]

  curr_sal_filename <- unlist(strsplit(sal_files[i], "/"))[12]
  sal_range <- unlist(strsplit(curr_sal_filename, "regrid_S_"))[2]

  if(do_range == temp_range & do_range == sal_range & temp_range == sal_range){
    #convert do_mmol_m3 to atm
    do_atm <- do_to_atm(do_file, temp_file, sal_file)

    #save
    writeCDF(do_atm, paste0("data/enviro/merc/regrid/o2_atm/", new_do_filename), overwrite = TRUE)
} else {paste0("ERROR DO: ", do_range, " T: ", temp_range, " S: ", sal_range)}

}
