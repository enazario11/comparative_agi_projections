#libraries
library(tidyverse)
library(here)
library(terra)
library(sf)

getvarGLORYS <- function(loc_dat, enviro_folder, depth, col_name){
  
  depth = paste0("_", depth, "m")

  if(col_name == "AGI"){
    enviro_files <- list.files(here(enviro_folder), full.names = TRUE, pattern = depth)
  } else {
    enviro_files <- list.files(here(paste0(enviro_folder, col_name)), full.names = TRUE, pattern = depth)
  }
  
  full_rast <- rast(enviro_files)
  loc_dat[[paste0(col_name, depth)]] <- NA

  for(i in 1:nrow(loc_dat)){
    #print(paste(i, "/", nrow(loc_dat)))
    date_lyr <- subset(full_rast, time(full_rast) == as.Date(loc_dat$date[i]))
    extract_value <- extract(date_lyr, vect(loc_dat$geometry[i]))[1,2]

    loc_dat[i, paste0(col_name, depth)] <- extract_value
  }

  loc_dat <- loc_dat %>% st_drop_geometry()
  return(loc_dat[,paste0(col_name, depth)])

}
