#libraries
library(tidyverse)
library(here)
library(terra)
library(sf)

getvarGLORYS <- function(loc_dat, enviro_folder, depth = "_0m", col_name){
  
  enviro_files <- list.files(here(enviro_folder), full.names = TRUE, pattern = "_0m")
  full_rast <- rast(enviro_files)
  loc_dat[[col_name]] <- NA

  for(i in 1:nrow(loc_dat)){
    print(paste(i, "/", nrow(loc_dat)))
    date_lyr <- subset(full_rast, time(full_rast) == as.Date(loc_dat$date[i]))
    extract_value <- extract(date_lyr, vect(loc_dat$geometry[i]))[1,2]

    loc_dat[i, col_name] <- extract_value
  }

  return(loc_dat)
}
