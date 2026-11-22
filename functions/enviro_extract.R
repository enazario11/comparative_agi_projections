#libraries
library(tidyverse)
library(here)
library(terra)
library(sf)
library(foreach)
library(doParallel)

getvarGLORYS_par <- function(loc_dat, enviro_folder, depth, col_name){
   
  depth <- paste0("_", depth, "m")
  out_name <- paste0(col_name, depth)
  cache_name <- paste0("full_rast_", out_name) 
  loc_dat[[out_name]] <- NA

  if(col_name == "AGI"){
    enviro_files <- list.files(here(enviro_folder), full.names = TRUE, pattern = depth)
  } else {
    enviro_files <- list.files(here(paste0(enviro_folder, col_name)), full.names = TRUE, pattern = depth)
  }
  
  doParallel::registerDoParallel(cl)

  enviro_col <- foreach(i=1:nrow(loc_dat), .combine = c, .export = c("enviro_files", "loc_dat", "cache_name")) %dopar% {
        
        if (i %% 500 == 0) {
          try(cat("Completed", i, cache_name, "\n",
              file = paste0("progress_", Sys.getpid(), ".txt"), append = TRUE),
              silent = TRUE)
        }

        if (!exists(cache_name, envir = .GlobalEnv, inherits = FALSE)) {
            assign(cache_name, terra::rast(enviro_files), envir = .GlobalEnv)
        }

        full_rast <- get(cache_name, envir = .GlobalEnv)

        date_lyr <- subset(full_rast, time(full_rast) == as.Date(loc_dat$date[i]))
        extract_value <- extract(date_lyr, vect(loc_dat$geometry[i]))[1,2]

        loc_dat[i, paste0(col_name, depth)] <- extract_value
  }
   
  loc_dat[[paste0(col_name, depth)]] <- enviro_col
  loc_dat <- loc_dat %>% st_drop_geometry()

  return(loc_dat[,paste0(col_name, depth)])
}
