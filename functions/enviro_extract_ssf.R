#libraries
library(tidyverse)
library(here)
library(terra)
library(sf)
library(foreach)
library(doParallel)
library(ISOweek)

getvarGLORYS_par_ssf <- function(loc_dat, enviro_folder, depth, col_name){
   
  depth <- paste0("_", depth, "m")
  out_name <- paste0(col_name, depth)
  cache_name <- paste0("full_rast_", out_name) 
  weekly_cache <- "weekly_rast"
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
          full_rast <- terra::rast(enviro_files)
          weekly_rast <- tapp(full_rast, "yearweek", mean, na.rm = TRUE)
          rast_iso <- sub("^yw_(\\d{4})(\\d{2})$", "\\1-W\\2-4", names(weekly_rast))
          time(weekly_rast) <- ISOweek2date(rast_iso)

          assign(cache_name, full_rast, envir = .GlobalEnv)
          assign(weekly_cache, weekly_rast, envir = .GlobalEnv)
        }

        weekly_rast <- get(weekly_cache, envir = .GlobalEnv)

        nearest_week <- which.min(abs(time(weekly_rast) - as.Date(loc_dat$date[i])))
    
        week_lyr <- weekly_rast[[unique(nearest_week)]]
        extract_value <- extract(week_lyr, vect(loc_dat$geometry[i]))[1,2]

        loc_dat[i, paste0(col_name, depth)] <- extract_value
  }
   
  loc_dat[[paste0(col_name, depth)]] <- enviro_col
  loc_dat <- loc_dat %>% st_drop_geometry()

  return(loc_dat[,paste0(col_name, depth)])
}


getbathy_par <- function(loc_dat, agg_factor){
   
  cache_name <- "full_rast_bathy"
  mean_cache <-"mean_bathy"
  sd_cache <- "sd_bathy"
  loc_dat[["bathy"]] <- NA
  loc_dat[["rugosity"]] <- NA

  enviro_files <- list.files(here("data/enviro/bathy"), full.names = TRUE, pattern = "sub_0m")
  
  doParallel::registerDoParallel(cl)

  enviro_col <- foreach(i=1:nrow(loc_dat), .combine = rbind, .export = c("enviro_files", "loc_dat", "agg_factor", "cache_name", "mean_cache", "sd_cache")) %dopar% {
        
        if (i %% 500 == 0) {
          try(cat("Completed", i, cache_name, "\n",
              file = paste0("progress_", Sys.getpid(), ".txt"), append = TRUE),
              silent = TRUE)
        }

        if (!exists(cache_name, envir = .GlobalEnv, inherits = FALSE)) {
          full_rast <- terra::rast(enviro_files)
          assign(cache_name, full_rast, envir = .GlobalEnv)
          assign(mean_cache, terra::aggregate(full_rast, fact = agg_factor, fun = mean, na.rm = TRUE), envir = .GlobalEnv)
          assign(sd_cache,   terra::aggregate(full_rast, fact = agg_factor, fun = sd,   na.rm = TRUE), envir = .GlobalEnv)
        }

        bathy_mean <- get(mean_cache, envir = .GlobalEnv)
        bathy_sd   <- get(sd_cache, envir = .GlobalEnv)

        extract_bathy <- extract(bathy_mean, vect(loc_dat$geometry[i]))[1,2] 
        extract_rug <- extract(bathy_sd, vect(loc_dat$geometry[i]))[1,2] 

        data.frame(bathy = extract_bathy, rugosity = extract_rug)
  }
   
  return(enviro_col)
}
