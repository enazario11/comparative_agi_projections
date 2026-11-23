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


getBathy <- function(nc,inpts,bathy_inpt,desired.resolution){
  
  if ((desired.resolution*10) %% 2 == 0) stop("Desired Resolution must be an odd number (e.g. 0.3 not 0.2)")
  nc.data <- nc_open(nc, write=FALSE)
  
  lat <- ncvar_get(nc.data,'latitude') #change to lat and lon here if I don't need to change the crs and res of the nc file
  lon <- ncvar_get(nc.data,'longitude')
  nrows <- length(lat); ncols <- length(lon)
  desired.resolution = 0.25
  desired.resolution = desired.resolution/2
  
  for (i in 1:nrow(inpts)){
    print(i)
    #print(i)
    # #Compute bathymetry: point extraction
    # c <- which.min(abs(lon-data$lon[i]))
    # r <- which.min(abs(lat-data$lat[i]))
    # data.var  <-  ncvar_get(nc.data,'altitude',start=c(c,r),count=c(1,1),verbose=FALSE)
    
    #Compute Bathymetry & rugosity
    c <- which.min(abs(lon-inpts$lon[i]))
    c_low <- which.min(abs(lon-(inpts$lon[i]-desired.resolution)))
    c_up <- which.min(abs(lon-(inpts$lon[i]+desired.resolution)))
    r <- which.min(abs(lat-inpts$lat[i]))
    r_low <- which.min(abs(lat-(inpts$lat[i]-desired.resolution)))
    r_up <- which.min(abs(lat-(inpts$lat[i]+desired.resolution)))
    numcols=abs(c_up-c_low)+1; numrows=abs(r_up-r_low)+1
    #change var below to 'elevation' if I don't need to change the res and CRS of the nc file
    data.var  <-  ncvar_get(nc.data, varid = bathy_inpt, start=c(c_low,r_low),count=c(numcols,numrows),verbose=FALSE)
    data.var <- data.var[data.var < 0]
    inpts$bathy[i] <- mean(data.var, na.rm=T)
    inpts$bathy_sd[i] <- sd(data.var, na.rm=T)
  }
  nc_close(nc.data)
  return(inpts)
}