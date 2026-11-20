#libraries
library(tidyverse)
library(here)
library(terra)

path_copernicusmarine <- "C:/Users/nazar/Documents/Projects/comparative_agi_projections/.venv/Scripts/copernicusmarine.exe"

ymin <- 0
ymax <- 55
xmin <- -170
xmax <- -100

depth_vect = c(0.51, 9.82, 19.43, 30.87, 41.18, 47.21, 61.11, 69.02, 77.61, 86.93, 97.04, 108.03, 120, 133.08, 147.41, 199.79, 300.89, 411.79, 508.64, 628.03, 697.26, 773.37, 856.68, 1045.85)
depth_vect_name <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 150, #remove 140 b/c would duplicate the 133.08 layer
             200, 300, 400, 500, 600, 700, 800, 900, 1000)
depth_tbl <- data.frame(
  level   = depth_vect,
  nominal = depth_vect_name
)
depth_vect_2d = 1
date_vect = seq(as.Date("2023-01-01"), as.Date("2025-12-31"), by = "6 months")

depth_date_vect = merge(depth_tbl, data.frame(date = date_vect), by = NULL)
depth_date_vect_2d = expand.grid(depth = depth_vect_2d, date = date_vect)

cmem_nc <- function(cmem_id, out_directory, out_file, cmem_var, start_date, end_date, lon_min = xmin, lat_min = ymin, lon_max = xmax, lat_max = ymax, min_depth, max_depth){

    command <- paste(
      shQuote(path_copernicusmarine),
      "subset",
      paste("--dataset-id", cmem_id),
      paste("--variable", cmem_var), 
      paste("--start-datetime", start_date),
      paste("--end-datetime", end_date),
      paste("--minimum-longitude", lon_min),
      paste("--maximum-longitude", lon_max),
      paste("--minimum-latitude", lat_min),
      paste("--maximum-latitude", lat_max),
      paste("--minimum-depth", min_depth),
      paste("--maximum-depth", max_depth),
      "-o", shQuote(out_directory), #save directory
      paste("--output-filename", out_file), #save filename within directory
      sep = " "
    )

  system(command)
}

#chlorophyll
for(i in 1:nrow(depth_date_vect_2d)){
  curr_dep_date <- depth_date_vect_2d[i, 1:2]
  min_depth = 0 #only getting surface values
  max_depth = curr_dep_date$depth
  start = curr_dep_date$date
  end = start %m+% months(6) - days(1)

  cmem_nc(cmem_id = "cmems_mod_glo_bgc_my_0.25deg_P1D-m",
          out_directory = here("data/enviro/cmem/chl_temp"), 
          out_file = paste0("chl_", start, "_", end, "_", "depth_", depth,"m", "_0.25_D.nc"),
          cmem_var = "chl",
          start_date = paste0(start,"T00:00:00"),
          end_date = paste0(end,"T00:00:00"), 
          min_depth = min_depth, 
          max_depth = max_depth) 
}

#load and clean rasters, regrid to 0.25 res
chl_files <- list.files(here("data/enviro/cmem/chl_temp/"), full.names = TRUE)
for(i in 1:length(chl_files)){
  curr_rast <- rast(chl_files[i]) 
  curr_filename <- unlist(strsplit(chl_files[i], "/"))[11]

  template_rast <- rast(
      crs = crs(curr_rast),
      extent = ext(-170, -100, 0, 55), 
      resolution = 0.25 
    )

  rast_clean <- terra::resample(curr_rast, template_rast)
  time(rast_clean) <- time(curr_rast)
  varnames(rast_clean) <- varnames(curr_rast)

  writeCDF(rast_clean, here(paste0("data/enviro/merc/regrid/chl/regrid_",curr_filename)), overwrite = TRUE)
}


#oxygen
for(i in 1:nrow(depth_date_vect)){
  curr_dep_date <- depth_date_vect[i, 1:3]
  min_depth = curr_dep_date$level - 0.2
  max_depth = curr_dep_date$level + 0.2
  depth_name = curr_dep_date$nominal
  start = curr_dep_date$date
  end = start %m+% months(6) - days(1)

  cmem_nc(cmem_id = "cmems_mod_glo_bgc_my_0.25deg_P1D-m",
          out_directory = here("data/enviro/cmem/o2_temp"), 
          out_file = paste0("o2_", start, "_", end, "_", "depth_", depth_name,"m", "_0.25_D.nc"),
          cmem_var = "o2",
          start_date = paste0(start,"T00:00:00"),
          end_date = paste0(end,"T00:00:00"), 
          min_depth = min_depth, 
          max_depth = max_depth) 
}

#load and clean rasters, regrid to 0.25 res
o2_files <- list.files(here("data/enviro/cmem/o2_temp/"), full.names = TRUE)
for(i in 1:length(o2_files)){
  curr_rast <- rast(o2_files[i]) 
  curr_filename <- unlist(strsplit(o2_files[i], "/"))[11]

  template_rast <- rast(
      crs = crs(curr_rast),
      extent = ext(-170, -100, 0, 55), 
      resolution = 0.25 
    )

  rast_clean <- terra::resample(curr_rast, template_rast)
  time(rast_clean) <- time(curr_rast)
  varnames(rast_clean) <- varnames(curr_rast)

  writeCDF(rast_clean, here(paste0("data/enviro/merc/regrid/o2/regrid_",curr_filename)), overwrite = TRUE)
}
