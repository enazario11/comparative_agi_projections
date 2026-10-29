
# libraries
library(tidyverse)
library(terra)
library(ncdf4)
library(here)
library(stars)
library(stringr)

stars_to_terra <- function(nc_files, nc_var_name, var_abbrv, file_2D = FALSE){
  
  stars_to_terra_all <- NULL

  template_rast <- rast(
    crs = 'EPSG:4326',
    extent = ext(-170, -100, 0, 55), 
    resolution = 0.25 
  )
    
  for(i in 1:length(nc_files)){
    # As stars object
    if(file_2D == FALSE ){
      merc_stars <- read_stars(nc_files[i]) %>% 
      st_set_dimensions("time_counter", as.POSIXct(st_get_dimension_values(., "time_counter"))) %>% 
      setNames(nc_var_name)
    } else if (file_2D == TRUE){
      merc_stars <- read_ncdf(nc_files[i], curvilinear = c("nav_lon", "nav_lat"), var = "somxlavt") %>% 
      st_set_dimensions("time_counter", as.POSIXct(st_get_dimension_values(., "time_counter")))
    }
    
    #template
    dest <- st_bbox(merc_stars) %>% 
      st_as_stars(dx = 0.25, dy = 0.25, crs = "EPSG:4326")
    
    #seq along each time step and warp each curv time layer to reg (aka template aka dest)
    #stupid st_warp doesn't work with time...maybe later check github comment
    if(file_2D == FALSE){
      merc_list <- lapply(seq_along(st_get_dimension_values(merc_stars, "time_counter")),
                        \(i) st_warp(merc_stars[1, , , 1, i], dest)) %>% suppressMessages()
    } else if(file_2D == TRUE){
      merc_list <- lapply(seq_along(st_get_dimension_values(merc_stars, "time_counter")),
                        \(i) st_warp(merc_stars[1, , , i], dest)) %>% suppressMessages()
    }
    
    merc_reg <- do.call(c, append(merc_list, list(along = "time_counter")))[drop = TRUE]
    
    #save as a .nc file 
    dat_year <- substr(as.character(st_get_dimension_values(merc_stars, which = "time_counter")[1]), 1, 4)
    dat_mo <- substr(as.character(st_get_dimension_values(merc_stars, which = "time_counter")[1]), 6, 7)

    out_name <- paste0("regrid_", strsplit(nc_files[i], split = "/")[[1]][12])
    write_mdim(merc_reg, here(paste0("data/enviro/merc/temporary/", var_abbrv,"/", out_name)))
    
    #read with terra 
    temp_terra <- rast(here(paste0("data/enviro/merc/temporary/", var_abbrv,"/", out_name)))
    
    #assign time 
    d1 <- seq(as.Date(paste0(dat_year,"-","01-01")), by = "day",length = 181)
    d1_L <- seq(as.Date(paste0(dat_year,"-","01-01")), by = "day", length = 182) #leap year number of days
    
    d2 <- seq(as.Date(paste0(dat_year,"-","07-01")), by = "day",length = 184)
    
    # add date info to raster files
    if(str_detect(dat_mo, "01") == TRUE & str_detect(dat_year, "2004|2008|2012|2016|2020|2024") == TRUE){
      time(temp_terra) <- d1_L
    } else if(str_detect(dat_mo, "01") == TRUE & str_detect(dat_year, "2004|2008|2012|2016|2020|2024") == FALSE) {
      time(temp_terra) <- d1
    } else if(str_detect(dat_mo, "07") == TRUE) {
      time(temp_terra) <- d2
    }
    
    #resample and save
    terra_resamp <- resample(temp_terra, template_rast)
    writeCDF(terra_resamp, here(paste0("data/enviro/merc/regrid/", var_abbrv,"/", out_name)))
  }
  
}
