# libraries
library(tidyverse)
library(terra)
library(here)
library(stars)
library(parallel)
source(here("functions/stars_to_terra.R"))

#regrid and save
#set up parallel computing
cl <- makeCluster(8)

clusterEvalQ(cl, {
  library(tidyverse)
  library(terra)
  library(ncdf4)
  library(here)
  library(stars)
  library(stringr)
})

clusterExport(cl, varlist = c("stars_to_terra"))

  #temperature
temp_files <- list.files(here("data/enviro/merc/raw/T"), full.names = TRUE)
parLapply(cl, temp_files, stars_to_terra, nc_var_name = "thetao", var_abbrv = "T")
stopCluster(cl)

  #salinity
temp_files <- list.files(here("data/enviro/merc/raw/S"), full.names = TRUE)
parLapply(cl, temp_files, stars_to_terra, nc_var_name = "vosaline", var_abbrv = "S")
stopCluster(cl)

  #2D
temp_files <- list.files(here("data/enviro/merc/raw/2D"), full.names = TRUE)
parLapply(cl, temp_files, stars_to_terra, nc_var_name = "somxlavt", var_abbrv = "MLD", file_2D = TRUE)
stopCluster(cl)

  #U, eastward surface velocity (to later calculate EKE) and wind stress
temp_files <- list.files(here("data/enviro/merc/raw/U"), full.names = TRUE)
parLapply(cl, temp_files, stars_to_terra, nc_var_name = "vozocrtx", var_abbrv = "su") 

parLapply(cl, temp_files, stars_to_terra, nc_var_name = "sozotaux", var_abbrv = "sustr", file_2D = TRUE)
stopCluster(cl)

  #V, northward surface velocity (to later calculate EKE) and wind stress
temp_files <- list.files(here("data/enviro/merc/raw/V"), full.names = TRUE)
parLapply(cl, temp_files, stars_to_terra, nc_var_name = "vomecrty", var_abbrv = "sv", file_2D = TRUE)

parLapply(cl, temp_files, stars_to_terra, nc_var_name = "sometauy", var_abbrv = "svstr", file_2D = TRUE)
stopCluster(cl)

  #oxygen
temp_files <- list.files(here("data/enviro/merc/raw/o2"), full.names = TRUE)
parLapply(cl, temp_files, stars_to_terra, nc_var_name = "o2", var_abbrv = "o2", file_2D = TRUE)
stopCluster(cl)

  #chlorophyll
temp_files <- list.files(here("data/enviro/merc/raw/chl"), full.names = TRUE)
parLapply(cl, temp_files, stars_to_terra, nc_var_name = "chl", var_abbrv = "chl", file_2D = TRUE)
stopCluster(cl)



