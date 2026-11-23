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
parLapply(cl, temp_files, stars_to_terra, nc_var_name = "vozocrtx", var_abbrv = "su", file_2D = TRUE) 

parLapply(cl, temp_files, stars_to_terra, nc_var_name = "sozotaux", var_abbrv = "sustr", file_2D = TRUE)
stopCluster(cl)

  #V, northward surface velocity (to later calculate EKE) and wind stress
temp_files <- list.files(here("data/enviro/merc/raw/V"), full.names = TRUE)
parLapply(cl, temp_files, stars_to_terra, nc_var_name = "vomecrty", var_abbrv = "sv", file_2D = TRUE)

parLapply(cl, temp_files, stars_to_terra, nc_var_name = "sometauy", var_abbrv = "svstr", file_2D = TRUE)
stopCluster(cl)

#rename files in a regrid folder to get rid of duplicate var name
curr_var <- "su"
folder_path <- paste0("data/enviro/merc/regrid/", curr_var)

old_files <- list.files(path = folder_path, pattern = paste0(curr_var, "_", "U", "_"), full.names = TRUE)
new_files <- file.path(
  dirname(old_files),
  sub(paste0(curr_var, "_", "U", "_"), paste0(curr_var, "_"), basename(old_files), fixed = TRUE)
)

data.frame(from = basename(old_files), to = basename(new_files))
file.rename(from = old_files, to = new_files)
