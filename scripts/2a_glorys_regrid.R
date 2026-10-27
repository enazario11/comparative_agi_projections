# libraries
library(tidyverse)
library(terra)
library(here)
library(stars)
library(parallel)
source(here("functions/stars_to_terra.R"))

#regrid and save
#set up parallel computing
cl <- makeCluster(10)

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
parLapply(cl, temp_files, stars_to_terra, nc_var_name = "thetao")
stopCluster(cl)
