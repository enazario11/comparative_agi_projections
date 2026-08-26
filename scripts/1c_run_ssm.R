#libraries
library(tidyverse)
library(here)
library(rnaturalearth)
library(tidyquant)
library(aniMotum)

### to do 
#1) confirm from hists the average time step between locations. Ideally step will be set to 24hr but could be longer depending on resolution.


### figure features ####
# custom theme
theme_custom <- function(){
 theme_tq() + 
 theme(axis.text = element_text(size = 12, color = "black"), 
       axis.title = element_text(size = 14, color = "black"), 
       legend.text = element_text(size = 12, color = "black"), 
       legend.title = element_text(size = 14, color = "black"), 
       plot.title = element_text(size = 14, color = "black"),
       strip.text = element_text(size = 12), 
       legend.position = "right")
}

# land file
land <- ne_countries(scale = "large", returnclass = "sf") %>% st_make_valid()
land <- st_transform(land, crs = 4326)
land_pac <- st_transform(land, crs = 3832)

### fit ssm ####
#skip albacore bc cleaning already happened -- ssm not needed
#blue sharks
blu <- readRDS(here("data/loc_data/processed/pre_ssm/blu_dat.rds")) %>% select(-c(diff)) %>% filter(lc != "P")
colnames(blu) <- c("id", "date", "lc", "lon", "lat", "sp")

blu_ssm <- fit_ssm(blu, 
                   spdf = FALSE, #turn off sda filter
                   date = "date", 
                   coord = c("lon", "lat"), 
                   model = "crw", 
                   time.step = 24) #control = ssm_control(verbose = 0)



#mako sharks


#swordfish 

