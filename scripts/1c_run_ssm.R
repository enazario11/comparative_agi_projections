#libraries
library(tidyverse)
library(here)
library(rnaturalearth)
library(tidyquant)

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
land <- ne_countries(scale = "large", returnclass = "sf")
land_pac <- st_transform(land, crs = 3832)

### fit ssm ####
#skip albacore bc cleaning already happened -- ssm not needed

#blue sharks


#mako sharks


#swordfish 
