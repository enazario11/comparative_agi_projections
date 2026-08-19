#libraries
library(tidyverse)
library(here)
library(argosfilter)
library(rnaturalearth)
library(tidyquant)

# TO DO
#1) SDA filter from freitas et al., 2008
# 2) segment loc gaps > 10 days into separate tracks. If tracks < 20 days excluded. 
# 3) apply SSM (next script)

#custom theme
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

#land file
land <- ne_countries(scale = "large", returnclass = "sf")


#argos loc filter 
#mako sharks
mako <- read.csv(here("data/loc_data/mako_tag/PacificArgos.csv")) 
mako <- mako %>% 
  mutate(sp = "Mako sharks", 
         date = as.Date(date, format = "%Y-%m-%d"))
colnames(mako) <- c("tag", "date", "lc", "lon", "lat", "sp")

mako_filter <- sdafilter(mako$lat, mako$lon, mako$date, mako$lc, vmax = 5) # 5 m/s to capture slightly higher swim speeds of study sp, other values set to default
mako_clean <- filter(mako, mako_filter == "not")

  #check tracks after filtering
ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    coord_sf(xlim = c(min(mako_clean$lon) - 2, max(mako_clean$lon) + 2),
      ylim = c(min(mako_clean$lat) - 2, max(mako_clean$lat) + 2),
      expand = FALSE) +
    geom_point(data = mako_clean, aes(lon, lat, color = tag)) + 
    theme_custom() + 
    theme(legend.position = "none")
