#libraries
library(tidyverse)
library(here)
library(rnaturalearth)
library(sf)
library(tidyquant)

#read in data
land <- ne_countries(scale = "large", returnclass = "sf")

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

alb <- readRDS(here("data/loc_data/alb_tag/validTagsLocns_NOAAonly_bathym_corrected.rds"))
blu <- read.csv(here("data/loc_data/blu_tag/blue shark 1 per day.csv"))
mako <- read.csv(here("data/loc_data/mako_tag/mako_spot_filtered_1_step_per_day.csv")) 
swo <- read.csv(here("data/loc_data/swo_tag/PIER_swordfish_tag_data_exp_EM.csv"))

#clean dfs
#albacore
alb <- alb %>% 
  mutate(lon = ifelse(lon360 > 180, lon360 - 360, lon360), 
         date = as.Date(dateRd, format = "%Y-%m-%d")) 

locs_alb <- alb %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

#blue sharks
blu <- blu %>% 
  mutate(date = as.Date(Date, format = "%m/%d/%y"), 
         id = paste0(ptt, "_", FL))

#mako 
mako <- mako %>% 
  mutate(date = as.Date(PDT.date, format = "%m/%d/%Y"), 
         id = paste0(PTT, "_", size))

#swordfish 
swo <- swo %>%
  mutate(date = as.Date(Date..dd.mm.yyyy., format = "%d/%m/%Y"), 
         Ptt = as.character(Ptt)) 

#summary stats of sp tag data
#albacore
table(alb$tag)
length(unique(alb$tag)) #21
min(alb$date) # 2003-08-24
max(alb$date) # 2013-08-23

alb %>% group_by(tag) %>% mutate(deploy_dur = max(date, na.rm = TRUE) - min(date, na.rm = TRUE)) %>% ungroup() %>% summarise(mean_deploy = mean(deploy_dur), sd_deploy = sd(deploy_dur))

ggplot() + 
  geom_sf(data = land_pac, fill = "grey85", color = "grey30", linewidth = 0.2) +
  geom_sf(data = locs_alb, aes(color = tag), size = 2) +
  coord_sf(xlim = c(-1500000, 12000000), ylim = c(2000000, 7000000)) +
  geom_point(data = alb, aes(lon, lat, color = tag)) + 
  theme_custom()

#blue sharks
table(blu$tag)
length(unique(blu$tag)) #60
min(blu$date) # 2004-11-07
max(blu$date) # 2012-11-21

blu %>% group_by(ptt) %>% mutate(deploy_dur = max(date, na.rm = TRUE) - min(date, na.rm = TRUE)) %>% filter(!is.infinite(deploy_dur)) %>% ungroup() %>% summarise(mean_deploy = mean(deploy_dur, na.rm = TRUE), sd_deploy = sd(deploy_dur, na.rm = TRUE))

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    coord_sf(xlim = c(min(blu$lon) - 2, max(blu$lon) + 2),
      ylim = c(min(blu$lat) - 2, max(blu$lat) + 2),
      expand = FALSE) +
    geom_point(data = blu, aes(lon, lat, color = tag)) + 
    theme_custom()

#makos
table(mako$tag)
length(unique(mako$tag)) # 121
min(mako$date) # 2003-06-25
max(mako$date) # 2016-02-03

mako %>% group_by(PTT, size) %>% mutate(deploy_dur = max(date, na.rm = TRUE) - min(date, na.rm = TRUE)) %>% filter(!is.infinite(deploy_dur)) %>% ungroup() %>% summarise(mean_deploy = mean(deploy_dur, na.rm = TRUE), sd_deploy = sd(deploy_dur, na.rm = TRUE))

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    coord_sf(xlim = c(min(mako$lon) - 2, max(mako$lon) + 2),
      ylim = c(min(mako$lat) - 2, max(mako$lat) + 2),
      expand = FALSE) +
    geom_point(data = mako, aes(lon, lat, color = tag)) + 
    theme_custom() + 
    theme(legend.position = "none")

#swordfish 
table(swo$tag) 
length(unique(swo$tag)) #57
min(swo$date) # 2017-10-25
max(swo$date) # 2026-01-06

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    coord_sf(xlim = c(min(swo$lon) - 2, max(swo$lon) + 2),
      ylim = c(min(swo$lat) - 2, max(swo$lat) + 2),
      expand = FALSE) +
    geom_point(data = swo, aes(lon, lat, color = tag)) + 
    theme_custom() + 
    theme(legend.position = "none")
