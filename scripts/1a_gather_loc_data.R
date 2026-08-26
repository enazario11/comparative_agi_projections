#libraries
library(tidyverse)
library(here)
library(rnaturalearth)
library(sf)
library(tidyquant)

#read in data
land <- ne_countries(scale = "large", returnclass = "sf")
land_pac <- st_transform(land, crs = 3832)

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

alb <- readRDS(here("data/loc_data/alb_tag/validTagsLocns_NOAAonly_bathym_corrected.rds")) %>% mutate(tag = as.character(tag), sp = "Albacore tuna")
blu <- read.csv(here("data/loc_data/blu_tag/blue shark 1 per day.csv"))
mako <- read.csv(here("data/loc_data/mako_tag/mako_spot_filtered_1_step_per_day.csv")) 
swo <- read.csv(here("data/loc_data/swo_tag/PIER_swordfish_tag_data_exp_EM.csv"))

#clean dfs
#albacore
alb <- alb %>% 
  mutate(lon = ifelse(lon360 > 180, lon360 - 360, lon360))

locs_alb <- alb %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |> 
  st_transform(crs = 3832)


#blue sharks
blu <- blu %>% 
  filter(PresAbs == 1) %>% 
  select(-c(X, iter, PresAbs, flag, id)) %>%
  mutate(sp = "Blue sharks")

colnames(blu) <- c("lon", "lat", "date", "tag", "sp")

blu <- blu %>% 
  mutate(date = as.Date(date, format = "%m/%d/%Y"), 
         tag = as.character(tag))

#mako 
mako <- mako %>% 
  mutate(sp = "Mako sharks", 
         date = as.Date(date, format = "%Y-%m-%d")) %>%
  select(-c(lc))
colnames(mako) <- c("tag", "date", "lon", "lat", "sp")

#swordfish 
swo <- swo %>%
  mutate(sp = "Swordfish", 
         date = as.Date(Date..dd.mm.yyyy., format = "%d/%m/%Y"), 
         Ptt = as.character(Ptt))  %>%
         select(c(Ptt, date, Loc.Class, Longitude, Latitude, sp))
colnames(swo) <- c("tag", "date", "lc", "lon", "lat", "sp")

#summary stats of sp tag data
#albacore
table(alb$tag)
length(unique(alb$tag)) #21
min(alb$date) # 2003-08-24
max(alb$date) # 2013-08-23

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
