#libraries
library(tidyverse)
library(here)
library(argosfilter)
library(rnaturalearth)
library(tidyquant)
library(sf)


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
land <- ne_countries(scale = "large", returnclass = "sf") %>% st_make_valid() %>% st_shift_longitude()
land <- st_transform(land, crs = 4326)

#clean data files, switch to 0-360 longitude (MOM6 format), and remove tracks < 20 days. 

#albacore
alb <- readRDS(here("data/loc_data/alb_tag/validTagsLocns_NOAAonly_bathym_corrected.rds")) %>% 
  mutate(tag = as.character(tag), sp = "Albacore tuna") %>%
  select(c("tag", "dateRd", "lon360", "lat", "sp"))
colnames(alb) <- c("tag", "date", "lon", "lat", "sp")

  #filter out locs on land
alb_filt <- alb %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

alb_intersect <- st_intersects(alb_filt, land, sparse = FALSE)
alb_filt2 <- st_difference(alb_filt, st_union(land)) %>% st_shift_longitude()

ggplot() + 
  geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
  geom_sf(data = alb_filt2, aes(color = tag), size = 2) +
  coord_sf(xlim = c(st_bbox(alb_filt2$geometry)[1], st_bbox(alb_filt2$geometry)[3]), ylim = c(st_bbox(alb_filt2$geometry)[2], st_bbox(alb_filt2$geometry)[4])) +
  theme_custom() +
  theme(legend.position = "none")

alb2 <- alb_filt2 %>% 
  mutate(lon = st_coordinates(.)[,1], 
         lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

  #remove tracks < 20 days
alb_clean <- alb2 %>%
  group_by(tag) %>%
  filter(n() > 20) %>%
  ungroup()

saveRDS(alb_clean, here("data/loc_data/processed/pre_ssm/alb_dat.rds"))

  #albacore summary stats
table(alb_clean$tag)
length(unique(alb_clean$tag)) #21
min(alb_clean$date) # 2003-08-24
max(alb_clean$date) # 2013-08-23

#blue sharks
blu <- read.csv(here("data/loc_data/blu_tag/blue shark 1 per day.csv"))
blu <- blu %>% 
  mutate(sp = "Blue sharks", 
         date2 = paste0(year, "-", mo, "-", day), 
         date2 = as.Date(date2, format = "%Y-%m-%d")) %>%
  select(c(Shark.ID, date2, lc, long, lat, sp)) %>%
  filter(lc != "P")
colnames(blu) <- c("tag", "date", "lc", "lon", "lat", "sp")

locs_blu <- blu %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_shift_longitude()

#check tracks
ggplot() + 
  geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
  geom_sf(data = locs_blu, aes(color = tag), size = 2) +
  coord_sf(xlim = c(st_bbox(locs_blu$geometry)[1], st_bbox(locs_blu$geometry)[3]), ylim = c(st_bbox(locs_blu$geometry)[2], st_bbox(locs_blu$geometry)[4])) +
  #geom_point(data = locs_blu, aes(lon, lat, color = tag)) + 
  theme_custom() + 
  theme(legend.position = "none")

blu2 <- locs_blu %>% 
  mutate(lon = st_coordinates(.)[,1], 
         lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

  #remove tracks < 20 days
blu_clean <- blu2 %>%
  group_by(tag) %>%
  filter(n() > 20) %>%
  ungroup()

saveRDS(blu_clean, here("data/loc_data/processed/pre_ssm/blu_dat.rds"))

  #blue shark summary stats
table(blu_clean$tag)
length(unique(blu_clean$tag)) #62
min(blu_clean$date) # 2002-06-28
max(blu_clean$date) # 2015-10-29

#mako sharks
mako <- read.csv(here("data/loc_data/mako_tag/mako_spot_filtered_1_step_per_day.csv"))
mako <- mako %>% 
  mutate(sp = "Mako sharks", 
         date = as.Date(date, format = "%m/%d/%Y"), 
         PTT = as.character(PTT)) %>% 
         select(c(PTT, date, lc, long, lat, sp)) %>%
  filter(lc != "D")
colnames(mako) <- c("tag", "date", "lc", "lon", "lat", "sp")
  
locs_mako <- mako %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_shift_longitude()

  #check tracks
ggplot() + 
  geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
  geom_sf(data = locs_mako, aes(color = tag), size = 2) +
  coord_sf(xlim = c(st_bbox(locs_mako$geometry)[1], st_bbox(locs_mako$geometry)[3]), ylim = c(st_bbox(locs_mako$geometry)[2], st_bbox(locs_mako$geometry)[4])) +
  #geom_point(data = locs_blu, aes(lon, lat, color = tag)) + 
  theme_custom() + 
  theme(legend.position = "none")

mako2 <- locs_mako %>% 
  mutate(lon = st_coordinates(.)[,1], 
         lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

  #remove tracks < 20 days
mako_clean <- mako2 %>%
  group_by(tag) %>%
  filter(n() > 20) %>%
  ungroup()

saveRDS(mako_clean, here("data/loc_data/processed/pre_ssm/mako_dat.rds"))

  #mako summary stats
table(mako_clean$tag)
length(unique(mako_clean$tag)) # 79
min(mako_clean$date) # 2003-06-26
max(mako_clean$date) # 2016-02-21

#swordfish
swo <- read.csv(here("data/loc_data/swo_tag/PIER_swordfish_tag_data_exp_EM.csv")) 
swo <- swo %>%
  mutate(sp = "Swordfish", 
         date = as.Date(Date..dd.mm.yyyy., format = "%d/%m/%Y"), 
         Ptt = as.character(Ptt))  %>%
         select(c(Ptt, date, Loc.Class, Longitude, Latitude, sp)) %>%
         filter(Loc.Class != "Deploy" & Loc.Class != "Recapture" & Loc.Class != "Mote" & Loc.Class != "" & Loc.Class != "Z")
colnames(swo) <- c("tag", "date", "lc", "lon", "lat", "sp")

locs_swo <- swo %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_shift_longitude()

  #check tracks
ggplot() + 
  geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
  geom_sf(data = locs_swo, aes(color = tag), size = 2) +
  coord_sf(xlim = c(st_bbox(locs_swo$geometry)[1], st_bbox(locs_swo$geometry)[3]), ylim = c(st_bbox(locs_swo$geometry)[2], st_bbox(locs_swo$geometry)[4])) +
  #geom_point(data = locs_blu, aes(lon, lat, color = tag)) + 
  theme_custom() + 
  theme(legend.position = "none")

swo2 <- locs_swo %>% 
  mutate(lon = st_coordinates(.)[,1], 
         lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

  #remove tracks < 20 days
swo_clean <- swo2 %>%
  group_by(tag) %>%
  filter(n() > 20) %>%
  ungroup()

saveRDS(swo_clean, here("data/loc_data/processed/pre_ssm/swo_dat.rds"))

  #swordfish summary stats
table(swo_clean$tag)
length(unique(swo_clean$tag)) # 28
min(swo_clean$date) # 2017-11-22
max(swo_clean$date) # 2025-12-19
