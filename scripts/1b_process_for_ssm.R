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
land <- ne_countries(scale = "large", returnclass = "sf") %>% st_transform(land, crs = 4326)
land_valid <- land %>% st_make_valid()

### clean data files, filter points off land, and remove tracks < 20 days ####
#albacore
alb <- readRDS(here("data/loc_data/alb_tag/validTagsLocns_NOAAonly_bathym_corrected.rds")) %>% 
  mutate(tag = as.character(tag), 
         lon = ifelse(lon360 > 180, lon360 - 360, lon360), 
         date = as.POSIXct(dateRd, format = "%Y-%m-%d")) %>%
  select(c(id = tag, date, lon, lat))

  #filter out locs on land
locs_alb <- alb %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

alb_filt <- locs_alb %>%
        mutate(land_intersect = rowSums(st_intersects(locs_alb, land_valid, sparse = FALSE)) > 0,
               omit_keep = if_else(land_intersect, "omit", "keep")) %>%  #keep the locations that do not intersect with land
        filter(omit_keep == "keep") %>%
        dplyr::select(-c("land_intersect", "omit_keep"))

ggplot() + 
  geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
  geom_sf(data = alb_filt, aes(color = id), size = 2) +
  #coord_sf(xlim = c(st_bbox(alb_filt2$geometry)[1], st_bbox(alb_filt2$geometry)[3]), ylim = c(st_bbox(alb_filt2$geometry)[2], st_bbox(alb_filt2$geometry)[4])) +
  theme_custom() +
  theme(legend.position = "none")

alb_filt2 <- alb_filt %>% 
  mutate(lon = st_coordinates(.)[,1], 
         lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

  #remove tracks < 20 days
alb_clean <- alb_filt2 %>%
  group_by(id) %>%
  filter(n() > 20) %>%
  ungroup()

saveRDS(alb_clean, here("data/loc_data/processed/pre_ssm/alb_dat.rds"))

  #albacore summary stats
alb_clean %>% group_by(id) %>% mutate(deploy_dur = max(date, na.rm = TRUE) - min(date, na.rm = TRUE)) %>% ungroup() %>% summarise(mean_deploy = mean(deploy_dur), sd_deploy = sd(deploy_dur))

table(alb_clean$id)
length(unique(alb_clean$id)) #21
min(alb_clean$date) # 2003-08-24
max(alb_clean$date) # 2013-08-23

#blue sharks
blu <- read.csv(here("data/loc_data/blu_tag/blue shark 1 per day.csv"))
blu <- blu %>% 
  mutate(date = as.POSIXct(Date, format = "%m/%d/%y"), 
         id = paste0(ptt, "_", FL)) %>%
  select(c(id, date, lc, lon = long, lat)) %>%
  filter(lc != "P")

locs_blu <- blu %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) 

blu_filt <- locs_blu %>%
        mutate(land_intersect = rowSums(st_intersects(locs_blu, land_valid, sparse = FALSE)) > 0,
               omit_keep = if_else(land_intersect, "omit", "keep")) %>%  #keep the locations that do not intersect with land
        filter(omit_keep == "keep") %>%
        dplyr::select(-c("land_intersect", "omit_keep"))

#check tracks
ggplot() + 
  geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
  geom_sf(data = blu_filt, aes(color = id), size = 2) +
  coord_sf(xlim = c(st_bbox(locs_blu$geometry)[1], st_bbox(locs_blu$geometry)[3]), ylim = c(st_bbox(locs_blu$geometry)[2], st_bbox(locs_blu$geometry)[4])) +
  #geom_point(data = locs_blu, aes(lon, lat, color = tag)) + 
  theme_custom() + 
  theme(legend.position = "none")

blu_filt2 <- blu_filt %>% 
  mutate(lon = st_coordinates(.)[,1], 
         lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

  #remove tracks < 20 days
blu_clean <- blu_filt2 %>%
  group_by(id) %>%
  filter(n() > 20) %>%
  ungroup() %>% 
  drop_na(date)

saveRDS(blu_clean, here("data/loc_data/processed/pre_ssm/blu_dat.rds"))

  #blue shark summary stats
blu_clean %>% group_by(id) %>% mutate(deploy_dur = max(date, na.rm = TRUE) - min(date, na.rm = TRUE)) %>% ungroup() %>% summarise(mean_deploy = mean(deploy_dur)/86400
, sd_deploy = sd(deploy_dur)/86400)

table(blu_clean$id)
length(unique(blu_clean$id)) # 61
min(blu_clean$date) # 2002-06-28
max(blu_clean$date) # 2015-10-29

#mako sharks
mako <- read.csv(here("data/loc_data/mako_tag/mako_spot_filtered_1_step_per_day.csv"))
mako <- mako %>% 
  mutate(date = as.POSIXct(PDT.date, format = "%m/%d/%Y"), 
         id = paste0(PTT, "_", size)) %>% 
         select(c(id, date, lc, lon = long, lat))
  
locs_mako <- mako %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) 

mako_filt <- locs_mako %>%
        mutate(land_intersect = rowSums(st_intersects(locs_mako, land_valid, sparse = FALSE)) > 0,
               omit_keep = if_else(land_intersect, "omit", "keep")) %>%  #keep the locations that do not intersect with land
        filter(omit_keep == "keep") %>%
        dplyr::select(-c("land_intersect", "omit_keep"))

  #check tracks
ggplot() + 
  geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
  geom_sf(data = mako_filt, aes(color = id), size = 2) +
  coord_sf(xlim = c(st_bbox(locs_mako$geometry)[1], st_bbox(locs_mako$geometry)[3]), ylim = c(st_bbox(locs_mako$geometry)[2], st_bbox(locs_mako$geometry)[4])) +
  #geom_point(data = locs_blu, aes(lon, lat, color = tag)) + 
  theme_custom() + 
  theme(legend.position = "none")

mako_filt2 <- mako_filt %>% 
  mutate(lon = st_coordinates(.)[,1], 
         lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

  #remove tracks < 20 days
mako_clean <- mako_filt2 %>%
  group_by(id) %>%
  filter(n() > 20) %>%
  ungroup() %>% 
  drop_na(date)

saveRDS(mako_clean, here("data/loc_data/processed/pre_ssm/mako_dat.rds"))

  #mako summary stats
mako_clean %>% group_by(id) %>% mutate(deploy_dur = max(date, na.rm = TRUE) - min(date, na.rm = TRUE)) %>% ungroup() %>% summarise(mean_deploy = mean(deploy_dur), sd_deploy = sd(deploy_dur))

table(mako_clean$id)
length(unique(mako_clean$id)) # 85
min(mako_clean$date) # 2003-06-26
max(mako_clean$date) # 2016-02-21

#swordfish
swo <- read.csv(here("data/loc_data/swo_tag/PIER_swordfish_tag_data_exp_EM.csv")) 
swo <- swo %>%
  mutate(sp = "Swordfish", 
         date = as.POSIXct(Date..dd.mm.yyyy., format = "%d/%m/%Y"), 
         Ptt = as.character(Ptt))  %>%
         select(c(id = Ptt, date, lc = Loc.Class, lon = Longitude, lat = Latitude)) 

locs_swo <- swo %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) 

swo_filt <- locs_swo %>%
        mutate(land_intersect = rowSums(st_intersects(locs_swo, land_valid, sparse = FALSE)) > 0,
               omit_keep = if_else(land_intersect, "omit", "keep")) %>%  #keep the locations that do not intersect with land
        filter(omit_keep == "keep") %>%
        dplyr::select(-c("land_intersect", "omit_keep"))

  #check tracks
ggplot() + 
  geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
  geom_sf(data = locs_swo, aes(color = id), size = 2) +
  coord_sf(xlim = c(st_bbox(locs_swo$geometry)[1], st_bbox(locs_swo$geometry)[3]), ylim = c(st_bbox(locs_swo$geometry)[2], st_bbox(locs_swo$geometry)[4])) +
  #geom_point(data = locs_blu, aes(lon, lat, color = tag)) + 
  theme_custom() + 
  theme(legend.position = "none")

swo_filt2 <- swo_filt %>% 
  mutate(lon = st_coordinates(.)[,1], 
         lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

  #remove tracks < 20 days
swo_clean <- swo_filt2 %>%
  group_by(id) %>%
  filter(n() > 20) %>%
  ungroup() %>%
  drop_na(date)

saveRDS(swo_clean, here("data/loc_data/processed/pre_ssm/swo_dat.rds"))

  #swordfish summary stats
swo_clean %>% group_by(id) %>% mutate(deploy_dur = max(date, na.rm = TRUE) - min(date, na.rm = TRUE)) %>% ungroup() %>% summarise(mean_deploy = mean(deploy_dur), sd_deploy = sd(deploy_dur))

table(swo_clean$id)
length(unique(swo_clean$id)) # 28
min(swo_clean$date) # 2017-11-22
max(swo_clean$date) # 2025-12-19
