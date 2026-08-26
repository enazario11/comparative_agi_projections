#libraries
library(tidyverse)
library(here)
library(argosfilter)
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
land <- ne_countries(scale = "large", returnclass = "sf") %>% st_make_valid()
land <- st_transform(land, crs = 4326)
land_pac <- st_transform(land, crs = 3832)

### argos SDA filter ####
group_sda_filter <- function (sp_df){

df_clean <- data.frame() #empty df for filtered tracks grouped by ptt
  
for(i in 1:length(unique(sp_df$tag))){
  curr_id = unique(sp_df$tag)[i]
  curr_df <- sp_df %>%
    filter(tag == curr_id) %>%
    arrange(date)

 if(nrow(curr_df) > 20){
  filter <- sdafilter(curr_df$lat, curr_df$lon, curr_df$date, curr_df$lc, vmax = 3) #consider setting vmax to 3 and defaults for dist/ang for all sp. following O'Toole et al., 2021
  curr_clean <- filter(curr_df, filter == "not")
 } else {
   curr_clean = data.frame()
 }
  
 df_clean <- rbind(df_clean, curr_clean)

 }

 return(df_clean)
}

#skip for albacore bc archival tags not satellite

#blue sharks
blu <- read.csv(here("data/loc_data/blu_tag/blue shark 1 per day.csv"))
blu <- blu %>% 
  mutate(sp = "Blue sharks", 
         date2 = paste0(year, "-", mo, "-", day), 
         date2 = as.Date(date2, format = "%Y-%m-%d")) %>%
  select(c(Shark.ID, date2, lc, long, lat, sp)) %>%
  filter(lc != "P")
colnames(blu) <- c("tag", "date", "lc", "lon", "lat", "sp")

blu_clean <- group_sda_filter(blu)
locs_blu <- blu_clean %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |> 
  st_transform(crs = 3832)

#check tracks after filtering
ggplot() + 
  geom_sf(data = land_pac, fill = "grey85", color = "grey30", linewidth = 0.2) +
  geom_sf(data = locs_blu, aes(color = tag), size = 2) +
  coord_sf(xlim = c(st_bbox(locs_blu$geometry)[1], st_bbox(locs_blu$geometry)[3]), ylim = c(st_bbox(locs_blu$geometry)[2], st_bbox(locs_blu$geometry)[4])) +
  geom_point(data = alb, aes(lon, lat, color = tag)) + 
  theme_custom() + 
  theme(legend.position = "none")

#mako sharks
mako <- read.csv(here("data/loc_data/mako_tag/mako_spot_filtered_1_step_per_day.csv"))
mako <- mako %>% 
  mutate(sp = "Mako sharks", 
         date = as.Date(date, format = "%m/%d/%Y"), 
         PTT = as.character(PTT)) %>% 
         select(c(PTT, date, lc, long, lat, sp)) %>%
  filter(lc != "D")

colnames(mako) <- c("tag", "date", "lc", "lon", "lat", "sp")

mako_clean <- group_sda_filter(mako)
  
  #check tracks after filtering
ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    coord_sf(xlim = c(min(mako_clean$lon) - 2, max(mako_clean$lon) + 2),
      ylim = c(min(mako_clean$lat) - 2, max(mako_clean$lat) + 2),
      expand = FALSE) +
    geom_point(data = mako_clean, aes(lon, lat, color = tag)) + 
    theme_custom() +
    theme(legend.position = "none")

#swordfish
swo <- read.csv(here("data/loc_data/swo_tag/PIER_swordfish_tag_data_exp_EM.csv")) 
swo <- swo %>%
  mutate(sp = "Swordfish", 
         date = as.Date(Date..dd.mm.yyyy., format = "%d/%m/%Y"), 
         Ptt = as.character(Ptt))  %>%
         select(c(Ptt, date, Loc.Class, Longitude, Latitude, sp)) %>%
         filter(Loc.Class != "Deploy" & Loc.Class != "Recapture" & Loc.Class != "Mote" & Loc.Class != "" & Loc.Class != "Z")
colnames(swo) <- c("tag", "date", "lc", "lon", "lat", "sp")

swo_clean <- group_sda_filter(swo)

  #check tracks after filtering
ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    coord_sf(xlim = c(min(swo_clean$lon) - 2, max(swo_clean$lon) + 2),
      ylim = c(min(swo_clean$lat) - 2, max(swo_clean$lat) + 2),
      expand = FALSE) +
    geom_point(data = swo_clean, aes(lon, lat, color = tag)) + 
    theme_custom() +
    theme(legend.position = "none")

### Segment tracks w/ gaps ####
segment_tracks <- function(sp_df){

segmented_df <- data.frame()
for(i in 1:length(unique(sp_df$tag))){
  curr_id = unique(sp_df$tag)[i]
  curr_df <- sp_df %>%
    filter(tag == curr_id) %>%
    arrange(date) %>%
    mutate(diff = date - lag(date), 
           new_segment = is.na(diff) | diff > 10, 
           seg_num = cumsum(new_segment), 
           tag = paste0(tag, letters[seg_num])) %>%
    group_by(tag) %>%
    filter(n() > 20) %>%
    ungroup() %>% 
    select(c(tag, date, lc, lon, lat, sp, diff))
   
   
  segmented_df <- rbind(segmented_df, curr_df)
 } #end tag id loop
return(segmented_df)
} #end function

#albacore
alb <- readRDS(here("data/loc_data/alb_tag/validTagsLocns_NOAAonly_bathym_corrected.rds")) %>% 
  mutate(tag = as.character(tag), sp = "Albacore tuna", lc = NA) %>%
  select(c("tag", "dateRd", "lc", "lon360", "lat", "sp"))

alb <- alb %>% 
  mutate(lon = ifelse(lon360 > 180, lon360 - 360, lon360)) %>%
  select(-c(lon360))

colnames(alb) <- c("tag", "date", "lc", "lat", "sp", "lon")

  #filter out locs on land
alb_filt <- alb %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

alb_intersect <- st_intersects(alb_filt, land, sparse = FALSE)
alb_filt2 <- st_difference(alb_filt, st_union(land))
alb <- alb_filt2 %>% mutate(lon = st_coordinates(.)[,1], 
                             lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

#run segment function
alb_segment <- segment_tracks(alb) 

locs_alb <- alb_segment %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = 3832)

ggplot() + 
  geom_sf(data = land_pac, fill = "grey85", color = "grey30", linewidth = 0.2) +
  geom_sf(data = locs_alb, aes(color = tag), size = 2) +
  coord_sf(xlim = c(st_bbox(locs_alb$geometry)[1], st_bbox(locs_alb$geometry)[3]), ylim = c(st_bbox(locs_alb$geometry)[2], st_bbox(locs_alb$geometry)[4])) +
  theme_custom() +
  theme(legend.position = "none")

saveRDS(alb_segment, here("data/loc_data/processed/pre_ssm/alb_dat.rds"))

  #albacore summary stats
table(alb_segment$tag)
length(unique(alb_segment$tag)) #23
min(alb_segment$date) # 2003-08-24
max(alb_segment$date) # 2013-08-23

#blue sharks
blu_segment <- segment_tracks(blu_clean)

locs_blu <- blu_segment %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |> 
  st_transform(crs = 3832)

#check tracks after filtering
ggplot() + 
  geom_sf(data = land_pac, fill = "grey85", color = "grey30", linewidth = 0.2) +
  geom_sf(data = locs_blu, aes(color = tag), size = 2) +
  coord_sf(xlim = c(st_bbox(locs_blu$geometry)[1], st_bbox(locs_blu$geometry)[3]), ylim = c(st_bbox(locs_blu$geometry)[2], st_bbox(locs_blu$geometry)[4])) +
  geom_point(data = alb, aes(lon, lat, color = tag)) + 
  theme_custom() + 
  theme(legend.position = "none")

saveRDS(blu_segment, here("data/loc_data/processed/pre_ssm/blu_dat.rds"))

  #blue shark summary stats
table(blu_segment$tag)
length(unique(blu_segment$tag)) #60
min(blu_segment$date) # 2002-07-01
max(blu_segment$date) # 2015-10-25

#mako sharks
mako_segment <- segment_tracks(mako_clean)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    coord_sf(xlim = c(min(mako_segment$lon) - 2, max(mako_segment$lon) + 2),
      ylim = c(min(mako_segment$lat) - 2, max(mako_segment$lat) + 2),
      expand = FALSE) +
    geom_point(data = mako_segment, aes(lon, lat, color = tag)) + 
    theme_custom() +
    theme(legend.position = "none")

saveRDS(mako_segment, here("data/loc_data/processed/pre_ssm/mako_dat.rds"))

  #mako summary stats
table(mako_segment$tag)
length(unique(mako_segment$tag)) # 115
min(mako_segment$date) # 2003-07-04
max(mako_segment$date) # 2016-02-03

#swordfish
swo_segment <- segment_tracks(swo_clean)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    coord_sf(xlim = c(min(swo_segment$lon) - 2, max(swo_segment$lon) + 2),
      ylim = c(min(swo_segment$lat) - 2, max(swo_segment$lat) + 2),
      expand = FALSE) +
    geom_point(data = swo_segment, aes(lon, lat, color = tag)) + 
    theme_custom() +
    theme(legend.position = "none")

saveRDS(swo_segment, here("data/loc_data/processed/pre_ssm/swo_dat.rds"))

  #swordfish summary stats
table(swo_segment$tag)
length(unique(swo_segment$tag)) # 29
min(swo_segment$date) # 2019-03-13
max(swo_segment$date) # 2025-11-01
