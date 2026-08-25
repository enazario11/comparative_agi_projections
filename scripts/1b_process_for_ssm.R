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
land <- ne_countries(scale = "large", returnclass = "sf")
land_pac <- st_transform(land, crs = 3832)

### argos SDA filter ####
group_sda_filter <- function (sp_df){

df_clean <- data.frame() #empty df for filtered tracks grouped by ptt
  
for(i in 1:length(unique(sp_df$tag))){
  curr_id = unique(sp_df$tag)[i]
  curr_df <- sp_df %>%
    filter(tag == curr_id) %>%
    arrange(date)

 if(nrow(curr_df) > 50){
  filter <- sdafilter(curr_df$lat, curr_df$lon, curr_df$date, curr_df$lc, vmax = 3) #consider setting vmax to 3 and defaults for dist/ang for all sp. following O'Toole et al., 2021
  curr_clean <- filter(curr_df, filter == "not")
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
  select(c(Shark.ID, date2, lc, long, lat, sp))
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
         select(c(PTT, date, lc, long, lat, sp))

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
alb <- readRDS(here("data/loc_data/alb_tag/allValidTagsLocns_NOAAonly.rds")) %>% 
  mutate(tag = as.character(tag), sp = "Albacore tuna", lc = NA) %>%
  select(c("tag", "dateRd", "lc", "lon360", "lat", "sp"))

alb <- alb %>% 
  mutate(lon = ifelse(lon360 > 180, lon360 - 360, lon360)) %>%
  select(-c(lon360))

colnames(alb) <- c("tag", "date", "lc", "lat", "sp", "lon")

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

#swordfish


