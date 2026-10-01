#libraries
library(tidyverse)
library(here)
library(sf)
library(terra)
library(rnaturalearth)
library(tidyquant)
library(aniMotum)
source(here("functions/keep_windows.R"))
source(here("functions/CRW_PA.R"))
set.seed(0902)

# land file
land <- ne_countries(scale = "large", returnclass = "sf") %>% st_make_valid() 
land <- st_transform(land, crs = 4326)

#load re-routed SSMs
blu <- readRDS("data/loc_data/processed/ssm/blu_ssm.rds")
mako <- readRDS("data/loc_data/processed/ssm/mako_ssm.rds")
swo <- readRDS("data/loc_data/processed/ssm/swo_ssm.rds")

#get min/max lats across species to filter domain down to 
#albacore
alb <- readRDS("data/loc_data/processed/pre_ssm/alb_dat.rds") %>%
  mutate(lon = ifelse(lon > 180, lon - 360, lon))
colnames(alb) <- c("id", "date", "sp", "lon", "lat")

min(alb$lat) - 2 #23.1
max(alb$lat) + 2 #54.6
min(alb$lon) - 2 #-182.0

#blue sharks
blu_ssm <- grab(blu, what = "predicted")
min(blu_ssm$lat) - 2 #2.2
max(blu_ssm$lat) + 2 #52.7
min(blu_ssm$lon) -2 #-163.1

#mako sharks
mako_ssm <- grab(mako, what = "predicted")
min(mako_ssm$lat) - 2 #0.8
max(mako_ssm$lat) + 2 #49.2
min(mako_ssm$lon) - 2 #-157.3

#swordfish
swo_ssm <- grab(swo, what = "predicted")
min(swo_ssm$lat) - 2 #-5.3
max(swo_ssm$lat) + 2 #51.1
min(swo_ssm$lon) - 2 #-168.4

### PA generation ####
#### albacore #####
  #have to use in house function bc tracks already regularized and light-level geolocations without error so cannot fit_ssm with animotum.
alb_pa <- data.frame()
alb <- alb %>% 
        group_by(id) %>% 
        arrange(date) %>% 
        ungroup() %>%
        mutate(time = paste('00', '00', sep = ":"),
               date_time = paste(date, time, sep = " "),
               dTime = as.POSIXct(strptime(as.character(date_time), "%Y-%m-%d %H:%M")),
               tagid = as.character(id), 
               long = as.numeric(lon), 
               lat = as.numeric(lat)) %>%
        select(c("tagid", "long", "lat", "dTime"))

#run PA generation in parallel
#create the cluster
n.cores <- parallel::detectCores() - 2
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

foreach(tagid = unique(alb$tagid)[unique(alb$tagid)>0], .packages = c("tidyverse", "adehabitatLT", "maps", "mapdata", "maptools", "sp", "raster")) %dopar% {
  #simulate CRWs -- takes
  alb_pa_r <- createCRW(alb, tagid, n.sim = 50)
  #out.csv2 = sprintf('%s/crw_sim_all_%s.csv', out.dir, tagid) #keeps all iterations in a csv file -- EN ADDED
  #write.csv(sim.alldata, file = out.csv2, row.names = F)
}

parallel::stopCluster(cl = my.cluster)

saveRDS(alb_pa_r, here("data/loc_data/processed/pa/alb_pa_routed.rds"))

#### SSM Species #####
#get bounding polygon for MOM6 domain 
ylims <- c(0, 55)
xlims <- c(-170, -100)
box_coords <- tibble(x = xlims, y = ylims) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  st_set_crs(st_crs(4326))

bounding_box <- st_bbox(box_coords) %>% st_as_sfc()

  #reproject to CRS that aligns with SSM output
bb_merc <- st_transform(bounding_box, crs = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs")
  
#read in continents polygon
land_merc = land %>%
  st_as_sfc() %>%
  st_transform(crs = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs")

land_subset <- st_intersection(land_merc, bb_merc)
plot(land_subset)

#crop where ROMS domain polygon and continents polygons intersect to get a final polygon of the CMEMS domain
grad_poly <- st_difference(bb_merc, land_subset)
plot(grad_poly)

df <- data.frame(id = seq(length(grad_poly)))
df$geometry <- grad_poly
grad_sf <- st_as_sf(df)

grad_spatVect <- vect(grad_poly)

# calculate gradient between CMEMS polygon (that fits within CMEMS dataset) and template that is 3x the CMEMS domain
  #domain that is x3 area of ROMS domain
CMEMS_large_rast <- rast(
  crs = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs",
  extent = ext(-18924.31, -9448.548, -2016.141, 10116.68), 
  resolution =  32.01272
)

#create 2D gradient
x <- rasterize(grad_spatVect, CMEMS_large_rast, fun = "mean") 

## generate gradient rasters
dist <- distance(x)
x1 <- terrain(dist, v = "slope", unit = "radians")
y1 <- terrain(dist, v = "aspect", unit = "radians")
grad.x <- -1 * x1 * cos(0.5 * pi - y1)
grad.y <- -1 * x1 * sin(0.5 * pi - y1)
grad <- c(grad.x, grad.y)

plot(grad)

#### blue sharks ######
blu_pa <- sim_fit(blu, what = "predicted", reps = 100, grad = grad)
blu_pa_filt <- sim_filter(blu_pa, keep = 0.25, flag = 1) #flag based on hazen et al., 2017 journal of applied ecology
blu_pa_r <- route_path(blu_pa_filt, centroids = TRUE)

plot(blu_pa_r[1,])
saveRDS(blu_pa_r, here("data/loc_data/processed/pa/blu_pa_routed.rds"))

#### mako sharks #####
mako_pa <- sim_fit(mako, what = "predicted", reps = 100, grad = grad)
mako_pa_filt <- sim_filter(mako_pa, keep = 0.25, flag = 1) #flag based on hazen et al., 2017 journal of applied ecology
mako_pa_r <- route_path(mako_pa_filt, centroids = TRUE)

plot(mako_pa_r[4,])
saveRDS(mako_pa_r, here("data/loc_data/processed/pa/mako_pa_routed.rds"))

#### swordfish#####
swo_pa <- sim_fit(swo, what = "predicted", reps = 100, grad = grad)
swo_pa_filt <- sim_filter(swo_pa, keep = 0.25, flag = 1) #flag based on hazen et al., 2017 journal of applied ecology
swo_pa_r <- route_path(swo_pa_filt, centroids = TRUE)

plot(swo_pa_r[1,])
saveRDS(swo_pa_r, here("data/loc_data/processed/pa/swo_pa_routed.rds"))

#### remove locations in gap windows, and tracks that have locations on land and outside of study domain ####
#combine rerouted presences and PAs into one df, then make locs a spatial feature
#albacore 
alb_locs <- readRDS("data/loc_data/processed/pre_ssm/alb_dat.rds") %>%
  mutate(lon = ifelse(lon > 180, lon - 360, lon), 
         loc_type = "presence", 
         rep = NA)
colnames(alb_locs) <- c("id", "date", "sp", "lon", "lat", "loc_type", "rep")
alb_pa <- readRDS(here("data/loc_data/processed/pa/alb_pa.rds")) %>%
  mutate(sp = "Albacore tuna")

alb <- rbind(alb_locs, alb_pa) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)

#blue sharks
blu_locs <- readRDS(here("data/loc_data/processed/ssm/blu_ssm.rds")) %>% 
  grab(what = "rerouted") %>%
  mutate(loc_type = "presence", 
         rep = NA) %>%
  select(-c("x", "y", "x.se", "y.se"))
blu_pa <- readRDS(here("data/loc_data/processed/pa/blu_pa_routed.rds")) %>% 
  unnest(sims) %>%
  mutate(loc_type = "absence") %>%
  select(-c("x", "y", "model"))

blu <- rbind(blu_locs, blu_pa) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)

#mako sharks
mako_locs <- readRDS(here("data/loc_data/processed/ssm/mako_ssm.rds")) %>% 
  grab(what = "rerouted") %>%
  mutate(loc_type = "presence", 
         rep = NA) %>%
  select(-c("x", "y", "x.se", "y.se"))
mako_pa <- readRDS(here("data/loc_data/processed/pa/mako_pa_routed.rds")) %>% 
  unnest(sims) %>%
  mutate(loc_type = "absence") %>%
  select(-c("x", "y", "model"))

mako <- rbind(mako_locs, mako_pa) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)

#swordfish
swo_locs <- readRDS(here("data/loc_data/processed/ssm/swo_ssm.rds")) %>% 
  grab(what = "rerouted") %>%
  mutate(loc_type = "presence", 
         rep = NA) %>%
  select(-c("x", "y", "x.se", "y.se"))
swo_pa <- readRDS(here("data/loc_data/processed/pa/swo_pa_routed.rds")) %>% 
  unnest(sims) %>%
  mutate(loc_type = "absence") %>%
  select(-c("x", "y", "model"))

swo <- rbind(swo_locs, swo_pa) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)

#### Gap windows #####
#albacore
alb_windows <- keep_windows(alb_locs)

alb_no_gaps <- alb %>%
  inner_join(bind_rows(alb_windows), 
             by = join_by(id, between(date, start_time, end_time)))

#blue sharks
blu_raw <- readRDS("data/loc_data/processed/pre_ssm/blu_dat.rds") %>% filter(lc != "P" & lc != "D")
colnames(blu_raw) <- c("id", "date", "lc", "sp", "lon", "lat")
blu_windows <- keep_windows(blu_raw)

blu_no_gaps <- blu %>%
  inner_join(bind_rows(blu_windows), 
             by = join_by(id, between(date, start_time, end_time)))

#mako sharks
mako_raw <- readRDS("data/loc_data/processed/pre_ssm/mako_dat.rds") %>% filter(lc != "P" & lc != "D")
colnames(mako_raw) <- c("id", "date", "lc", "sp", "lon", "lat")
mako_windows <- keep_windows(mako_raw)

mako_no_gaps <- mako %>%
  inner_join(bind_rows(mako_windows), 
             by = join_by(id, between(date, start_time, end_time)))

#swordfish
swo_raw <- readRDS("data/loc_data/processed/pre_ssm/swo_dat.rds") %>% filter(lc != "P" & lc != "D")
colnames(swo_raw) <- c("id", "date", "lc", "sp", "lon", "lat")
swo_windows <- keep_windows(swo_raw)

swo_no_gaps <- swo %>%
  inner_join(bind_rows(swo_windows), 
             by = join_by(id, between(date, start_time, end_time)))

#### Land and domain filter #####
land_union <- st_union(land) #speeds up st_filter

land_dom_filt <- function(sp_dat){
  pa_dat <- sp_dat %>% filter(loc_type == "absence")
  loc_dat <- sp_dat %>% filter(loc_type == "presence")

  #filter tracks with locations outside of domain     
  bbox <- st_as_sfc(st_bbox(c(xmin = -170, ymin = 0, xmax = -100, ymax = 55), crs = 4326))
                   
  pa_dat_dom <- pa_dat %>%
    group_by(id, rep) %>%
    mutate(domain_intersect = any(st_intersects(geometry, bbox, sparse = FALSE)), 
           omit_keep = if_else(domain_intersect, "keep", "omit")) %>%
    filter(!any(omit_keep == "omit")) %>%
    select(-c("domain_intersect", "omit_keep"))

  #filter tracks that overlap with land
  pa_dat_land <- pa_dat_dom %>%
    mutate(land_intersect = any(st_intersects(geometry, land_union, sparse = FALSE)), 
           omit_keep = if_else(land_intersect, "keep", "omit")) %>%
    filter(!any(omit_keep == "omit")) %>%
    select(-c("land_intersect", "omit_keep"))

  all_dat <- rbind(loc_dat, pa_dat_land)

  return(all_dat)
}


#albacore 
alb_filter <- land_dom_filt(alb_no_gaps)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = alb_filter, aes(color = loc_type)) + 
    coord_sf(xlim = c(-170, -100),
      ylim = c(0, 55),
      expand = FALSE) +
    theme_bw() 

#blue sharks
blu_filter <- land_dom_filt(blu_no_gaps)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = blu_filter, aes(color = loc_type)) + 
    coord_sf(xlim = c(-170, -100),
      ylim = c(0, 55),
      expand = FALSE) +
    theme_bw() 

#mako sharks
mako_filter <- land_dom_filt(mako_no_gaps)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = mako_filter, aes(color = loc_type)) + 
    coord_sf(xlim = c(-170, -100),
      ylim = c(0, 55),
      expand = FALSE) +
    theme_bw() 

#swordfish
swo_filter <- land_dom_filt(swo_no_gaps)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = swo_filter, aes(color = loc_type)) + 
    coord_sf(xlim = c(-170, -100),
      ylim = c(0, 55),
      expand = FALSE) +
    theme_bw() 

#### randomly sample PAs to get 1:1 with presences #####
pa_ratio <- function(sp_dat, ratio = 1){

  pres_abs_df <- data.frame()

for(i in 1:length(unique(sp_dat$id))){

  curr_id = unique(sp_dat$id)[i]

  print(curr_id)
  temp_dat = sp_dat %>% filter(id == curr_id)
  num_presence = temp_dat %>% filter(loc_type == "presence") %>% summarise(n = n())
  num_presence = num_presence$n

  pres_df = temp_dat %>% filter(loc_type == "presence")
  abs_df = temp_dat %>% filter(loc_type == "absence")
  abs_df_sub = abs_df[sample(nrow(abs_df), size = num_presence*ratio, replace = FALSE), ]

  pres_abs_temp <- rbind(pres_df, abs_df_sub)

  pres_abs_df <- rbind(pres_abs_df, pres_abs_temp)

} #end for loop
  return(pres_abs_df)
} #end function

#albacore 
alb_pres_abs <- pa_ratio(alb_filter)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = test, aes(color = loc_type), size = 2, alpha = 0.8) + 
    coord_sf(xlim = c(-175, -98),
      ylim = c(-5, 60),
      expand = FALSE) +
    theme_bw() +
  scale_color_manual(values = c("#77ABD9", "dodgerblue4"))

saveRDS(alb_pres_abs, here("data/loc_data/processed/pres_abs/alb_pres_abs.rds"))

#blue sharks
blu_pres_abs <- pa_ratio(blu_filter)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = blu_pres_abs, aes(color = loc_type), size = 2, alpha = 0.8) + 
    coord_sf(xlim = c(-175, -98),
      ylim = c(-5, 60),
      expand = FALSE) +
    theme_bw() +
  scale_color_manual(values = c("#77ABD9", "dodgerblue4"))

saveRDS(blu_pres_abs, here("data/loc_data/processed/pres_abs/blu_pres_abs.rds"))

#mako sharks
mako_pres_abs <- pa_ratio(mako_filter)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = mako_pres_abs, aes(color = loc_type), size = 2, alpha = 0.8) + 
    coord_sf(xlim = c(-175, -98),
      ylim = c(-5, 60),
      expand = FALSE) +
    theme_bw() +
  scale_color_manual(values = c("#77ABD9", "dodgerblue4"))

saveRDS(mako_pres_abs, here("data/loc_data/processed/pres_abs/mako_pres_abs.rds"))

#swordfish
swo_pres_abs <- pa_ratio(swo_filter)

ggplot() + 
    geom_sf(data = land, fill = "grey85", color = "grey30", linewidth = 0.2) +
    geom_sf(data = swo_pres_abs, aes(color = loc_type), size = 2, alpha = 0.8) + 
    coord_sf(xlim = c(-175, -98),
      ylim = c(-5, 60),
      expand = FALSE) +
    theme_bw() +
  scale_color_manual(values = c("#77ABD9", "dodgerblue4"))

saveRDS(swo_pres_abs, here("data/loc_data/processed/pres_abs/swo_pres_abs.rds"))


