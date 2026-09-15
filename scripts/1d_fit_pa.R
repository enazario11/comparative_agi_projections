#libraries
library(tidyverse)
library(here)
library(sf)
library(terra)
library(rnaturalearth)
library(tidyquant)
library(aniMotum)

#TO DO 
#1) apply keep window filter to PA locs

#load re-routed SSMs
blu <- readRDS("data/loc_data/processed/ssm/blu_ssm.rds")
mako <- readRDS("data/loc_data/processed/ssm/mako_ssm.rds")
swo <- readRDS("data/loc_data/processed/ssm/swo_ssm.rds")

#get min/max lats across species to filter domain down to 
#albacore
alb <- readRDS("data/loc_data/processed/pre_ssm/alb_dat.rds") %>%
  mutate(lon = ifelse(lon > 180, lon - 360, lon))
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

### create distance gradient for PA generation ####
#get bounding polygon for MOM6 domain 
ylims <- c(0, 55)
xlims <- c(-170, -110)
box_coords <- tibble(x = xlims, y = ylims) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  st_set_crs(st_crs(4326))

domain <- st_bbox(box_coords) %>% st_as_sfc()

land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
land <- st_transform(land, crs = 4326)
land_vect <- vect(land)

land_subset <- st_intersection(land, domain)
plot(land_subset)

grad_poly <- st_difference(domain, land_subset)
plot(grad_poly)

df <- data.frame(id = seq(length(grad_poly)))
df$geometry <- grad_poly
grad_sf <- st_as_sf(df)

grad_spatVect <- vect(grad_poly)

#create shape that's larger than the domain
large_rast <- rast(
  crs = "EPSG:4326", 
  extent = ext(-180, -80, -10, 70), 
  resolution = 0.25
)

#create 2D gradient
#create 2D gradient
x <- rasterize(grad_spatVect, large_rast, fun = "mean") 
dist <- distance(x)

# create gradient rasters
x1 <- terrain(dist, v = "slope", unit = "radians")
y1 <- terrain(dist, v = "aspect", unit = "radians")
grad.x <- -1 * x1 * cos(0.5 * pi - y1)
grad.y <- -1 * x1 * sin(0.5 * pi - y1)
grad <- c(grad.x, grad.y)

plot(grad)
