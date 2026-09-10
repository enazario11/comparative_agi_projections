#libraries
library(tidyverse)
library(here)
library(sf)
library(rnaturalearth)
library(tidyquant)
library(aniMotum)

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

### fit ssm ####
#skip albacore bc cleaning already happened -- ssm not needed
#blue sharks
blu <- readRDS(here("data/loc_data/processed/pre_ssm/blu_dat.rds")) %>% select(-c(diff)) %>% filter(lc != "P")
colnames(blu) <- c("id", "date", "lc", "lon", "lat", "sp")

blu %>%
  group_by(id) %>%
  arrange(id, date) %>%
  mutate(diff = date - lag(date), 
         diff_hours = as.numeric(diff, units = "hours")) %>%
  summarise(med_diff = median(diff_hours, na.rm = T)) %>%
  ungroup() %>%
  summarise(all_mean = mean(med_diff)) #average time step btwn positions for all tracks is 42 hours

blu_ssm <- fit_ssm(blu, 
                   spdf = FALSE, #turn off sda filter
                   date = "date", 
                   coord = c("lon", "lat"), 
                   model = "crw", 
                   time.step = 42) 

#look at outputs
summary(blu_ssm)
plot(blu_ssm[1:4,], what = "predicted", type = 1, pages = 1)
plot(blu_ssm[2,], what = "predicted", type = 2)

blu_ssm_r <- route_path(blu_ssm, what = "predicted")
aniMotum::map(blu_ssm_r, what = "predicted")|aniMotum::map(blu_ssm_r, what = "rerouted")  
saveRDS(blu_ssm_r, "data/loc_data/processed/ssm/blu_ssm.rds")

#check residuals
resid_blu <- osar(blu_ssm_r)

plot(resid_blu, type = "qq", pages = 0)
plot(resid_blu, type = "acf", pages = 0)
plot(resid_blu, type = "ts", pages = 0)

#mako sharks
mako <- readRDS(here("data/loc_data/processed/pre_ssm/mako_dat.rds")) %>% select(-c(diff))
colnames(mako) <- c("id", "date", "lc", "lon", "lat", "sp")

mako %>%
  group_by(id) %>%
  arrange(id, date) %>%
  mutate(diff = date - lag(date), 
         diff_hours = as.numeric(diff, units = "hours")) %>%
  summarise(med_diff = median(diff_hours, na.rm = T)) %>%
  ungroup() %>%
  summarise(all_mean = mean(med_diff)) #average time step btwn positions for all tracks is 45 hours

mako_ssm <- fit_ssm(mako, 
                   spdf = FALSE, #turn off sda filter
                   date = "date", 
                   coord = c("lon", "lat"), 
                   model = "crw", 
                   time.step = 45) 

#look at outputs
summary(mako_ssm)
plot(mako_ssm[1:4,], what = "predicted", type = 1, pages = 1)
plot(mako_ssm[2,], what = "predicted", type = 2)

mako_ssm_r <- route_path(mako_ssm, what = "predicted")
aniMotum::map(mako_ssm_r, what = "predicted")|aniMotum::map(mako_ssm_r, what = "rerouted")  
saveRDS(mako_ssm_r, "data/loc_data/processed/ssm/mako_ssm.rds")

#check residuals
resid_mako <- osar(mako_ssm_r)

plot(resid_mako, type = "qq", pages = 0)
plot(resid_mako, type = "acf", pages = 0)
plot(resid_mako, type = "ts", pages = 0)

#swordfish 
swo <- readRDS(here("data/loc_data/processed/pre_ssm/swo_dat.rds")) %>% select(-c(diff))
colnames(swo) <- c("id", "date", "lc", "lon", "lat", "sp")

swo %>%
  group_by(id) %>%
  arrange(id, date) %>%
  mutate(diff = date - lag(date), 
         diff_hours = as.numeric(diff, units = "hours")) %>%
  summarise(med_diff = median(diff_hours, na.rm = T)) %>%
  ungroup() %>%
  summarise(all_mean = mean(med_diff)) #average time step btwn positions for all tracks is 45 hours

swo_ssm <- fit_ssm(swo, 
                   spdf = FALSE, #turn off sda filter
                   date = "date", 
                   coord = c("lon", "lat"), 
                   model = "crw", 
                   time.step = 45) 

#look at outputs
summary(swo_ssm)
plot(swo_ssm[1:4,], what = "predicted", type = 1, pages = 1)
plot(swo_ssm[2,], what = "predicted", type = 2)

swo_ssm_r <- route_path(swo_ssm, what = "predicted")
aniMotum::map(swo_ssm_r, what = "predicted")|aniMotum::map(swo_ssm_r, what = "rerouted")  
saveRDS(swo_ssm_r, "data/loc_data/processed/ssm/swo_ssm.rds")

#check residuals
resid_swo <- osar(swo_ssm_r)

plot(resid_swo, type = "qq", pages = 0)
plot(resid_swo, type = "acf", pages = 0)
plot(resid_swo, type = "ts", pages = 0)
