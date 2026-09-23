#libraries
library(tidyverse)
library(here)
library(sf)
library(rnaturalearth)
library(tidyquant)
library(aniMotum)
source(here("functions/keep_windows.R"))
set.seed(0902)

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

### load sp data
blu <- readRDS("data/loc_data/processed/pre_ssm/blu_dat.rds") %>% filter(lc != "P" & lc != "D")
colnames(blu) <- c("id", "date", "lc", "sp", "lon", "lat")

mako <- readRDS("data/loc_data/processed/pre_ssm/mako_dat.rds")
colnames(mako) <- c("id", "date", "lc", "sp", "lon", "lat")

swo <- readRDS("data/loc_data/processed/pre_ssm/swo_dat.rds")
colnames(swo) <- c("id", "date", "lc", "sp", "lon", "lat")

#calculate keep windows
blu_windows <- keep_windows(blu)
mako_windows <- keep_windows(mako)
swo_windows <- keep_windows(swo)

### fit ssm ####
#skip albacore bc cleaning already happened -- ssm not needed
#blue sharks
blu %>%
  group_by(id) %>%
  arrange(id, date) %>%
  mutate(diff = date - lag(date), 
         diff_hours = as.numeric(diff, units = "hours")) %>%
  summarise(med_diff = median(diff_hours, na.rm = T)) %>%
  ungroup() %>%
  summarise(all_mean = mean(med_diff)) #average time step btwn positions for all tracks is 43 hours

blu_ssm <- fit_ssm(blu, 
                   vmax = 3, #Poisson et al., 2024 Fish. Res. 
                   date = "date", 
                   coord = c("lon", "lat"), 
                   model = "crw", 
                   time.step = 43) 

#look at outputs
plot(blu_ssm[1:4,], what = "predicted", type = 1, pages = 1)
plot(blu_ssm[2,], what = "predicted", type = 2)

blu_ssm_r <- route_path(blu_ssm, what = "predicted")
aniMotum::map(blu_ssm_r, what = "predicted")|aniMotum::map(blu_ssm_r, what = "rerouted")  
saveRDS(blu_ssm_r, "data/loc_data/processed/ssm/blu_ssm.rds")

#apply keep windows filter
blu_ssm_df <- grab(blu_ssm_r, what = "rerouted")

blu_ssm_clean <- blu_ssm_df %>%
  inner_join(bind_rows(blu_windows), 
             by = join_by(id, between(date, start_time, end_time)))

ggplot() +
  geom_path(data = blu_ssm_clean, aes(lon, lat, color = id)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  coord_equal()

saveRDS(blu_ssm_clean, "data/loc_data/processed/ssm_df/blu_ssm_df.rds")

#check residuals
resid_blu <- osar(blu_ssm_r)

plot(resid_blu, type = "qq", pages = 0)
plot(resid_blu, type = "acf", pages = 0)
plot(resid_blu, type = "ts", pages = 0)

#mako sharks
mako %>%
  group_by(id) %>%
  arrange(id, date) %>%
  mutate(diff = date - lag(date), 
         diff_hours = as.numeric(diff, units = "hours")) %>%
  summarise(med_diff = median(diff_hours, na.rm = T)) %>%
  ungroup() %>%
  summarise(all_mean = mean(med_diff)) #average time step btwn positions for all tracks is 45 hours

mako_ssm <- fit_ssm(mako, 
                   vmax = 4.5, #Byrne et al., 2024 Divers. Distrib.
                   date = "date", 
                   coord = c("lon", "lat"), 
                   model = "rw", 
                   time.step = 45, 
                   control=ssm_control(verbose=1)) 

#look at outputs
summary(mako_ssm)
plot(mako_ssm[1:4,], what = "predicted", type = 1, pages = 1)
plot(mako_ssm[2,], what = "predicted", type = 2)

mako_ssm_r <- route_path(mako_ssm, what = "predicted")
aniMotum::map(mako_ssm_r, what = "predicted")|aniMotum::map(mako_ssm_r, what = "rerouted")  
saveRDS(mako_ssm_r, "data/loc_data/processed/ssm/mako_ssm.rds")

#apply keep windows filter
mako_ssm_df <- grab(mako_ssm_r, what = "rerouted")

mako_ssm_clean <- mako_ssm_df %>%
  inner_join(bind_rows(mako_windows), 
             by = join_by(id, between(date, start_time, end_time)))

ggplot() +
  geom_path(data = mako_ssm_clean, aes(lon, lat, color = id)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  coord_equal()

saveRDS(mako_ssm_clean, "data/loc_data/processed/ssm_df/mako_ssm_df.rds")

#check residuals
resid_mako <- osar(mako_ssm_r)

plot(resid_mako, type = "qq", pages = 0)
plot(resid_mako, type = "acf", pages = 0)
plot(resid_mako, type = "ts", pages = 0)

#swordfish 
swo %>%
  group_by(id) %>%
  arrange(id, date) %>%
  mutate(diff = date - lag(date), 
         diff_hours = as.numeric(diff, units = "hours")) %>%
  summarise(med_diff = median(diff_hours, na.rm = T)) %>%
  ungroup() %>%
  summarise(all_mean = mean(med_diff)) #average time step btwn positions for all tracks is 26 hours

swo_ssm <- fit_ssm(swo, 
                   vmax = 3, #specifics not available, used guidelines set by O'Toole et al., 2021 Methods Ecol Evol
                   date = "date", 
                   coord = c("lon", "lat"), 
                   model = "crw", 
                   time.step = 29) 

#look at outputs
summary(swo_ssm)
plot(swo_ssm[1:4,], what = "predicted", type = 1, pages = 1)
plot(swo_ssm[1,], what = "predicted", type = 2)

swo_ssm_r <- route_path(swo_ssm, what = "predicted")
aniMotum::map(swo_ssm_r, what = "predicted")|aniMotum::map(swo_ssm_r, what = "rerouted")  
saveRDS(swo_ssm_r, "data/loc_data/processed/ssm/swo_ssm.rds")

#apply keep windows filter
swo_ssm_df <- grab(swo_ssm_r, what = "rerouted")

swo_ssm_clean <- swo_ssm_df %>%
  inner_join(bind_rows(swo_windows), 
             by = join_by(id, between(date, start_time, end_time)))

ggplot() +
  geom_path(data = swo_ssm_clean, aes(lon, lat, color = id)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  coord_equal()

saveRDS(swo_ssm_clean, "data/loc_data/processed/ssm_df/swo_ssm_df.rds")

#check residuals
resid_swo <- osar(swo_ssm_r)

plot(resid_swo, type = "qq", pages = 0)
plot(resid_swo, type = "acf", pages = 0)
plot(resid_swo, type = "ts", pages = 0)
