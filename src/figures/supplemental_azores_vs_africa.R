# Load libraries
library(raster)
library(tidyverse)
library(rfishbase)
library(cowplot)
library(sf)
library(tmap)
library(here)

# Load data
list_files<- list.files(file.path(here::here(), "data-updated/model-data/inputs/all-rfmo-models/"), pattern = "1x1_count_kwh.csv", full.names = T)

count_1x1_gfw <- NULL

for(file in list_files) {
  temp <- read.csv(file)
  count_1x1_gfw <- bind_rows(count_1x1_gfw, temp)
} 

t2 <- count_1x1_obs %>% filter(latitude >= 33.5204 & latitude <=43.0676 & longitude >=-35.5793 & longitude <= -20.8553)

t2x <- t2 %>% 
    dplyr::select(-colnames(t2 %>% select(-effort_units) %>% # remove columns where bycatch associated effort sums to 0
    dplyr::select_if(grepl("effort", colnames(.))) %>% 
    dplyr::select_if(colSums(., na.rm = T) == 0)))

t2x2 <- t2x %>% select(latitude, longitude, year, matches("target_effort|bycatch_total_effort")) %>%
  distinct_all() %>% select(-latitude, -longitude, -year) %>%  
  summarise_all(sum, na.rm = T) %>% t() %>% as.data.frame() %>% 
  arrange(desc(V1))

t2 <- t2 %>% 
  group_by(latitude, longitude, species_commonname) %>% 
  summarise(total_catch = sum(catch), 
            total_target_effort = sum(target_effort), 
            total_bycatch_effort = sum(bycatch_total_effort)) %>% 
  ungroup() %>% 
  st_as_sf(., coords = c("longitude", "latitude"), crs = 4326)
  
  
ggplot() + 
  geom_sf(data = t2, mapping = aes(color = total_catch)) + 
  facet_wrap(vars(species_commonname)) + 
  coord_sf() + 
  ggtitle("Total shark catch")

ggplot() + 
  geom_sf(data = t2 %>% filter(species_commonname == "BLUE SHARK"), mapping = aes(color = total_target_effort)) + 
  coord_sf()+ 
  ggtitle("Total effort associated with target catch")

ggplot() + 
  geom_sf(data = t2 %>% filter(species_commonname == "BLUE SHARK"), mapping = aes(color = total_bycatch_effort)) + 
  coord_sf()+ 
  ggtitle("Total effort associated with bycatch")

t2b <- count_1x1_obs %>% 
  filter(latitude <= 5.217866 & latitude >= -27.414690 & longitude >=  -15.407754 & longitude <=  22.560994)

t2b <- t2b %>% 
  group_by(latitude, longitude, species_commonname) %>% 
  summarise(total_catch = sum(catch), 
            total_target_effort = sum(target_effort), 
            total_bycatch_effort = sum(bycatch_total_effort)) %>% 
  ungroup() %>% 
  st_as_sf(., coords = c("longitude", "latitude"), crs = 4326)


ggplot() + 
  geom_sf(data = t2b, mapping = aes(color = total_catch)) + 
  facet_wrap(vars(species_commonname)) + 
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
            mapping = aes(x=x, y=y), fill = "gray48", color = "gray48") + 
  coord_sf(xlim = c(-15.407754, 22.560994), 
           ylim = c(5.217866, -27.414690)) + 
  ggtitle("Total shark catch")

ggplot() + 
  geom_sf(data = t2b %>% filter(species_commonname == "BLUE SHARK"), mapping = aes(color = total_target_effort)) + 
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
            mapping = aes(x=x, y=y), fill = "gray48", color = "gray48") + 
  coord_sf(xlim = c(-15.407754, 22.560994), 
           ylim = c(5.217866, -27.414690)) + 
  ggtitle("Total effort associated with target catch")

ggplot() + 
  geom_sf(data = t2b %>% filter(species_commonname == "BLUE SHARK"), mapping = aes(color = total_bycatch_effort)) + 
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
            mapping = aes(x=x, y=y), fill = "gray48", color = "gray48") + 
  coord_sf(xlim = c(-15.407754, 22.560994), 
           ylim = c(5.217866, -27.414690)) + 
  ggtitle("Total effort associated with bycatch")
