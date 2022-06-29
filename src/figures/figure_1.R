# Figure 1. multi-panel figure with A) rfmo catch, B) RFMO cpue and C) rfmo effort

# Load libraries
library(tidyverse)
library(sf)
library(tmap)
library(here)

# Load plotting defaults
source(file.path(here::here(), "src/figures/plot_defaults.R"))

# Load data - use cleaned data at the 1x1 resolution using count (not mt converted to count)
list_files <- list.files(file.path(here::here(), "data/model-data/inputs/all-rfmo-models/"), 
                         pattern = "1x1_count_hooks", full.names = TRUE)

all_dat <- NULL
for(file in list_files) { 
  all_dat <- bind_rows(all_dat, read.csv(file))
}

# Calculate mean totals per year
rfmo_totals <- all_dat %>% 
  group_by(latitude, longitude, year, rfmo) %>% 
  summarise(total_catch = sum(catch, na.rm = T), 
            total_bycatch_effort = mean(bycatch_total_effort, na.rm = T), # repeated across individuals
            total_cpue = total_catch/total_bycatch_effort) %>% 
  ungroup() %>%
  group_by(latitude, longitude, rfmo) %>% 
  summarise(mean_total_catch = mean(total_catch, na.rm = T), 
            mean_bycatch_effort = mean(total_bycatch_effort, na.rm = T), 
            mean_cpue = mean(total_cpue, na.rm = T))%>% 
  ungroup() %>% 
  mutate(latitude_temp = latitude, 
         longitude_temp = longitude) %>% 
  st_as_sf(., coords = c("longitude_temp", "latitude_temp"), crs = 4326)

# Scale by rfmo
rfmo_totals_scaled <- rfmo_totals %>% 
  group_by(rfmo) %>% 
  mutate(rfmo_max_catch = max(mean_total_catch, na.rm = T), 
         rfmo_max_bycatch_effort = max(mean_bycatch_effort, na.rm = T),  
         rfmo_max_cpue = max(mean_cpue, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(mean_total_catch_scaled = mean_total_catch/rfmo_max_catch, 
         mean_bycatch_effort_scaled = mean_bycatch_effort/rfmo_max_bycatch_effort, 
         mean_cpue_scaled = mean_cpue/rfmo_max_cpue)


# Rasterize based on the center of each cell (a little annoying)

for(layer in c("mean_total_catch_scaled", "mean_bycatch_effort_scaled", "mean_cpue_scaled")) { 
  rfmo_catch_raster_1 <- rfmo_totals_scaled %>% 
    filter(latitude%%1 == 0 & longitude%%1 == 0) %>% 
    rasterize(., whole_numbers, field = layer, fun = sum, background = 0)
  
  rfmo_catch_raster_2 <- rfmo_totals_scaled %>% 
    filter(latitude%%1 != 0 & longitude%%1 == 0) %>% 
    rasterize(., whole_numbers_lon, field = layer, fun = sum, background = 0) %>% 
    resample(., whole_numbers)
  
  # not present in the data
  # rfmo_catch_raster_3 <- rfmo_totals_scaled %>% 
  #   filter(latitude%%1 == 0 & longitude%%1 != 0) %>% 
  #   rasterize(., whole_numbers_lat, field = layer, fun = sum) %>% 
  #   resample(., whole_numbers)
  
  rfmo_catch_raster_4 <- rfmo_totals_scaled %>% 
    filter(latitude%%1 != 0 & longitude%%1 != 0) %>% 
    rasterize(., fraction_numbers, field = layer, fun = sum, background = 0) %>% 
    resample(., whole_numbers)
  
  rfmo_catch_raster <- rfmo_catch_raster_1 + rfmo_catch_raster_2 + rfmo_catch_raster_4
  
  #rfmo_catch_raster <- raster::projectRaster(rfmo_catch_raster, basemap_raster)
  
  rfmo_catch_raster[rfmo_catch_raster == 0] <- NA
  
  rfmo_catch_raster <- raster::as.data.frame(rfmo_catch_raster, xy = TRUE) %>% 
    mutate(layer = ifelse(layer == 0, NA, layer))
  
  assign(layer, rfmo_catch_raster)
} 

# Figure 1a: catch
fig_1a <- ggplot() + 
  geom_tile(mean_total_catch_scaled, 
              mapping = aes(x=x, y=y, fill=layer)) + 
  scale_fill_distiller("", 
                       palette = "RdBu", na.value = NA) + 
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)), 
              mapping = aes(x=x, y=y), fill = "black") + 
  coord_sf() + 
  custom_theme 

# Figure 1b. cpue
fig_1b <- ggplot() + 
  geom_raster(mean_cpue, 
              mapping = aes(x=x, y=y, fill=layer)) + 
  scale_fill_distiller("", 
                       palette = "RdBu", na.value = NA, 
                       breaks = c(min(mean_cpue$layer, na.rm = T), max(mean_cpue$layer, na.rm = T)), 
                       labels = c("Low", "High")) + 
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)), 
              mapping = aes(x=x, y=y), fill = "black") + 
  coord_sf() + 
  custom_theme 

# Figure 1b. cpue
fig_1c <- ggplot() + 
  geom_raster(mean_bycatch_effort, 
              mapping = aes(x=x, y=y, fill=layer)) + 
  scale_fill_distiller("", 
                       palette = "RdBu", na.value = NA, trans=
                       breaks = c(min(mean_bycatch_effort$layer, na.rm = T), max(mean_bycatch_effort$layer, na.rm = T)), 
                       labels = c("Low", "High")) + 
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)), 
            mapping = aes(x=x, y=y), fill = "black") + 
  coord_sf() + 
  custom_theme 
