# Visualize confidence in predictions

# Load libraries
library(raster)
library(tidyverse)
library(cowplot)
library(sf)
library(tmap)
library(here)

# Load plotting defaults
source(file.path(here::here(), "src/figures/plot_defaults.R"))

# Load data - use output data from the best model and convert to 1x1 if necessary
list_files <- list.files(file.path(here::here(), "data-updated/model-data/outputs/all-rfmo-models/quantiles"), 
                         pattern = ".rds", full.names = TRUE)

all_dat <- NULL
for(file in list_files) { 
  
  temp_rds <- readRDS(file)
  
  res <- unique(temp_rds$final_predict$spatial_notes)
  
  temp <- temp_rds$confidence_class_cell %>% 
    mutate(spatial_notes = res) %>% 
    separate(indID, sep = "[|]", into = c("latitude", "longitude")) %>% 
    mutate(latitude = as.numeric(latitude), 
           longitude = as.numeric(longitude)) %>% 
    mutate(lat = latitude, 
           lon = longitude)
  
  if(unique(temp$spatial_notes) == "center of 5x5 cell") { 
    
    # Rescale to 1x1 degree 
    locations_y <- temp %>%
      filter(spatial_notes == "center of 5x5 cell") %>%
      select(latitude) %>%
      distinct_all() %>%
      mutate(latitude_lag0 = latitude,
             latitude_lag1 = latitude-1,
             latitude_lag2 = latitude-2,
             latitude_plus1 = latitude+1,
             latitude_plus2 = latitude+2) %>%
      pivot_longer(latitude_lag0:latitude_plus2, values_to = "latitude_rescaled") %>%
      select(-name)
    
    locations_x <- temp %>%
      filter(spatial_notes == "center of 5x5 cell") %>%
      select(longitude) %>%
      distinct_all() %>%
      mutate(longitude_lag0 = longitude,
             longitude_lag1 = longitude-1,
             longitude_lag2 = longitude-2,
             longitude_plus1 = longitude+1,
             longitude_plus2 = longitude+2) %>%
      pivot_longer(longitude_lag0:longitude_plus2, values_to = "longitude_rescaled") %>%
      select(-name)
    
    # New dataset
    temp <- temp %>%
      filter(spatial_notes == "center of 5x5 cell") %>%
      left_join(locations_y) %>%
      left_join(locations_x) %>%
      mutate(longitude_rescaled = ifelse(is.na(longitude_rescaled), longitude, longitude_rescaled),
             latitude_rescaled = ifelse(is.na(latitude_rescaled), latitude, latitude_rescaled)) %>% 
      mutate(spatial_notes = "center of 1x1 cell") %>%
      group_by(latitude_rescaled, longitude_rescaled, spatial_notes) %>% 
      summarise(confidence = mean(confidence, na.rm = T), 
                .pred_class = mean(.pred_class, na.rm = T)) %>% 
      ungroup() %>% 
      rename(latitude = latitude_rescaled, 
             longitude = longitude_rescaled) %>% 
      mutate(lat = latitude, 
             lon = longitude)
  }
  
  all_dat <- all_dat %>% rbind(temp)
}

all_dat <- all_dat %>% 
  st_as_sf(., coords = c("lon", "lat"), crs = 4326)

raster_stack <- stack()

raster_1 <- all_dat %>% 
  filter(latitude%%1 == 0 & longitude%%1 == 0) 

if(nrow(raster_1) > 0) { 
  raster_1 <- raster_1 %>% 
    rasterize(., whole_numbers, field = "confidence", fun = mean, background = NA)
  
  raster_stack <- stack(raster_stack, raster_1)
}

raster_2 <- all_dat %>% 
  filter(latitude%%1 != 0 & longitude%%1 == 0) 

if(nrow(raster_2) > 0) { 
  raster_2 <- raster_2 %>% 
    rasterize(., whole_numbers, field = "confidence", fun = mean, background = NA)
  
  raster_stack <- stack(raster_stack, raster_2)
}

raster_3 <- all_dat %>% 
  filter(latitude%%1 == 0 & longitude%%1 != 0)  

if(nrow(raster_3) > 0) { 
  raster_3 <- raster_3 %>% 
    rasterize(., whole_numbers, field = "confidence", fun = mean, background = NA)
  
  raster_stack <- stack(raster_stack, raster_3)
}

raster_4 <- all_dat %>% 
  filter(latitude%%1 != 0 & longitude%%1 != 0)  

if(nrow(raster_4) > 0) { 
  raster_4 <- raster_4 %>% 
    rasterize(., whole_numbers, field = "confidence", fun = mean, background = NA)
  
  raster_stack <- stack(raster_stack, raster_4)
}

raster_comb <- calc(raster_stack, mean, na.rm = T)

raster_comb_df <- raster::as.data.frame(raster_comb, xy = TRUE)

ggplot() + 
  geom_tile(raster_comb_df, 
            mapping = aes(x=x, y=y, fill=layer, color = layer)) + 
  scale_fill_distiller("Confidence in Classification\nPrediction", palette = "RdYlBu", na.value = NA, 
                       breaks = c(floor(min(raster_comb_df$layer, na.rm = T)/0.01)*0.01, ceiling(max(raster_comb_df$layer, na.rm = T)/0.01)*0.01), 
                       limits = c(floor(min(raster_comb_df$layer, na.rm = T)/0.01)*0.01, ceiling(max(raster_comb_df$layer, na.rm = T)/0.01)*0.01),
                       guide = guide_colorbar(title.vjust = 0.8)) + 
  scale_color_distiller("", palette = "RdYlBu", na.value = NA, 
                        breaks = c(floor(min(raster_comb_df$layer, na.rm = T)/0.01)*0.01, ceiling(max(raster_comb_df$layer, na.rm = T)/0.01)*0.01),
                        limits = c(floor(min(raster_comb_df$layer, na.rm = T)/0.01)*0.01, ceiling(max(raster_comb_df$layer, na.rm = T)/0.01)*0.01),
                        guide = "none") + 
  geom_sf(data = wcpfc_boundary, fill = NA, color = "black") +
  geom_sf(data = iotc_boundary, fill = NA, color = "black") +
  geom_sf(data = iccat_boundary, fill = NA, color = "black") +
  geom_sf(data = iattc_boundary, fill = NA, color = "black") +
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
            mapping = aes(x=x, y=y), fill = "black", color = "black") +
  coord_sf() + 
  custom_theme + 
  theme(legend.position = "bottom", 
        text = element_text(size = 18)) 
