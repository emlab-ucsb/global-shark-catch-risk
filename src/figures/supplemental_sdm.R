# Species distribution model maps

# Load libraries
library(raster)
library(tidyverse)
library(cowplot)
library(sf)
library(here)

# Load plotting defaults
source(file.path(here::here(), "src/figures/plot_defaults.R"))

# Load sdm data
sdms <- read.csv(file.path(here::here(), "data-updated/iucn-sdm-data/outputs/global_combined_iucn_sdm_1x1.csv"))

# Load species lists
species_list <- read.csv(file.path(here::here(), "data-updated/rfmo-data/outputs/all_data.csv")) %>% 
  filter(gear_group == "longline") %>% 
  select(species_sciname, species_commonname) %>% 
  distinct_all() 
  
# Create plots for each species that we have data for
sdms <- sdms %>% 
  filter(species_sciname %in% unique(species_list$species_sciname)) %>% 
  left_join(species_list)

for(layer in unique(sdms$species_commonname)) { 
  
  dat_temp <- sdms %>% 
    filter(grepl(layer, species_commonname))
  
  if(nrow(dat_temp) == 0) { next }
  
  dat_temp <- dat_temp %>% 
    group_by(latitude, longitude) %>% 
    summarise(sdm = sum(sdm, na.rm = T)) %>% 
    ungroup() %>%
    mutate(longitude_orig = longitude, 
           latitude_orig = latitude) %>% 
    st_as_sf(., coords = c("longitude_orig", "latitude_orig"), crs = 4326)
  
  raster_stack <- stack()
  
  raster_1 <- dat_temp %>% 
    filter(latitude%%1 == 0 & longitude%%1 == 0) 
  
  if(nrow(raster_1) > 0) { 
    raster_1 <- raster_1 %>% 
      rasterize(., whole_numbers, field = "sdm", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_1)
  }
  
  raster_2 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 == 0) 
  
  if(nrow(raster_2) > 0) { 
    raster_2 <- raster_2 %>% 
      rasterize(., whole_numbers, field = "sdm", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_2)
  }
  
  raster_3 <- dat_temp %>% 
    filter(latitude%%1 == 0 & longitude%%1 != 0)  
  
  if(nrow(raster_3) > 0) { 
    raster_3 <- raster_3 %>% 
      rasterize(., whole_numbers, field = "sdm", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_3)
  }
  
  raster_4 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 != 0)  
  
  if(nrow(raster_4) > 0) { 
    raster_4 <- raster_4 %>% 
      rasterize(., whole_numbers, field = "sdm", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_4)
  }
  
  if(nlayers(raster_stack) > 1) { 
    raster_comb <- calc(raster_stack, mean, na.rm = T) 
  } else { 
    raster_comb <- raster_stack }
  
  raster_comb <- raster::as.data.frame(raster_comb, xy = TRUE) 
  
  # Some species are all 0s
  if(min(raster_comb$layer, na.rm = T) == mean(raster_comb$layer, na.rm = T)) {next}
  
  fig_supp <- ggplot() + 
    geom_tile(raster_comb, 
              mapping = aes(x=x, y=y, fill=layer, color = layer)) + 
    scale_fill_distiller("Probability of Occurrence", palette = "RdYlBu", na.value = NA, 
                         limits = c(0,1), 
                         breaks = c(0, 0.5, 1),
                         guide = guide_colorbar(title.vjust = 0.8)) + 
    scale_color_distiller("", palette = "RdYlBu", na.value = NA, 
                          limits = c(0,1), guide = "none") + 
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
  
  # Save
  ggsave(paste0(here::here("figures/supplemental"), 
                "/sdm_", gsub(" ", "_", 
                                   gsub("[(]|[)]", "", str_to_lower(layer))), 
                ".png"), 
         width = 7, height = 4, units = "in", dpi = 600, bg = "white")
} 
