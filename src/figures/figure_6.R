# Figure 6 - 6 panel plot of hawaiian islands and interaction risk for C longimanus, S. lewini, A. pelagicus, I. paucus, A. superciliosus, S. zygaena

# Some important notes
# 1) Data we are using are 1x1 resolution (with 5x5 resolution evenly redistributed to 1x1 cells)
#    and count for bycatch catch units and metric tonnes for target catch units. All effort is
#    hooks from total target effort
# 2) Some RFMOs report their data at different spatial resolutions (1x1 degrees with degrees
#    centered around whole numbers [e.g., 150] vs centered around half numbers [e.g., 150.5]).
#    To combat this, we first rasterized data in groups using "like" spatial resolutions. We 
#    then re-sampled to a common, arbitrarily chosen CRS using nearest neighbor methods. This
#    method reduced artefacts of slight 0.5 degree shifts among RFMO reporting. If RFMO reporting
#    areas overlapped, we took the mean of groupings across overlapping regions.

# Load libraries
library(raster)
library(tidyverse)
library(rfishbase)
library(cowplot)
library(sf)
library(tmap)
library(here)

# Load plotting defaults
source(file.path(here::here(), "src/figures/plot_defaults.R"))

# Load data - use cleaned data at the 1x1 resolution using count (not mt converted to count)
list_files <- list.files(file.path(here::here(), "data-updated/model-data/outputs/all-rfmo-models"), 
                         pattern = "_untuned_final_predict", full.names = TRUE)

all_dat <- NULL
for(file in list_files) { 
  temp <- read.csv(file) %>% 
    select(rfmo, year, latitude, longitude, species_commonname, species_sciname, spatial_notes, .final_pred)
  
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
             latitude_rescaled = ifelse(is.na(latitude_rescaled), latitude, latitude_rescaled),
             .final_pred= .final_pred/25) %>% 
      mutate(spatial_notes = "center of 1x1 cell") %>%
      group_by(rfmo, year, latitude_rescaled, longitude_rescaled, species_commonname, species_sciname, spatial_notes) %>% 
      summarise(.final_pred = sum(.final_pred, na.rm = T)) %>% 
      ungroup() %>% 
      rename(latitude = latitude_rescaled, 
             longitude = longitude_rescaled)
  }
  
  all_dat <- all_dat %>% rbind(temp)
}

# Rasterize based on the center of each cell (a little annoying)
species_relevant <- c("CARCHARHINUS LONGIMANUS", "SPHYRNA LEWINI", 
                      "ALOPIAS PELAGICUS", "ISURUS PAUCUS", 
                      "ALOPIAS SUPERCILIOSUS", "SPHYRNA ZYGAENA")

for(layer in species_relevant) { 
  
  dat_temp <- all_dat %>% 
      filter(grepl(layer, species_sciname) & rfmo == "WCPFC")
  
  dat_temp <- dat_temp %>% 
    group_by(rfmo, year, latitude, longitude) %>% 
    summarise(total_pred = sum(.final_pred, na.rm = T)) %>% 
    ungroup() %>% 
    group_by(rfmo, latitude, longitude) %>% 
    summarise(mean_total_pred = mean(total_pred, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(longitude_orig = longitude, 
           latitude_orig = latitude) %>% 
    st_as_sf(., coords = c("longitude_orig", "latitude_orig"), crs = 4326)
  
  raster_stack <- stack()
  
  raster_1 <- dat_temp %>% 
    filter(latitude%%1 == 0 & longitude%%1 == 0) 
  
  if(nrow(raster_1) > 0) { 
    raster_1 <- raster_1 %>% 
      rasterize(., whole_numbers, field = "mean_total_pred", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_1)
  }
  
  raster_2 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 == 0) 
  
  if(nrow(raster_2) > 0) { 
    raster_2 <- raster_2 %>% 
      rasterize(., whole_numbers, field = "mean_total_pred", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_2)
  }
  
  raster_3 <- dat_temp %>% 
    filter(latitude%%1 == 0 & longitude%%1 != 0)  
  
  if(nrow(raster_3) > 0) { 
    raster_3 <- raster_3 %>% 
      rasterize(., whole_numbers, field = "mean_total_pred", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_3)
  }
  
  raster_4 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 != 0)  
  
  if(nrow(raster_4) > 0) { 
    raster_4 <- raster_4 %>% 
      rasterize(., whole_numbers, field = "mean_total_pred", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_4)
  }
  
  if(nlayers(raster_stack) > 1) { 
    raster_comb <- calc(raster_stack, mean, na.rm = T) 
  } else { 
    raster_comb <- raster_stack }
  
  raster_comb <- raster::as.data.frame(raster_comb, xy = TRUE) 
  
  assign(str_to_lower(gsub(" ", "_", layer)), raster_comb)
} 

# Loop for each species
for(iter in 1:length(species_relevant)) { 
  
  temp <- eval(as.name(str_to_lower(gsub(" ", "_", species_relevant[iter]))))
  
  temp_plot <- ggplot() + 
    geom_tile(temp, 
              mapping = aes(x=x, y=y, fill=layer, color = layer)) + 
    scale_fill_distiller("", palette = "RdYlBu", na.value = NA, 
                         breaks = c(min(temp$layer, na.rm = T), max(temp$layer, na.rm = T)), 
                         labels = c("Low", "High"), 
                         guide = guide_colorbar(title.vjust = 0.8)) + 
    scale_color_distiller("", palette = "RdYlBu", na.value = NA, 
                          breaks = c(min(temp$layer, na.rm = T), max(temp$layer, na.rm = T)), guide = "none") + 
    geom_sf(data = wcpfc_boundary, fill = NA, color = "black") +
    geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
              mapping = aes(x=x, y=y), fill = "black", color = "black") +
    coord_sf(ylim = c(5, 40), xlim = c(-179, -149.5), expand = 0) + 
    custom_theme + 
    theme(legend.position = "bottom",
          text = element_text(size = 18))
  
  if(iter == 1) { 
    legend <- get_legend(temp_plot)
  }
  
  temp_plot <- temp_plot + 
    theme(legend.position = "none")
  
  assign(paste0(str_to_lower(gsub(" ", "_", species_relevant[iter])), "_plot"), temp_plot)
}

# Final plot
final_plot <- ggdraw() + 
  draw_plot(carcharhinus_longimanus_plot, 0, 0.67, 0.5, 0.3) + 
  draw_plot(sphyrna_lewini_plot, 0.5, 0.67, 0.5, 0.3) + 
  draw_plot(alopias_pelagicus_plot, 0, 0.37, 0.5, 0.3) +
  draw_plot(isurus_paucus_plot, 0.5, 0.37, 0.5, 0.3) + 
  draw_plot(alopias_superciliosus_plot, 0, 0.07, 0.5, 0.3) +
  draw_plot(sphyrna_zygaena_plot, 0.5, 0.07, 0.5, 0.3) +
  draw_plot(legend, 0, 0, 1, 0.07) +
  draw_plot_label(label = LETTERS[1:6], x = c(0, 0.5, 0, 0.5, 0, 0.5), 
                  y = c(0.985, 0.985, 0.685, 0.685, 0.385, 0.385), 
                  hjust = 0, size = 30)

# Save
ggsave(here::here("figures/final/figure_6.png"), final_plot,
       width = 5.6, height = 14, units = "in", dpi = 600, bg = "white")
