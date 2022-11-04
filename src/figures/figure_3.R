# Figure 3. Hotspots of shark interactions with longline fisheries for with 
# (A) predicted spatial risk map for sharks globally, (B) endangered sharks, 
# (C) blue sharks, (D) shortfin mako

# Some important notes
# 1) Data we are using are 1x1 resolution (with 5x5 resolution evenly redistributed to 1x1 cells)
#    and count for shark catch units and metric tonnes for target catch units. All effort is
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
                         pattern = "_untuned_final_predict.csv", full.names = TRUE)

all_dat <- NULL
for(file in list_files) { 
  temp <- read.csv(file) %>% 
    select(rfmo, year, latitude, longitude, species_commonname, species_sciname, spatial_notes, .final_pred, catch)
  
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
             .final_pred= .final_pred/25, 
             catch = catch/25) %>% 
      mutate(spatial_notes = "center of 1x1 cell") %>%
      group_by(rfmo, year, latitude_rescaled, longitude_rescaled, species_commonname, species_sciname, spatial_notes) %>% 
      summarise(.final_pred = sum(.final_pred, na.rm = T), 
                catch = sum(catch, na.rm = T)) %>% 
      ungroup() %>% 
      rename(latitude = latitude_rescaled, 
             longitude = longitude_rescaled)
  }
  
  all_dat <- all_dat %>% rbind(temp)
}

# Get list of species and their IUCN code
# Had to look up by hand because rfishbase is not up to date
species_listing <- data.frame(
  species_sciname = c("Alopias pelagicus",
                      "Alopias superciliosus",
                      "Alopias vulpinus",
                      "Carcharhinus falciformis",
                      "Carcharhinus limbatus",
                      "Carcharhinus longimanus",
                      "Isurus oxyrinchus",
                      "Isurus paucus",
                      "Lamna nasus",
                      "Prionace glauca",
                      "Rhincodon typus",
                      "Sphyrna lewini",
                      "Sphyrna mokarran",
                      "Sphyrna zygaena"), 
  IUCN_Code = c("EN",
                "VU",
                "VU",
                "VU",
                "VU",
                "CR",
                "EN",
                "EN",
                "VU",
                "NT",
                "EN",
                "CR",
                "CR",
                "VU")) %>% 
  filter(IUCN_Code %in% c("EN", "VU", "CR")) %>% 
  mutate(species_sciname = str_to_upper(species_sciname))

# Rasterize based on the center of each cell (a little annoying)
for(layer in c("global", "endangered", "BLUE SHARK", "SHORTFIN MAKO SHARK")) { 
  
  if(layer == "global") { 
    dat_temp <- all_dat
  } else if(layer == "endangered") {
    dat_temp <- all_dat %>% 
      filter(species_sciname %in% unique(species_listing$species_sciname)) 
  } else {
    dat_temp <- all_dat %>% 
      filter(grepl(layer, species_commonname))}
  
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
  
  # # Scale by RFMO? 
  # rfmo_scaled <- NULL
  # 
  # for(rfmos in unique(dat_temp$rfmo)) { 
  #   dat_temp_rfmo <- dat_temp %>% 
  #     filter(rfmo == rfmos)
  #   
  #   quant_breaks_pred <- c(-0.1, unique(as.numeric(quantile(dat_temp_rfmo$mean_total_pred,
  #                                                            probs = seq(0, 1, 0.1), na.rm = T))))
  #   
  #   dat_temp_rfmo <- dat_temp_rfmo %>% 
  #     mutate(mean_total_pred_scaled = as.numeric(paste(cut(mean_total_pred, 
  #                                                           breaks = quant_breaks_pred, 
  #                                                           labels = c((13-length(quant_breaks_pred)):11)))))
  #   
  #   rfmo_scaled <- rfmo_scaled %>% 
  #     bind_rows(dat_temp_rfmo)
  # }
  
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
 
  raster_comb <- calc(raster_stack, mean, na.rm = T)
  
  raster_comb <- raster::as.data.frame(raster_comb, xy = TRUE) 
  
  assign(str_to_lower(gsub(" ", "_", layer)), raster_comb)
} 

# Figure 3a - all sharks
fig_3a <- ggplot() + 
  geom_tile(global, 
            mapping = aes(x=x, y=y, fill=layer, color = layer)) + 
  scale_fill_distiller("Catch Risk", palette = "RdYlBu", na.value = NA, 
                       breaks = c(min(global$layer, na.rm = T), max(global$layer, na.rm = T)), 
                       labels = c("Low", "High"), 
                       guide = guide_colorbar(title.vjust = 0.8)) + 
  scale_color_distiller("Catch Risk", palette = "RdYlBu", na.value = NA, 
                       breaks = c(min(global$layer, na.rm = T), max(global$layer, na.rm = T)), guide = "none") + 
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

legend <- get_legend(fig_3a)

fig_3a <- fig_3a + 
  theme(legend.position = "none")

# Figure 3b - endangered sharks
fig_3b <- ggplot() + 
  geom_tile(endangered, 
            mapping = aes(x=x, y=y, fill=layer, color = layer)) + 
  scale_fill_distiller("", palette = "RdYlBu", na.value = NA, 
                       breaks = c(min(endangered$layer, na.rm = T), max(endangered$layer, na.rm = T)), 
                       labels = c("Low", "High"), 
                       guide = guide_colorbar(title.vjust = 0.8)) + 
  scale_color_distiller("", palette = "RdYlBu", na.value = NA, 
                       breaks = c(min(endangered$layer, na.rm = T), max(endangered$layer, na.rm = T)), guide = "none") + 
  geom_sf(data = wcpfc_boundary, fill = NA, color = "black") +
  geom_sf(data = iotc_boundary, fill = NA, color = "black") +
  geom_sf(data = iccat_boundary, fill = NA, color = "black") +
  geom_sf(data = iattc_boundary, fill = NA, color = "black") +
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
            mapping = aes(x=x, y=y), fill = "black", color = "black") +
  coord_sf() + 
  custom_theme + 
  theme(legend.position = "none") 

# Figure 3c - blue shark
fig_3c <- ggplot() + 
  geom_tile(blue_shark, 
            mapping = aes(x=x, y=y, fill=layer, color = layer)) + 
  scale_fill_distiller("", palette = "RdYlBu", na.value = NA, 
                       breaks = c(min(blue_shark$layer, na.rm = T), max(blue_shark$layer, na.rm = T)), 
                       labels = c("Low", "High"), 
                       guide = guide_colorbar(title.vjust = 0.8)) + 
  scale_color_distiller("", palette = "RdYlBu", na.value = NA, 
                       breaks = c(min(blue_shark$layer, na.rm = T), max(blue_shark$layer, na.rm = T)), guide = "none") + 
  geom_sf(data = wcpfc_boundary, fill = NA, color = "black") +
  geom_sf(data = iotc_boundary, fill = NA, color = "black") +
  geom_sf(data = iccat_boundary, fill = NA, color = "black") +
  geom_sf(data = iattc_boundary, fill = NA, color = "black") +
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
            mapping = aes(x=x, y=y), fill = "black", color = "black") +
  coord_sf() + 
  custom_theme + 
  theme(legend.position = "none") 

# Figure 3d - shortfin mako shark
fig_3d <- ggplot() + 
  geom_tile(shortfin_mako_shark, 
            mapping = aes(x=x, y=y, fill=layer, color = layer)) + 
  scale_fill_distiller("", palette = "RdYlBu", na.value = NA, 
                       breaks = c(min(shortfin_mako_shark$layer, na.rm = T), max(shortfin_mako_shark$layer, na.rm = T)), 
                       labels = c("Low", "High"), 
                       guide = guide_colorbar(title.vjust = 0.8)) + 
  scale_color_distiller("", palette = "RdYlBu", na.value = NA, 
                       breaks = c(min(shortfin_mako_shark$layer, na.rm = T), max(shortfin_mako_shark$layer, na.rm = T)), 
                                  guide = "none") +
  geom_sf(data = wcpfc_boundary, fill = NA, color = "black") +
  geom_sf(data = iotc_boundary, fill = NA, color = "black") +
  geom_sf(data = iccat_boundary, fill = NA, color = "black") +
  geom_sf(data = iattc_boundary, fill = NA, color = "black") +
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
            mapping = aes(x=x, y=y), fill = "black", color = "black") +
  coord_sf() + 
  custom_theme + 
  theme(legend.position = "none") 

# Final plot - option 1 
final_plot <- ggdraw() + 
  draw_plot(fig_3a, 0, 0.55, 0.5, 0.45) + 
  draw_plot(fig_3b, 0.5, 0.55, 0.5, 0.45) + 
  draw_plot(fig_3c, 0, 0.1, 0.5, 0.45) + 
  draw_plot(fig_3d, 0.5, 0.1, 0.5, 0.45) + 
  draw_plot(legend, 0, 0, 1, 0.1) +
  draw_plot_label(label = LETTERS[1:4], x = c(0, 0.5, 0, 0.5), y = c(1,1,0.5,0.5), 
                  hjust = 0, size = 30)

# Save
ggsave(here::here("figures/final/figure_3_option1.png"), final_plot,
       width = 14, height = 8, units = "in", dpi = 600, bg = "white")

# Final plot - Option 2 just all sharks and vulnerable
final_plot <- ggdraw() + 
  draw_plot(fig_3a, 0, 0.1, 0.5, 0.9) + 
  draw_plot(fig_3b, 0.5, 0.1, 0.5, 0.9) + 
  draw_plot(legend, 0, 0, 1, 0.14) +
  draw_plot_label(label = LETTERS[1:2], x = c(0, 0.5), y = 1, 
                  hjust = 0, size = 30)

# Save
ggsave(here::here("figures/final/figure_3_option2.png"), final_plot,
       width = 14, height = 4, units = "in", dpi = 600, bg = "white")

# Final plot, but vertical
final_plot <- ggdraw() + 
  draw_plot(fig_3a, 0, 0.53, 1, 0.45) + 
  draw_plot(fig_3b, 0, 0.09, 1, 0.45) + 
  draw_plot(legend, 0, 0, 1, 0.14) +
  draw_plot_label(label = paste0(LETTERS[1:2], ") ", c("All shark species", 
                                                       "Threatened shark species")), x = c(0,0), y = c(1,0.56), 
                  hjust = 0, size = 30)

# Save
ggsave(here::here("figures/final/figure_3_vertical_revised.png"), final_plot,
       width = 10, height = 11, units = "in", dpi = 600, bg = "white")

# Run similar plots but for every species independently for the supplemental
for(layer in unique(all_dat$species_commonname)) { 
  
  dat_temp <- all_dat %>% 
    filter(grepl(layer, species_commonname))
  
  if(nrow(dat_temp) == 0) { next }
  
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
  
  # Some species are all 0s
  if(min(raster_comb$layer, na.rm = T) == max(raster_comb$layer, na.rm = T)) {next}
  
  fig_supp <- ggplot() + 
    geom_tile(raster_comb, 
              mapping = aes(x=x, y=y, fill=layer, color = layer)) + 
    scale_fill_distiller("Catch Risk", palette = "RdYlBu", na.value = NA, 
                         breaks = c(min(raster_comb$layer, na.rm = T), max(raster_comb$layer, na.rm = T)), 
                         labels = c("Low", "High"), 
                         guide = guide_colorbar(title.vjust = 0.8)) + 
    scale_color_distiller("", palette = "RdYlBu", na.value = NA, 
                         breaks = c(min(raster_comb$layer, na.rm = T), max(raster_comb$layer, na.rm = T)), guide = "none") + 
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
                "/figure_3_", gsub(" ", "_", 
                                   gsub("[(]|[)]", "", str_to_lower(layer))), 
                ".png"), 
         width = 7, height = 4, units = "in", dpi = 600, bg = "white")
} 

# One more supplemental plot - for observed vs expected
dat_temp <- all_dat
  
dat_temp <- dat_temp %>% 
  group_by(rfmo, year, latitude, longitude) %>% 
  summarise(total_obs = sum(catch, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(rfmo, latitude, longitude) %>% 
  summarise(mean_total_obs = mean(total_obs, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(longitude_orig = longitude, 
         latitude_orig = latitude) %>% 
  st_as_sf(., coords = c("longitude_orig", "latitude_orig"), crs = 4326)

raster_stack <- stack()

raster_1 <- dat_temp %>% 
  filter(latitude%%1 == 0 & longitude%%1 == 0) 

if(nrow(raster_1) > 0) { 
  raster_1 <- raster_1 %>% 
    rasterize(., whole_numbers, 
              field = "mean_total_obs", fun = mean, background = NA)
  
  raster_stack <- stack(raster_stack, raster_1)
}

raster_2 <- dat_temp %>% 
  filter(latitude%%1 != 0 & longitude%%1 == 0) 

if(nrow(raster_2) > 0) { 
  raster_2 <- raster_2 %>% 
    rasterize(., whole_numbers, 
              field = "mean_total_obs", fun = mean, background = NA)
  
  raster_stack <- stack(raster_stack, raster_2)
}

raster_3 <- dat_temp %>% 
  filter(latitude%%1 == 0 & longitude%%1 != 0)  

if(nrow(raster_3) > 0) { 
  raster_3 <- raster_3 %>% 
    rasterize(., whole_numbers, 
              field = "mean_total_obs", fun = mean, background = NA)
  
  raster_stack <- stack(raster_stack, raster_3)
}

raster_4 <- dat_temp %>% 
  filter(latitude%%1 != 0 & longitude%%1 != 0)  

if(nrow(raster_4) > 0) { 
  raster_4 <- raster_4 %>% 
    rasterize(., whole_numbers, 
              field = "mean_total_obs", fun = mean, background = NA)
  
  raster_stack <- stack(raster_stack, raster_4)
}

raster_comb <- calc(raster_stack, mean, na.rm = T)

raster_comb <- raster::as.data.frame(raster_comb, xy = TRUE) 

fig_supp_expected <- ggplot() + 
  geom_tile(raster_comb, 
            mapping = aes(x=x, y=y, fill=layer, color = layer)) + 
  scale_fill_distiller("", palette = "RdYlBu", na.value = NA, 
                       breaks = c(min(raster_comb$layer, na.rm = T), max(raster_comb$layer, na.rm = T)), 
                       labels = c("Low", "High"), 
                       guide = guide_colorbar(title.vjust = 0.8)) + 
  scale_color_distiller("", palette = "RdYlBu", na.value = NA, 
                        breaks = c(min(raster_comb$layer, na.rm = T), max(raster_comb$layer, na.rm = T)), guide = "none") + 
  geom_sf(data = wcpfc_boundary, fill = NA, color = "black") +
  geom_sf(data = iotc_boundary, fill = NA, color = "black") +
  geom_sf(data = iccat_boundary, fill = NA, color = "black") +
  geom_sf(data = iattc_boundary, fill = NA, color = "black") +
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
            mapping = aes(x=x, y=y), fill = "black", color = "black") +
  coord_sf() + 
  custom_theme + 
  theme(legend.position = "none") 

# Final plot
final_plot <- ggdraw() + 
  draw_plot(fig_supp_expected, 0, 0.08, 0.5, 0.88) + 
  draw_plot(fig_3a, 0.5, 0.08, 0.5, 0.88) + 
  draw_plot(legend, 0, 0, 1, 0.14) +
  draw_plot_label(label = paste0(LETTERS[1:2], ") ", c("Reported catch", "Predicted catch risk")), x = c(0, 0.5), y = 1, 
                  hjust = 0, size = 30)

# Save
ggsave(here::here("figures/supplemental/observed_vs_predicted_revised.png"), final_plot,
       width = 14, height = 4.2, units = "in", dpi = 600, bg = "white")
