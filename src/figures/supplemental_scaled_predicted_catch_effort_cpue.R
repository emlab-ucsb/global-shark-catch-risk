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
    select(rfmo, year, latitude, longitude, species_commonname, species_sciname, spatial_notes, .final_pred, catch, 
           matches("target_effort$|bycatch_total_effort$")) %>% 
    rename(effort = matches("target_effort$|bycatch_total_effort$"))
  
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
             catch = catch/25, 
             effort = effort/25) %>% 
      mutate(spatial_notes = "center of 1x1 cell") %>%
      group_by(rfmo, year, latitude_rescaled, longitude_rescaled, species_commonname, species_sciname, spatial_notes) %>% 
      summarise(.final_pred = sum(.final_pred, na.rm = T), 
                catch = sum(catch, na.rm = T), 
                effort = sum(effort, na.rm = T)) %>% 
      ungroup() %>% 
      rename(latitude = latitude_rescaled, 
             longitude = longitude_rescaled)
  }
  
  all_dat <- all_dat %>% rbind(temp)
}

# Calculate mean totals per year
rfmo_totals <- all_dat %>% 
  group_by(latitude, longitude, year, rfmo) %>% 
  summarise(.final_pred = sum(.final_pred, na.rm = T), 
            effort = mean(effort, na.rm = T), # repeated across individuals
            cpue = .final_pred/effort) %>% 
  ungroup() %>%
  group_by(latitude, longitude, rfmo) %>% 
  summarise(.final_pred = mean(.final_pred, na.rm = T), 
            effort = mean(effort, na.rm = T), 
            cpue = mean(cpue, na.rm = T))%>% 
  ungroup() %>% 
  mutate(latitude_temp = latitude, 
         longitude_temp = longitude) %>% 
  st_as_sf(., coords = c("longitude_temp", "latitude_temp"), crs = 4326) 

# Scale by rfmo
rfmo_totals_scaled <- NULL
for(rfmos in unique(rfmo_totals$rfmo)) { 
  temp_dat <- rfmo_totals %>% 
    filter(rfmo == rfmos)
  
  quant_breaks_catch <- c(-0.1, unique(as.numeric(quantile(temp_dat$.final_pred, probs = seq(0, 1, 0.1), na.rm = T))))
  quant_breaks_effort <- c(-0.1,unique(as.numeric(quantile(temp_dat$effort, probs = seq(0, 1, 0.1), na.rm = T))))
  quant_breaks_cpue <- c(-0.1, unique(as.numeric(quantile(temp_dat$cpue, probs = seq(0, 1, 0.1), na.rm = T))))
  
  temp_dat <- temp_dat %>% 
    mutate(mean_total_catch_scaled = as.numeric(paste(cut(.final_pred, 
                                                          breaks = quant_breaks_catch, 
                                                          labels = c((13-length(quant_breaks_catch)):11)))), 
           mean_bycatch_effort_scaled = as.numeric(paste(cut(effort, 
                                                             breaks = quant_breaks_effort, 
                                                             labels = c((13-length(quant_breaks_effort)):11)))), 
           mean_cpue_scaled = as.numeric(paste(cut(cpue, 
                                                   breaks = quant_breaks_cpue, 
                                                   labels = c((13-length(quant_breaks_cpue)):11)))))
  
  rfmo_totals_scaled <- rfmo_totals_scaled %>% 
    bind_rows(temp_dat)
  
}

# Rasterize based on the center of each cell (a little annoying)
for(layer in c("mean_total_catch_scaled", "mean_bycatch_effort_scaled", "mean_cpue_scaled")) { 
  rfmo_catch_raster_1 <- rfmo_totals_scaled %>% 
    filter(latitude%%1 == 0 & longitude%%1 == 0) %>% 
    rasterize(., whole_numbers, field = layer, fun = mean, background = NA)
  
  rfmo_catch_raster_2 <- rfmo_totals_scaled %>% 
    filter(latitude%%1 != 0 & longitude%%1 == 0) %>% 
    rasterize(., whole_numbers_lon, field = layer, fun = mean, background = NA) %>% 
    resample(., whole_numbers, method = "ngb")
  
  # not present in the data
  # rfmo_catch_raster_3 <- rfmo_totals_scaled %>% 
  #   filter(latitude%%1 == 0 & longitude%%1 != 0) %>% 
  #   rasterize(., whole_numbers_lat, field = layer, fun = sum) %>% 
  #   resample(., whole_numbers)
  
  rfmo_catch_raster_4 <- rfmo_totals_scaled %>% 
    filter(latitude%%1 != 0 & longitude%%1 != 0) %>% 
    rasterize(., fraction_numbers, field = layer, fun = mean, background = NA) %>% 
    resample(., whole_numbers, method = "ngb")
  
  rfmo_catch_raster <- calc(stack(rfmo_catch_raster_1, rfmo_catch_raster_2, rfmo_catch_raster_4), 
                            mean, na.rm = T)
  
  rfmo_catch_raster[rfmo_catch_raster == 0] <- NA
  
  rfmo_catch_raster <- raster::as.data.frame(rfmo_catch_raster, xy = TRUE) %>% 
    mutate(layer = ifelse(layer == 0, NA, layer))
  
  assign(layer, rfmo_catch_raster)
} 

# Figure 1a: catch
fig_1a <- ggplot() + 
  geom_tile(mean_total_catch_scaled, 
            mapping = aes(x=x, y=y, fill=layer)) + 
  scale_fill_distiller("", palette = "RdYlBu", na.value = NA, 
                       breaks = c(min(mean_total_catch_scaled$layer, na.rm = T), max(mean_total_catch_scaled$layer, na.rm = T)), 
                       labels = c("Low", "High"), 
                       guide = guide_colorbar(title.vjust = 0.8)) + 
  geom_sf(data = wcpfc_boundary, fill = NA, color = "black") + 
  geom_sf(data = iotc_boundary, fill = NA, color = "black") + 
  geom_sf(data = iccat_boundary, fill = NA, color = "black") + 
  geom_sf(data = iattc_boundary, fill = NA, color = "black") + 
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
            mapping = aes(x=x, y=y), fill = "black", color = "black") +
  coord_sf() + 
  custom_theme + 
  theme(legend.position = "bottom") 

legend <- get_legend(fig_1a)

fig_1a <- fig_1a + 
  theme(legend.position = "none")

# Figure 1b. effort
fig_1b <- ggplot() + 
  geom_raster(mean_bycatch_effort_scaled, 
              mapping = aes(x=x, y=y, fill=layer)) + 
  scale_fill_distiller("", 
                       palette = "RdYlBu", na.value = NA,
                       breaks = c(min(mean_bycatch_effort_scaled$layer, na.rm = T), max(mean_bycatch_effort_scaled$layer, na.rm = T)), 
                       labels = c("Low", "High")) + 
  geom_sf(data = wcpfc_boundary, fill = NA, color = "black") + 
  geom_sf(data = iotc_boundary, fill = NA, color = "black") + 
  geom_sf(data = iccat_boundary, fill = NA, color = "black") + 
  geom_sf(data = iattc_boundary, fill = NA, color = "black") + 
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)), 
            mapping = aes(x=x, y=y), fill = "black", color = "black") + 
  coord_sf() + 
  custom_theme + 
  theme(legend.position = "none")

# Figure 1c. cpue
fig_1c <- ggplot() + 
  geom_raster(mean_cpue_scaled, 
              mapping = aes(x=x, y=y, fill=layer)) + 
  scale_fill_distiller("", 
                       palette = "RdYlBu", na.value = NA, 
                       breaks = c(min(mean_cpue_scaled$layer, na.rm = T), max(mean_cpue_scaled$layer, na.rm = T)), 
                       labels = c("Low", "High")) + 
  geom_sf(data = wcpfc_boundary, fill = NA, color = "black") + 
  geom_sf(data = iotc_boundary, fill = NA, color = "black") + 
  geom_sf(data = iccat_boundary, fill = NA, color = "black") + 
  geom_sf(data = iattc_boundary, fill = NA, color = "black") + 
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)), 
            mapping = aes(x=x, y=y), fill = "black", color = "black") + 
  coord_sf() + 
  custom_theme + 
  theme(legend.position = "none")

# Put them all together
final_plot <- ggdraw() + 
  draw_plot(fig_1a, 0, 0.1, 0.33, 0.9) + 
  draw_plot(fig_1b, 0.33, 0.1, 0.33, 0.9) + 
  draw_plot(fig_1c, 0.66, 0.1, 0.33, 0.9) + 
  draw_plot(legend, 0, 0, 1, 0.22) + 
  draw_plot_label(label = c("A", "B", "C"), 
                  x = c(0, 0.33, 0.66), y = 1, hjust = 0)

# Save
ggsave(here::here("figures/supplemental/rfmo_scaled_predict.png"), final_plot,
       width = 10, height = 2.5, units = "in", dpi = 600, bg = "white")
