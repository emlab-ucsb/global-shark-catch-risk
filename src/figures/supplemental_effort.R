# Figures for effort data used

# Load libraries
library(raster)
library(tidyverse)
library(cowplot)
library(sf)
library(here)

# Load plotting defaults
source(file.path(here::here(), "src/figures/plot_defaults.R"))

# Load data
gfw <-  read.csv(file.path(here::here(), "data-updated/global-fishing-watch/outputs/binned_global_gfw_1x1.csv"))

rfmo_data <- list.files(file.path(here::here(), "data-updated/model-data/inputs/all-rfmo-models/"), 
                        pattern = "_ll_data_1x1_count_hooks.csv", full.names = TRUE)

gfw <- gfw %>% 
  filter(gear_group == "longline") %>%
  group_by(year, latitude, longitude) %>% 
  summarise(total_fishing_kwh = sum(total_fishing_kwh)+0.001) %>% 
  ungroup() 

rfmo <- NULL
for(file in rfmo_data) { 
  rfmo <- rfmo %>% bind_rows(read.csv(file))
}

rfmo <- rfmo %>%
  select(year, latitude, longitude, bycatch_total_effort, target_effort) %>% 
  distinct_all()

# GFW effort
for(layer in unique(gfw$year)) { 
  
  dat_temp <- gfw %>% 
    filter(year == layer)
  
  if(nrow(dat_temp) == 0) { next }
  
  dat_temp <- dat_temp %>% 
    mutate(longitude_orig = longitude, 
           latitude_orig = latitude) %>% 
    st_as_sf(., coords = c("longitude_orig", "latitude_orig"), crs = 4326)
  
  raster_stack <- stack()
  
  raster_1 <- dat_temp %>% 
    filter(latitude%%1 == 0 & longitude%%1 == 0) 
  
  if(nrow(raster_1) > 0) { 
    raster_1 <- raster_1 %>% 
      rasterize(., whole_numbers, field = "total_fishing_kwh", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_1)
  }
  
  raster_2 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 == 0) 
  
  if(nrow(raster_2) > 0) { 
    raster_2 <- raster_2 %>% 
      rasterize(., whole_numbers, field = "total_fishing_kwh", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_2)
  }
  
  raster_3 <- dat_temp %>% 
    filter(latitude%%1 == 0 & longitude%%1 != 0)  
  
  if(nrow(raster_3) > 0) { 
    raster_3 <- raster_3 %>% 
      rasterize(., whole_numbers, field = "total_fishing_kwh", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_3)
  }
  
  raster_4 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 != 0)  
  
  if(nrow(raster_4) > 0) { 
    raster_4 <- raster_4 %>% 
      rasterize(., whole_numbers, field = "total_fishing_kwh", fun = mean, background = NA)
    
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
    scale_fill_distiller("Fishing Effort (kwh)", palette = "RdYlBu", na.value = NA, 
                         labels = scales::comma_format(),
                         guide = guide_colorbar(title.vjust = 0.8)) + 
    scale_color_distiller("", palette = "RdYlBu",   
                         na.value = NA, guide = "none") + 
    geom_sf(data = wcpfc_boundary, fill = NA, color = "black") +
    geom_sf(data = iotc_boundary, fill = NA, color = "black") +
    geom_sf(data = iccat_boundary, fill = NA, color = "black") +
    geom_sf(data = iattc_boundary, fill = NA, color = "black") +
    geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
              mapping = aes(x=x, y=y), fill = "black", color = "black") +
    coord_sf() + 
    custom_theme + 
    theme(legend.position = "bottom", 
          text = element_text(size = 18), 
          legend.text = element_text(angle = -90, hjust = 0.5, vjust = 0.5)) 
  
  if(layer == 2012) { 
    legend <- get_legend(fig_supp)}
  
  # fig_supp <- fig_supp + theme(legend.position = "none")
  
  assign(paste0("gfw_mean_", layer), fig_supp)
} 

gfw_mean_final_plot <- ggdraw() + 
  draw_plot(gfw_mean_2012, 0, 0.66, 0.33, 0.33) + 
  draw_plot(gfw_mean_2013, 0.33, 0.66, 0.33, 0.33) +  
  draw_plot(gfw_mean_2014, 0.66, 0.66, 0.33, 0.33) +  
  draw_plot(gfw_mean_2015, 0, 0.33, 0.33, 0.33) + 
  draw_plot(gfw_mean_2016, 0.33, 0.33, 0.33, 0.33) +  
  draw_plot(gfw_mean_2017, 0.66, 0.33, 0.33, 0.33) +  
  draw_plot(gfw_mean_2018, 0, 0, 0.33, 0.33) + 
  draw_plot(gfw_mean_2019, 0.33, 0, 0.33, 0.33) + 
  draw_plot(gfw_mean_2020, 0.66, 0, 0.33, 0.33) + 
  draw_plot_label(label = paste0(LETTERS[1:9], ") ", 2012:2020), x = rep(c(0, 0.33, 0.66), 3), 
                  y = rep(c(1, 0.66, 0.33), each = 3))

ggsave(file.path(here::here(), "figures/supplemental/mean_gfw_kwh.png"), plot = gfw_mean_final_plot,
       bg = "white", dpi = 600, width = 20, height = 12.5)

# Target effort
for(layer in unique(rfmo$year)) { 
  
  dat_temp <- rfmo %>% 
    filter(year == layer)
  
  if(nrow(dat_temp) == 0) { next }
  
  dat_temp <- dat_temp %>% 
    mutate(longitude_orig = longitude, 
           latitude_orig = latitude) %>% 
    st_as_sf(., coords = c("longitude_orig", "latitude_orig"), crs = 4326)
  
  raster_stack <- stack()
  
  raster_1 <- dat_temp %>% 
    filter(latitude%%1 == 0 & longitude%%1 == 0) 
  
  if(nrow(raster_1) > 0) { 
    raster_1 <- raster_1 %>% 
      rasterize(., whole_numbers, field = "target_effort", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_1)
  }
  
  raster_2 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 == 0) 
  
  if(nrow(raster_2) > 0) { 
    raster_2 <- raster_2 %>% 
      rasterize(., whole_numbers, field = "target_effort", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_2)
  }
  
  raster_3 <- dat_temp %>% 
    filter(latitude%%1 == 0 & longitude%%1 != 0)  
  
  if(nrow(raster_3) > 0) { 
    raster_3 <- raster_3 %>% 
      rasterize(., whole_numbers, field = "target_effort", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_3)
  }
  
  raster_4 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 != 0)  
  
  if(nrow(raster_4) > 0) { 
    raster_4 <- raster_4 %>% 
      rasterize(., whole_numbers, field = "target_effort", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_4)
  }
  
  if(nlayers(raster_stack) > 1) { 
    raster_comb <- calc(raster_stack, mean, na.rm = T) 
  } else { 
    raster_comb <- raster_stack }
  
  #raster_comb <- log(raster_comb)
  raster_comb <- raster::as.data.frame(raster_comb, xy = TRUE) 
  
  # Some species are all 0s
  if(min(raster_comb$layer, na.rm = T) == mean(raster_comb$layer, na.rm = T)) {next}
  
  fig_supp <- ggplot() + 
    geom_tile(raster_comb, 
              mapping = aes(x=x, y=y, fill=layer, color = layer)) + 
    scale_fill_distiller("Fishing Effort (hooks)", palette = "RdYlBu", na.value = NA, 
                         labels = scales::label_comma(),
                         guide = guide_colorbar(title.vjust = 0.8)) + 
    scale_color_distiller("", palette = "RdYlBu",  
                          na.value = NA, guide = "none") + 
    geom_sf(data = wcpfc_boundary, fill = NA, color = "black") +
    geom_sf(data = iotc_boundary, fill = NA, color = "black") +
    geom_sf(data = iccat_boundary, fill = NA, color = "black") +
    geom_sf(data = iattc_boundary, fill = NA, color = "black") +
    geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
              mapping = aes(x=x, y=y), fill = "black", color = "black") +
    coord_sf() + 
    custom_theme + 
    theme(legend.position = "bottom", 
          text = element_text(size = 18), 
          legend.text = element_text(angle = -90, hjust = 0.5, vjust = 0.5)) 
  
  if(layer == 2012) { 
    legend <- get_legend(fig_supp)}
  
  #fig_supp <- fig_supp + theme(legend.position = "none")
  
  assign(paste0("rfmo_mean_", layer), fig_supp)
} 

rfmo_mean_final_plot <- ggdraw() + 
  draw_plot(rfmo_mean_2012, 0, 0.66, 0.33, 0.33) + 
  draw_plot(rfmo_mean_2013, 0.33, 0.66, 0.33, 0.33) +  
  draw_plot(rfmo_mean_2014, 0.66, 0.66, 0.33, 0.33) +  
  draw_plot(rfmo_mean_2015, 0, 0.33, 0.33, 0.33) + 
  draw_plot(rfmo_mean_2016, 0.33, 0.33, 0.33, 0.33) +  
  draw_plot(rfmo_mean_2017, 0.66, 0.33, 0.33, 0.33) +  
  draw_plot(rfmo_mean_2018, 0, 0, 0.33, 0.33) + 
  draw_plot(rfmo_mean_2019, 0.33, 0, 0.33, 0.33) + 
  draw_plot(rfmo_mean_2020, 0.66, 0, 0.33, 0.33) + 
  #draw_plot(legend, 0, 0, 1, 0.1) + 
  draw_plot_label(label = paste0(LETTERS[1:9], ") ", 2012:2020), x = rep(c(0, 0.33, 0.66), 3), 
                  y = rep(c(1, 0.66, 0.33), each = 3))

ggsave(file.path(here::here(), "figures/supplemental/mean_rfmo_target_effort.png"), plot = rfmo_mean_final_plot,
       bg = "white", dpi = 600, width = 20, height = 12)

# Bycatch total effort
for(layer in unique(rfmo$year)) { 
  
  dat_temp <- rfmo %>% 
    filter(year == layer)
  
  if(nrow(dat_temp) == 0) { next }
  
  dat_temp <- dat_temp %>% 
    mutate(longitude_orig = longitude, 
           latitude_orig = latitude) %>% 
    st_as_sf(., coords = c("longitude_orig", "latitude_orig"), crs = 4326)
  
  raster_stack <- stack()
  
  raster_1 <- dat_temp %>% 
    filter(latitude%%1 == 0 & longitude%%1 == 0) 
  
  if(nrow(raster_1) > 0) { 
    raster_1 <- raster_1 %>% 
      rasterize(., whole_numbers, field = "bycatch_total_effort", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_1)
  }
  
  raster_2 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 == 0) 
  
  if(nrow(raster_2) > 0) { 
    raster_2 <- raster_2 %>% 
      rasterize(., whole_numbers, field = "bycatch_total_effort", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_2)
  }
  
  raster_3 <- dat_temp %>% 
    filter(latitude%%1 == 0 & longitude%%1 != 0)  
  
  if(nrow(raster_3) > 0) { 
    raster_3 <- raster_3 %>% 
      rasterize(., whole_numbers, field = "bycatch_total_effort", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_3)
  }
  
  raster_4 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 != 0)  
  
  if(nrow(raster_4) > 0) { 
    raster_4 <- raster_4 %>% 
      rasterize(., whole_numbers, field = "bycatch_total_effort", fun = mean, background = NA)
    
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
    scale_fill_distiller("Fishing Effort (hooks)", palette = "RdYlBu", na.value = NA, 
                         labels = scales::label_comma(),
                         guide = guide_colorbar(title.vjust = 0.8)) + 
    scale_color_distiller("", palette = "RdYlBu",  
                          na.value = NA, guide = "none") + 
    geom_sf(data = wcpfc_boundary, fill = NA, color = "black") +
    geom_sf(data = iotc_boundary, fill = NA, color = "black") +
    geom_sf(data = iccat_boundary, fill = NA, color = "black") +
    geom_sf(data = iattc_boundary, fill = NA, color = "black") +
    geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
              mapping = aes(x=x, y=y), fill = "black", color = "black") +
    coord_sf() + 
    custom_theme + 
    theme(legend.position = "bottom", 
          text = element_text(size = 18), 
          legend.text = element_text(angle = -90, hjust = 0.5, vjust = 0.5)) 
  
  if(layer == 2012) { 
    legend <- get_legend(fig_supp)}
  
  #fig_supp <- fig_supp + theme(legend.position = "none")
  
  assign(paste0("rfmo_mean_", layer), fig_supp)
} 

rfmo_mean_final_plot <- ggdraw() + 
  draw_plot(rfmo_mean_2012, 0, 0.66, 0.33, 0.33) + 
  draw_plot(rfmo_mean_2013, 0.33, 0.66, 0.33, 0.33) +  
  draw_plot(rfmo_mean_2014, 0.66, 0.66, 0.33, 0.33) +  
  draw_plot(rfmo_mean_2015, 0, 0.33, 0.33, 0.33) + 
  draw_plot(rfmo_mean_2016, 0.33, 0.33, 0.33, 0.33) +  
  draw_plot(rfmo_mean_2017, 0.66, 0.33, 0.33, 0.33) +  
  draw_plot(rfmo_mean_2018, 0, 0, 0.33, 0.33) + 
  draw_plot(rfmo_mean_2019, 0.33, 0, 0.33, 0.33) + 
  draw_plot(rfmo_mean_2020, 0.66, 0, 0.33, 0.33) + 
  #draw_plot(legend, 0, 0, 1, 0.1) + 
  draw_plot_label(label = paste0(LETTERS[1:9], ") ", 2012:2020), x = rep(c(0, 0.33, 0.66), 3), 
                  y = rep(c(1, 0.66, 0.33), each = 3))

ggsave(file.path(here::here(), "figures/supplemental/mean_rfmo_bycatch_total_effort.png"), plot = rfmo_mean_final_plot,
       bg = "white", dpi = 600, width = 20, height = 12)

