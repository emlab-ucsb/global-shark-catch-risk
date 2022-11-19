# Figures for environmentals used in modeling (sst, chl-a, ssh)

# Load libraries
library(raster)
library(tidyverse)
library(cowplot)
library(sf)
library(here)

# Load plotting defaults
source(file.path(here::here(), "src/figures/plot_defaults.R"))

# Load data
sst <- read.csv(file.path(here::here(), "data-updated/sea-surface-temperature/outputs/binned_global_sst_1x1.csv"))
chla <- read.csv(file.path(here::here(), "data-updated/chlorophyll-a/outputs/binned_global_chla_1x1.csv"))
ssh <- read.csv(file.path(here::here(), "data-updated/sea-surface-height/outputs/binned_global_ssh_1x1.csv"))

# Paneled plots for each year using mean and cv

## SST - Mean

max_mean <- 0
min_mean <- 100
mean_list <- list()

for(layer in unique(sst$year)) { 
  
  dat_temp <- sst %>% 
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
      rasterize(., whole_numbers, field = "mean_sst", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_1)
  }
  
  raster_2 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 == 0) 
  
  if(nrow(raster_2) > 0) { 
    raster_2 <- raster_2 %>% 
      rasterize(., whole_numbers, field = "mean_sst", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_2)
  }
  
  raster_3 <- dat_temp %>% 
    filter(latitude%%1 == 0 & longitude%%1 != 0)  
  
  if(nrow(raster_3) > 0) { 
    raster_3 <- raster_3 %>% 
      rasterize(., whole_numbers, field = "mean_sst", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_3)
  }
  
  raster_4 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 != 0)  
  
  if(nrow(raster_4) > 0) { 
    raster_4 <- raster_4 %>% 
      rasterize(., whole_numbers, field = "mean_sst", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_4)
  }
  
  if(nlayers(raster_stack) > 1) { 
    raster_comb <- calc(raster_stack, mean, na.rm = T) 
  } else { 
    raster_comb <- raster_stack }
  
  raster_comb <- raster::as.data.frame(raster_comb, xy = TRUE) 
  
  # Adjust max/min
  if(max(raster_comb$layer, na.rm = T) > max_mean) { 
    max_mean <- max(raster_comb$layer, na.rm = T)}
  if(min(raster_comb$layer, na.rm = T) < min_mean) { 
    min_mean <- min(raster_comb$layer, na.rm = T)}
  
  mean_list <- append(mean_list, list(raster_comb))
} 
  
names(mean_list) <- 2012:2020
  
for(layer in names(mean_list)) { 
  # Some species are all 0s
  if(min(mean_list[[layer]]$layer, na.rm = T) == mean(mean_list[[layer]]$layer, na.rm = T)) {next}
  
  fig_supp <- ggplot() + 
    geom_tile(mean_list[[layer]], 
              mapping = aes(x=x, y=y, fill=layer, color = layer)) + 
    scale_fill_distiller("Mean SST (°C)", palette = "RdYlBu", na.value = NA, 
                         limits = c(0, max(sst$mean_sst, na.rm = T)),  
                         labels = scales::label_comma(accuracy = 1), 
                         guide = guide_colorbar(title.vjust = 0.8)) + 
    scale_color_distiller("", palette = "RdYlBu", limits = c(0, max(sst$mean_sst, na.rm = T)), na.value = NA, guide = "none") + 
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
  
  if(layer == 2012) { 
    legend <- get_legend(fig_supp)}
  
  fig_supp <- fig_supp + theme(legend.position = "none")
  
  assign(paste0("sst_mean_", layer), fig_supp)
} 

sst_mean_final_plot <- ggdraw() + 
  draw_plot(sst_mean_2012, 0, 0.7, 0.33, 0.3) + 
  draw_plot(sst_mean_2013, 0.33, 0.7, 0.33, 0.3) +  
  draw_plot(sst_mean_2014, 0.66, 0.7, 0.33, 0.3) +  
  draw_plot(sst_mean_2015, 0, 0.4, 0.33, 0.3) + 
  draw_plot(sst_mean_2016, 0.33, 0.4, 0.33, 0.3) +  
  draw_plot(sst_mean_2017, 0.66, 0.4, 0.33, 0.3) +  
  draw_plot(sst_mean_2018, 0, 0.1, 0.33, 0.3) + 
  draw_plot(sst_mean_2019, 0.33, 0.1, 0.33, 0.3) + 
  draw_plot(sst_mean_2020, 0.66, 0.1, 0.33, 0.3) + 
  draw_plot(legend, 0, 0, 1, 0.1) + 
  draw_plot_label(label = paste0(LETTERS[1:9], ") ", 2012:2020), x = rep(c(0, 0.33, 0.66), 3), 
                  y = rep(c(1, 0.7, 0.4), each = 3))

ggsave(file.path(here::here(), "figures/supplemental/mean_sst.png"), plot = sst_mean_final_plot,
       bg = "white", dpi = 600, width = 20, height = 12)

## Chl-a - Mean 

max_mean <- 0
min_mean <- 100
mean_list <- list()

for(layer in unique(chla$year)) { 
  
  dat_temp <- chla %>% 
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
      rasterize(., whole_numbers, field = "mean_chla", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_1)
  }
  
  raster_2 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 == 0) 
  
  if(nrow(raster_2) > 0) { 
    raster_2 <- raster_2 %>% 
      rasterize(., whole_numbers, field = "mean_chla", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_2)
  }
  
  raster_3 <- dat_temp %>% 
    filter(latitude%%1 == 0 & longitude%%1 != 0)  
  
  if(nrow(raster_3) > 0) { 
    raster_3 <- raster_3 %>% 
      rasterize(., whole_numbers, field = "mean_chla", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_3)
  }
  
  raster_4 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 != 0)  
  
  if(nrow(raster_4) > 0) { 
    raster_4 <- raster_4 %>% 
      rasterize(., whole_numbers, field = "mean_chla", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_4)
  }
  
  if(nlayers(raster_stack) > 1) { 
    raster_comb <- calc(raster_stack, mean, na.rm = T) 
  } else { 
    raster_comb <- raster_stack }
  
  raster_comb <- log(raster_comb)
  raster_comb <- raster::as.data.frame(raster_comb, xy = TRUE) 
  
  # Adjust max/min
  if(max(raster_comb$layer, na.rm = T) > max_mean) { 
    max_mean <- max(raster_comb$layer, na.rm = T)}
  if(min(raster_comb$layer, na.rm = T) < min_mean) { 
    min_mean <- min(raster_comb$layer, na.rm = T)}
  
  mean_list <- append(mean_list, list(raster_comb))
} 

names(mean_list) <- 2012:2020

for(layer in names(mean_list)) {
  # Some species are all 0s
  if(min(mean_list[[layer]]$layer, na.rm = T) == mean(mean_list[[layer]]$layer, na.rm = T)) {next}
  
  fig_supp <- ggplot() + 
    geom_tile(mean_list[[layer]], 
              mapping = aes(x=x, y=y, fill=layer, color = layer)) + 
    scale_fill_distiller("Mean Chl-A (milligram m-3)", palette = "RdYlBu", na.value = NA, 
                         limits = c(min_mean, max_mean), 
                         labels = scales::label_comma(accuracy = 1), 
                         guide = guide_colorbar(title.vjust = 0.8)) + 
    scale_color_distiller("", palette = "RdYlBu", 
                          limits = c(min_mean, max_mean),
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
          text = element_text(size = 18)) 
  
  if(layer == 2012) { 
    legend <- get_legend(fig_supp)}
  
  fig_supp <- fig_supp + theme(legend.position = "none")
  
  assign(paste0("chla_mean_", layer), fig_supp)
} 

chla_mean_final_plot <- ggdraw() + 
  draw_plot(chla_mean_2012, 0, 0.7, 0.33, 0.3) + 
  draw_plot(chla_mean_2013, 0.33, 0.7, 0.33, 0.3) +  
  draw_plot(chla_mean_2014, 0.66, 0.7, 0.33, 0.3) +  
  draw_plot(chla_mean_2015, 0, 0.4, 0.33, 0.3) + 
  draw_plot(chla_mean_2016, 0.33, 0.4, 0.33, 0.3) +  
  draw_plot(chla_mean_2017, 0.66, 0.4, 0.33, 0.3) +  
  draw_plot(chla_mean_2018, 0, 0.1, 0.33, 0.3) + 
  draw_plot(chla_mean_2019, 0.33, 0.1, 0.33, 0.3) + 
  draw_plot(chla_mean_2020, 0.66, 0.1, 0.33, 0.3) + 
  draw_plot(legend, 0, 0, 1, 0.1) + 
  draw_plot_label(label = paste0(LETTERS[1:9], ") ", 2012:2020), x = rep(c(0, 0.33, 0.66), 3), 
                  y = rep(c(1, 0.7, 0.4), each = 3))

ggsave(file.path(here::here(), "figures/supplemental/mean_chla.png"), plot = chla_mean_final_plot,
       bg = "white", dpi = 600, width = 20, height = 12)

## SSH - Mean 

max_mean <- 0
min_mean <- 100
mean_list <- list()

for(layer in unique(ssh$year)) { 
  
  dat_temp <- ssh %>% 
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
      rasterize(., whole_numbers, field = "mean_ssh", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_1)
  }
  
  raster_2 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 == 0) 
  
  if(nrow(raster_2) > 0) { 
    raster_2 <- raster_2 %>% 
      rasterize(., whole_numbers, field = "mean_ssh", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_2)
  }
  
  raster_3 <- dat_temp %>% 
    filter(latitude%%1 == 0 & longitude%%1 != 0)  
  
  if(nrow(raster_3) > 0) { 
    raster_3 <- raster_3 %>% 
      rasterize(., whole_numbers, field = "mean_ssh", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_3)
  }
  
  raster_4 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 != 0)  
  
  if(nrow(raster_4) > 0) { 
    raster_4 <- raster_4 %>% 
      rasterize(., whole_numbers, field = "mean_ssh", fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_4)
  }
  
  if(nlayers(raster_stack) > 1) { 
    raster_comb <- calc(raster_stack, mean, na.rm = T) 
  } else { 
    raster_comb <- raster_stack }
  
  raster_comb <- raster::as.data.frame(raster_comb, xy = TRUE) 
  
  # Adjust max/min
  if(max(raster_comb$layer, na.rm = T) > max_mean) { 
    max_mean <- max(raster_comb$layer, na.rm = T)}
  if(min(raster_comb$layer, na.rm = T) < min_mean) { 
    min_mean <- min(raster_comb$layer, na.rm = T)}
  
  mean_list <- append(mean_list, list(raster_comb))
} 

names(mean_list) <- 2012:2019

for(layer in names(mean_list)) {
  
  # Some species are all 0s
  if(min(mean_list[[layer]]$layer, na.rm = T) == mean(mean_list[[layer]]$layer, na.rm = T)) {next}
  
  fig_supp <- ggplot() + 
    geom_tile(mean_list[[layer]], 
              mapping = aes(x=x, y=y, fill=layer, color = layer)) + 
    scale_fill_distiller("Mean SSH (m)", palette = "RdYlBu", na.value = NA, 
                         limits = c(min_mean, max_mean), 
                         labels = scales::label_comma(accuracy = 0.01),
                         guide = guide_colorbar(title.vjust = 0.8)) + 
    scale_color_distiller("", palette = "RdYlBu", limits = c(min_mean, max_mean), na.value = NA, guide = "none") + 
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
          legend.text = element_text(angle = -45, hjust = 0.5, vjust = 0.5)) 
  
  if(layer == 2012) { 
    legend <- get_legend(fig_supp)}
  
  fig_supp <- fig_supp + theme(legend.position = "none")
  
  assign(paste0("ssh_mean_", layer), fig_supp)
} 

ssh_mean_final_plot <- ggdraw() + 
  draw_plot(ssh_mean_2012, 0, 0.7, 0.33, 0.3) + 
  draw_plot(ssh_mean_2013, 0.33, 0.7, 0.33, 0.3) +  
  draw_plot(ssh_mean_2014, 0.66, 0.7, 0.33, 0.3) +  
  draw_plot(ssh_mean_2015, 0, 0.4, 0.33, 0.3) + 
  draw_plot(ssh_mean_2016, 0.33, 0.4, 0.33, 0.3) +  
  draw_plot(ssh_mean_2017, 0.66, 0.4, 0.33, 0.3) +  
  draw_plot(ssh_mean_2018, 0, 0.1, 0.33, 0.3) + 
  draw_plot(ssh_mean_2019, 0.33, 0.1, 0.33, 0.3) + 
  draw_plot(legend, 0, 0, 1, 0.1) + 
  draw_plot_label(label = paste0(LETTERS[1:8], ") ", 2012:2019), x = c(0, 0.33, 0.66, 0, 0.33, 0.66, 0, 0.33), 
                  y = c(1, 1, 1, 0.7, 0.7, 0.7, 0.4, 0.4))

ggsave(file.path(here::here(), "figures/supplemental/mean_ssh.png"), plot = ssh_mean_final_plot,
       bg = "white", dpi = 600, width = 20, height = 12)
  