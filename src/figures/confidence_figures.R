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
file_list <- list.files(file.path(here::here(), "data-updated/model-data/outputs/all-rfmo-models/quantiles"), 
                         pattern = ".rds", full.names = TRUE)

combine_files <- function(list_files, filter = NULL) { 
  all_dat <- NULL
  for(file in list_files) { 
  
    temp_rds <- readRDS(file)
    
    res <- unique(temp_rds$final_predict$spatial_notes)
    
    temp <- temp_rds$cell_class_confidence %>% 
      mutate(spatial_notes = res) %>% 
      separate(indID, sep = "[|]", into = c("latitude", "longitude")) %>% 
      mutate(latitude = as.numeric(latitude), 
             longitude = as.numeric(longitude)) %>% 
      mutate(lat = latitude, 
             lon = longitude)
    
    if(!is.null(filter)) { 
      if(filter == 1) { 
        temp <- temp %>% 
          filter(.pred_class == 1)
      } else if(filter == 0) { 
        temp <- temp %>% 
          filter(.pred_class == 0)
      }
    }
  
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
  return(all_dat)
} 

spatial_confidence_plots <- function(file, scale = NULL) { 

  all_dat <- file %>% 
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
  
  plot <- ggplot() + 
    geom_tile(raster_comb_df, 
              mapping = aes(x=x, y=y, fill=layer, color = layer)) + 
    geom_sf(data = wcpfc_boundary, fill = NA, color = "black") +
    geom_sf(data = iotc_boundary, fill = NA, color = "black") +
    geom_sf(data = iccat_boundary, fill = NA, color = "black") +
    geom_sf(data = iattc_boundary, fill = NA, color = "black") +
    geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
              mapping = aes(x=x, y=y), fill = "black", color = "black") +
    coord_sf() + 
    custom_theme + 
    theme(legend.position = "bottom") 
  
  if(is.null(scale)) { 
    plot <- plot + 
      scale_fill_distiller("Confidence in Classification\nPrediction", palette = "RdYlBu", na.value = NA, 
                           breaks = c(floor(min(raster_comb_df$layer, na.rm = T)/0.01)*0.01, ceiling(max(raster_comb_df$layer, na.rm = T)/0.01)*0.01), 
                           limits = c(floor(min(raster_comb_df$layer, na.rm = T)/0.01)*0.01, ceiling(max(raster_comb_df$layer, na.rm = T)/0.01)*0.01),
                           guide = guide_colorbar(title.vjust = 0.8)) + 
      scale_color_distiller("", palette = "RdYlBu", na.value = NA, 
                            breaks = c(floor(min(raster_comb_df$layer, na.rm = T)/0.01)*0.01, ceiling(max(raster_comb_df$layer, na.rm = T)/0.01)*0.01),
                            limits = c(floor(min(raster_comb_df$layer, na.rm = T)/0.01)*0.01, ceiling(max(raster_comb_df$layer, na.rm = T)/0.01)*0.01),
                            guide = "none") 
  } else { 
      plot <- plot + 
        scale_fill_distiller("Confidence in Classification\nPrediction", palette = "RdYlBu", na.value = NA, 
                             breaks = c((scale[1]/0.01)*0.01, (scale[2]/0.01)*0.01), 
                             limits = c((scale[1]/0.01)*0.01, (scale[2]/0.01)*0.01),
                             guide = guide_colorbar(title.vjust = 0.8)) + 
        scale_color_distiller("", palette = "RdYlBu", na.value = NA, 
                              breaks = c((scale[1]/0.01)*0.01, (scale[2]/0.01)*0.01), 
                              limits = c((scale[1]/0.01)*0.01, (scale[2]/0.01)*0.01),
                              guide = "none") 
      }
  
  return(plot)
} 

# Save plots for all classifications, just sharks present, and just sharks absent
# Play around with different scales 

all_classes <- combine_files(file_list)
plot1a <- spatial_confidence_plots(all_classes)
plot1b <- spatial_confidence_plots(all_classes, scale = c(0,1))

pres_classes <- combine_files(file_list, filter = 1)
plot2a <- spatial_confidence_plots(pres_classes)
plot2b <- spatial_confidence_plots(pres_classes, scale = c(0,1))

abs_classes <- combine_files(file_list, filter = 0)
plot3a <- spatial_confidence_plots(abs_classes)
plot3b <- spatial_confidence_plots(abs_classes, scale = c(0,1))

final <- ggdraw() + 
  draw_plot(plot1a, 0, 0.5, 0.33, 0.5) + 
  draw_plot(plot2a, 0.33, 0.5, 0.33, 0.5) + 
  draw_plot(plot3a, 0.66, 0.5, 0.33, 0.5) + 
  draw_plot(plot1b, 0, 0, 0.33, 0.5) + 
  draw_plot(plot2b, 0.33, 0, 0.33, 0.5) + 
  draw_plot(plot3b, 0.66, 0, 0.33, 0.5) + 
  draw_plot_label(label = c("All Classes", "Shark Presence Only", "Shark Absence Only"), 
                  x = c(0+0.165, 0.33+0.165, 0.66+0.165), 
                  y = 1, hjust = 0.5)

ggsave("/Users/echelleburns/Desktop/longline_classification_confidence.png", final,
               height = 6, width = 12, units = "in", dpi = 600, bg = "white")

# Save plots for distribution of predictions
predictions <- purrr::map_df(file_list, 
                          ~{readRDS(.x)$id_predictions})

beta_dist <- purrr::map_df(file_list, 
                           ~{readRDS(.x)$beta_parameters})
  
predictions_cell <- predictions %>% 
  filter(run == "cell")

beta_dist_cell <- beta_dist %>% 
  filter(run == "cell")

plot_distributions <- function(predictions_data, beta_data) { 
  beta_plot <- ggplot(data = predictions_data %>% 
                   filter(indID %in% unique(beta_data$indID))) + 
    geom_histogram(mapping = aes(x = predictions, y = ..density..), binwidth = 0.1, 
                   fill = "navy", alpha = 0.7, color = "white") + 
    scale_y_continuous(name = "Density", limits = c(0, 5)) + 
    xlab("Classification Prediction") + 
    theme_classic()
  
  for(i in 1:nrow(beta_data)){
    beta_plot <- beta_plot +
    stat_function(fun = dbeta, args = c(beta_data$shape1[i], 
                                        beta_data$shape2[i]),
                  alpha = 0.1)
  }
  return(beta_plot)
} 

# Where mean prediction for that ID is 1
beta_dist_cell_false <- beta_dist_cell %>% 
  filter(lower_tail == FALSE) 

# Where mean prediction for that ID is 0
beta_dist_cell_true <- beta_dist_cell %>% 
  filter(lower_tail == TRUE) 

beta_1 <- plot_distributions(predictions_cell, beta_dist_cell_false)
beta_0 <- plot_distributions(predictions_cell, beta_dist_cell_true[sample(1:nrow(beta_dist_cell_true), 600),])

final <- ggdraw() + 
  draw_plot(beta_1, 0, 0, 0.5, 1) + 
  draw_plot(beta_0, 0.5, 0, 0.5, 1) +
  draw_plot_label(label = c("Predicted Shark Presence", "Predicted Shark Absence"), 
                  x = c(.25, 0.75), 
                  y = 1, hjust = 0.5)

ggsave("/Users/echelleburns/Desktop/longline_classification_distributions.png", final,
       height = 6, width = 12, units = "in", dpi = 600, bg = "white")



