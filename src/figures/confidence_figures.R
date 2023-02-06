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

combine_files <- function(list_files, dataset, filter = NULL) { 
  all_dat <- NULL
  for(file in list_files) { 
  
    temp_rds <- readRDS(file)
    
    res <- unique(temp_rds$final_predict$spatial_notes)
    
    temp <- temp_rds[[dataset]] %>% 
      mutate(spatial_notes = res) 
    
    if(dataset == "cell_class_confidence") { 
      temp <- temp %>% 
        separate(indID, sep = "[|]", into = c("latitude", "longitude"))
    }
    
    temp <- temp %>% 
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
        mutate(spatial_notes = "center of 1x1 cell") 
      
      if(dataset == "cell_class_confidence") { 
        temp <- temp %>%
          group_by(latitude_rescaled, longitude_rescaled, spatial_notes) %>% 
          summarise(confidence = mean(confidence, na.rm = T), 
                    .pred_class = mean(.pred_class, na.rm = T)) %>% 
          ungroup() 
      } else if (dataset == "final_predict") { 
        temp <- temp %>%
          group_by(year, species_commonname, species_sciname, latitude_rescaled, longitude_rescaled, spatial_notes) %>% 
          summarise(.pred_upper = .pred_upper/25,
                    .pred_lower = .pred_lower/25, 
                    .final_pred = .final_pred/25) %>% 
          ungroup() 
      } 
      
      temp <- temp %>% 
        rename(latitude = latitude_rescaled, 
               longitude = longitude_rescaled) 
    } 
    
    if(dataset == "cell_class_confidence") { 
      temp <- temp %>% 
        select(latitude, longitude, confidence, .pred_class)
    } else if(dataset == "final_predict") { 
      temp <- temp %>% 
        select(year, latitude, longitude, species_commonname, species_sciname,
               .pred_upper, .final_pred, .pred_lower)
    }
    
  all_dat <- all_dat %>% 
    rbind(temp %>%
            mutate(lat = latitude,
                   lon = longitude))
  } 
  return(all_dat)
} 

spatial_confidence_plots <- function(file, field = "confidence", legendtitle = "Confidence in Classification\nPrediction", scale = NULL) { 

  all_dat <- file %>% 
    st_as_sf(., coords = c("lon", "lat"), crs = 4326)
  
  raster_stack <- stack()
  
  raster_1 <- all_dat %>% 
    filter(latitude%%1 == 0 & longitude%%1 == 0) 
  
  if(nrow(raster_1) > 0) { 
    raster_1 <- raster_1 %>% 
      rasterize(., whole_numbers, field = field, fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_1)
  }
  
  raster_2 <- all_dat %>% 
    filter(latitude%%1 != 0 & longitude%%1 == 0) 
  
  if(nrow(raster_2) > 0) { 
    raster_2 <- raster_2 %>% 
      rasterize(., whole_numbers, field = field, fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_2)
  }
  
  raster_3 <- all_dat %>% 
    filter(latitude%%1 == 0 & longitude%%1 != 0)  
  
  if(nrow(raster_3) > 0) { 
    raster_3 <- raster_3 %>% 
      rasterize(., whole_numbers, field = field, fun = mean, background = NA)
    
    raster_stack <- stack(raster_stack, raster_3)
  }
  
  raster_4 <- all_dat %>% 
    filter(latitude%%1 != 0 & longitude%%1 != 0)  
  
  if(nrow(raster_4) > 0) { 
    raster_4 <- raster_4 %>% 
      rasterize(., whole_numbers, field = field, fun = mean, background = NA)
    
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
      scale_fill_distiller(legendtitle, palette = "RdYlBu", na.value = NA, 
                           breaks = c(floor(min(raster_comb_df$layer, na.rm = T)/0.01)*0.01, ceiling(max(raster_comb_df$layer, na.rm = T)/0.01)*0.01), 
                           limits = c(floor(min(raster_comb_df$layer, na.rm = T)/0.01)*0.01, ceiling(max(raster_comb_df$layer, na.rm = T)/0.01)*0.01),
                           guide = guide_colorbar(title.vjust = 0.8)) + 
      scale_color_distiller("", palette = "RdYlBu", na.value = NA, 
                            breaks = c(floor(min(raster_comb_df$layer, na.rm = T)/0.01)*0.01, ceiling(max(raster_comb_df$layer, na.rm = T)/0.01)*0.01),
                            limits = c(floor(min(raster_comb_df$layer, na.rm = T)/0.01)*0.01, ceiling(max(raster_comb_df$layer, na.rm = T)/0.01)*0.01),
                            guide = "none") 
  } else { 
      plot <- plot + 
        scale_fill_distiller(legendtitle, palette = "RdYlBu", na.value = NA, 
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

### 
# Individual Sharks
### 

# Save plots for all classifications, just sharks present, and just sharks absent
# Play around with different scales 

all_classes <- combine_files(file_list, dataset = "cell_class_confidence")
plot1a <- spatial_confidence_plots(all_classes)
plot1b <- spatial_confidence_plots(all_classes, scale = c(0,1))

pres_classes <- combine_files(file_list, dataset = "cell_class_confidence", filter = 1)
plot2a <- spatial_confidence_plots(pres_classes)
plot2b <- spatial_confidence_plots(pres_classes, scale = c(0,1))

abs_classes <- combine_files(file_list, dataset = "cell_class_confidence", filter = 0)
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

ggsave(file.path(here::here(), "figures", "supplemental", "longline_classification_confidence.png"), 
       final, height = 6, width = 12, units = "in", dpi = 600, bg = "white")

# Save plots for distribution of predictions
predictions <- purrr::map_df(file_list, 
                          ~{readRDS(.x)$id_predictions})

beta_dist <- purrr::map_df(file_list, 
                           ~{readRDS(.x)$beta_parameters})
  
## By cell
predictions_cell <- predictions %>% 
  filter(run == "cell")

beta_dist_cell <- beta_dist %>% 
  filter(run == "cell")

# Where mean prediction for that ID is 1
beta_dist_cell_false <- beta_dist_cell %>% 
  filter(lower_tail == FALSE) 

# Where mean prediction for that ID is 0
beta_dist_cell_true <- beta_dist_cell %>% 
  filter(lower_tail == TRUE) 

if(nrow(beta_dist_cell_false) > 0 & nrow(beta_dist_cell_true) > 0) { 
  if(nrow(beta_dist_cell_false) > 600) { 
    beta_dist_cell_false <- beta_dist_cell_false[sample(1:nrow(beta_dist_cell_false), 600),]
  }
  if(nrow(beta_dist_cell_true) > 600) { 
    beta_dist_cell_true <- beta_dist_cell_true[sample(1:nrow(beta_dist_cell_true), 600),]
  }
  
  beta_1 <- plot_distributions(predictions_cell, beta_dist_cell_false)
  beta_0 <- plot_distributions(predictions_cell, beta_dist_cell_true)
  
  final <- ggdraw() + 
    draw_plot(beta_1, 0, 0, 0.5, 1) + 
    draw_plot(beta_0, 0.5, 0, 0.5, 1) +
    draw_plot_label(label = c("Predicted Shark Presence", "Predicted Shark Absence"), 
                    x = c(.25, 0.75), 
                    y = 1, hjust = 0.5)
  
  ggsave(file.path(here::here(), "figures", "supplemental", "longline_classification_distributions_cell.png"),
         final, height = 6, width = 12, units = "in", dpi = 600, bg = "white")
} 

## By species
predictions_spp <- predictions %>% 
  filter(run == "spp")

beta_dist_spp <- beta_dist %>% 
  filter(run == "spp")

# Where mean prediction for that ID is 1
beta_dist_spp_false <- beta_dist_spp %>% 
  filter(lower_tail == FALSE) 

# Where mean prediction for that ID is 0
beta_dist_spp_true <- beta_dist_spp %>% 
  filter(lower_tail == TRUE) 

if(nrow(beta_dist_spp_false) > 0 & nrow(beta_dist_spp_true) > 0) { 
  if(nrow(beta_dist_spp_false) > 600) { 
    beta_dist_spp_false <- beta_dist_spp_false[sample(1:nrow(beta_dist_spp_false), 600),]
  }
  if(nrow(beta_dist_spp_true) > 600) { 
    beta_dist_spp_true <- beta_dist_spp_true[sample(1:nrow(beta_dist_spp_true), 600),]
  }
  
  beta_1 <- plot_distributions(predictions_spp, beta_dist_spp_false)
  beta_0 <- plot_distributions(predictions_spp, beta_dist_spp_true)
  
  final <- ggdraw() + 
    draw_plot(beta_1, 0, 0, 0.5, 1) + 
    draw_plot(beta_0, 0.5, 0, 0.5, 1) +
    draw_plot_label(label = c("Predicted Shark Presence", "Predicted Shark Absence"), 
                    x = c(.25, 0.75), 
                    y = 1, hjust = 0.5)
  
  ggsave(file.path(here::here(), "figures", "supplemental", "longline_classification_distributions_spp.png"),
         final, height = 6, width = 12, units = "in", dpi = 600, bg = "white")
} 

# High and low values from the regression
all_pred <- combine_files(file_list, dataset = "final_predict") %>% 
  mutate(.pred_lower = ifelse(.final_pred == 0, 0, .pred_lower), 
         .pred_upper = ifelse(.final_pred == 0, 0, .pred_upper)) %>% 
  group_by(latitude, longitude, year, lat, lon) %>% 
  summarise(.final_pred = sum(.final_pred, na.rm = T), 
            .pred_lower = sum(.pred_lower, na.rm = T), 
            .pred_upper = sum(.pred_upper, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(latitude, longitude, lat, lon) %>% 
  summarise(.final_pred = mean(.final_pred, na.rm = T), 
            .pred_lower = mean(.pred_lower, na.rm = T), 
            .pred_upper = mean(.pred_upper, na.rm = T))

plot1a <- spatial_confidence_plots(all_pred, field = ".pred_lower", legendtitle = "Predicted catch (lower bound)")
plot1b <- spatial_confidence_plots(all_pred, field = ".final_pred", legendtitle = "Predicted catch")
plot1c <- spatial_confidence_plots(all_pred, field = ".pred_upper", legendtitle = "Predicted catch (upper bound)")

final <- ggdraw() + 
  draw_plot(plot1a, 0, 0, 0.33, 1) + 
  draw_plot(plot1b, 0.33, 0, 0.33, 1) + 
  draw_plot(plot1c, 0.66, 0, 0.33, 1) + 
  draw_plot_label(label = c("Lower Bound", "Predicted Catch", "Upper Bound"), 
                  x = c(0+0.165, 0.33+0.165, 0.66+0.165), 
                  y = 1, hjust = 0.5)

ggsave(file.path(here::here(), "figures", "supplemental", "longline_regression.png"), 
       final, height = 3, width = 12, units = "in", dpi = 600, bg = "white")

# Save tables for outputs as well
file_list <- list.files(file.path(here::here(), "data-updated", "model-data", "outputs", "all-rfmo-models", "quantiles"), 
                        pattern = ".rds", full.names = TRUE)

spp_conf <- NULL
for(file in file_list) { 
  temp_rds <- readRDS(file)$species_class_confidence
  spp_conf <- bind_rows(spp_conf, 
                        temp_rds %>%
                          mutate(rfmo = gsub("[^IATTC|IOTC|WCPFC|ICCAT]+", "", file)) %>% 
                          mutate(rfmo = gsub("\\<C", "", rfmo)))
} 

spp_conf <- spp_conf %>% 
  # pivot_wider(names_from = .pred_class, values_from = confidence) %>% 
  # pivot_longer(`0`:`1`, names_to = ".pred_class", values_to = "confidence") %>% 
  pivot_wider(names_from = "rfmo", values_from = "confidence") %>% 
  arrange(indID, .pred_class)

write.csv(spp_conf, 
          file.path(here::here(), "tables", "supplemental", "longline_classification_spp.csv"), 
          row.names =  FALSE)

###
# All sharks - not species specific
###
file_list <- list.files(file.path(here::here(), "data-updated/model-data/outputs/all-rfmo-models/all-sharks-model"), 
                        pattern = ".rds", full.names = TRUE)


# Save plots for all classifications, just sharks present, and just sharks absent
# Play around with different scales 

all_classes <- combine_files(file_list, dataset = "cell_class_confidence")
plot1a <- spatial_confidence_plots(all_classes)
plot1b <- spatial_confidence_plots(all_classes, scale = c(0,1))

pres_classes <- combine_files(file_list, dataset = "cell_class_confidence", filter = 1)
plot2a <- spatial_confidence_plots(pres_classes)
plot2b <- spatial_confidence_plots(pres_classes, scale = c(0,1))

abs_classes <- combine_files(file_list, dataset = "cell_class_confidence", filter = 0)
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

ggsave(file.path(here::here(), "figures", "supplemental", "longline_classification_confidence_justsharks.png"), 
       final, height = 6, width = 12, units = "in", dpi = 600, bg = "white")

# Save plots for distribution of predictions
predictions <- purrr::map_df(file_list, 
                             ~{readRDS(.x)$id_predictions})

beta_dist <- purrr::map_df(file_list, 
                           ~{readRDS(.x)$beta_parameters})

## By cell
predictions_cell <- predictions %>% 
  filter(run == "cell")

beta_dist_cell <- beta_dist %>% 
  filter(run == "cell")

# Where mean prediction for that ID is 1
beta_dist_cell_false <- beta_dist_cell %>% 
  filter(lower_tail == FALSE) 

# Where mean prediction for that ID is 0
beta_dist_cell_true <- beta_dist_cell %>% 
  filter(lower_tail == TRUE) 

if(nrow(beta_dist_cell_false) > 0 & nrow(beta_dist_cell_true) > 0) { 
  if(nrow(beta_dist_cell_false) > 600) { 
    beta_dist_cell_false <- beta_dist_cell_false[sample(1:nrow(beta_dist_cell_false), 600),]
  }
  if(nrow(beta_dist_cell_true) > 600) { 
    beta_dist_cell_true <- beta_dist_cell_true[sample(1:nrow(beta_dist_cell_true), 600),]
  }
  
  beta_1 <- plot_distributions(predictions_cell, beta_dist_cell_false)
  beta_0 <- plot_distributions(predictions_cell, beta_dist_cell_true)
  
  final <- ggdraw() + 
    draw_plot(beta_1, 0, 0, 0.5, 1) + 
    draw_plot(beta_0, 0.5, 0, 0.5, 1) +
    draw_plot_label(label = c("Predicted Shark Presence", "Predicted Shark Absence"), 
                    x = c(.25, 0.75), 
                    y = 1, hjust = 0.5)
  
  ggsave(file.path(here::here(), "figures", "supplemental", "longline_classification_distributions_cell_justsharks.png"),
         final, height = 6, width = 12, units = "in", dpi = 600, bg = "white")
} 

## By species
predictions_spp <- predictions %>% 
  filter(run == "spp")

beta_dist_spp <- beta_dist %>% 
  filter(run == "spp")

# Where mean prediction for that ID is 1
beta_dist_spp_false <- beta_dist_spp %>% 
  filter(lower_tail == FALSE) 

# Where mean prediction for that ID is 0
beta_dist_spp_true <- beta_dist_spp %>% 
  filter(lower_tail == TRUE) 

if(nrow(beta_dist_spp_false) > 0 & nrow(beta_dist_spp_true) > 0) { 
  if(nrow(beta_dist_spp_false) > 600) { 
    beta_dist_spp_false <- beta_dist_spp_false[sample(1:nrow(beta_dist_spp_false), 600),]
  }
  if(nrow(beta_dist_spp_true) > 600) { 
    beta_dist_spp_true <- beta_dist_spp_true[sample(1:nrow(beta_dist_spp_true), 600),]
  }

  beta_1 <- plot_distributions(predictions_spp, beta_dist_spp_false)
  beta_0 <- plot_distributions(predictions_spp, beta_dist_spp_true)

  final <- ggdraw() + 
    draw_plot(beta_1, 0, 0, 0.5, 1) + 
    draw_plot(beta_0, 0.5, 0, 0.5, 1) +
    draw_plot_label(label = c("Predicted Shark Presence", "Predicted Shark Absence"), 
                    x = c(.25, 0.75), 
                    y = 1, hjust = 0.5)
  
  ggsave(file.path(here::here(), "figures", "supplemental", "longline_classification_distributions_spp_justsharks.png"),
         final, height = 6, width = 12, units = "in", dpi = 600, bg = "white")
  } 


# High and low values from the regression
all_pred <- combine_files(file_list, dataset = "final_predict") %>% 
  mutate(.pred_lower = ifelse(.final_pred == 0, 0, .pred_lower), 
         .pred_upper = ifelse(.final_pred == 0, 0, .pred_upper)) %>% 
  group_by(latitude, longitude, year, lat, lon) %>% 
  summarise(.final_pred = sum(.final_pred, na.rm = T), 
            .pred_lower = sum(.pred_lower, na.rm = T), 
            .pred_upper = sum(.pred_upper, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(latitude, longitude, lat, lon) %>% 
  summarise(.final_pred = mean(.final_pred, na.rm = T), 
            .pred_lower = mean(.pred_lower, na.rm = T), 
            .pred_upper = mean(.pred_upper, na.rm = T))

plot1a <- spatial_confidence_plots(all_pred, field = ".pred_lower", legendtitle = "Predicted catch (lower bound)")
plot1b <- spatial_confidence_plots(all_pred, field = ".final_pred", legendtitle = "Predicted catch")
plot1c <- spatial_confidence_plots(all_pred, field = ".pred_upper", legendtitle = "Predicted catch (upper bound)")

final <- ggdraw() + 
  draw_plot(plot1a, 0, 0, 0.33, 1) + 
  draw_plot(plot1b, 0.33, 0, 0.33, 1) + 
  draw_plot(plot1c, 0.66, 0, 0.33, 1) + 
  draw_plot_label(label = c("Lower Bound", "Predicted Catch", "Upper Bound"), 
                  x = c(0+0.165, 0.33+0.165, 0.66+0.165), 
                  y = 1, hjust = 0.5)

ggsave(file.path(here::here(), "figures", "supplemental", "longline_regression_justsharks.png"), 
       final, height = 3, width = 12, units = "in", dpi = 600, bg = "white")

# Save tables for outputs as well
file_list <- list.files(file.path(here::here(), "data-updated", "model-data", "outputs", "all-rfmo-models", "all-sharks-model"), 
                        pattern = ".rds", full.names = TRUE)

spp_conf <- NULL
for(file in file_list) { 
  temp_rds <- readRDS(file)$species_class_confidence
  spp_conf <- bind_rows(spp_conf, 
                        temp_rds %>%
                          mutate(rfmo = gsub("[^IATTC|IOTC|WCPFC|ICCAT]+", "", file)) %>% 
                          mutate(rfmo = gsub("\\<C", "", rfmo)))
} 

# if(length(unique(spp_conf$.pred_class)) > 1) { 
#   spp_conf <- spp_conf %>% 
#     pivot_wider(names_from = .pred_class, values_from = confidence) %>% 
#     pivot_longer(`0`:`1`, names_to = ".pred_class", values_to = "confidence") 
# }

spp_conf <- spp_conf %>% 
  pivot_wider(names_from = "rfmo", values_from = "confidence") %>% 
  arrange(indID, .pred_class)

write.csv(spp_conf, 
          file.path(here::here(), "tables", "supplemental", "longline_classification_spp_justsharks.csv"), 
          row.names =  FALSE)
