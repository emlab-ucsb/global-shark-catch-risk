# Figure 7 - might get moved around a little

# Calculate temporal autocorrelation 

# Load libraries
library(raster)
library(tidyverse)
library(rfishbase)
library(cowplot)
library(sf)
library(tmap)
library(here)
library(ape)

# Load plotting defaults
source(file.path(here::here(), "src/figures/plot_defaults.R"))

# Load data - use cleaned data at the 1x1 resolution using count (not mt converted to count)
list_files <- list.files(file.path(here::here(), "data-updated/model-data/outputs/all-rfmo-models"), 
                         pattern = "_untuned_final_predict", full.names = TRUE)

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
  
  remove(temp)
}

all_dat <- all_dat %>% mutate(cell_id = paste(longitude, latitude, rfmo, sep = "|"))

# Save output
write.csv(all_dat, file.path(here::here(), "data-updated", "model-data", "outputs", "predicted_catch_yearly.csv"), row.names = F)


# # Calculate Moran's I for each year - using this example (https://stats.oarc.ucla.edu/r/faq/how-can-i-calculate-morans-i-in-r/)
# 
# results <- list() 
# 
# for(yr in unique(all_dat$year)) { 
#   
#   temp <- all_dat %>% 
#     filter(year == yr)
#   
#   for(rfmos in unique(temp$rfmo)) { 
#   
#     temp <- temp %>% 
#       filter(rfmo == rfmos) %>% 
#       mutate(.final_pred = sum(.final_pred, na.rm = T))
#     
#     # Step 1: Generate matrix oc inverse distance weights
#     temp_dists <- as.matrix(dist(cbind(temp$longitude, temp$latitude)))
#     temp_dists_inv <- 1/temp_dists
#     diag(temp_dists_inv) <- 0
#     
#     # Step 2: Calculate Moran's I
#     temp_moran <- Moran.I(temp$.final_pred, temp_dists_inv)
#     
#     # Save results
#     results <- append(results, temp_moran)
#   }
# }
# 
# # Try again... 
# results <- NULL
# 
# for(spp in unique(all_dat$species_commonname)) { 
# all_dat_2 <- all_dat %>% 
#   filter(species_commonname == spp) %>% 
#   group_by(latitude, longitude, year) %>% 
#   summarise(total_pred = sum(.final_pred, na.rm = TRUE)) %>% 
#   ungroup() 
# 
# if(sum(all_dat_2$total_pred) == 0) { next } 
# 
# all_dat_2 <- all_dat_2 %>% 
#   arrange(year) %>% 
#   pivot_wider(names_from = year, values_from = total_pred, values_fill = 0)
# 
# for(i in c(4:ncol(all_dat_2))) { 
#   temp_results <- cor.test(all_dat_2[[4]], all_dat_2[[i]])
#   
#   results <- bind_rows(results, 
#                        data.frame("species_commonname" = spp, 
#                                   "year" = colnames(all_dat_2)[i], 
#                                   "p_val" = temp_results$p.value, 
#                                   "coef" = temp_results$estimate))
# }
# } 
# 
# results <- results %>% 
#   filter(p_val < 0.05) %>% 
#   select(-p_val) %>% 
#   pivot_wider(names_from = year, values_from = coef) %>% 
#   mutate(`2012` = ifelse(!is.na(`2013`), "baseline", NA), 
#          `2013` = ifelse(is.na(`2013`), "baseline", `2013`)) %>% 
#   relocate(`2012`, .before = `2013`)
# 
# write.csv(results, file.path(here::here(), "tables", "supplemental", "yearly_predict_correlation.csv"), 
#           row.names = F)

# Run autoregressive (AR) model
# An AR is an arima(1,0,0) model (https://financetrain.com/estimating-autoregressive-ar-model-in-r)
all_dat_pivot <- all_dat %>% 
  mutate(cell_id = paste(longitude, latitude, rfmo, sep = "|")) %>% 
  group_by(year, latitude, longitude, rfmo, cell_id) %>% 
  summarise(.final_pred = sum(.final_pred, na.rm = T)) %>% 
  ungroup() %>% 
  pivot_wider(values_from = .final_pred, names_from = year) %>% 
  pivot_longer(-c(cell_id, latitude, longitude, rfmo), names_to = "year", values_to = ".final_pred") %>% 
  arrange(year, longitude, latitude)

ar_results <- NULL 
for(cell in unique(all_dat_pivot$cell_id)) { 
  temp <- all_dat_pivot %>% filter(cell_id == cell)
  
  min_year <- min(temp$year[which(!is.na(temp$.final_pred))])
  max_year <- max(temp$year[which(!is.na(temp$.final_pred))])
  
  temp <- temp %>% 
    filter(year >= min_year & year <= max_year)
  
  if(nrow(temp) < 2) {next}
  
  temp_test <- temp %>% filter(!is.na(.final_pred))
  
  if(temp_test$year[1] == min_year & temp_test$year[2] == max_year) { next }
  
  if(length(unique(temp$.final_pred[!is.na(temp$.final_pred)])) == 1) { next }
  
  ts_temp <- ts(temp$.final_pred, start = min_year, end = max_year, frequency = 1)
  
  arima_temp <- arima(ts_temp[!is.na(ts_temp)], order = c(1,0,0), method = "ML")
  
  ar_results <- bind_rows(ar_results, 
                          data.frame("cell_id" = cell, 
                                     "year_start" = min_year, 
                                     "year_end" = max_year,
                                     "n_nas" = length(ts_temp[is.na(ts_temp)]),
                                     "ar1" = arima_temp$coef[1], 
                                     "intercept" = arima_temp$coef[2], 
                                     "sigma2" = arima_temp$sigma2, 
                                     "log_likelihood" = arima_temp$loglik, 
                                     "aic" = arima_temp$aic))
  }

write.csv(ar_results,  file.path(here::here(), "tables", "supplemental", "ar_results.csv"), 
                    row.names = F)

ar_results <- read.csv(file.path(here::here(), "tables", "supplemental", "ar_results.csv"))

ar_results <- ar_results %>% 
  separate(cell_id, into = c("longitude", "latitude", "rfmo"), sep = "[|]", remove = FALSE)

ggplot() + 
  geom_boxplot(data = ar_results, mapping = aes(x = rfmo, y = ar1))

# Calculate cells that are high risk (0.9), low risk (0.1), and intermediate
summarized_dat <- all_dat %>% 
  group_by(latitude, longitude, rfmo) %>% 
  summarise(.final_pred = sum(.final_pred, na.rm = T)) %>% 
  ungroup() 

thresholds <- NULL
for(rfmos in unique(summarized_dat$rfmo)) { 
  temp <- summarized_dat %>% filter(rfmo == rfmos) %>% filter(.final_pred > 0)
  
  thresholds <- thresholds %>% 
    bind_rows(data.frame("rfmo" = rfmos, 
                         "low_threshold" = as.numeric(quantile(temp$.final_pred, 0.1)), 
                         "high_threshold" = as.numeric(quantile(temp$.final_pred, 0.9))))
}

summarized_dat <- summarized_dat %>% 
  left_join(thresholds) %>% 
  group_by(rfmo, latitude, longitude) %>% 
  mutate(risk_level = case_when(.final_pred >= high_threshold ~ "high", 
                                .final_pred <= low_threshold ~ "low", 
                                TRUE ~ "intermediate")) %>% 
  ungroup() %>% 
  select(latitude, longitude, rfmo, risk_level)

ar_results <- ar_results %>% 
  mutate(latitude = as.numeric(latitude), 
         longitude = as.numeric(longitude)) %>% 
  left_join(summarized_dat)


ggplot() + 
  geom_boxplot(data = ar_results, mapping = aes(x = rfmo, y = ar1, fill = risk_level), 
               alpha = 0.5) + 
  xlab("") + 
  ylab("AR slope (AR1)") + 
  scale_fill_manual("Risk Level", values = c("high" = "firebrick1", 
                                             "intermediate" = "gold", 
                                             "low" = "olivedrab")) + 
  theme_classic()

ggsave(file.path(here::here(), "figures", "supplemental", "ar_results.png"), bg = "white")


## Clustering plot for high vs low vs intermediate cells
# Create yearly thresholds 
# Calculate cells that are high risk (0.9), low risk (0.1), and intermediate
summarized_dat_yr <- all_dat %>% 
  group_by(latitude, longitude, year) %>% 
  summarise(.final_pred = sum(.final_pred, na.rm = T)) %>% 
  ungroup() 

thresholds <- NULL
for(yr in unique(summarized_dat_yr$year)) { 
  temp <- summarized_dat_yr %>% filter(year == yr) %>% filter(.final_pred > 0)
  
  thresholds <- thresholds %>% 
    bind_rows(data.frame("year" = yr, 
                         "low_threshold" = as.numeric(quantile(temp$.final_pred, 0.1)), 
                         "high_threshold" = as.numeric(quantile(temp$.final_pred, 0.9))))
}

summarized_dat_yr <- summarized_dat_yr %>% 
  left_join(thresholds) %>% 
  group_by(year, latitude, longitude) %>% 
  mutate(risk_level = case_when(.final_pred >= high_threshold ~ "high", 
                                .final_pred <= low_threshold ~ "low", 
                                TRUE ~ "intermediate")) %>% 
  ungroup() %>% 
  select(latitude, longitude, year, risk_level)

all_dat2 <- all_dat %>% 
  group_by(rfmo, year, latitude, longitude) %>% 
  summarise(.final_pred = sum(.final_pred, na.rm = T), 
            catch = sum(catch, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(summarized_dat_yr)

for(yr in unique(all_dat2$year)) { 
  
  dat_temp <- all_dat2 %>% 
      filter(year == yr)
  
  dat_temp <- dat_temp %>%
    mutate(risk_level = case_when(risk_level == "low" ~ 1, 
                                  risk_level == "intermediate" ~ 2, 
                                  risk_level == "high" ~ 3), 
           longitude_orig = longitude, 
           latitude_orig = latitude) %>% 
    st_as_sf(., coords = c("longitude_orig", "latitude_orig"), crs = 4326)
  
  raster_stack <- stack()
  
  raster_1 <- dat_temp %>% 
    filter(latitude%%1 == 0 & longitude%%1 == 0) 
  
  if(nrow(raster_1) > 0) { 
    raster_1 <- raster_1 %>% 
      rasterize(., whole_numbers, field = "risk_level", fun = max, background = NA)
    
    raster_stack <- stack(raster_stack, raster_1)
  }
  
  raster_2 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 == 0) 
  
  if(nrow(raster_2) > 0) { 
    raster_2 <- raster_2 %>% 
      rasterize(., whole_numbers, field = "risk_level", fun = max, background = NA)
    
    raster_stack <- stack(raster_stack, raster_2)
  }
  
  raster_3 <- dat_temp %>% 
    filter(latitude%%1 == 0 & longitude%%1 != 0)  
  
  if(nrow(raster_3) > 0) { 
    raster_3 <- raster_3 %>% 
      rasterize(., whole_numbers, field = "risk_level", fun = max, background = NA)
    
    raster_stack <- stack(raster_stack, raster_3)
  }
  
  raster_4 <- dat_temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 != 0)  
  
  if(nrow(raster_4) > 0) { 
    raster_4 <- raster_4 %>% 
      rasterize(., whole_numbers, field = "risk_level", fun = max, background = NA)
    
    raster_stack <- stack(raster_stack, raster_4)
  }
  
  raster_comb <- calc(raster_stack, max, na.rm = T)
  
  raster_comb <- raster::as.data.frame(raster_comb, xy = TRUE) %>% 
    mutate(layer = as.factor(layer))
  
  plot_t <- ggplot() + 
    geom_tile(raster_comb, 
              mapping = aes(x=x, y=y, fill=layer, color = layer)) + 
    scale_fill_manual("Risk Level", values = c("khaki", "orange", "darkred"), na.value = NA, 
                         breaks = c(1,2,3), 
                         labels = c("Low", "Intermediate", "High"), 
                         guide = guide_legend(title.vjust = 0.8)) + 
    scale_color_manual("", values = c("khaki", "orange", "darkred"), na.value = NA, 
                          breaks = c(1,2,3), guide = "none") + 
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
  
  if(yr == 2012) { 
    legend_year <- get_legend(plot_t)
  }
  
  plot_t <- plot_t + 
    theme(legend.position = "none")
  
  assign(paste0("yearly_", yr, "_plot"), plot_t)
  
  # assign(paste0("yearly_", yr), raster_comb)
} 

# 3x3 plot

final_plot <- 
  ggdraw() + 
  draw_plot(yearly_2012_plot, 0, 0.7, 0.33, 0.3) + 
  draw_plot(yearly_2013_plot, 0.33, 0.7, 0.33, 0.3) + 
  draw_plot(yearly_2014_plot, 0.66, 0.7, 0.33, 0.3) + 
  draw_plot(yearly_2015_plot, 0, 0.4, 0.33, 0.3) + 
  draw_plot(yearly_2016_plot, 0.33, 0.4, 0.33, 0.3) + 
  draw_plot(yearly_2017_plot, 0.66, 0.4, 0.33, 0.3) + 
  draw_plot(yearly_2018_plot, 0, 0.1, 0.33, 0.32) + 
  draw_plot(yearly_2019_plot, 0.33, 0.1, 0.33, 0.3) + 
  draw_plot(yearly_2020_plot, 0.66, 0.1, 0.33, 0.3) + 
  draw_plot(legend_year, 0, 0, 1, 0.1) + 
  draw_plot_label(label = paste0(LETTERS[1:9], ") ", 2012:2020), x = rep(c(0, 0.33, 0.66), 3), 
                  y = rep(c(1, 0.7, 0.4), each = 3))

# Save
ggsave(file.path(here::here(), "figures", "supplemental", "yearly_hotspots.png"), 
       final_plot, bg = "white", dpi = 600, width = 20, height = 12)
  
