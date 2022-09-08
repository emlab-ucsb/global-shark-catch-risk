# Figure 1. multi-panel figure with A) rfmo catch, B) RFMO cpue and C) rfmo effort

# Some important notes: 
# 1) We use the mean annual catch, effort, CPUE for each RFMO.
# 2) Data we are using are 1x1 resolution (with 5x5 resolution evenly redistributed to 1x1 cells)
#    and count for catch units
# 3) In each RFMO, we grouped things based on quantiles from 0-1 using 0.1 increments. We added
#    -0.1 to the beginning of each grouping so that we could include 0 values (instead of them
#    being treated as NA in the final grouping). We kept the unique values for the quantiles for
#    each metric. If the RFMO and metric had the full number of available quantiles (11), we 
#    assigned each group a value between 1 (smallest values) to 11 (highest values). If the 
#    RFMOs had quantiles that repeated (0, 0, 0, 0.03, [...]), we just took the unique
#    increments (0, 0.03, [...]) and adjusted the group assignment so that the minimum was
#    the (13 - number of unique quantiles) and the maxium was 11. For example, if the RFMO
#    had only 8 unique quantiles, we would begin group counts at (13-8) = 5 and end at 11 so 
#    that 7 groups of data had "names". 
# 4) Some RFMOs report their data at different spatial resolutions (1x1 degrees with degrees
#    centered around whole numbers [e.g., 150] vs centered around half numbers [e.g., 150.5]).
#    To combat this, we first rasterized data in groups using "like" spatial resolutions. We 
#    then re-sampled to a common, arbitrarily chosen CRS using nearest neighbor methods. This
#    method reduced artefacts of slight 0.5 degree shifts among RFMO reporting. If RFMO reporting
#    areas overlapped, we took the mean of groupings across overlapping regions.

# Load libraries
library(raster)
library(tidyverse)
library(cowplot)
library(sf)
library(tmap)
library(here)

# Load plotting defaults
source(file.path(here::here(), "src/figures/plot_defaults.R"))

# Load data - use cleaned data at the 1x1 resolution using count (not mt converted to count)
list_files <- list.files(file.path(here::here(), "data-updated/model-data/inputs/all-rfmo-models"), 
                         pattern = "1x1_count_hooks", full.names = TRUE)

all_dat <- NULL
for(file in list_files) { 
  all_dat <- bind_rows(all_dat, read.csv(file))
}

# Calculate mean totals per year
rfmo_totals <- all_dat %>% 
  group_by(latitude, longitude, year, rfmo) %>% 
  summarise(total_catch = sum(catch, na.rm = T), 
            total_bycatch_effort = mean(bycatch_total_effort, na.rm = T), # repeated across individuals
            total_cpue = total_catch/total_bycatch_effort) %>% 
  ungroup() %>%
  group_by(latitude, longitude, rfmo) %>% 
  summarise(mean_total_catch = mean(total_catch, na.rm = T), 
            mean_bycatch_effort = mean(total_bycatch_effort, na.rm = T), 
            mean_cpue = mean(total_cpue, na.rm = T))%>% 
  ungroup() %>% 
  mutate(latitude_temp = latitude, 
         longitude_temp = longitude) %>% 
  st_as_sf(., coords = c("longitude_temp", "latitude_temp"), crs = 4326) 

# Scale by rfmo
rfmo_totals_scaled <- NULL
for(rfmos in unique(rfmo_totals$rfmo)) { 
  temp_dat <- rfmo_totals %>% 
    filter(rfmo == rfmos)
  
  quant_breaks_catch <- c(-0.1, unique(as.numeric(quantile(temp_dat$mean_total_catch, probs = seq(0, 1, 0.1), na.rm = T))))
  quant_breaks_effort <- c(-0.1,unique(as.numeric(quantile(temp_dat$mean_bycatch_effort, probs = seq(0, 1, 0.1), na.rm = T))))
  quant_breaks_cpue <- c(-0.1, unique(as.numeric(quantile(temp_dat$mean_cpue, probs = seq(0, 1, 0.1), na.rm = T))))
  
  temp_dat <- temp_dat %>% 
    mutate(mean_total_catch_scaled = as.numeric(paste(cut(mean_total_catch, 
                                         breaks = quant_breaks_catch, 
                                         labels = c((13-length(quant_breaks_catch)):11)))), 
           mean_bycatch_effort_scaled = as.numeric(paste(cut(mean_bycatch_effort, 
                                            breaks = quant_breaks_effort, 
                                            labels = c((13-length(quant_breaks_effort)):11)))), 
           mean_cpue_scaled = as.numeric(paste(cut(mean_cpue, 
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
ggsave(here::here("figures/final/figure_1.png"), final_plot,
       width = 10, height = 2.5, units = "in", dpi = 600, bg = "white")

# Calculate how many cells overlap
overlapping_cells <- mean_total_catch_scaled %>% 
  filter(layer == max(layer, na.rm = T)) %>% 
  mutate(name = "catch") %>% 
  bind_rows(mean_bycatch_effort_scaled %>% 
              filter(layer == max(layer, na.rm = T)) %>% 
              mutate(name = "effort")) %>% 
  bind_rows(mean_cpue_scaled %>% 
              filter(layer == max(layer, na.rm = T)) %>% 
              mutate(name = "cpue")) %>% 
  select(-layer) %>% 
  mutate(value = 1) %>% 
  pivot_wider(names_from = name, values_from = value, values_fill = 0) %>% 
  st_as_sf(., coords = c("x", "y"), crs = 4326)

overlapping_proportions <- NULL

for(layer in c("wcpfc_boundary", "iccat_boundary", "iotc_boundary", "iattc_boundary")) { 
  
  overlapping_proportions <- overlapping_proportions %>% 
    bind_rows(st_intersection(overlapping_cells, eval(as.name(layer))) %>% 
      as.data.frame() %>% 
      select(-geometry) %>% 
      mutate(combination = str_squish(paste0(ifelse(catch == 1, " catch ", ""), 
                                             ifelse(effort == 1, " effort ", ""),
                                             ifelse(cpue == 1, " cpue ", ""))), 
             combination = case_when(combination %in% c("catch", "effort", "cpue") ~ "no overlap", 
             TRUE ~ combination)) %>%
      group_by(combination) %>% 
      summarise(n = n()) %>% 
      ungroup() %>% 
      mutate(perc_overlap = n/sum(n)*100, 
             rfmo = gsub("_boundary", "", layer)))
}

overlapping_proportions <- overlapping_proportions %>% 
  select(-n) %>% 
  pivot_wider(names_from = combination, values_from = perc_overlap)

write.csv(overlapping_proportions, 
          here::here("tables/supplemental/overlapping_catch_effort_cpue.csv"), 
          row.names = FALSE)

# Figure of areas with high CPUE - for discussion/results
ggplot() + 
  geom_tile(mean_cpue_scaled %>% 
            filter(layer == max(layer, na.rm = T)), 
            mapping = aes(x=x, y=y), fill = "blue") + 
  geom_sf(data = wcpfc_boundary, fill = NA, color = "black") + 
  geom_sf(data = iotc_boundary, fill = NA, color = "black") + 
  geom_sf(data = iccat_boundary, fill = NA, color = "black") + 
  geom_sf(data = iattc_boundary, fill = NA, color = "black") + 
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)), 
            mapping = aes(x=x, y=y), fill = "black", color = "black") + 
  coord_sf() + 
  custom_theme + 
  theme(legend.position = "none")

# Grab species catches within those areas
rfmo_spp_totals <- all_dat %>% 
  group_by(latitude, longitude, year, rfmo, species_commonname) %>% 
  summarise(total_catch = sum(catch, na.rm = T), 
            total_bycatch_effort = mean(bycatch_total_effort, na.rm = T), # repeated across individuals
            total_cpue = total_catch/total_bycatch_effort) %>% 
  ungroup() %>%
  group_by(latitude, longitude, rfmo, species_commonname) %>% 
  summarise(mean_total_catch = mean(total_catch, na.rm = T), 
            mean_bycatch_effort = mean(total_bycatch_effort, na.rm = T), 
            mean_cpue = mean(total_cpue, na.rm = T))%>% 
  ungroup() %>% 
  mutate(longitude_orig = longitude, 
         latitude_orig = latitude) %>% 
  st_as_sf(., coords = c("longitude", "latitude"), crs = 4326) 

overlapped_raster <- rasterize(overlapping_cells %>% filter(cpue == 1), whole_numbers, field = "cpue", method = "ngb")
overlapped_raster <- rasterToPolygons(overlapped_raster)
overlapped_raster <- st_as_sf(overlapped_raster)

rfmo_spp_totals_overlap <- lengths(st_intersects(rfmo_spp_totals,
                                                 overlapped_raster))  > 0

rfmo_spp_totals <- rfmo_spp_totals[rfmo_spp_totals_overlap, ]

ggplot() + 
  geom_sf(data = rfmo_spp_totals %>% 
            group_by(rfmo, longitude_orig, latitude_orig) %>% 
            slice_max(mean_cpue, n = 1, with_ties = FALSE) %>% 
            ungroup(), aes(color = species_commonname), size = 0.2) +
  geom_sf(data = wcpfc_boundary, fill = NA, color = "black") + 
  geom_sf(data = iotc_boundary, fill = NA, color = "black") + 
  geom_sf(data = iccat_boundary, fill = NA, color = "black") + 
  geom_sf(data = iattc_boundary, fill = NA, color = "black") + 
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)), 
            mapping = aes(x=x, y=y), fill = "black", color = "black") + 
  coord_sf() + 
  custom_theme + 
  theme(legend.position = "bottom")
