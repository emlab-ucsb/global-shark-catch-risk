# Figure 5 - Map of hotspots (or top 10% of risk areas) for all threatened species 
# (VU, EN, CR from IUCN red list), color coded by species

# Some important notes:
# 1) Data we are using are 1x1 resolution (with 5x5 resolution evenly redistributed to 1x1 cells)
#    and count for bycatch catch units and metric tonnes for target catch units. All effort is
#    hooks from total target effort
# 2) Some RFMOs report their data at different spatial resolutions (1x1 degrees with degrees
#    centered around whole numbers [e.g., 150] vs centered around half numbers [e.g., 150.5]).
#    To combat this, we first rasterized data in groups using "like" spatial resolutions. We
#    then re-sampled to a common, arbitrarily chosen CRS using nearest neighbor methods. This
#    method reduced artefacts of slight 0.5 degree shifts among RFMO reporting. If RFMO reporting
#    areas overlapped, we took the mean of groupings across overlapping regions.
# 3) I calculated the top 10% of cells for each species using only cells in which
#    mean annual predicted catch > 0

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
      group_by(rfmo, year, latitude_rescaled, longitude_rescaled, species_commonname, 
               species_sciname, spatial_notes) %>% 
      summarise(.final_pred = sum(.final_pred, na.rm = T)) %>% 
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
                  "VU")) 

# Remove species groups with no IUCN code
all_dat <- all_dat %>% 
  filter(str_to_sentence(species_sciname) %in% species_listing$species_sciname)

# Map of hotspots for threatened species... probably one per IUCN listing
# Rasterize based on each species first
species_left <- NULL
for(spp in unique(all_dat$species_sciname)) { 
  
  dat_temp <- all_dat %>% 
    filter(grepl(spp, species_sciname))
  
  dat_temp <- dat_temp %>% 
    group_by(latitude, longitude, species_commonname, species_sciname) %>% 
    summarise(mean_total_pred = mean(.final_pred)) %>% 
    ungroup() %>% 
    filter(mean_total_pred > 0)
  
  if(nrow(dat_temp) == 0) { next }
  
  top_ten_q <- as.numeric(quantile(dat_temp$mean_total_pred, 0.9))
  
  dat_temp <- dat_temp %>% 
    filter(mean_total_pred >= top_ten_q) %>% 
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
    raster_comb <- raster_stack
  }
  
  raster_comb <- raster::as.data.frame(raster_comb, xy = TRUE) %>% 
    mutate(species_sciname = spp)
  
  assign(str_to_lower(gsub(" ", "_", spp)), raster_comb)
  
  species_left <- c(species_left, str_to_lower(gsub(" ", "_", spp)))
} 

# Divvy up into different vulnerability listings
endangered <- map(.x = species_left[which(species_left %in% str_to_lower(gsub(" ", "_", species_listing$species_sciname[which(species_listing$IUCN_Code == "EN")])))], 
                  .f = ~{ list(eval(as.name(.x)))})
endangered <- bind_rows(endangered)

critically_endangered <- map(.x = species_left[which(species_left %in% str_to_lower(gsub(" ", "_", species_listing$species_sciname[which(species_listing$IUCN_Code == "CR")])))], 
                             .f = ~{ list(eval(as.name(.x)))})
critically_endangered <- bind_rows(critically_endangered)

vulnerable <- map(.x = species_left[which(species_left %in% str_to_lower(gsub(" ", "_", species_listing$species_sciname[which(species_listing$IUCN_Code == "VU")])))], 
                  .f = ~{ list(eval(as.name(.x)))})
vulnerable <- bind_rows(vulnerable)

near_threatened <- map(.x = species_left[which(species_left %in% str_to_lower(gsub(" ", "_", species_listing$species_sciname[which(species_listing$IUCN_Code == "NT")])))], 
                  .f = ~{ list(eval(as.name(.x)))})
near_threatened <- bind_rows(near_threatened)

# Create color schemes - each species a different color
color_palette <- c(
  # critically endangered
  "CARCHARHINUS LONGIMANUS" = "#C5000BFF", 
  "SPHYRNA LEWINI" = "#FFD320FF", 
  "SPHYRNA MOKARRAN" = "black", # no data
  # endangered
  "ALOPIAS PELAGICUS" = "#FF950EFF", 
  "ISURUS OXYRINCHUS" = "darkorange4",
  "ISURUS PAUCUS" = "orchid3", 
  # vulnerable
  "ALOPIAS SUPERCILIOSUS" = "limegreen",
  "ALOPIAS VULPINUS" = "blue",
  "CARCHARHINUS FALCIFORMIS" = "gray48", 
  "CARCHARHINUS LIMBATUS" = "wheat2", 
  "LAMNA NASUS" = "lightskyblue",
  "SPHYRNA ZYGAENA" = "violetred4",  
  # near threatened
  "PRIONACE GLAUCA" = "navy"
  )

# paletteer::paletteer_d("ggthemes::calc")

# Plot 5a - critically endangered
fig_5a <-
  ggplot() +
  geom_tile(critically_endangered %>%
              filter(!is.na(layer)) %>% 
              arrange(species_sciname) %>% 
              mutate(species_sciname = factor(species_sciname)),
            mapping = aes(x=x, y=y, fill=species_sciname),
            height = 1, width = 1, color = NA) +
  scale_fill_manual(name = "", values = sort(color_palette[unique(critically_endangered$species_sciname)]),
                    breaks = sort(names(color_palette[unique(critically_endangered$species_sciname)])),
                    labels = sort(str_to_sentence(names(color_palette[unique(critically_endangered$species_sciname)])))) +
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic"), 
        text = element_text(size = 18))

legend_5a <- get_legend(fig_5a)

fig_5a <- fig_5a + 
  geom_sf(data = wcpfc_boundary, fill = NA, color = "black") +
  geom_sf(data = iotc_boundary, fill = NA, color = "black") +
  geom_sf(data = iccat_boundary, fill = NA, color = "black") +
  geom_sf(data = iattc_boundary, fill = NA, color = "black") +
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
            mapping = aes(x=x, y=y), fill = "black", color = "black") +
  coord_sf() +
  custom_theme + 
  theme(legend.position = "none")

# Plot 5b - endangered
fig_5b <-
  ggplot() +
  geom_tile(endangered %>%
              filter(!is.na(layer)) %>% 
              arrange(species_sciname) %>% 
              mutate(species_sciname = factor(species_sciname)),
            mapping = aes(x=x, y=y, fill=species_sciname),
            height = 1, width = 1, color = NA) +
  scale_fill_manual(name = "", values = sort(color_palette[unique(endangered$species_sciname)]),
                    breaks = sort(names(color_palette[unique(endangered$species_sciname)])),
                    labels = sort(str_to_sentence(names(color_palette[unique(endangered$species_sciname)])))) +
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic"), 
        text = element_text(size = 18))

legend_5b <- get_legend(fig_5b)

fig_5b <- fig_5b + 
  ggpattern::geom_tile_pattern(endangered %>%
                               filter(!is.na(layer)) %>%
                               select(-layer) %>%
                               group_by(x, y) %>%
                               mutate(n = n()) %>%
                               ungroup() %>%
                               filter(n > 1) %>% 
                               pivot_wider(names_from = species_sciname, values_from = n) %>% 
                               mutate(`ISURUS PAUCUS` = ifelse(!is.na(`ISURUS PAUCUS`), "ISURUS PAUCUS", NA),  
                                      `ALOPIAS PELAGICUS` = ifelse(!is.na(`ALOPIAS PELAGICUS`), "ALOPIAS PELAGICUS", NA)) %>% 
                               rowwise() %>% 
                               rename(species_1 = `ISURUS PAUCUS`, species_2 = `ALOPIAS PELAGICUS`),
                               mapping = aes(x=x, y=y, pattern_colour = species_1, pattern_fill = species_1, fill = species_2),
                               pattern = "stripe", color = NA, height = 1, width = 1, 
                               pattern_spacing = 0.01, 
                               pattern_density = 0.2, 
                               pattern_size = 0.1) +
  ggpattern::scale_pattern_colour_manual(name = "", values = sort(color_palette[unique(endangered$species_sciname)]),
                                         breaks = sort(names(color_palette[unique(endangered$species_sciname)])),
                                         labels = sort(str_to_sentence(names(color_palette[unique(endangered$species_sciname)])))) + 
  ggpattern::scale_pattern_fill_manual(name = "", values = sort(color_palette[unique(endangered$species_sciname)]),
                                       breaks = sort(names(color_palette[unique(endangered$species_sciname)])),
                                       labels = sort(str_to_sentence(names(color_palette[unique(endangered$species_sciname)])))) +
  geom_sf(data = wcpfc_boundary, fill = NA, color = "black") +
  geom_sf(data = iotc_boundary, fill = NA, color = "black") +
  geom_sf(data = iccat_boundary, fill = NA, color = "black") +
  geom_sf(data = iattc_boundary, fill = NA, color = "black") +
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
            mapping = aes(x=x, y=y), fill = "black", color = "black") +
  coord_sf() +
  custom_theme + 
  theme(legend.position = "none")

# Plot 5c - vulnerable
fig_5c <-
  ggplot() +
  geom_tile(vulnerable %>%
              filter(!is.na(layer)) %>% 
              arrange(species_sciname) %>% 
              mutate(species_sciname = factor(species_sciname)),
            mapping = aes(x=x, y=y, fill=species_sciname),
            height = 1, width = 1, color = NA) +
  scale_fill_manual(name = "", values = sort(color_palette[unique(vulnerable$species_sciname)]),
                    breaks = sort(names(color_palette[unique(vulnerable$species_sciname)])),
                    labels = sort(str_to_sentence(names(color_palette[unique(vulnerable$species_sciname)]))),  
                    guide = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic"), 
        text = element_text(size = 18))

legend_5c <- get_legend(fig_5c)

fig_5c <- fig_5c + 
  ggpattern::geom_tile_pattern(vulnerable %>%
                                 filter(!is.na(layer)) %>%
                                 select(-layer) %>%
                                 group_by(x, y) %>%
                                 mutate(n = n()) %>%
                                 ungroup() %>%
                                 filter(n > 1) %>% 
                                 pivot_wider(names_from = species_sciname, values_from = n) %>% 
                                 mutate(`CARCHARHINUS FALCIFORMIS` = ifelse(!is.na(`CARCHARHINUS FALCIFORMIS`), "CARCHARHINUS FALCIFORMIS", NA),  
                                        `CARCHARHINUS LIMBATUS` = ifelse(!is.na(`CARCHARHINUS LIMBATUS`), "CARCHARHINUS LIMBATUS", NA)) %>% 
                                 rowwise() %>% 
                                 rename(species_2 = `CARCHARHINUS LIMBATUS`) %>% 
                                 rename(species_1 = `CARCHARHINUS FALCIFORMIS`),
                               mapping = aes(x=x, y=y, pattern_colour = species_1, pattern_fill = species_1, fill = species_2),
                               pattern = "stripe", color = NA, height = 1, width = 1, 
                               pattern_spacing = 0.01, 
                               pattern_density = 0.2, 
                               pattern_size = 0.1) +
  ggpattern::scale_pattern_colour_manual(name = "", values = sort(color_palette[unique(vulnerable$species_sciname)]),
                                         breaks = sort(names(color_palette[unique(vulnerable$species_sciname)])),
                                         labels = sort(str_to_sentence(names(color_palette[unique(vulnerable$species_sciname)])))) + 
  ggpattern::scale_pattern_fill_manual(name = "", values = sort(color_palette[unique(vulnerable$species_sciname)]),
                                       breaks = sort(names(color_palette[unique(vulnerable$species_sciname)])),
                                       labels = sort(str_to_sentence(names(color_palette[unique(vulnerable$species_sciname)])))) +
  geom_sf(data = wcpfc_boundary, fill = NA, color = "black") +
  geom_sf(data = iotc_boundary, fill = NA, color = "black") +
  geom_sf(data = iccat_boundary, fill = NA, color = "black") +
  geom_sf(data = iattc_boundary, fill = NA, color = "black") +
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
            mapping = aes(x=x, y=y), fill = "black", color = "black") +
  coord_sf() +
  custom_theme + 
  theme(legend.position = "none")

# Plot 5d - near threatened
fig_5d <-
  ggplot() +
  geom_tile(near_threatened %>%
              filter(!is.na(layer)) %>% 
              arrange(species_sciname) %>% 
              mutate(species_sciname = factor(species_sciname)),
            mapping = aes(x=x, y=y, fill=species_sciname),
            height = 1, width = 1, color = NA) +
  scale_fill_manual(name = "", values = sort(color_palette[unique(near_threatened$species_sciname)]),
                    breaks = sort(names(color_palette[unique(near_threatened$species_sciname)])),
                    labels = sort(str_to_sentence(names(color_palette[unique(near_threatened$species_sciname)])))) +
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic"), 
        text = element_text(size = 18))
  
legend_5d <- get_legend(fig_5d)

fig_5d <- fig_5d + 
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
  draw_plot(fig_5a, 0, 0.57, 0.5, 0.38) + 
  draw_plot(legend_5a, 0, 0.5, 0.5, 0.07) + 
  draw_plot(fig_5b, 0.5, 0.57, 0.5, 0.38) + 
  draw_plot(legend_5b, 0.5, 0.5, 0.5, 0.07) + 
  draw_plot(fig_5c, 0, 0.07, 0.5, 0.38) + 
  draw_plot(legend_5c, 0, 0, 0.5, 0.07) + 
  draw_plot(fig_5d, 0.5, 0.07, 0.5, 0.38) + 
  draw_plot(legend_5d, 0.5, 0, 0.5, 0.07) + 
  draw_plot_label(label = c("Critically Endangered", "Endangered", 
                            "Vulnerable", "Near Threatened"),
                  x = c(0.25, 0.75, 0.25, 0.75), y = c(1,1,0.5,0.5), 
                  hjust = 0.5, size = 30)



# Save
ggsave(here::here("figures/final/figure_5.png"), final_plot,
       width = 14, height = 8, units = "in", dpi = 600, bg = "white")

