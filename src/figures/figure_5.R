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
species_listing <- stocks(str_to_sentence(unique(all_dat$species_sciname)), 
                          fields = c("Species", "IUCN_Code", "IUCN_DateAssessed")) %>% 
  filter(!is.na(IUCN_Code)) %>% 
  group_by(Species) %>% 
  slice_max(order_by = IUCN_DateAssessed, n = 1) %>% 
  ungroup() 

# Remove species groups with no IUCN code
all_dat <- all_dat %>% 
  filter(str_to_sentence(species_sciname) %in% species_listing$Species)

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
endangered <- map(.x = species_left[which(species_left %in% str_to_lower(gsub(" ", "_", species_listing$Species[which(species_listing$IUCN_Code == "EN")])))], 
                  .f = ~{ list(eval(as.name(.x)))})
endangered <- bind_rows(endangered)

critically_endangered <- map(.x = species_left[which(species_left %in% str_to_lower(gsub(" ", "_", species_listing$Species[which(species_listing$IUCN_Code == "CR")])))], 
                             .f = ~{ list(eval(as.name(.x)))})
critically_endangered <- bind_rows(critically_endangered)

vulnerable <- map(.x = species_left[which(species_left %in% str_to_lower(gsub(" ", "_", species_listing$Species[which(species_listing$IUCN_Code == "VU")])))], 
                  .f = ~{ list(eval(as.name(.x)))})
vulnerable <- bind_rows(vulnerable)

not_threatened <- map(.x = species_left[which(species_left %in% str_to_lower(gsub(" ", "_", species_listing$Species[which(species_listing$IUCN_Code == "NT")])))], 
                  .f = ~{ list(eval(as.name(.x)))})
not_threatened <- bind_rows(not_threatened)

# Create color schemes - each species a different color
color_palette <- c(
  # critically endangered
  "CARCHARHINUS LONGIMANUS" = "#800000", # maroon
  "SPHYRNA LEWINI" = "#dcbeff", # lavender
  # endangered
  "ISURUS PAUCUS" = "#808000", # olive 
  "ALOPIAS PELAGICUS" = "darkorchid4", # purple 
  # vulnerable
  "CARCHARHINUS FALCIFORMIS" = "gray48", # match fig 2
  "LAMNA NASUS" = "#4363d8", # blue 
  "ALOPIAS SUPERCILIOSUS" = "#f032e6", # magenta 
  "SPHYRNA ZYGAENA" = "#469990", # teal 
  "ALOPIAS VULPINUS" = "#fabed4", # pink 
  # not threatened
  "PRIONACE GLAUCA" = "navy", # to match fig 2 
  "ISURUS OXYRINCHUS" = "darkorange4"# match fig 2
  )

# Plot 5a - critically endangered
# fig_5a <- 
  ggplot() + 
  geom_tile(vulnerable %>% 
              filter(!is.na(layer)) %>% 
              mutate(species_sciname = str_to_sentence(species_sciname)), 
            mapping = aes(x=x, y=y, fill=species_sciname), 
            alpha = 0.7, size = 1) + 
  # scale_fill_manual(name = "", values = color_palette,
  #                   labels = str_to_sentence(names(color_palette)), 
  #                   drop = TRUE) + 
  # scale_fill_viridis_d(name = "") + 
  paletteer::scale_fill_paletteer_d(palette = "ggthemes::calc") + 
  geom_sf(data = wcpfc_boundary, fill = NA, color = "black") +
  geom_sf(data = iotc_boundary, fill = NA, color = "black") +
  geom_sf(data = iccat_boundary, fill = NA, color = "black") +
  geom_sf(data = iattc_boundary, fill = NA, color = "black") +
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
            mapping = aes(x=x, y=y), fill = "black") +
  coord_sf() + 
  custom_theme + 
  theme(legend.position = "bottom", 
        legend.text = element_text(face = "italic"))


