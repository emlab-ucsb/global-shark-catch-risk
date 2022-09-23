# Calculate total catch for each RFMO, in a 1x1 cell
# 1) latitude, longitude, year, gear, shark catch (count)
# 2) latitude, longitude, year, gear, shark catch (converted to mt)
# 3) latitude, longitude, year, gear, tuna catch (mt)

# Load Libraries
library(tidyverse)
library(readxl)
library(googledrive)
library(sf)
library(rgdal)
library(here)

# Load data
all_data <- read.csv(file.path(here::here(), "data-updated/rfmo-data/outputs/all_data.csv")) %>% 
  mutate(sources = ifelse(sources == "ublicLLSharkNum", "PublicLLSharkNum", sources))

# Sharks to mt conversion
mt_count <- read.csv(file.path(here::here(), "data-updated/species-information/spp_weight_count_conversion.csv")) %>%
  mutate(scientific_name = str_to_upper(scientific_name))

# Add sharks nei to the mix for sharks to mt conversion
mt_count <- mt_count %>%
  bind_rows(data.frame(scientific_name = "SHARKS NEI",
                       weight_mt = mean(mt_count$weight_mt, na.rm = T))) %>%
  dplyr::select(scientific_name, weight_mt)

# Calculate RFMO totals
rfmo_sharks_count <- NULL
rfmo_sharks_mt <- NULL
rfmo_tunas_mt <- NULL 

for(rfmos in unique(all_data$rfmo)) { 
  
  temp <- all_data %>% 
    filter(rfmo == rfmos &
             spatial_notes %in% c("center of 5x5 cell", "center of 1x1 cell") &
             gear_group %in% c("longline", "purse seine") &
             catch_units %in% c("metric tonnes", "count") & 
             was_generated == "no" &
             year >= 2012)
  
  if(rfmos == "WCPFC") { 
    temp <- temp %>% 
      select(rfmo, year, country, latitude, longitude, gear_group, time_period, effort, effort_units, spatial_notes, 
             species_commonname, species_sciname, catch, catch_units, species_group, sources) %>% 
      distinct_all() %>% 
      group_by(rfmo, year, country, latitude, longitude, gear_group, time_period, effort_units,
               spatial_notes, species_commonname, species_sciname, catch_units, species_group, sources) %>% 
      summarise(catch = sum(catch, na.rm = T), 
                effort = sum(effort, na.rm = T)) %>% 
      ungroup() 
  } else {
    temp <- temp %>%
      select(rfmo, year, country, latitude, longitude, gear_group, time_period, effort, effort_units, spatial_notes, 
             species_commonname, species_sciname, catch, catch_units, species_group, sources) %>% 
      distinct_all() %>% 
      group_by(rfmo, year, country, latitude, longitude, gear_group, time_period, effort_units,
               spatial_notes, species_commonname, species_sciname, catch_units, species_group, sources) %>% 
      summarise(catch = sum(catch, na.rm = T), 
                effort = sum(effort, na.rm = T)) %>% 
      ungroup() %>% 
      mutate(was_generated = "yes", 
             aggregation_period = time_period, 
             time_period = "yearly") %>% 
      filter(aggregation_period == "monthly")
  } 
  
  short_data <- temp
  
  # Scale to lowest possible resolution
  short_data <- short_data %>% 
    mutate(species_sciname = ifelse(is.na(species_sciname), species_commonname, species_sciname)) %>%
    group_by(rfmo, country, year, latitude, longitude, gear_group, effort_units, species_commonname, spatial_notes,
             species_sciname, catch_units, species_group, sources) %>% 
    summarise(catch = sum(catch, na.rm = T), 
              effort = sum(effort, na.rm = T)) %>% 
    ungroup() 
  
  # short_data <- short_data %>% 
  #   left_join(species_groups %>% 
  #               rename(species_group = group, 
  #                      species_commonname = english_name) %>% 
  #               mutate(species_commonname = str_to_upper(species_commonname)) %>% 
  #               select(-species_commonname),
  #             by = c("species_sciname" = "scientific_name")) %>% 
  #   select(rfmo, year, country, latitude, longitude, gear_group, effort, effort_units, species_group, spatial_notes,
  #          species_commonname, species_sciname, catch, catch_units) %>%
  #   mutate(species_sciname = case_when(species_sciname == "TETRAPTURUS AUDAX" ~ "KAJIKIA AUDAX",
  #                                      species_sciname == "MAKAIRA INDICA" ~ "ISTIOMPAX INDICA", 
  #                                      species_sciname == "TETRAPTURUS ALBIDUS" ~ "KAJIKIA ALBIDA",
  #                                      species_sciname == "TETRAPTURUS PFLUEGERI" ~ "TETRAPTURUS PFLUEGERI", 
  #                                      TRUE ~ species_sciname))
  
  # how much data contributed to "nei"s
  short_data <- short_data %>% 
    mutate(species_group = case_when(is.na(species_group) & 
                                       grepl("NEI|OTHER", species_commonname) & 
                                       grepl("SHARK", species_commonname) ~ "sharks and rays", 
                                     species_commonname == "SILKY OR BLACKTIP SHARK" ~ "sharks and rays",
                                     is.na(species_group) & 
                                       grepl("NEI|OTHER", species_commonname) & 
                                       grepl("MARLIN|BILLFISH|TUNA-LIKE", species_commonname) ~ "tuna-like species", 
                                     species_commonname == "FRIGATE AND BULLET TUNAS" ~ "tunas",
                                     is.na(species_group) & 
                                       grepl("NEI|OTHER", species_commonname) & 
                                       grepl("TUNA", species_commonname) ~ "tunas",
                                     TRUE ~ species_group))
  
  # Now standardize nei species if can...
  short_data <- short_data %>%
    mutate(species_resolution = case_when(grepl("NEI|OTHER|SPP", species_commonname) ~ "NEI",
                                          grepl("NEI|OTHER|SPP", species_sciname) ~ "NEI",
                                          species_commonname == "SILKY OR BLACKTIP SHARK" ~ "NEI",
                                          species_commonname == "FRIGATE AND BULLET TUNAS" ~ "NEI",
                                          !grepl(" ", species_sciname) ~ "NEI",
                                          TRUE ~ "SPECIES SPECIFIC"),
           species_commonname = case_when(species_commonname == "OTHER SHARKS" ~ "SHARKS NEI",
                                          species_commonname == "OTHER TUNA" ~ "TUNAS NEI",
                                          TRUE ~ species_commonname),
           species_sciname = case_when(species_commonname == "TUNAS NEI" ~ "THUNNUS",
                                       species_commonname == "SHARKS NEI" ~ "SHARKS NEI", 
                                       is.na(species_sciname) ~ species_commonname,
                                       TRUE ~ gsub(" SPP|[.]", "", species_sciname))) %>% 
    filter(species_group %in% c("tunas", "tuna-like species", "sharks and rays"))
  
  # Rescale to 1x1 degree 
  locations_y <- short_data %>%
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
  
  locations_x <- short_data %>%
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
  
  short_data <- short_data %>%
    select(rfmo, year, latitude, longitude, gear_group, species_group, spatial_notes, 
           species_commonname, species_sciname, catch, catch_units, sources)
  
  # New dataset
  short_data_5x5 <- short_data %>%
    filter(spatial_notes == "center of 5x5 cell") %>%
    left_join(locations_y) %>%
    left_join(locations_x) %>%
    mutate(longitude_rescaled = ifelse(is.na(longitude_rescaled), longitude, longitude_rescaled),
           latitude_rescaled = ifelse(is.na(latitude_rescaled), latitude, latitude_rescaled),
           catch = catch/25) %>%
    rename(latitude_orig = latitude,
           longitude_orig = longitude,
           latitude = latitude_rescaled,
           longitude = longitude_rescaled) %>%
    mutate(spatial_notes = "center of 1x1 cell") %>%
    select(colnames(short_data)) %>% 
    group_by_at(c(colnames(short_data)[colnames(short_data) != "catch"])) %>% 
    summarise_all(sum) %>% 
    ungroup()
  
  short_data_1x1 <- short_data %>% 
    filter(spatial_notes != "center of 5x5 cell") %>% 
    bind_rows(short_data_5x5) %>% 
    group_by_at(c(colnames(short_data)[colnames(short_data) != "catch"])) %>% 
    summarise_all(sum) %>% 
    ungroup()
  
  # Save sharks by count
  sharks_count_temp <- short_data_1x1 %>% 
    filter(catch_units == "count" & species_group == "sharks and rays") %>% 
    group_by(rfmo, latitude, longitude, year, gear_group, sources) %>% 
    summarise(total_catch = sum(catch, na.rm = T)) %>% 
    ungroup()
  
  rfmo_sharks_count <- rfmo_sharks_count %>% 
    bind_rows(sharks_count_temp)
  
  # Run conversions to mt
  short_data_shark_mt <- short_data_1x1 %>%
    filter(catch_units == "count" & species_group == "sharks and rays") %>% 
    left_join(mt_count, by = c("species_sciname" = "scientific_name")) %>%
    mutate(catch = catch*weight_mt,
           catch_units = "metric tonnes") %>%
    select(colnames(short_data))
  
  # Save sharks by metric tonnes
  sharks_mt_temp <- short_data_1x1 %>% 
    filter(catch_units == "metric tonnes" & species_group == "sharks and rays") %>% 
    bind_rows(short_data_shark_mt) %>% 
    group_by(rfmo, latitude, longitude, year, gear_group, sources) %>% 
    summarise(total_catch = sum(catch, na.rm = T)) %>% 
    ungroup()
  
  rfmo_sharks_mt <- rfmo_sharks_mt %>% 
    bind_rows(sharks_mt_temp)
  
  # Save tunas and tuna-like species in metric tonnes
  tunas_mt_temp <- short_data_1x1 %>% 
    filter(catch_units == "metric tonnes" & species_group %in% c("tunas", "tuna-like species")) %>% 
    group_by(rfmo, latitude, longitude, year, gear_group, sources) %>% 
    summarise(total_catch = sum(catch, na.rm = T)) %>% 
    ungroup()
  
  rfmo_tunas_mt <- rfmo_tunas_mt %>% 
    bind_rows(tunas_mt_temp)
} 

# Save outputs
write.csv(rfmo_sharks_count, 
          file.path(here::here(), "data-updated/rfmo-data/outputs/total_shark_catch_count.csv"), 
          row.names = F)

write.csv(rfmo_sharks_mt, 
          file.path(here::here(), "data-updated/rfmo-data/outputs/total_shark_catch_mt.csv"), 
          row.names = F)

write.csv(rfmo_tunas_mt, 
          file.path(here::here(), "data-updated/rfmo-data/outputs/total_tunas_catch_mt.csv"), 
          row.names = F)
