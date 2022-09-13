# Table 1 - common name, scientific name, rfmo, catch units, spatial resolution for raw data

# Load libraries
library(raster)
library(tidyverse)
library(cowplot)
library(sf)
library(tmap)
library(here)

# Load raw data
all_data <- read.csv(file.path(here::here(), "data-updated/rfmo-data/outputs/all_data.csv"))

# Grab the data used in the final models
all_data <- all_data %>% 
  filter(spatial_notes %in% c("center of 5x5 cell", "center of 1x1 cell") &
           species_group == "sharks and rays" & 
           gear_group == "longline" & 
           catch_units %in% c("metric tonnes", "count") & 
           was_generated == "no" &
           year <= 2020 & sources != "IOTC-2020-WPEB16-DATA12_CE") # used for PS not LL

# Create table 1
table_1 <- all_data %>% 
  group_by(species_commonname, species_sciname) %>% 
  arrange(rfmo, catch_units, spatial_notes) %>% 
  summarise(rfmo = paste(unique(rfmo), collapse = ", "), 
            catch_units = paste(unique(catch_units), collapse = ", "), 
            spatial_notes = gsub("center of | cell", "", paste(unique(spatial_notes), collapse = ", "))) %>% 
  ungroup() %>% 
  mutate(species_commonname = str_to_sentence(species_commonname), 
         species_sciname = str_to_sentence(species_sciname)) %>% 
  arrange(species_sciname) 

# Save
write.csv(table_1, file.path(here::here(), "tables/table_1.csv"), row.names = FALSE)
