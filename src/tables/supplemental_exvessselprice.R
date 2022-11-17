# Table for ex-vessel prices

# Load Libraries
library(tidyverse)
library(here)

min_yr = 2012 

# Ex-vessel price data
prices <- read.csv(file.path(here::here(), "data-updated/ex-vessel-prices/exvessel_price_database_1976_2019.csv")) %>% 
  distinct_all()

# Load species lists
species_list <- read.csv(file.path(here::here(), "data-updated/rfmo-data/outputs/all_data.csv")) %>% 
  filter(gear_group == "longline") %>% 
  select(species_sciname, species_commonname) %>% 
  distinct_all() 

# Fix the prices data
prices <- prices %>% 
  mutate(ASFIS_species = str_to_upper(ASFIS_species), 
         scientific_name = str_to_upper(scientific_name), 
         ASFIS_species = case_when(ASFIS_species == "INDO-PACIF. BOTTLENOSE DOLPHIN" ~  "INDO-PACIFIC BOTTLENOSE DOLPHIN", 
                                   ASFIS_species == "SHORTFIN MAKO" ~ "SHORTFIN MAKO SHARK", 
                                   ASFIS_species == "ALBACORE" ~ "ALBACORE TUNA", 
                                   TRUE ~ ASFIS_species), 
         ISSCAAP_group = case_when(ISSCAAP_group == "Sharks, rays, chimaeras" ~ "sharks and rays", 
                                   ISSCAAP_group == "Miscellaneous aquatic mammals" ~ "marine mammals", 
                                   ISSCAAP_group == "Turtles" ~ "turtles", 
                                   TRUE ~ paste(ISSCAAP_group))) %>% 
  dplyr::rename(year = Year, 
                species_commonname = ASFIS_species, 
                species_group = ISSCAAP_group)

###
# median prices by species group
###
spp_prices <- prices %>% 
  mutate(scientific_name = gsub(" SPP", "", scientific_name)) %>%
  filter(scientific_name %in% unique(species_list$species_sciname)) %>% 
  group_by(year) %>% 
  summarise(median_price_group = median(exvessel, na.rm = T)) %>% 
  ungroup() %>% 
  filter(year >= min_yr) 

# adjust for 2020
spp_prices <- spp_prices %>% 
  bind_rows(spp_prices %>% 
              filter(year == 2019) %>% 
              mutate(year = 2020, 
                     median_price_group = median_price_group * 1.012)) %>% 
  bind_rows(spp_prices %>% 
              filter(year == 2019) %>% 
              mutate(year = 2021, 
                     median_price_group = median_price_group * 1.06))

###
# target prices for individuals
###
spp_prices_ind <- prices %>% 
  mutate(scientific_name = gsub(" SPP", "", scientific_name)) %>% 
  bind_rows(prices %>% filter(grepl("SPHYRNA", scientific_name)) %>% 
              group_by(year) %>% 
              summarise(exvessel = median(exvessel, na.rm = T)) %>% 
              ungroup() %>% 
              mutate(scientific_name = "SPHYRNA")) %>% 
  filter(scientific_name %in% unique(species_list$species_sciname)) %>% 
  group_by(year, scientific_name) %>% 
  summarise(median_price_species = median(exvessel, na.rm = T)) %>% 
  ungroup() %>% 
  filter(year >= min_yr) 

# Adjust for 2020
spp_prices_ind <- spp_prices_ind  %>% 
  bind_rows(spp_prices_ind %>% 
              filter(year == 2019) %>% 
              mutate(year = 2020, 
                     median_price_species = median_price_species * 1.012)) %>% 
  bind_rows(spp_prices_ind %>% 
              filter(year == 2019) %>% 
              mutate(year = 2021, 
                     median_price_species = median_price_species * 1.06))

# Add species that aren't there
spp_prices_ind <- spp_prices_ind  %>% 
  bind_rows(spp_prices %>% 
              mutate(scientific_name = "SHARKS NEI") %>% 
              dplyr::rename(median_price_species = median_price_group) %>% 
              select(year, scientific_name, median_price_species)) %>% 
  bind_rows(spp_prices %>% 
              mutate(scientific_name = "RHINCODON TYPUS") %>% # don't have any data for whaleshark so we will just use all
              dplyr::rename(median_price_species = median_price_group) %>% 
              select(year, scientific_name, median_price_species))

# Save output
write.csv(spp_prices %>% rename(median_price_sharks = median_price_group), file.path(here::here(), "tables/supplemental/sharks_exvessel_price.csv"), 
          row.names = F)

write.csv(spp_prices_ind, file.path(here::here(), "tables/supplemental/species_exvessel_price.csv"), 
          row.names = F)
