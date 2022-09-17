# Supplemental table for raw data of perc total catch (total) per RFMO (and total)

# Load libraries
library(tidyverse)
library(here)

# Load data
all_data <- read.csv(file.path(here::here(), "data-updated", "rfmo-data", "outputs", "all_data.csv"))

# Calculate totals for shark species
all_data <- all_data %>% 
  filter(species_group == "sharks and rays" & 
           spatial_notes %in% c("center of 5x5 cell", "center of 1x1 cell") &
           gear_group == "longline" & 
           catch_units %in% c("metric tonnes", "count") & 
           was_generated == "no" & 
           sources != "IOTC-2020-WPEB16-DATA12_CE") # only used for ps 

totals_table <- all_data %>% 
  group_by(rfmo, species_commonname, species_sciname, catch_units) %>% 
  summarise(sum_catch = sum(catch, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(rfmo, catch_units) %>% 
  mutate(total_catch = sum(sum_catch, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(perc_total_catch = round(sum_catch/total_catch*100/0.01)*0.01) %>% 
  bind_rows(all_data %>% 
              group_by(species_commonname, species_sciname, catch_units) %>% 
              summarise(sum_catch = sum(catch, na.rm = T)) %>% 
              ungroup() %>% 
              group_by(catch_units) %>% 
              mutate(total_catch = sum(sum_catch, na.rm = T)) %>% 
              ungroup() %>% 
              mutate(perc_total_catch = round(sum_catch/total_catch*100/0.01)*0.01, 
                     rfmo = "Global Total")) %>% 
  select(- sum_catch, - total_catch) %>% 
  pivot_wider(names_from = rfmo, values_from = perc_total_catch, values_fill = 0) 

write.csv(totals_table, file.path(here::here(), "tables", "supplemental", "raw_catch_rfmo_species.csv"), 
          row.names = FALSE)
  