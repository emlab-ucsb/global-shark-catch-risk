# Figure 4 - Some way to visualize species level information in terms of total number of cells at high risk 
# (top 10% per species) to show if fishing interactions are diffuse or concentrated 

# Load libraries
library(tidyverse)
library(cowplot)
library(sf)
library(tmap)
library(here)

# Load plotting defaults
source(file.path(here::here(), "src/figures/plot_defaults.R"))

# Load data - use cleaned data at the 1x1 resolution using count (not mt converted to count)
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

# Calculate top 10 cells of data for each species
top_ten_dat <- NULL
for(spp in unique(all_dat$species_commonname)) { 
  for(rfmos in unique(all_dat$rfmo)) { 
  temp <- all_dat %>% 
    filter(species_commonname == spp & rfmo == rfmos)
  
  top_ten_q <- as.numeric(quantile(temp$.final_pred, 0.9))
  
  temp_ten <- all_dat %>% 
    filter(.final_pred >= top_ten_q)
  
  top_ten_dat <- top_ten_dat %>% 
    bind_rows(temp_ten)
  }
} 

# Stacked bar plot for how many cells in each species for each rfmo
ggplot() + 
  geom_bar(data = top_ten_dat %>% 
             filter(species_sciname != "SHARKS NEI") %>% 
             mutate(species_sciname = str_to_sentence(species_sciname)), 
           aes(y = species_sciname), fill = "white", color = "white") + 
  geom_bar(data = top_ten_dat %>% 
             filter(species_sciname != "SHARKS NEI") %>% 
             mutate(species_sciname = str_to_sentence(species_sciname)), 
           aes(y = species_sciname, fill = rfmo), alpha = 0.6, color = "white") + 
  scale_x_continuous(name = "Number of High-Risk Cells", labels = scales::comma, 
                     expand = c(0,0)) + 
  ylab("") + 
  scale_fill_manual(name = "", values = c("IATTC" = "darkorchid4", 
                                          "ICCAT" = "darkolivegreen4", 
                                          "IOTC" = "darkorange4", 
                                          "WCPFC" = "dodgerblue4")) + 
  theme_classic() + 
  theme(panel.grid.major.x = element_line(color = "gray58"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(face = "italic"))


