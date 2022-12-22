# Figure 4 - Some way to visualize species level information in terms of total number of cells at high risk 
# (top 10% per species) to show if fishing interactions are diffuse or concentrated 

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
# 3) The top 10% of catch for each species in each RFMO is going to be 10% of the total cells
#    within that RFMO, Since we have a lot of 0s in the data, even for predicted catch, I
#    calculated the top 10% of cells for each species in each RFMO using only cells in which
#    mean annual predicted catch > 0. Since a large proportion of the data could contain
#    0s, this returns a >10% of the total cells in the RFMO in every case and lets us 
#    see how much of each RFMO is actually high-risk.

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
                         pattern = "_untuned_final_predict.csv", full.names = TRUE)

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
    filter(species_commonname == spp & rfmo == rfmos) %>% 
    group_by(latitude, longitude, species_commonname, species_sciname, rfmo) %>% 
    summarise(mean_pred = mean(.final_pred)) %>% 
    ungroup() %>% 
    filter(mean_pred > 0)
  
  top_ten_q <- as.numeric(quantile(temp$mean_pred, 0.9))
  
  temp_ten <- temp %>% 
    filter(mean_pred >= top_ten_q) 
  
  top_ten_dat <- top_ten_dat %>% 
    bind_rows(temp_ten)
  }
} 

# Calculate number of cells compared to all cells present
plot_data <- top_ten_dat %>% 
  filter(species_sciname != "SHARKS NEI") %>% 
  mutate(species_sciname = str_to_sentence(species_sciname)) %>% 
         # species_commonname = gsub(" shark$| shark [(]vulpinus[)]", "", species_commonname), 
         # species_commonname = gsub("sharks,porbeagles", "sharks, porbeagles", species_commonname)) %>% 
  group_by(rfmo, species_sciname) %>% 
  summarise(n_cells = n()) %>% 
  ungroup() %>% 
  left_join(all_dat %>% 
              rowwise() %>% 
              mutate(cell_id = paste0(latitude, "|", longitude)) %>% 
              ungroup() %>% 
              select(rfmo, cell_id) %>% 
              distinct_all() %>% 
              group_by(rfmo) %>% 
              summarise(total_cells = n()) %>% 
              ungroup()
              ) %>% 
  mutate(n_cells_perc = n_cells/total_cells*100)

## Plot a: map of rfmo locations and colors
fig_4a <- ggplot() + 
  geom_sf(data = wcpfc_boundary, fill = "dodgerblue4", alpha = 0.6, color = "black") + 
  geom_sf(data = iotc_boundary, fill = "darkorange4", alpha = 0.6, color = "black") + 
  geom_sf(data = iccat_boundary, fill = "darkolivegreen4", alpha = 0.6, color = "black") + 
  geom_sf(data = iattc_boundary, fill = "darkorchid4", alpha = 0.6, color = "black") + 
  geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
            mapping = aes(x=x, y=y), fill = "black",  color = "black") +
  coord_sf() + 
  custom_theme 

## Plot b: stacked bar plot for how many cells in each species for each rfmo
fig_4b <- ggplot() + 
  geom_col(data = plot_data, 
           aes(y = species_sciname, x = n_cells_perc), fill = "white", color = "white") + 
  geom_col(data = plot_data, 
           aes(y = species_sciname, x = n_cells_perc, fill = rfmo), alpha = 0.6, color = "white") + 
  scale_x_continuous(name = "% of cells in each tRFMO considered high-risk", labels = scales::comma, 
                     expand = c(0,0)) + 
  ylab("") + 
  scale_fill_manual(name = "", values = c("IATTC" = "darkorchid4", 
                                          "ICCAT" = "darkolivegreen4", 
                                          "IOTC" = "darkorange4", 
                                          "WCPFC" = "dodgerblue4")) + 
  theme_classic() + 
  theme(panel.grid.major.x = element_line(color = "gray58"), 
        axis.text.y = element_text(face = "italic"),
        legend.position = "bottom", 
        text = element_text(size = 18))

legend <- get_legend(fig_4b)

fig_4b <- fig_4b + 
  theme(legend.position = "none")

# Final plot
final_plot <- ggdraw() + 
  draw_plot(fig_4a, 0, 0.1, 0.4, 0.9) + 
  draw_plot(fig_4b, 0.4, 0, 0.6, 0.91) + 
  draw_plot(legend, 0, 0, 0.4, 0.1) + 
  draw_plot_label(label = c("A) tRFMO regions", 
                            "B) High catch risk cells"), x = c(0, 0.4), y = c(1,1), 
                  hjust = 0, size = 30)

# Save
ggsave(here::here("figures/final/figure_4_revised.png"), final_plot,
       width = 15, height = 5, units = "in", dpi = 600, bg = "white")

