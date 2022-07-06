# Figure 2. multipanel figure with highest catch of (3) tuna species, 
# cpue of the tuna species in that area and cpue of shark species in that
# area

# Some important notes: 
# 1) WCPFC does not contain 2018 data for target catch. It looks like we knew this already though. 
#    Trying to re-download and add the data back would likely require us to re-run everything
# 2) Areas that are selected are those that are greater than or equal to the 90th quantile for 
#    mean annual catch of each species in each RFMO. Each RFMO should be represented in the
#    output map.

# Load libraries
library(tidyverse)
library(cowplot)
library(sf)
library(tmap)
library(here)

# Load plotting defaults
source(file.path(here::here(), "src/figures/plot_defaults.R"))

# Load data - use cleaned data at the 1x1 resolution using count (not mt converted to count)
list_files <- list.files(file.path(here::here(), "data/model-data/inputs/all-rfmo-models/"), 
                         pattern = "1x1_tuna_hooks", full.names = TRUE)

top_ten <- NULL
top_ten_yearly <- NULL

for(file in list_files) { 
  temp <- read.csv(file)
  
  # Get mean annual catch and effort
  temp_sub <- temp %>%
    group_by(latitude, longitude, species_commonname, species_group, effort_units, catch_units, rfmo) %>%
    summarise(mean_annual_catch = mean(catch, na.rm = T),
              mean_annual_target_effort = mean(target_effort, na.rm = T)) %>%
    ungroup() %>%
    filter(mean_annual_catch > 0)

  # Areas in each RFMO with the highest 10% of tuna catch
  top_ten_tuna <- temp_sub %>%
    filter(species_group == "tunas" & catch_units == "metric tonnes") %>%
    group_by(species_commonname) %>%
    mutate(catch_quantile_90 = quantile(mean_annual_catch, probs = 0.9, na.rm = T)) %>%
    ungroup() %>%
    filter(mean_annual_catch >= catch_quantile_90) %>%
    mutate(latlon_group = paste(latitude, longitude, sep = "|"))

  # The same areas but with shark catch
  top_ten_sharks <- temp_sub %>%
    mutate(latlon_group = paste(latitude, longitude, sep = "|")) %>%
    filter(species_group == "sharks and rays" & catch_units == "count" &
             latlon_group %in% top_ten_tuna$latlon_group)

  top_ten <- top_ten %>%
    bind_rows(top_ten_tuna) %>%
    bind_rows(top_ten_sharks)

  # Calculate yearly means and standard deviation of hotspot areas
  temp_yearly <- temp %>%
    mutate(latlon_group = paste(latitude, longitude, sep = "|")) %>% 
    filter(latlon_group %in% top_ten_tuna$latlon_group & 
             target_effort > 0) %>% 
    mutate(cpue = catch/target_effort) %>% 
    select(year, rfmo, species_commonname, species_group, 
             effort_units, catch_units, cpue, catch, latlon_group) 

  top_ten_yearly <- top_ten_yearly %>% 
    bind_rows(temp_yearly)
  
}

# Species selection (which tuna species have the highest overall catch)
top_species <- top_ten %>% 
  filter(species_group == "tunas" & catch_units == "metric tonnes") %>% 
  group_by(species_commonname) %>% 
  summarise(tot_catch = sum(mean_annual_catch, na.rm = T)) %>% 
  ungroup() %>% 
  arrange(desc(tot_catch)) %>% 
  slice_head(n = 3) 

# Grab data... 
top_three <- top_ten %>% 
  filter(species_commonname %in% top_species$species_commonname | species_group == "sharks and rays") 

# Rasterize for tunas
tuna_rasters <- stack() 

for(spp in unique(top_species$species_commonname)) { 
  temp <- top_three %>% 
    filter(species_commonname == spp) %>% 
    mutate(longitude_orig = longitude, 
           latitude_orig = latitude) %>%
    st_as_sf(., coords = c("longitude_orig", "latitude_orig"), crs = 4326) %>% 
    mutate(mean_annual_catch = 1)
  
  rfmo_catch_raster_1 <- temp %>% 
    filter(latitude%%1 == 0 & longitude%%1 == 0) %>% 
    rasterize(., whole_numbers, field = "mean_annual_catch", fun = mean, background = NA)
  
  rfmo_catch_raster_2 <- temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 == 0) %>% 
    rasterize(., whole_numbers_lon, field = "mean_annual_catch", fun = mean, background = NA) %>% 
    resample(., whole_numbers, method = "ngb")
  
  # not in the date
  # rfmo_catch_raster_3 <- temp %>%
  #   filter(latitude%%1 == 0 & longitude%%1 != 0) %>%
  #   rasterize(., whole_numbers_lat, field = "mean_annual_catch", fun = mean, background = NA) %>%
  #   resample(., whole_numbers)

  rfmo_catch_raster_4 <- temp %>% 
    filter(latitude%%1 != 0 & longitude%%1 != 0) %>% 
    rasterize(., fraction_numbers, field = "mean_annual_catch", fun = mean, background = NA) %>% 
    resample(., whole_numbers, method = "ngb")
  
  rfmo_catch_raster <- calc(stack(rfmo_catch_raster_1, rfmo_catch_raster_2, rfmo_catch_raster_4), 
                            mean, na.rm = T)
  
  tuna_rasters <- stack(tuna_rasters, rfmo_catch_raster)
}

names(tuna_rasters) <- str_to_lower(gsub(" ", "_", unique(top_species$species_commonname)))

tuna_raster_df <- as.data.frame(tuna_rasters, xy = TRUE) %>%
  mutate_at(vars(bigeye_tuna, albacore_tuna, yellowfin_tuna), as.factor)

# Plot A, D, and G (tuna location plots)
location_plots <- list()

for(name in names(tuna_rasters)) { 
  temp_plot <- ggplot() + 
    geom_raster(tuna_raster_df, 
              mapping = aes_string(x="x", y="y", fill = name)) + 
    scale_fill_manual("", values = c("1" = "darkolivegreen4", "NaN" = NA), na.value = NA) + 
    geom_sf(data = wcpfc_boundary, fill = NA, color = "black") +
    geom_sf(data = iotc_boundary, fill = NA, color = "black") +
    geom_sf(data = iccat_boundary, fill = NA, color = "black") +
    geom_sf(data = iattc_boundary, fill = NA, color = "black") +
    geom_tile(basemap_df %>% filter(!is.na(land_low_res_moll)),
              mapping = aes(x=x, y=y), fill = "black") +
    coord_sf() + 
    custom_theme + 
    theme(legend.position = "none") 
  
  location_plots[[name]] <- temp_plot
}


# Plot B, E, H (lines for tuna - cpue per year)
top_ten_yearly_mean <- top_ten_yearly %>% 
  group_by(year, species_commonname, species_group,
           effort_units, catch_units) %>% 
  summarise(mean_cpue = mean(cpue, na.rm = T), 
            sd_cpue = sd(cpue, na.rm = T))

tuna_plots <- list()

for(name in names(tuna_rasters)) { 
  temp_plot <- ggplot(data = top_ten_yearly_mean %>% 
                        filter(species_commonname == str_to_upper(gsub("_", " ", name))) %>% 
                        filter(catch_units == "metric tonnes"), 
                      aes(x = year, y = mean_cpue)) + 
    geom_point(color = "darkolivegreen4") + 
    geom_errorbar(aes(ymin = mean_cpue-sd_cpue, ymax = mean_cpue+sd_cpue), width = 0.2, 
                  color = "darkolivegreen4") +
    geom_line(color = "darkolivegreen4") + 
    scale_x_continuous(breaks = 2012:2018) +
    xlab("") + 
    ylab("Mean CPUE (mt/hook)") + 
    theme_classic() + 
    theme(text = element_text(size = 18))
  
  tuna_plots[[name]] <- temp_plot
} 

# Plot C, F, I (lines for shark species - cpue per year)
top_sharks <- top_ten %>% 
  filter(species_group == "sharks and rays" & catch_units == "count" & !grepl("NEI", species_commonname)) %>% 
  group_by(species_commonname) %>% 
  summarise(tot_catch = sum(mean_annual_catch, na.rm = T)) %>% 
  ungroup() %>% 
  arrange(desc(tot_catch)) %>% 
  slice_head(n = 3)

shark_plots <- list()

for(name in names(tuna_rasters)) { 
  
  relevant_locations <- tuna_raster_df %>%
    select(x,y,as.name(name)) %>% 
    rename(val = as.name(name)) %>% 
    filter(val == 1) %>% 
    mutate(latlon_group = paste(y,x,sep = "|")) %>% 
    distinct_all()
  
  sharks_temp <- top_ten_yearly %>% 
    filter(latlon_group %in% relevant_locations$latlon_group & 
             species_commonname %in% top_sharks$species_commonname & 
             catch_units == "count") %>% 
    group_by(year, species_commonname, species_group,
             effort_units, catch_units) %>% 
    summarise(mean_cpue = mean(cpue, na.rm = T), 
              sd_cpue = sd(cpue, na.rm = T)) 
  
  temp_plot <- ggplot(data = sharks_temp %>% 
                        mutate(species_commonname = str_to_lower(species_commonname)), 
                      aes(x = year, y = mean_cpue, color = species_commonname)) + 
    geom_point() + 
    geom_errorbar(aes(ymin = mean_cpue-sd_cpue, ymax = mean_cpue+sd_cpue), width = 0.2) +
    geom_line() + 
    facet_wrap(vars(species_commonname), ncol = 1, scales = "free_y") + 
    scale_color_manual(values = c("blue shark" = "navy", "shortfin mako shark" = "darkorange4", 
                                  "silky shark" = "gray48")) + 
    scale_x_continuous(breaks = 2012:2018) + 
    xlab("") + 
    ylab("Mean CPUE (count/hook)") + 
    theme_classic() + 
    theme(text = element_text(size = 18), 
          legend.position = "none")
  
  shark_plots[[name]] <- temp_plot
}

# Final plot
final_plot <- ggdraw() + 
  draw_plot(location_plots[[1]], 0, 0.66, 0.33, 0.33) + 
  draw_plot(tuna_plots[[1]], 0.33, 0.66, 0.33, 0.33) + 
  draw_plot(shark_plots[[1]], 0.66, 0.66, 0.33, 0.33) + 
  draw_plot(location_plots[[2]], 0, 0.33, 0.33, 0.33) + 
  draw_plot(tuna_plots[[2]], 0.33, 0.33, 0.33, 0.33) + 
  draw_plot(shark_plots[[2]], 0.66, 0.33, 0.33, 0.33) + 
  draw_plot(location_plots[[3]], 0, 0.0, 0.33, 0.33) + 
  draw_plot(tuna_plots[[3]], 0.33, 0.0, 0.33, 0.33) + 
  draw_plot(shark_plots[[3]], 0.66, 0.0, 0.33, 0.33) + 
  draw_plot_label(label = LETTERS[1:9], x = rep(c(0, 0.33, 0.66), 3), 
                  y = rep(c(1, 0.66, 0.33), each = 3), hjust = 0, size = 30)

# Save
ggsave(here::here("figures/final/figure_2.png"), final_plot,
       width = 16, height = 14, units = "in", dpi = 600, bg = "white")
