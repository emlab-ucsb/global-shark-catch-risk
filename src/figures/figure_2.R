# Figure 2. multipanel figure with highest catch of (3) tuna species, 
# cpue of the tuna species in that area and cpue of shark species in that
# area

# Some important notes: 
# 1) Data we are using are 1x1 resolution (with 5x5 resolution evenly redistributed to 1x1 cells)
#    and count for bycatch catch units and metric tonnes for target catch units. All effort is
#    hooks from total target effort
# 2) WCPFC does not contain 2018 data for target catch. It looks like we knew this already though. 
#    Trying to re-download and add the data back would likely require us to re-run everything
# 3) Areas that are selected are those that are greater than or equal to the 90th quantile for 
#    mean annual catch of each species in each RFMO. Each RFMO should be represented in the
#    output map.
# 4) Some RFMOs report their data at different spatial resolutions (1x1 degrees with degrees
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

# Load data - use cleaned data at the 1x1 resolution using count (not mt converted to count)
list_files <- list.files(file.path(here::here(), "data-updated/model-data/inputs/all-rfmo-models"), 
                         pattern = "1x1_tuna_hooks", full.names = TRUE)

top_ten <- NULL
top_ten_yearly <- NULL

for(file in list_files) { 
  temp <- read.csv(file)
  
  if(grepl("WCPFC", file)) { # no shark data for 2012 for WCPFC 
    temp <- temp %>% 
      filter(year != 2012)
    }
  
  # Get mean annual catch and effort
  temp_sub <- temp %>%
    group_by(latitude, longitude, species_sciname, species_commonname, species_group, effort_units, catch_units, rfmo) %>%
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
    select(year, rfmo, species_commonname, species_sciname, species_group, 
             effort_units, catch_units, cpue, catch, latlon_group) 

  top_ten_yearly <- top_ten_yearly %>% 
    bind_rows(temp_yearly)
  
}

# Species selection (which tuna species have the highest overall catch)
top_species <- top_ten %>% 
  filter(species_group == "tunas" & catch_units == "metric tonnes") %>% 
  group_by(species_commonname, species_sciname) %>% 
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
              mapping = aes(x=x, y=y), fill = "black", color = "black") +
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
                        filter(catch_units == "metric tonnes") %>% 
                        mutate(mean_cpue = mean_cpue*10000, 
                               sd_cpue = sd_cpue*10000), 
                      aes(x = year, y = mean_cpue)) + 
    geom_point(color = "darkolivegreen4") + 
    geom_errorbar(aes(ymin = max(0, mean_cpue-sd_cpue, na.rm = T), ymax = mean_cpue+sd_cpue), width = 0.2, 
                  color = "darkolivegreen4") +
    geom_line(color = "darkolivegreen4") + 
    scale_y_continuous(limits = c(0, 0.000175*10000)) +
    scale_x_continuous(breaks = 2012:2020) +
    xlab("") + 
    ylab("Mean CPUE (mt/10,000 hooks)") + 
    theme_classic() + 
    theme(text = element_text(size = 18))
  
  tuna_plots[[name]] <- temp_plot
} 

# Plot C, F, I (lines for shark species - all sharks and endangered sharks)
top_sharks <- top_ten %>% 
  filter(species_group == "sharks and rays" & catch_units == "count")

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
                "VU")) %>% 
  filter(IUCN_Code %in% c("EN", "VU", "CR")) %>% 
  mutate(species_sciname = str_to_upper(species_sciname))

top_sharks <- species_listing %>% 
  filter(IUCN_Code %in% c("EN", "VU", "CR")) %>% 
  mutate(species_sciname = str_to_upper(Species))
  
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
             species_sciname %in% top_sharks$species_sciname & 
             catch_units == "count") %>% 
    group_by(year, effort_units, catch_units) %>% 
    summarise(mean_cpue = mean(cpue, na.rm = T), 
              sd_cpue = sd(cpue, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(species_cat = "Vulnerable Shark Species") %>% 
    bind_rows(top_ten_yearly %>% 
                filter(latlon_group %in% relevant_locations$latlon_group & 
                         catch_units == "count") %>% 
                group_by(year, effort_units, catch_units) %>% 
                summarise(mean_cpue = mean(cpue, na.rm = T), 
                          sd_cpue = sd(cpue, na.rm = T)) %>% 
                ungroup() %>% 
                mutate(species_cat = "All Reported Shark Species")) %>%
    bind_rows(data.frame(mean_cpue = c(1.3e-05, 0.0035),
                         species_cat = c("Vulnerable Shark Species",
                                         "All Reported Shark Species")))
  
  temp_plot <- ggplot(data = sharks_temp %>% 
                        mutate(mean_cpue = mean_cpue*10000, 
                               sd_cpue = sd_cpue*10000), 
                      aes(x = year, y = mean_cpue, color = species_cat)) + 
    geom_point() + 
    geom_line() + 
    geom_errorbar(aes(ymin = max(0, mean_cpue-sd_cpue, na.rm = T), ymax = mean_cpue+sd_cpue), width = 0.2) +
    facet_wrap(vars(species_cat), ncol = 1, scales = "free_y") + 
    scale_color_manual(values = c("Vulnerable Shark Species" = "navy", 
                                  "All Reported Shark Species" = "darkorange4")) + 
    scale_x_continuous(breaks = 2012:2020) + 
    xlab("") + 
    ylab("Mean CPUE (count/10,000 hooks)") + 
    theme_classic() + 
    theme(text = element_text(size = 18), 
          legend.position = "none")
  
  shark_plots[[name]] <- temp_plot
}

# Final plot
final_plot <- ggdraw() + 
  draw_plot(location_plots[[1]], 0, 0.66, 0.33, 0.32) + 
  draw_plot(tuna_plots[[1]], 0.33, 0.66, 0.33, 0.31) + 
  draw_plot(shark_plots[[1]], 0.66, 0.66, 0.33, 0.31) + 
  draw_plot(location_plots[[2]], 0, 0.33, 0.33, 0.32) + 
  draw_plot(tuna_plots[[2]], 0.33, 0.33, 0.33, 0.31) + 
  draw_plot(shark_plots[[2]], 0.66, 0.33, 0.33, 0.31) + 
  draw_plot(location_plots[[3]], 0, 0.0, 0.33, 0.32) + 
  draw_plot(tuna_plots[[3]], 0.33, 0.0, 0.33, 0.31) + 
  draw_plot(shark_plots[[3]], 0.66, 0.0, 0.33, 0.31) + 
  draw_plot_label(label = c("A) High bigeye tuna catch",
                            "B) Bigeye tuna CPUE", 
                            "C) Shark CPUE", 
                            "D) High albacore tuna catch", 
                            "E) Albacore tuna CPUE", 
                            "F) Shark CPUE", 
                            "G) High yellowfin tuna catch", 
                            "H) Yellowfin tuna CPUE", 
                            "I) Shark CPUE"), 
                  x = rep(c(0.0, 0.33, 0.66), 3), 
                  y = rep(c(0.985, 0.66, 0.33), each = 3), hjust = 0, vjust = 0.5,
                  size = 30)

# Save
ggsave(here::here("figures/final/figure_2.png"), final_plot,
       width = 20.5, height = 14, units = "in", dpi = 600, bg = "white")
