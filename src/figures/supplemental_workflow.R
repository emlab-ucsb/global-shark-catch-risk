# Supplementary figures - workflow
# A series of figures/tables that show with full transparency the impacts of our assumptions

# For each RFMO - how many cells have fishing effort vs shark catch in the raw data, 
# compared to how many cells have fishing effort vs shark catch in the predicted data

# For each RFMO - how many cells are considered high-risk in the raw data vs how many are predicted 
# to be high risk in the model

# Percent carryover of zero vs non-zero catch cells

# Change in catch based on mt to count conversions?? 

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
list_files <- list.files(file.path(here::here(), "data-updated/model-data/outputs/all-rfmo-models"), 
                         pattern = "_untuned_final_predict", full.names = TRUE)

all_dat <- NULL
for(file in list_files) { 
  temp <- read.csv(file) %>% 
    select(rfmo, year, latitude, longitude, species_commonname, species_sciname, spatial_notes, .final_pred, 
           catch, matches("target_effort$|bycatch_total_effort$")) %>% 
    rename(effort = matches("target_effort$|bycatch_total_effort$"))
  
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
             .final_pred= .final_pred/25, 
             catch = catch/25, 
             effort = effort/25) %>% 
      mutate(spatial_notes = "center of 1x1 cell") %>%
      group_by(rfmo, year, latitude_rescaled, longitude_rescaled, species_commonname, species_sciname, spatial_notes) %>% 
      summarise(.final_pred = sum(.final_pred, na.rm = T), 
                catch = sum(catch, na.rm = T), 
                effort = sum(effort, na.rm = T)) %>% 
      ungroup() %>% 
      rename(latitude = latitude_rescaled, 
             longitude = longitude_rescaled)
  }
  
  all_dat <- all_dat %>% rbind(temp)
}

# For each RFMO, the values of cells with effort that have 0 shark catch
rfmo_effort_catch <- all_dat %>% 
  group_by(rfmo, year, latitude, longitude) %>% 
  summarise(.final_pred = sum(.final_pred, na.rm = T), # summarise so all species added together
            catch = sum(catch, na.rm = T), 
            effort = mean(effort, na.rm = T)) %>% 
  group_by(rfmo, latitude, longitude) %>%
  summarise(.final_pred = mean(.final_pred, na.rm = T), # grab total mean
            catch = mean(catch, na.rm = T),
            effort = mean(effort, na.rm = T)) %>%
  filter(effort > 0) %>% 
  mutate(catch_class_orig = ifelse(catch == 0, 0, 1), 
         catch_class_final = ifelse(.final_pred == 0, 0, 1)) 

risk_thresholds <- NULL
for(rfmos in unique(all_dat$rfmo)) { 
  
  temp <- rfmo_effort_catch %>% filter(rfmo == rfmos)
  risk_thresholds <- bind_rows(risk_thresholds, 
                               data.frame("rfmo" = rfmos, 
                                          "original_threshold" = as.numeric(quantile(temp$catch[which(temp$catch > 0)], 0.9)), 
                                          "predicted_threshold" = as.numeric(quantile(temp$.final_pred[which(temp$.final_pred > 0)], 0.9))))
}

rfmo_effort_catch <- rfmo_effort_catch %>% 
  left_join(risk_thresholds) %>% 
  group_by(rfmo) %>% 
  mutate(total_rows = n(), 
         risk_level_orig = case_when(catch >= original_threshold ~ "high", 
                                     catch == 0 ~ "not included",  
                                     TRUE ~ "low"), 
         risk_level_final = case_when(.final_pred >= predicted_threshold ~ "high",
                                      .final_pred == 0 ~ "not included", 
                                      TRUE ~ "low")) %>% 
  ungroup()

# Part 1a - how many cells have > 0 catch in the raw data? 
raw_sharks <- rfmo_effort_catch %>% 
  group_by(rfmo, catch_class_orig, total_rows) %>% 
  summarise(n = n()) %>% 
  ungroup() #%>%
  # pivot_wider(names_from = catch_class_orig, values_from = n) %>%
  # mutate(perc_total = `1`/total_rows*100)

# Part 2a - How many cells have > 0 catch in the raw data and are considered high-risk? 
raw_risk <- rfmo_effort_catch %>% 
  group_by(rfmo, risk_level_orig, total_rows) %>% 
  summarise(n = n()) %>% 
  ungroup() #%>% 
  # pivot_wider(names_from = risk_level_orig, values_from = n) %>%
  # mutate(perc_total = high/total_rows*100)

# Part 1b - how many cells have > 0 catch in the predicted data? 
pred_sharks <- rfmo_effort_catch %>% 
  group_by(rfmo, catch_class_final, total_rows) %>% 
  summarise(n = n()) %>% 
  ungroup() #%>% 
  # pivot_wider(names_from = catch_class_final, values_from = n) %>%
  # mutate(perc_total = `1`/total_rows*100)

# Part 2a - How many cells have > 0 catch in the predicted data and are considered high-risk? 
pred_risk <- rfmo_effort_catch %>% 
  group_by(rfmo, risk_level_final, total_rows) %>% 
  summarise(n = n()) %>% 
  ungroup() # %>% 
  # pivot_wider(names_from = risk_level_final, values_from = n) %>%
  # mutate(perc_total = high/total_rows*100)

# Draw the plots... 
# panel a - raw data with shark catch > 0
fig_a <- ggplot() + 
  geom_col(data = raw_sharks %>% mutate(catch_class_orig = factor(catch_class_orig, levels = c(1, 0))), 
           mapping = aes(x = rfmo, y = n, fill = catch_class_orig)) + 
  scale_fill_manual(name = "", values = c("0" = "gray", "1" = "navy"), 
                    labels = c("0" = "No sharks caught", "1" = "Sharks caught"), 
                    guide = guide_legend(ncol = 1)) + 
  scale_y_continuous("Number of cells", expand = c(0,0)) + 
  scale_x_discrete("", expand = c(0,0)) + 
  theme_classic() + 
  theme(text = element_text(size = 18), 
        legend.position = "bottom")

legend_a <- get_legend(fig_a)

fig_a <- fig_a + 
  theme(legend.position = "none")

# panel b - raw data with high risk cells
fig_b <- ggplot() + 
  geom_col(data = raw_risk, 
           mapping = aes(x = rfmo, y = n, fill = risk_level_orig)) + 
  scale_fill_manual(name = "", values = c("high" = "firebrick4", "low" = "darkorange3", "not included" = "gray"), 
                    labels = c("high" = "High-risk cells", "low" = "Low-risk cells", "not included" = "Cells with no sharks caught"), 
                    guide = guide_legend(ncol = 1)) + 
  scale_y_continuous("Number of cells", expand = c(0,0)) + 
  scale_x_discrete("", expand = c(0,0)) + 
  theme_classic() + 
  theme(text = element_text(size = 18), 
        legend.position = "bottom")

legend_b <- get_legend(fig_b) 

fig_b <- fig_b + 
  theme(legend.position = "none")

# panel c - predicted data with shark catch > 0
fig_c <- ggplot() + 
  geom_col(data = pred_sharks %>% mutate(catch_class_final = factor(catch_class_final, levels = c(1, 0))), 
           mapping = aes(x = rfmo, y = n, fill = catch_class_final)) + 
  scale_fill_manual(name = "", values = c("0" = "gray", "1" = "navy"), labels = c("0" = "No sharks caught", "1" = "Sharks caught")) + 
  scale_y_continuous("Number of cells", expand = c(0,0)) + 
  scale_x_discrete("", expand = c(0,0)) + 
  theme_classic() + 
  theme(text = element_text(size = 18), 
        legend.position = "none")

# panel d - predicted data with high risk cells
fig_d <- ggplot() + 
  geom_col(data = pred_risk, 
           mapping = aes(x = rfmo, y = n, fill = risk_level_final)) + 
  scale_fill_manual(name = "", values = c("high" = "firebrick4", "low" = "darkorange3", "not included" = "gray"), 
                    labels = c("high" = "High-risk cells", "low" = "Low-risk cells", "not included" = "Cells with no sharks caught")) + 
  scale_y_continuous("Number of cells", expand = c(0,0)) + 
  scale_x_discrete("", expand = c(0,0)) + 
  theme_classic() + 
  theme(text = element_text(size = 18), 
        legend.position = "none")

# combine
final_plot <- ggdraw() + 
  draw_plot(fig_a, 0, 0.5, 0.5, 0.45) + 
  draw_plot(fig_b, 0.5, 0.5, 0.5, 0.45) + 
  draw_plot(fig_c, 0, 0.06, 0.5, 0.45) + 
  draw_plot(fig_d, 0.5, 0.06, 0.5, 0.45) + 
  draw_plot(legend_a, 0, 0, 0.5, 0.1) + 
  draw_plot(legend_b, 0.55, 0, 0.5, 0.1) + 
  draw_plot_label(label = LETTERS[1:4], x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.55, 0.55), 
             hjust = 0, size = 30)

# Save
ggsave(here::here("figures/supplemental/workflow_zero_cells.png"), final_plot,
       width = 12, height = 7, units = "in", dpi = 600, bg = "white")

