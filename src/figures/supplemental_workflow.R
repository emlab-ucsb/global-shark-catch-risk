# Supplementary figures - workflow
# A series of figures/tables that show with full transparency the impacts of our assumptions

# Change in catch based on mt to count conversions?? Panel A) total reported count by rfmo, B) total reported count (and mt converted to count for RFMOs that we used that dataset for) by RFMO, C) total predicted count by RFMO

# For each RFMO - how many cells have fishing effort vs shark catch in the raw data, 
# compared to how many cells have fishing effort vs shark catch in the predicted data

# For each RFMO - how many cells are considered high-risk in the raw data vs how many are predicted 
# to be high risk in the model

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

# Load data - raw datasets
# 1x1 data - count only
list_files<- list.files(file.path(here::here(), "data-updated/model-data/inputs/all-rfmo-models/"), pattern = "1x1_count_hooks.csv", full.names = T)

count_1x1_obs <- NULL 

for(file in list_files) {
  temp <- read.csv(file)
  count_1x1_obs <- bind_rows(count_1x1_obs, temp)
} 

# 1x1 data - mt and count
list_files<- list.files(file.path(here::here(), "data-updated/model-data/inputs/all-rfmo-models/"), pattern = "1x1_mt_to_count_hooks.csv", full.names = T)

mt_1x1_obs <- NULL 

for(file in list_files) {
  temp <- read.csv(file)
  mt_1x1_obs <- bind_rows(mt_1x1_obs, temp)
} 

# Load final predicted data
list_files <- list.files(file.path(here::here(), "data-updated/model-data/outputs/all-rfmo-models"), 
                         pattern = "_untuned_final_predict.csv", full.names = TRUE)

all_dat <- NULL
for(file in list_files) { 
  temp <- read.csv(file) %>% 
    select(rfmo, year, latitude, longitude, species_commonname, species_sciname, spatial_notes, .final_pred, 
           catch, original_effort) 
  
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
             original_effort = original_effort/25) %>% 
      mutate(spatial_notes = "center of 1x1 cell") %>%
      group_by(rfmo, year, latitude_rescaled, longitude_rescaled, species_commonname, species_sciname, spatial_notes) %>% 
      summarise(.final_pred = sum(.final_pred, na.rm = T), 
                catch = sum(catch, na.rm = T), 
                original_effort = sum(original_effort, na.rm = T)) %>% 
      ungroup() %>% 
      rename(latitude = latitude_rescaled, 
             longitude = longitude_rescaled)
  }
  
  all_dat <- all_dat %>% rbind(temp)
}

###
# First workflow figure - Change in catch based on mt to count conversions and predictions vs reported
### 

# Calculate totals by rfmo and species
count_1x1_rfmo <- count_1x1_obs %>% 
  group_by(rfmo, year) %>%
  summarise(total_catch = sum(catch, na.rm = T)) %>% 
  ungroup()

mt_1x1_rfmo <- mt_1x1_obs %>% 
  group_by(rfmo, year) %>%
  summarise(total_catch = sum(catch, na.rm = T)) %>% 
  ungroup()

# Create bargraphs for each 
options(scipen = 1000000)

count_1x1_rfmo_fig <- ggplot() + 
  geom_col(data = count_1x1_rfmo, mapping = aes(x = year, y = total_catch, fill = rfmo), 
           color = "white", alpha = 0.6) + 
  scale_fill_manual(name = "", values = c("IATTC" = "darkorchid4", 
                                          "ICCAT" = "darkolivegreen4", 
                                          "IOTC" = "darkorange4", 
                                          "WCPFC" = "dodgerblue4")) +
  scale_x_continuous(name = "", breaks = 2012:2020, expand = c(0,0)) +
  scale_y_continuous(name = "Total Reported Catch (count)", expand = c(0,0),limits = c(0, 850000), labels = scales::comma) + 
  theme_classic() + 
  theme(legend.position = "bottom", 
        text = element_text(size = 18), 
        axis.text.x = element_text(angle = 90))

rfmo_legend <- get_legend(count_1x1_rfmo_fig)

count_1x1_rfmo_fig <- count_1x1_rfmo_fig + 
  theme(legend.position = "none")

mt_1x1_rfmo_fig <- ggplot() + 
  geom_col(data = mt_1x1_rfmo %>% filter(rfmo %in% c("IATTC", "ICCAT")) %>% # only rfmos we used this for
             bind_rows(count_1x1_rfmo %>% filter(rfmo %in% c("WCPFC", "IOTC"))), mapping = aes(x = year, y = total_catch, fill = rfmo), 
           color = "white", alpha = 0.6) + 
  scale_fill_manual(name = "", values = c("IATTC" = "darkorchid4", 
                                          "ICCAT" = "darkolivegreen4", 
                                          "IOTC" = "darkorange4", 
                                          "WCPFC" = "dodgerblue4")) +
  scale_x_continuous(name = "", breaks = 2012:2020, expand = c(0,0)) +
  scale_y_continuous(name = "Total Reported Catch (count)", limits = c(0, 850000), expand = c(0,0), labels = scales::comma) + 
  theme_classic() + 
  theme(legend.position = "none", 
        text = element_text(size = 18), 
        axis.text.x = element_text(angle = 90))

# barplot for predicted data by rfmo
pred_rfmo_fig <- ggplot() + 
  geom_col(data = all_dat %>% 
             group_by(rfmo, year) %>%
             summarise(total_catch = sum(catch, na.rm = T)) %>% 
             ungroup(), 
           mapping = aes(x = year, y = total_catch, fill = rfmo), 
           color = "white", alpha = 0.6) + 
  scale_fill_manual(name = "", values = c("IATTC" = "darkorchid4", 
                                          "ICCAT" = "darkolivegreen4", 
                                          "IOTC" = "darkorange4", 
                                          "WCPFC" = "dodgerblue4")) +
  scale_x_continuous(name = "", breaks = 2012:2020, expand = c(0,0)) +
  scale_y_continuous(name = "Total Predicted Catch (count)", expand = c(0,0), limits = c(0, 850000), labels = scales::comma) + 
  theme_classic() + 
  theme(legend.position = "none", 
        text = element_text(size = 18), 
        axis.text.x = element_text(angle = 90))

# Save figure
final_plot <- ggdraw() + 
  draw_plot(count_1x1_rfmo_fig, 0, 0.05, 0.33, 0.8) + 
  draw_plot(mt_1x1_rfmo_fig, 0.33, 0.05, 0.33, 0.8) + 
  draw_plot(pred_rfmo_fig, 0.66, 0.05, 0.33, 0.8) + 
  draw_plot(rfmo_legend, 0, 0, 1, 0.1) + 
  draw_plot_label(label = paste0(LETTERS[1:3], ") ", 
                                 c("Reported shark catch", 
                                   "Processed shark catch", 
                                   "Predicted shark catch")), x = c(0, 0.33, 0.66), y = 1, hjust = 0, size = 28)

# Save
ggsave(here::here("figures/supplemental/workflow_mt_to_count_revised.png"), final_plot,
       width = 15, height = 5, units = "in", dpi = 600, bg = "white")
  
###
# Second and third workflow figure:the values of cells with original_effort that have 0 shark catch and those
# that have high-risk vs not
###

rfmo_effort_catch <- all_dat %>% 
  group_by(rfmo, year, latitude, longitude) %>% 
  summarise(.final_pred = sum(.final_pred, na.rm = T), # summarise so all species added together
            catch = sum(catch, na.rm = T), 
            original_effort = mean(original_effort, na.rm = T)) %>% 
  group_by(rfmo, latitude, longitude) %>%
  summarise(.final_pred = mean(.final_pred, na.rm = T), # grab total mean
            catch = mean(catch, na.rm = T),
            original_effort = mean(original_effort, na.rm = T)) %>%
  filter(original_effort > 0) %>% 
  mutate(catch_class_orig = ifelse(catch == 0, 0, 1), 
         catch_class_final = ifelse(.final_pred == 0, 0, 1)) 

risk_thresholds <- NULL
for(rfmos in unique(all_dat$rfmo)) { 
  
  temp <- rfmo_effort_catch %>% filter(rfmo == rfmos)
  risk_thresholds <- bind_rows(risk_thresholds, 
                               data.frame("rfmo" = rfmos, 
                                          "original_threshold" = as.numeric(quantile(temp$catch[which(temp$catch > 0)], 0.9)), 
                                          "original_low_threshold" = as.numeric(quantile(temp$catch[which(temp$catch > 0)], 0.1)), 
                                          "predicted_threshold" = as.numeric(quantile(temp$.final_pred[which(temp$.final_pred > 0)], 0.9)), 
                                          "predicted_low_threshold" = as.numeric(quantile(temp$.final_pred[which(temp$.final_pred > 0)], 0.1))))
}

rfmo_effort_catch <- rfmo_effort_catch %>% 
  left_join(risk_thresholds) %>% 
  group_by(rfmo) %>% 
  mutate(total_rows = n(), 
         risk_level_orig = case_when(catch == 0 ~ "not included", 
                                     catch >= original_threshold ~ "high", 
                                     catch <= original_low_threshold ~ "low",
                                     TRUE ~ "intermediate"), 
         risk_level_final = case_when(.final_pred == 0 ~ "not included", 
                                      .final_pred >= predicted_threshold ~ "high",
                                      .final_pred <= predicted_low_threshold ~ "low",
                                      TRUE ~ "intermediate")) %>% 
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
  scale_fill_manual(name = "", values = c("high" = "darkred", "intermediate" = "orange", "low" = "khaki", "not included" = "gray"), 
                    breaks = c("high", "intermediate", "low", "not included"),
                    labels = c("high" = "High risk cells", "intermediate" = "Intermediate risk cells", "low" = "Low risk cells", "not included" = "Cells with no sharks caught"), 
                    guide = guide_legend(ncol = 2)) + 
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
  scale_fill_manual(name = "", values = c("high" = "darkred", "intermediate" = "orange", "low" = "khaki", "not included" = "gray"), 
                    breaks = c("high", "intermediate", "low", "not included"),
                    labels = c("high" = "High risk cells", "intermediate" = "Intermediate risk cells", "low" = "Low risk cells", "not included" = "Cells with no sharks caught")) + 
  scale_y_continuous("Number of cells", expand = c(0,0)) + 
  scale_x_discrete("", expand = c(0,0)) + 
  theme_classic() + 
  theme(text = element_text(size = 18), 
        legend.position = "none")

# combine
final_plot <- ggdraw() + 
  draw_plot(fig_a, 0, 0.5, 0.5, 0.43) + 
  draw_plot(fig_b, 0.5, 0.5, 0.5, 0.43) + 
  draw_plot(fig_c, 0, 0.06, 0.5, 0.43) + 
  draw_plot(fig_d, 0.5, 0.06, 0.5, 0.43) + 
  draw_plot(legend_a, 0, 0, 0.5, 0.1) + 
  draw_plot(legend_b, 0.52, 0, 0.5, 0.1) + 
  draw_plot_label(label = paste0(LETTERS[1:4], ") ", 
                                 c("Reported shark catch", 
                                   "Reported shark catch", 
                                   "Predicted shark catch", 
                                   "Predicted shark catch")), x = c(0, 0.5, 0, 0.5), 
                  y = c(1, 1, 0.55, 0.55), hjust = 0, size = 28)

# Save
ggsave(here::here("figures/supplemental/workflow_zero_cells_revised.png"), final_plot,
       width = 12, height = 7, units = "in", dpi = 600, bg = "white")


## Check if the reason we're seeing a decrease in spatial footprint of shark catch (why are there more 0s 
# in the predicted data for some RFMOs compared to the raw data?)

# Select rows with shark catch that are predicted as 0s by the model
cells_0 <- all_dat %>% 
  group_by(rfmo, latitude, longitude) %>% 
  summarise(.final_pred = sum(.final_pred, na.rm = T), # summarise so all species added together
            catch = sum(catch, na.rm = T), 
            original_effort = mean(original_effort, na.rm = T)) %>% 
  ungroup() %>% 
  filter(original_effort > 0 & .final_pred == 0 & catch > 0)

# cells_0 <- rfmo_effort_catch %>% 
#   filter(.final_pred == 0 & catch > 0) #%>% 
  # select(latitude, longitude, rfmo) %>% 
  # distinct_all() %>% 
  # mutate(cell_id = paste(longitude, latitude, rfmo, sep = "|")) 
  
# sus <- all_dat %>% 
#   mutate(cell_id = paste(longitude, latitude, rfmo, sep = "|")) %>% 
#   filter(cell_id %in% cells_0$cell_id & catch > 0)

sus_n <- cells_0 %>% 
  group_by(rfmo) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  left_join(cells_0 %>% 
              group_by(rfmo) %>% 
              summarise(total_catch = round(sum(catch)/0.01)*0.01, 
                        mean_catch = round(mean(catch)/0.01)*0.01) %>% 
              ungroup())

plot_labels <- paste0(unique(cells_0$rfmo), " (n = ", sus_n$n, ")\nmean reported = ", sus_n$mean_catch, " sharks\ntotal reported = ", 
                      sus_n$total_catch, " sharks")
names(plot_labels) <- unique(cells_0$rfmo)

temp_plot <- ggplot() + 
  geom_histogram(cells_0, mapping = aes(x = catch), fill = "navy") + 
  facet_wrap(vars(rfmo), scales = "free", labeller = labeller(rfmo = plot_labels)) + 
  ylab("Number of cells") + 
  xlab("Catch reported by tRFMOs (count)") + 
  theme_classic()

ggsave(file.path(here::here(), "figures/supplemental/predicted_anomalies.png"), 
              temp_plot, dpi = 600, bg = "white", height = 4, width = 7)


# for(rfmos in unique(sus$rfmo)) { 

  # temp_plot <- ggplot() +
  #   geom_boxplot(data = sus %>% filter(rfmo == rfmos), mapping = aes(x = cell_id, y = catch)) +
  #   facet_wrap(vars(rfmo), scales = "free") +
  #   theme(axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank())
  # 
  # assign(paste0(rfmos, "_plot"), temp_plot)
  # 
  # ggsave(file.path(here::here(), paste0("figures/supplemental/predicted_0_anomalies_", str_to_lower(rfmos), "_boxplot.png")),
  #        temp_plot, dpi = 600, bg = "white", height = 7, width = 20)
  # 
  # temp_plot2 <- ggplot() +
  #   geom_point(data = sus %>% filter(rfmo == rfmos), mapping = aes(x = cell_id, y = catch),
  #              position = position_jitter(0.1,0.1), alpha = 0.2) +
  #   facet_wrap(vars(rfmo), scales = "free") +
  #   theme(axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank())
  # 
  # ggsave(file.path(here::here(), paste0("figures/supplemental/predicted_0_anomalies_", str_to_lower(rfmos), "_scatterplot.png")),
  #        temp_plot2, dpi = 600, bg = "white", height = 7, width = 20)
  # 
  # temp_plot3 <- ggplot() +
  #   geom_col(data = sus %>% filter(rfmo == rfmos) %>% group_by(cell_id, rfmo) %>% summarise(catch = sum(catch)) %>% ungroup(),
  #              mapping = aes(x = cell_id, y = catch)) +
  #   facet_wrap(vars(rfmo), scales = "free") +
  #   theme(axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank())
  # 
  # ggsave(file.path(here::here(), paste0("figures/supplemental/predicted_0_anomalies_", str_to_lower(rfmos), "_total_catch.png")),
  #        temp_plot3, dpi = 600, bg = "white", height = 7, width = 20)
  # 
  # temp_plot4 <- ggplot() +
  #   geom_histogram(sus %>% filter(rfmo == rfmos), 
  #                  mapping = aes(x = catch), bins = max(5, max((sus %>% filter(rfmo == rfmos))$catch)/10))
  # 
  # ggsave(file.path(here::here(), paste0("figures/supplemental/predicted_0_anomalies_", str_to_lower(rfmos), "_histogram.png")), 
  #        temp_plot4, dpi = 600, bg = "white", height = 7, width = 20)
  
# } 

