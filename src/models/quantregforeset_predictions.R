# General
library(raster)
library(sf)
library(fasterize)
library(scales)
library(knitr)
library(RColorBrewer)
library(readxl)
library(googledrive)
library(cowplot)
library(tidyverse)
library(sp)
library(rgdal)
library(geosphere)
library(paletteer)

# Neg Binomial Libraries
library(pscl)

# ML Libraries
library(tidymodels)
library(readr)
library(broom.mixed)
library(MultivariateRandomForest)
library(vroom)
library(quantregForest)

options(scipen = 1000000)

train_data_class_orig <- read.csv("training_dataset.csv")

test_data_class_orig <- read.csv("testing_dataset.csv")

train_data_reg_orig <- train_data_class_orig %>% 
  filter(pres_abs == "present")

test_data_reg_orig <- test_data_class_orig %>% 
  filter(pres_abs == "present")

train_reg <- train_data_reg_orig %>% 
  dplyr::select(catch, latitude, longitude, matches("species_commonname$|sdm$|mean_sst$|mean_chla$|median_price_species$|bycatch_total_effort$")) %>% 
  distinct_all()

final_test <- test_data_class_orig %>% 
  dplyr::select(pres_abs, catch, matches("species_commonname$|sdm$|mean_sst$|mean_chla$|median_price_species$|bycatch_total_effort$"), 
                latitude, longitude, year) %>% 
  distinct_all()

# Run model
# Need to remove missing data in columns 
regforest <- quantregForest(x = train_reg %>% na.omit() %>% dplyr::select(-catch), 
                          y = (train_reg %>% na.omit())$catch, 
                          nthreads = 8)

predict <- predict(regforest, final_test %>% na.omit()) %>% as.data.frame()

# Grab results from tidymodels
results <- read.csv("quantiles_test_tidymodels.csv")
results <- results %>% 
  dplyr::select(-.pred_upper, -.pred_lower, -.pred_mid) %>% 
  na.omit()

# Add the intervals here
results <- results %>% 
  bind_cols(predict %>% 
              rename(.pred_lower = "quantile= 0.1", 
                     .pred_upper = "quantile= 0.9", 
                     .pred_mid = "quantile= 0.5")) %>% 
  relocate(c(".pred_lower", ".pred_upper", ".pred_mid"), .before = .pred)

# How many rows are not within the interval 
nrow(results %>% filter(!(.pred > .pred_lower & .pred < .pred_upper)))/nrow(results)*100
# about 25.9%

# Save
write.csv(results, "quantregForest_teset.csv", row.names = F)