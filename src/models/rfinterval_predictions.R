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

getwd()

# Neg Binomial Libraries
library(pscl)

# ML Libraries
library(tidymodels)
library(readr)
library(broom.mixed)
library(MultivariateRandomForest)
library(vroom)
library(rfinterval)

options(scipen = 1000000)

train_data_class_orig <- read.csv("training_dataset.csv")

test_data_class_orig <- read.csv("testing_dataset.csv")

train_data_reg_orig <- train_data_class_orig %>% 
  filter(pres_abs == "present")

test_data_reg_orig <- test_data_class_orig%>% 
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
predict <- rfinterval(formula = catch ~ ., train_data = train_reg %>% na.omit(), test_data = final_test %>% na.omit())

# Grab results from tidymodels
results <- read.csv("quantiles_test_tidymodels.csv")
results <- results %>% 
  dplyr::select(-.pred_upper, -.pred_lower, -.pred_mid) %>% 
  na.omit()

# Add the intervals here
results2 <- results %>% 
  bind_cols(predict$quantreg_interval %>% 
              rename(.pred_lower = lower, 
                     .pred_upper = upper)) %>% 
  relocate(c(".pred_lower", ".pred_upper"), .before = .pred)

# How many rows are not within the interval 
nrow(results2 %>% filter(!(.pred > .pred_lower & .pred < .pred_upper)))/nrow(results2)*100
# about 6.15%

# Save
write.csv(results2, "rfinterval_quant_test.csv", row.names = F)

# Do again for oob
# Add the intervals here
results2 <- results %>% 
  bind_cols(predict$oob_interval %>% 
              rename(.pred_lower = lower, 
                     .pred_upper = upper)) %>% 
  relocate(c(".pred_lower", ".pred_upper"), .before = .pred)

# How many rows are not within the interval 
nrow(results2 %>% filter(!(.pred > .pred_lower & .pred < .pred_upper)))/nrow(results2)*100
# about 7.29%

write.csv(results2, "rfinterval_oob_test.csv", row.names = F)
