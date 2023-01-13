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

set.seed(1234)
reg_recipe <- recipe(catch ~ ., data = train_reg) %>%
  step_zv(all_predictors(), - latitude, - longitude) %>%
  step_center(all_numeric(), -all_outcomes(), - latitude, - longitude) %>% 
  step_scale(all_numeric(), -all_outcomes(), - latitude, - longitude) %>%
  step_unknown(all_nominal(), -all_outcomes()) %>%
  step_impute_knn(all_numeric(), -all_outcomes(), neighbors = 3, impute_with = imp_vars(latitude, longitude)) %>%
  update_role(latitude, longitude, new_role = "none") %>% 
  step_dummy(all_nominal(), -all_outcomes()) 
  
###
# Regression Model
###
reg_model <- rand_forest(trees = 500, 
                         mtry = 5, 
                         min_n = 5) %>% 
  set_engine("ranger", seed = 1234, quantreg = TRUE,
             importance = "impurity") %>% 
  set_mode("regression") %>%
  translate()
  
# Set Workflow
reg_wflow <- 
  workflow() %>% 
  add_model(reg_model) %>% 
  add_recipe(reg_recipe)
  
reg_fit <- reg_wflow %>% 
  fit(data = train_reg)

rf_wf <- workflows::workflow() %>% 
  add_model(reg_model) %>% 
  add_recipe(reg_recipe) %>% 
  fit(train_reg)

preds_bind <- function(data_fit, lower = 0, upper = 1){
  predict(
    rf_wf$fit$fit$fit, 
    workflows::pull_workflow_prepped_recipe(rf_wf) %>% bake(data_fit),
    type = "quantiles",
    quantiles = c(lower, upper, 0.5)
   ) %>% 
    with(predictions) %>%
    as_tibble() %>%
    set_names(paste0(".pred", c("_lower", "_upper", "_mid"))) %>%
    bind_cols(data_fit)
  }

# Quantile prediction
rf_preds_test <- preds_bind(final_test)
  
# Final prediction
final_testing_reg <- predict(rf_wf, final_test) %>% 
  bind_cols(final_test)

# Put them together
final <- full_join(final_testing_reg, rf_preds_test) %>% 
  relocate(c(".pred_lower", ".pred_upper", ".pred_mid"), .before = ".pred")

# How many rows are not within the interval 
nrow(final %>% filter(!(.pred > .pred_lower & .pred < .pred_upper)))/nrow(final)*100
# about 1.84%

# Save
write.csv(final, "quantiles_test_tidymodels.csv", row.names = F)
