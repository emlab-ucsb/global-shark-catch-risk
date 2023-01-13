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

# Load data
train_data_class_orig <- read.csv("training_dataset.csv") %>% 
  mutate(pres_abs_num = ifelse(pres_abs == "absent", 0, 1))

test_data_class_orig <- read.csv("testing_dataset.csv") %>% 
  mutate(pres_abs_num = ifelse(pres_abs == "absent", 0, 1))

train_data_reg_orig <- train_data_class_orig %>% 
  filter(pres_abs == "present")

test_data_reg_orig <- test_data_class_orig %>% 
  filter(pres_abs == "present")

train_class <- train_data_class_orig %>% 
  dplyr::select(pres_abs, latitude, longitude, matches("species_commonname$|sdm$|mean_sst$|mean_chla$|median_price_species$|bycatch_total_effort$")) %>% 
  distinct_all()

final_test <- test_data_class_orig %>% 
  dplyr::select(pres_abs, catch, matches("species_commonname$|sdm$|mean_sst$|mean_chla$|median_price_species$|bycatch_total_effort$"), 
                latitude, longitude, year) %>% 
  distinct_all()

# Start Recipes
set.seed(1234)
class_recipe <- recipe(pres_abs ~ ., data = train_class) %>%
  themis::step_upsample(pres_abs, over_ratio = 1) %>% 
  step_zv(all_predictors(), - latitude, - longitude) %>%
  step_center(all_numeric(), -all_outcomes(), - latitude, - longitude) %>% 
  step_scale(all_numeric(), -all_outcomes(), - latitude, - longitude) %>%
  step_unknown(all_nominal(), -all_outcomes()) %>%
  step_impute_knn(all_numeric(), -all_outcomes(), neighbors = 3, impute_with = imp_vars(latitude, longitude)) %>%
  update_role(latitude, longitude, new_role = "none") %>% 
  step_dummy(all_nominal(), -all_outcomes())

###
# Classification Model
###
class_model <- rand_forest(trees = 500, 
                           mtry = 5, 
                           min_n = 5) %>% 
  set_engine("ranger", seed = 1234, 
             importance = "impurity") %>% 
  set_mode("classification") %>% 
  translate()

# Set Workflow
class_wflow <- 
  workflow() %>% 
  add_model(class_model) %>% 
  add_recipe(class_recipe)

class_fit <- class_wflow %>% 
  fit(data = train_class)

final_testing_class <- predict(class_fit, final_test, type = "prob") %>% 
  bind_cols(predict(class_fit, final_test) %>% 
              mutate(.pred_class = ifelse(.pred_class == "absent", 0, 1))) %>%
  bind_cols(final_test) %>% 
  mutate(indID = year)

final_testing_class_averaged <- final_testing_class %>% 
  group_by(indID) %>% 
  summarise(.pred_class = (mean(.pred_class))) %>% 
  ungroup()

# Get pred interval
data <- final_testing_class

pred_conf <- final_testing_class_averaged |>
  
dplyr::mutate(confidence = furrr::future_map_dbl(.data$indID, function(x){
  
  line_classif <- which(.data$indID == x) # in averaged data frame
  predictions <- data$.pred_present[which(data$indID == x)]
  
  if (length(predictions) > 1 && (all(predictions == 1) || all(predictions == 0))){
    conf <- 1
  }else{
    # beta fitting
    beta_par <- EnvStats::ebeta(predictions, method = "mle")$parameters
    
    if (.data$.pred_class[line_classif] == 1){
      conf <- stats::pbeta(q = 0.5,
                           shape1 = beta_par[1], shape2 = beta_par[2], lower.tail = FALSE)
      
    }else{
      conf <- stats::pbeta(q = 0.5,
                           shape1 = beta_par[1], shape2 = beta_par[2], lower.tail = TRUE)
    }
    
  }
  
}, .options = furrr::furrr_options(seed = TRUE)))

