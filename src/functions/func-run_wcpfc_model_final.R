## Run WCPFC-specific model scenario - Final (pre tuning)

#' Function that runs a wide variety of model scenarios for the WCPFC region
#' @param data dataframe of catch values
#' @param save_loc output file loc to save r2
#' 

run_wcpfc_model_final <- function(data, save_loc){
  

  # Filter to only include relevant data
  prep_ll <- data %>% 
    filter(gear == "longline" & 
             catch_units == "count" & 
             effort_units == "hooks") %>% 
    mutate_if(is_character, as.factor) %>% 
    arrange(year, latitude, longitude) %>% 
    fill(mean_sst:cv_ssh, .direction = "downup") %>% 
    mutate(pres_abs = factor(ifelse(catch > 0, "present", "absent")),
           zone = paste(latitude, longitude, sep = "|"))  
  
  prep_ll <- prep_ll %>% 
    filter(target_sources == "self-reported catch and effort by flag + year")
  
  prep_ll <- prep_ll %>% 
    dplyr::select(-colnames(prep_ll %>% # remove columns where target effort or catch sums to 0
                              dplyr::select_if(grepl("target_effort|target_catch", colnames(.))) %>% 
                              dplyr::select_if(colSums(.) == 0))) %>% 
    filter(!is.na(target_effort))
  
  # Add spatial groups
  unique_pts <- prep_ll %>% 
    dplyr::select(latitude, longitude) %>% 
    distinct_all()
  
  # Convert to Spatial Points
  xy <- SpatialPointsDataFrame(
    matrix(c(unique_pts$longitude, unique_pts$latitude), ncol = 2), 
    data.frame(ID = seq(1:nrow(unique_pts))), 
    proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  
  # Calculate distances between points
  mdist <- distm(xy)
  
  # Complete heirarchical clustering based on dissimilarities
  hc <- hclust(as.dist(mdist), method = "complete")
  
  # Use clustering to group locations into 5 or 24 groups
  xy$clust <- cutree(hc, k = 5)
  
  unique_pts$location_cluster <- xy$clust
  
  # Add groups back to data
  prep_ll <- prep_ll %>% 
    left_join(unique_pts)
  
  # Clean column names 
  prep_ll <- janitor::clean_names(prep_ll)
  
  ### 
  # Run Models
  ###
  
  split_data <- train_test_split(data = prep_ll, 
                                 loc_group = "location_cluster", seed = 1234, prop = 0.25)
  
  # Assign training and testing data
  train_data_class_orig <- split_data$training
  
  test_data_class_orig <- split_data$testing
  
  train_data_reg_orig <- split_data$training %>% 
    filter(pres_abs == "present")
  
  test_data_reg_orig <- split_data$testing %>% 
    filter(pres_abs == "present")
  
  # Only run the rest if there is data present
  if(nrow(train_data_class_orig) > 0 & 
     nrow(test_data_class_orig) > 0 & 
     nrow(train_data_reg_orig) > 0 & 
     nrow(test_data_reg_orig) > 0 ) {
    
    # Grab data
    train_class <- train_data_class_orig %>% 
      dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, 
                    colnames(train_data_class_orig)[grepl("target_effort", colnames(train_data_class_orig)) &
                                                      grepl("longline", colnames(train_data_class_orig))]) %>% 
      distinct_all()
    
    train_reg <- train_data_reg_orig %>% 
      dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla,
                    colnames(train_data_reg_orig)[grepl("target_effort", colnames(train_data_reg_orig)) &
                                                    grepl("longline", colnames(train_data_reg_orig))]) %>% 
      distinct_all()
    
    final_test <- test_data_class_orig %>% 
      dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla,
                    colnames(train_data_class_orig)[grepl("target_effort", colnames(train_data_class_orig)) & 
                                                      grepl("longline", colnames(train_data_class_orig))], 
                            latitude, longitude, year) %>% 
      distinct_all()
            
     
     # Start Recipes
     class_recipe <- recipe(pres_abs ~ ., data = train_class) %>%
       themis::step_upsample(pres_abs, over_ratio = 1) %>% 
       step_center(all_numeric(), -all_outcomes()) %>% 
       step_scale(all_numeric(), -all_outcomes()) %>%
       step_unknown(all_nominal(), -all_outcomes()) %>%
       step_dummy(all_nominal(), -all_outcomes()) %>% 
       step_zv(all_predictors())
     
     reg_recipe <- recipe(catch ~ ., data = train_reg) %>% 
       step_center(all_numeric(), -all_outcomes()) %>% 
       step_scale(all_numeric(), -all_outcomes()) %>%
       step_unknown(all_nominal(), -all_outcomes()) %>%
       step_dummy(all_nominal(), -all_outcomes()) %>% 
       step_zv(all_predictors())
            
     ###
     # Classification Model
     ###
     ntrees = 100 
     
     class_model <- rand_forest(trees = ntrees) %>% 
       set_engine("ranger", 
                  importance = "impurity") %>% 
       set_mode("classification") %>% 
       translate()
     
     # Set Workflow
     class_wflow <- 
       workflow() %>% 
       add_model(class_model) %>% 
       add_recipe(class_recipe)
     
     # Fit Model
     class_fit <- class_wflow %>% 
       fit(data = train_class)
     
     # Return feature importance
     importance_class <- class_fit %>% 
       extract_fit_parsnip() %>% 
       vip::vi()
     
     ###
     # Regression Model
     ###
     reg_model <- rand_forest(trees = ntrees) %>% 
       set_engine("ranger", 
                  importance = "impurity") %>% 
       set_mode("regression") %>% 
       translate()
     
     # Set Workflow
     reg_wflow <- 
       workflow() %>% 
       add_model(reg_model) %>% 
       add_recipe(reg_recipe)
     
     # Fit Model
     reg_fit <- reg_wflow %>% 
       fit(data = train_reg)
     
     # Return feature importance
     importance_reg <- reg_fit %>% 
       extract_fit_parsnip() %>% 
       vip::vi()
     
     ###
     # Feature Importance
     ###
     total_importance <- importance_class %>% 
       mutate(component = "classification") %>% 
       bind_rows(importance_reg %>% 
                   mutate(component = "regression"))
     
     ###
     # Predict Testing Data
     ###
     final_class <- predict(class_fit, final_test) %>% 
       mutate(.pred_class = ifelse(.pred_class == "absent", 0, 1)) %>% 
       bind_cols(final_test)
     
     final_reg <- predict(reg_fit, final_test) %>% 
       bind_cols(final_test)
     
     final_predict <- final_class %>% 
       left_join(final_reg) %>% 
       mutate(.final_pred = .pred_class*.pred) 
     
     final_metrics <- metrics(data = final_predict, truth = catch, estimate = .final_pred)
     
     ### 
     # Predict Globally
     ###
     
     prep_pred <- prep_ll %>% 
       dplyr::select(year, latitude, longitude, pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla,
                     colnames(prep_ll)[grepl("target_effort", colnames(prep_ll)) &
                                               grepl("longline", colnames(prep_ll))]) %>% 
       group_by(latitude, longitude) %>% 
       mutate_at(vars(mean_chla:target_effort_portugal_longline), mean, na.rm = T) %>% 
       mutate_at(vars(mean_chla:target_effort_portugal_longline), as.integer) %>% 
       distinct_all()
     
     pred_class <- predict(class_fit, prep_pred) %>% 
       mutate(.pred_class = ifelse(.pred_class == "absent", 0, 1)) %>% 
       bind_cols(prep_pred)
     
     pred_reg <- predict(reg_fit, prep_pred) %>% 
       bind_cols(prep_pred)
     
     full_predict <- pred_class %>% 
       left_join(pred_reg) %>% 
       mutate(.final_pred = .pred_class*.pred)
     
     full_metrics <- metrics(data = full_predict, truth = catch, estimate = .final_pred)
     
     output_fit <- list("class_fit", class_fit, 
                        "reg_fit" = reg_fit, 
                        "class_predict" = final_class, 
                        "reg_predict" = final_reg, 
                        "final_predict" = final_predict, 
                        "final_metrics" = final_metrics,
                        "full_predict" = full_predict, 
                        "full_metrcis" = full_metrics)
     
     saveRDS(output_fit, paste0(save_loc, "_wcpfc_final_untuned_model.rds"))
     
     write.csv(full_predict, paste0(save_loc, "wcpfc_full_predict.csv"), row.names = FALSE)
  } 
} 
