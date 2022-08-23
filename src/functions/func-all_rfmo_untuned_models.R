## Run RFMO specific model scenario without tuning

#' Function that runs a wide variety of model scenarios
#' @param data dataframe of catch values
#' @param save_loc output file loc to save r2
#' @param rfmos which rfmo to choose
#' @param effort_source which effort source to choose
#' @param logtrans whether to log transform the catch variable
#' @param mtry_class the original model's mtry for classification
#' @param mtry_reg the original model's mtry for regression
#' @param min_n_class the original model's min_n for classification
#' @param min_n_reg the original model's min_n for classification

all_rfmo_untuned_models <- function(data, save_loc, rfmos, effort_source, logtrans = FALSE, classification_variables, regression_variables, mtry_class, mtry_reg, min_n_class, min_n_reg){

  # Filter to only include relevant data
  prep_ll <- data %>% 
    filter(gear_group == "longline" & 
             catch_units == "count" & 
             effort_units == "hooks") %>% 
    mutate_if(is.integer, as.numeric) %>% 
    mutate_if(is_character, as.factor) %>% 
    mutate_if(is.integer, as.numeric) %>% 
    arrange(year, latitude, longitude) %>% 
    mutate(pres_abs = factor(ifelse(catch > 0, "present", "absent")),
           zone = paste(latitude, longitude, sep = "|"))  
  
  if(effort_source == "effort reported with target catch (flag)") { 
    prep_ll <- prep_ll %>% 
      filter(rfmo == rfmos) %>% 
      dplyr::select_if(!grepl("bycatch_total_effort", colnames(.)))
    
    prep_ll <- prep_ll %>% 
      dplyr::select(-colnames(prep_ll %>% # remove columns where target effort sums to 0
                                dplyr::select_if(grepl("target_effort", colnames(.))) %>% 
                                dplyr::select_if(colSums(., na.rm = TRUE) == 0))) 
  } 
  
  if(effort_source == "effort reported with bycatch (flag)") { 
    prep_ll <- prep_ll %>% 
      filter(rfmo == rfmos) %>% 
      dplyr::select_if(!grepl("target_effort", colnames(.)))
    
    prep_ll <- prep_ll %>% 
      dplyr::select(-colnames(prep_ll %>% # remove columns where bycatch associated effort sums to 0
                                dplyr::select_if(grepl("bycatch_total_effort", colnames(.))) %>% 
                                dplyr::select_if(colSums(., na.rm = T) == 0)))
    
    if(length(colnames(prep_ll)[grepl("bycatch_total_effort", colnames(prep_ll))]) == 0) { 
      prep_ll <- data.frame()
    }
  } 
  
  if(effort_source == "gfw effort") { 
    # Subset by RFMO
    prep_ll <- prep_ll %>% 
      filter(rfmo == rfmos) 
  }
  
  if(nrow(prep_ll) > 0) {
    
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
      dplyr::select(pres_abs, latitude, longitude, matches(classification_variables)) %>% 
      distinct_all()
    
    train_reg <- train_data_reg_orig %>% 
      dplyr::select(catch, latitude, longitude, matches(regression_variables)) %>% 
      distinct_all()
    
    final_test <- test_data_class_orig %>% 
      dplyr::select(pres_abs, catch, matches(paste0(classification_variables, regression_variables, sep = "|")), 
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
    
    set.seed(1234)
    reg_recipe <- recipe(catch ~ ., data = train_reg) %>%
      step_zv(all_predictors(), - latitude, - longitude) %>%
      step_center(all_numeric(), -all_outcomes(), - latitude, - longitude) %>% 
      step_scale(all_numeric(), -all_outcomes(), - latitude, - longitude) %>%
      step_unknown(all_nominal(), -all_outcomes()) %>%
      step_impute_knn(all_numeric(), -all_outcomes(), neighbors = 3, impute_with = imp_vars(latitude, longitude)) %>%
      update_role(latitude, longitude, new_role = "none") %>% 
      step_dummy(all_nominal(), -all_outcomes()) 
    
    if(logtrans == TRUE) { 
      reg_recipe <- reg_recipe %>% 
        step_log(catch, offset = 1, skip = TRUE)
    }
            
     ###
     # Classification Model
     ###
    class_model <- rand_forest(trees = 500, 
                               mtry = mtry_class, 
                               min_n = min_n_class) %>% 
      set_engine("ranger", seed = 1234) %>% 
      set_mode("classification") %>% 
      translate()
    
    # Set Workflow
    class_wflow <- 
      workflow() %>% 
      add_model(class_model) %>% 
      add_recipe(class_recipe)
    
    class_fit <- class_wflow %>% 
      fit(data = train_class)
    
    ###
    # Regression Model
    ###
    reg_model <- rand_forest(trees = 500, 
                             mtry = mtry_reg, 
                             min_n = min_n_reg) %>% 
      set_engine("ranger", seed = 1234) %>% 
      set_mode("regression") %>% 
      translate()
    
    # Set Workflow
    reg_wflow <- 
      workflow() %>% 
      add_model(reg_model) %>% 
      add_recipe(reg_recipe)
    
    reg_fit <- reg_wflow %>% 
      fit(data = train_reg)
    
    # Predict Final Testing
    final_testing_class <- predict(class_fit, final_test) %>% 
      mutate(.pred_class = ifelse(.pred_class == "absent", 0, 1)) %>% 
      bind_cols(final_test)
    
    final_testing_reg <- predict(reg_fit, final_test) %>% 
      bind_cols(final_test)
    
    test_predict <- final_testing_class %>% 
      left_join(final_testing_reg) %>% 
      mutate(.final_pred = .pred_class*.pred) 
    
    # Save outputs
    test_metrics <- metrics(truth = catch, estimate = .final_pred, data = test_predict)
    
    test_metrics <- test_metrics %>% 
      dplyr::select(-.estimator) %>% 
      mutate(rfmo = rfmos) %>% 
      pivot_wider(names_from = .metric, values_from = .estimate)
    
    # Predict on full dataset
    prep_pred <- prep_ll %>% select(-effort_units) %>% 
      group_by(latitude, longitude) %>% 
      mutate_at(colnames(.)[grepl("sst|chla|ssh|effort|kwh", colnames(.))], mean, na.rm = T) %>% 
      ungroup() %>% 
      dplyr::select(pres_abs, catch, 
                    matches(paste0(classification_variables, 
                                   regression_variables, sep = "|"))) %>% 
      distinct_all()
    
    pred_class <- predict(class_fit, prep_pred) %>% 
      mutate(.pred_class = ifelse(.pred_class == "absent", 0, 1)) %>% 
      bind_cols(prep_pred)
    
    pred_reg <- predict(reg_fit, prep_pred) %>% 
      bind_cols(prep_pred)
    
    final_predict <- pred_class %>% 
      left_join(pred_reg) %>% 
      mutate(.final_pred = .pred_class*.pred)
    
    # Save outputs
    final_metrics <- metrics(truth = catch, estimate = .final_pred, data = final_predict)
    
    final_metrics <- final_metrics %>% 
      dplyr::select(-.estimator) %>% 
      mutate(rfmo = rfmos) %>% 
      pivot_wider(names_from = .metric, values_from = .estimate)
    
    output_fit <- list("class_fit" = class_fit, 
                       "reg_fit" = reg_fit, 
                       "test_predict" = test_predict, 
                       "test_metrics" = test_metrics,
                       "final_predict" = final_predict,
                       "final_metrics" = final_metrics)
    
    saveRDS(output_fit, paste0(save_loc, rfmos, "_ll_untuned_model.rds"))
    write.csv(final_predict, paste0(save_loc, rfmos, "_ll_untuned_final_predict.csv"), row.names = F)
    
    return(test_metrics)
    } 
  }
} 
