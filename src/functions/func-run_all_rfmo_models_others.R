## Run rfmospecific model scenario - Environmentals and Other Options

#' Function that runs a wide variety of model scenarios for the WCPFC region
#' @param data dataframe of catch values
#' @param save_loc output file loc to save r2
#' @param rfmos which rfmo to choose
#' @param effort_source which effort source to choose

run_all_rfmo_models_others <- function(data, save_loc, rfmos, effort_source){
  
  # Save final metrics output
  final_metrics_total <- NULL
  
  # Filter to only include relevant data
  prep_ll <- data %>% 
    filter(gear_group == "longline" & 
             catch_units == "count" & 
             effort_units == "hooks") %>% 
    mutate_if(is_character, as.factor) %>% 
    arrange(year, latitude, longitude) %>% 
    mutate(pres_abs = factor(ifelse(catch > 0, "present", "absent")),
           zone = paste(latitude, longitude, sep = "|"))  
  
  if(effort_source == "effort reported with target catch (flag)") { 
    prep_ll <- prep_ll %>% 
      filter(rfmo == rfmos) %>% 
      dplyr::select_if(!grepl("bycatch_total_effort", colnames(.)))
    
    prep_ll <- prep_ll %>% 
    dplyr::select(-colnames(prep_ll %>% # remove columns where target effort or catch sums to 0
                              dplyr::select_if(grepl("target_effort|target_catch", colnames(.))) %>% 
                              dplyr::select_if(colSums(., na.rm = TRUE) == 0))) %>% 
    filter(!is.na(target_effort))
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
      filter(rfmo == rfmos) %>% 
      filter(!is.na(total_fishing_kwh)) 
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
    
    # Run series of models using several combinations of parameters
    # 1) using sst mean only & chl-a mean only, using sst mean & chl-a mean & ssh mean, using sst mean and cv & chl-a mean and cv, using sst mean and cv & chl-a mean and cv & ssh mean and cv
    # 2) excluding price data, using price data (species groups), using price data (species)
    # 3) log transform the catch variable
    
    for(enviro_options in c("meancv", "mean")) { # use just the mean or mean and cv
      for(ssh_options in c(TRUE, FALSE)) { # use SSH data? 
        for(price_options in c("noprice", "speciesprice", "groupprice")) { # include price options
          for(catch_options in c("none", "log")) { # do any transforming to the catch variable
            # Grab data
            train_class <- train_data_class_orig %>% 
              dplyr::select(pres_abs, sdm, species_commonname, mean_sst, mean_chla, mean_ssh, cv_sst, cv_chla, cv_ssh, 
                            median_price_group, median_price_species, latitude, longitude,
                            colnames(train_data_class_orig)[grepl("target_effort", colnames(train_data_class_orig)) &
                                                              grepl("longline", colnames(train_data_class_orig))]) %>% 
              distinct_all()
            
            train_reg <- train_data_reg_orig %>% 
              dplyr::select(catch, sdm, species_commonname, mean_sst, mean_chla, mean_ssh, cv_sst, cv_chla, cv_ssh, latitude, longitude,
                            median_price_group, median_price_species, colnames(train_data_reg_orig)[grepl("target_effort", colnames(train_data_reg_orig)) &
                                                                                                        grepl("longline", colnames(train_data_reg_orig))]) %>% 
              distinct_all()
            
            final_test <- test_data_class_orig %>% 
              dplyr::select(pres_abs, catch, sdm, species_commonname, mean_sst, mean_chla, mean_ssh, cv_sst, cv_chla, cv_ssh, 
                            median_price_group, median_price_species, 
                            colnames(train_data_class_orig)[grepl("target_effort", colnames(train_data_class_orig)) & grepl("longline", colnames(train_data_class_orig))], 
                            latitude, longitude, year) %>% 
              distinct_all()
            
            if(ssh_options == FALSE) { 
              train_class <- train_class %>% 
                dplyr::select(-mean_ssh, -cv_ssh) %>% 
                distinct_all()
              
              train_reg <- train_reg %>% 
                dplyr::select(-mean_ssh, -cv_ssh) %>% 
                distinct_all()
              
              final_test <- final_test %>% 
                dplyr::select(-mean_ssh, -cv_ssh) %>% 
                distinct_all()
            }
            
            if(enviro_options == "mean" & ssh_options == TRUE) { 
              train_class <- train_class %>% 
                dplyr::select(-cv_sst, -cv_chla, -cv_ssh) %>% 
                distinct_all()
              
              train_reg <- train_reg %>% 
                dplyr::select(-cv_sst, -cv_chla, -cv_ssh) %>% 
                distinct_all()
              
              final_test <- final_test %>% 
                dplyr::select(-cv_sst, -cv_chla, -cv_ssh) %>% 
                distinct_all()
            }
            
            if(enviro_options == "mean" & ssh_options  == FALSE) { 
              train_class <- train_class %>% 
                dplyr::select(-cv_sst, -cv_chla) %>% 
                distinct_all()
              
              train_reg <- train_reg %>% 
                dplyr::select(-cv_sst, -cv_chla) %>% 
                distinct_all()
              
              final_test <- final_test %>% 
                dplyr::select(-cv_sst, -cv_chla) %>% 
                distinct_all()
            }
            
            if(price_options == "noprice") { 
              train_class <- train_class %>% 
                dplyr::select(-median_price_species, -median_price_group) %>% 
                distinct_all()
              
              train_reg <- train_reg %>% 
                dplyr::select(-median_price_species, -median_price_group) %>% 
                distinct_all()
              
              final_test <- final_test %>% 
                dplyr::select(-median_price_species, -median_price_group) %>% 
                distinct_all()
            }
            
            if(price_options == "speciesprice") { 
              train_class <- train_class %>% 
                dplyr::select(-median_price_group) %>% 
                distinct_all()
              
              train_reg <- train_reg %>% 
                dplyr::select(-median_price_group) %>% 
                distinct_all()
              
              final_test <- final_test %>% 
                dplyr::select(-median_price_group) %>% 
                distinct_all()            
            }
            
            if(price_options == "groupprice") { 
              train_class <- train_class %>% 
                dplyr::select(-median_price_species) %>% 
                distinct_all()
              
              train_reg <- train_reg %>% 
                dplyr::select(-median_price_species) %>% 
                distinct_all()
              
              final_test <- final_test %>% 
                dplyr::select(-median_price_species) %>% 
                distinct_all()              
            }
            
            # Start Recipes
            class_recipe <- recipe(pres_abs ~ ., data = train_class) %>%
              themis::step_upsample(pres_abs, over_ratio = 1) %>% 
              step_center(all_numeric(), -all_outcomes(), - latitude, - longitude) %>% 
              step_scale(all_numeric(), -all_outcomes(), - latitude, - longitude) %>%
              step_unknown(all_nominal(), -all_outcomes()) %>%
              step_zv(all_predictors(), - latitude, - longitude) %>%
              step_impute_knn(all_numeric(), -all_outcomes(), neighbors = 3, impute_with = imp_vars(latitude, longitude)) %>%
              update_role(latitude, longitude, new_role = "none") %>% 
              step_dummy(all_nominal(), -all_outcomes())
            
            reg_recipe <- recipe(catch ~ ., data = train_reg) %>%
              step_center(all_numeric(), -all_outcomes(), - latitude, - longitude) %>% 
              step_scale(all_numeric(), -all_outcomes(), - latitude, - longitude) %>%
              step_unknown(all_nominal(), -all_outcomes()) %>%
              step_zv(all_predictors(), - latitude, - longitude) %>%
              step_impute_knn(all_numeric(), -all_outcomes(), neighbors = 3, impute_with = imp_vars(latitude, longitude)) %>%
              update_role(latitude, longitude, new_role = "none") %>% 
              step_dummy(all_nominal(), -all_outcomes()) 
            
            if(catch_options == "log") { 
              reg_recipe <- reg_recipe %>% 
                step_log(catch, offset = 1, skip = TRUE)
            }
            
            ###
            # Classification Model
            ###
            ntrees = 500 
            
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
            
            final_metrics <- metrics(truth = catch, estimate = .final_pred, data = final_predict)
            
            final_metrics <- final_metrics %>% 
              dplyr::select(-.estimator) %>% 
              mutate(environmental_value = enviro_options,
                     include_ssh = ssh_options, 
                     price = price_options,
                     catch_transformation = catch_options) %>% 
              pivot_wider(names_from = .metric, values_from = .estimate)
            
            final_metrics_total <- final_metrics_total %>% 
              bind_rows(final_metrics)
          }
        }
      }
    }
  }
  } 
  write.csv(final_metrics_total, paste0(save_loc, paste0(rfmos,"_models_others_results.csv")), row.names = FALSE)
} 
