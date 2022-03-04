## Run WCPFC-specific model scenario - Effort Option

#' Function that runs a wide variety of model scenarios for the WCPFC region
#' @param data dataframe of catch values
#' @param save_loc output file loc to save r2

run_wcpfc_models_effort <- function(data, save_loc){
  
  final_metrics_total <- NULL
  
  data <- data %>% 
    filter(gear == "longline" & 
             catch_units == "count" & 
             effort_units == "hooks") %>% 
    mutate_if(is_character, as.factor) %>% 
    arrange(year, latitude, longitude) %>% 
    fill(mean_sst:cv_ssh, .direction = "downup") %>% 
    mutate(pres_abs = factor(ifelse(catch > 0, "present", "absent")),
           zone = paste(latitude, longitude, sep = "|"))  
  
  for(effort_source in c("self-reported catch and effort by month", 
                         "self-reported catch and effort by flag + quarter", 
                         "self-reported catch and effort by flag + year", 
                         "gfw")){ 
    
    ###
    # Select the data and to the spatial splitting for train/test
    ###
    
    if(effort_source == "gfw") { 
      prep_ll <- data %>%
        mutate(target_effort = total_fishing_kwh) %>% 
        filter(!is.na(target_effort))
    } else { 
      prep_ll <- data %>% 
        filter(target_sources == effort_source)
      
      prep_ll <- prep_ll %>% 
        dplyr::select(-colnames(prep_ll %>% # remove columns where target effort or catch sums to 0
                                  dplyr::select_if(grepl("target_effort|target_catch", colnames(.))) %>% 
                                  dplyr::select_if(colSums(.) == 0))) %>% 
        filter(!is.na(target_effort))
    }
    
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
      
      # Run series of models: 
      # 1) bycatch sdms
      # 2) target effort
      # 3) bycatch sdms, species name, sst mean, chla mean, total effort
      # 4) target catch in mt (only data from flagquarter or flagyear)
      # 5) bycatch sdms, species name, sst mean, chla mean, total effort by country (only data from flagquarter or flagyear)
      # 6) bycatch sdms, species name, sst mean, chla mean, total effort by country, target catch by species in mt (only data from flagquarter or flagyear)
      model_runs <- data.frame(model = c("model 1", "model 2", "model 3", "model 4", "model 5", "model 6"), 
                               parameters = c("bycatch sdms", 
                                              "target effort",
                                              "bycatch sdms, species name, sst mean, chla mean, total effort", 
                                              "target catch in mt (only data from flagquarter or flagyear)", 
                                              "bycatch sdms, species name, sst mean, chla mean, total effort by country  (only data from flagquarter or flagyear)", 
                                              "bycatch sdms, species name, sst mean, chla mean, total effort by country, target catch by species in mt (only data from flagquarter or flagyear)"), 
                               effort_sources = c(rep("self-reported catch and effort by month", 3), 
                                                  rep(NA, 3), 
                                                  rep("self-reported catch and effort by flag + quarter", 6), 
                                                  rep("self-reported catch and effort by flag + year", 6), 
                                                  rep("gfw", 3), rep(NA, 3)))
      
      models_to_run <- model_runs %>% 
        filter(effort_sources == effort_source) 
      
      for(models in models_to_run$model) { 
        
        if(models == "model 1") { # bycatch sdms
          train_class <- train_data_class_orig %>%  
            dplyr::select(pres_abs, sdm_probability, species_commonname) %>% 
            distinct_all() %>% 
            dplyr::select(-species_commonname)
          
          train_reg <- train_data_reg_orig %>%  
            dplyr::select(catch, sdm_probability, species_commonname) %>% 
            distinct_all() %>% 
            dplyr::select(-species_commonname)
          
          final_test <- test_data_class_orig %>% 
            dplyr::select(pres_abs, catch, sdm_probability, species_commonname, latitude, longitude, year) %>% 
            distinct_all() 
        }
        
        if(models == "model 2") { # target effort 
          train_class <- train_data_class_orig %>%  
            dplyr::select(pres_abs, target_effort, species_commonname) %>% 
            distinct_all() %>% 
            dplyr::select(-species_commonname)
          
          train_reg <- train_data_reg_orig %>%  
            dplyr::select(catch, target_effort, species_commonname) %>% 
            distinct_all() %>% 
            dplyr::select(-species_commonname)
          
          final_test <- test_data_class_orig %>% 
            dplyr::select(pres_abs, catch, target_effort, species_commonname, latitude, longitude, year) %>% 
            distinct_all() 
          
        }
        
        if(models == "model 3") { # bycatch sdms, species name, sst mean, chla mean, total effort
          train_class <- train_data_class_orig %>%  
            dplyr::select(pres_abs, target_effort, sdm_probability, species_commonname, mean_sst, mean_chla) %>% 
            distinct_all() 
          
          train_reg <- train_data_reg_orig %>%  
            dplyr::select(catch, target_effort, sdm_probability, species_commonname, mean_sst, mean_chla) %>% 
            distinct_all()
          
          final_test <- test_data_class_orig %>% 
            dplyr::select(pres_abs, catch, target_effort, sdm_probability, species_commonname, mean_sst, 
                          mean_chla, latitude, longitude, year) %>% 
            distinct_all()
          
        }
        
        if(models == "model 4") { # target catch in mt (only data from flagquarter or flagyear)
          train_class <- train_data_class_orig %>%  
            dplyr::select(pres_abs, target_catch_tonnes, species_commonname) %>% 
            distinct_all() %>% 
            dplyr::select(-species_commonname)
          
          train_reg <- train_data_reg_orig %>%  
            dplyr::select(catch, target_catch_tonnes, species_commonname) %>% 
            distinct_all() %>% 
            dplyr::select(-species_commonname)
          
          final_test <- test_data_class_orig %>% 
            dplyr::select(pres_abs, catch, target_catch_tonnes, species_commonname, latitude, longitude, year) %>% 
            distinct_all() 
          
        }
        
        if(models == "model 5") { # bycatch sdms, species name, sst mean, chla mean, total effort by country (only data from flagquarter or flagyear)
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
            dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, latitude, longitude, year, 
                          colnames(train_data_class_orig)[grepl("target_effort", colnames(train_data_class_orig)) &
                                                            grepl("longline", colnames(train_data_class_orig))]) %>% 
            distinct_all() 
          
        }
        
        if(models == "model 6") { # bycatch sdms, species name, sst mean, chla mean, total effort by country, target catch by species in mt (only data from flagquarter or flagyear)
          train_class <- train_data_class_orig %>%  
            dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, 
                          colnames(train_data_class_orig)[grepl("target_effort", colnames(train_data_class_orig)) &
                                                            grepl("longline", colnames(train_data_class_orig))], 
                          colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) &
                                                            grepl("tonnes", colnames(train_data_class_orig)) & 
                                                            grepl("longline", colnames(train_data_class_orig)) & 
                                                            !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
            distinct_all() 
          
          train_reg <- train_data_reg_orig %>%  
            dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, 
                          colnames(train_data_reg_orig)[grepl("target_effort", colnames(train_data_reg_orig)) &
                                                          grepl("longline", colnames(train_data_reg_orig))], 
                          colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) &
                                                          grepl("tonnes", colnames(train_data_reg_orig)) & 
                                                          grepl("longline", colnames(train_data_reg_orig)) & 
                                                          !grepl("target_catch_tonnes", colnames(train_data_reg_orig))]) %>% 
            distinct_all() 
          
          final_test <- test_data_class_orig %>% 
            dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, latitude, longitude, year, 
                          colnames(train_data_class_orig)[grepl("target_effort", colnames(train_data_class_orig)) &
                                                            grepl("longline", colnames(train_data_class_orig))], 
                          colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) &
                                                            grepl("tonnes", colnames(train_data_class_orig)) & 
                                                            grepl("longline", colnames(train_data_class_orig)) & 
                                                            !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
            distinct_all() 
          
        } 
        
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
        
        final_metrics <- metrics(truth = catch, estimate = .final_pred, data = final_predict)
        
        final_metrics <- final_metrics %>% 
          dplyr::select(-.estimator) %>% 
          mutate(model = models, 
                 parameters = models_to_run$parameters[which(models_to_run$model == models)], 
                 effort_source = effort_source) %>% 
          pivot_wider(names_from = .metric, values_from = .estimate)
        
        final_metrics_total <- final_metrics_total %>% 
          bind_rows(final_metrics)
      }
    }
  }
  
  write.csv(final_metrics_total, paste0(save_loc, "wcpfc_models_effort_results.csv"), row.names = FALSE)
}