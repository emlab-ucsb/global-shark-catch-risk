## Run All ML models - testing different efforts (reported with target species by flag, reported with bycatch by flag, reported by GFW) for each RFMO (without tuning) 

#' Function that runs a wide variety of model scenarios for the WCPFC region
#' @param data_gfw dataframe effort by gfw
#' @param data_effort dataframe of catch values and effort by flag
#' @param save_loc output file loc to save r2
#' 

all_rfmo_effort_models <- function(data_gfw, data_effort, save_loc){
  
  # Fix data
  data_gfw <- data_gfw %>% 
    filter(gear_group == "longline") %>% 
    mutate_if(is_character, as.factor) %>% 
    arrange(year, latitude, longitude) %>% 
    fill(mean_sst:cv_chla, .direction = "downup") %>% 
    fill(sdm, .direction = "downup") %>% 
    mutate(pres_abs = factor(ifelse(catch > 0, "present", "absent")),
           zone = paste(latitude, longitude, sep = "|"))  
  
  data_effort <- data_effort %>%
    filter(gear_group == "longline") %>% 
    mutate_if(is_character, as.factor) %>% 
    arrange(year, latitude, longitude) %>% 
    fill(mean_sst:cv_chla, .direction = "downup") %>% 
    fill(sdm, .direction = "downup") %>% 
    mutate(pres_abs = factor(ifelse(catch > 0, "present", "absent")),
           zone = paste(latitude, longitude, sep = "|"))  
  
  # Save final metrics output
  final_metrics_total <- NULL
  
  for(rfmos in unique(c(data_gfw$rfmo, data_effort$rfmo))) { 
    
    for(effort_source in c("effort reported with target catch (flag)", "effort reported with bycatch (flag)", "gfw effort")) { 
      
      prep_ll <- data.frame() 
      
      if(effort_source == "effort reported with target catch (flag)") {
        
        # Subset by RFMO
        prep_ll <- data_effort %>% 
          filter(rfmo == rfmos) %>% 
          dplyr::select_if(!grepl("bycatch_total_effort", colnames(.)))
          
        # Remove columns where all target effort sums to 0
        prep_ll <- prep_ll %>% 
            dplyr::select(-colnames(prep_ll %>% # remove columns where target effort sums to 0
                                      dplyr::select_if(grepl("target_effort", colnames(.))) %>% 
                                      dplyr::select_if(colSums(., na.rm = T) == 0))) 
      }
      
      if(effort_source == "effort reported with bycatch (flag)") {
        
        # Subset by RFMO
        prep_ll <- data_effort %>% 
          filter(rfmo == rfmos) %>% 
          dplyr::select_if(!grepl("target_effort", colnames(.)))
        
        # Remove columns where all bycatch associated effort sums to 0
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
        prep_ll <- data_gfw %>% 
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
      
      # Run series of models: 
      # 1) self-reported effort associated with target catch by flag
      # 2) self-reported effort associated with bycatch by flag
      # 3) effort reported by gfw
    
      model_runs <- data.frame(model = c("model 1", "model 2", "model 3"), 
                               parameters = c("effort reported with target catch (flag)", "effort reported with bycatch (flag)", "gfw effort"))
      
      models = model_runs$model[which(model_runs$parameters == effort_source)]
        
        if(models == "model 1") { # effort reported with target catch (flag)
          train_class <- train_data_class_orig %>%  
            dplyr::select(pres_abs, sdm, species_commonname, mean_sst, mean_chla, 
                          colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))]) %>% 
            distinct_all() 
          
          train_reg <- train_data_reg_orig %>%  
            dplyr::select(catch, sdm, species_commonname, mean_sst, mean_chla, 
                          colnames(train_data_reg_orig)[grepl("target_effort_", colnames(train_data_reg_orig))]) %>% 
            distinct_all() 
          
          final_test <- test_data_class_orig %>% 
            dplyr::select(pres_abs, catch, sdm, species_commonname, mean_sst, mean_chla, 
                          colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))]) %>% 
            distinct_all() 
        }
        
        if(models == "model 2") { # effort reported with bycatch (flag)
          train_class <- train_data_class_orig %>%  
            dplyr::select(pres_abs, sdm, species_commonname, mean_sst, mean_chla, 
                          colnames(train_data_class_orig)[grepl("bycatch_total_effort_", colnames(train_data_class_orig))]) %>% 
            distinct_all() 
          
          train_reg <- train_data_reg_orig %>%  
            dplyr::select(catch, sdm, species_commonname, mean_sst, mean_chla, 
                          colnames(train_data_reg_orig)[grepl("bycatch_total_effort_", colnames(train_data_reg_orig))]) %>% 
            distinct_all() 
          
          final_test <- test_data_class_orig %>% 
            dplyr::select(pres_abs, catch, sdm, species_commonname, mean_sst, mean_chla, 
                          colnames(train_data_class_orig)[grepl("bycatch_total_effort_", colnames(train_data_class_orig))]) %>% 
            distinct_all() 
          
        }
        
        if(models == "model 3") { # gfw_effort
          train_class <- train_data_class_orig %>%  
            dplyr::select(pres_abs, sdm, species_commonname, mean_sst, mean_chla, total_fishing_kwh) %>% 
            distinct_all() 
          
          train_reg <- train_data_reg_orig %>%  
            dplyr::select(catch, sdm, species_commonname, mean_sst, mean_chla, total_fishing_kwh) %>% 
            distinct_all()
          
          final_test <- test_data_class_orig %>% 
            dplyr::select(pres_abs, catch, sdm, species_commonname, mean_sst, 
                          mean_chla, total_fishing_kwh, latitude, longitude, year) %>% 
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
                 effort_source = effort_source, 
                 rfmo = rfmos) %>% 
          pivot_wider(names_from = .metric, values_from = .estimate)
        
        final_metrics_total <- final_metrics_total %>% 
          bind_rows(final_metrics)
    } 
    }
    }
  }
  write.csv(final_metrics_total, paste0(save_loc, "all_rfmos_effort_results.csv"), row.names = FALSE)
  
  return(final_metrics_total)
} 
