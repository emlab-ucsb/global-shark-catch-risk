## Run All ML models - one for each RFMO 

#' Function that runs a wide variety of model scenarios for the WCPFC region
#' @param data dataframe of catch values
#' @param save_loc output file loc to save r2
#' @param trees a numeric value for number of trees to use in the modeling framework
#' 

train_all_models <- function(data, save_loc, trees = 1000){
  
  for(rfmos in unique(data$rfmo)) { 
    
    # Prep data
    prep_ll <- data %>% 
      filter(rfmos == rfmo) %>% 
      mutate_if(is_character, as.factor) %>% 
      arrange(year, latitude, longitude) %>% 
      fill(mean_sst:cv_chla, .direction = "downup") %>% 
      fill(sdm, .direction = "downup") %>% 
      mutate(pres_abs = factor(ifelse(catch > 0, "present", "absent")),
             zone = paste(latitude, longitude, sep = "|"))  %>% 
      filter(!is.na(total_fishing_kwh)) 
    
    if(nrow(prep_ll) == 0) { 
      next 
    }
    
    # Rename for consistency
    prep_ll <- prep_ll %>% 
      rename(source_effort = total_fishing_kwh)
    
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
    
    # Use clustering to group locations into 5 groups
    xy$clust <- cutree(hc, k = 5)
    
    unique_pts$location_cluster <- xy$clust
    
    # Add groups back to data
    prep_ll <- prep_ll %>% 
      left_join(unique_pts)
    
    # Clean column names 
    prep_ll <- janitor::clean_names(prep_ll)
    
    ###
    # Run Model
    ###
    
    # Split data into testing and training
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
        dplyr::select(pres_abs, sdm, species_commonname, mean_sst, mean_chla, cv_sst, cv_chla, source_effort) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm, species_commonname, mean_sst, mean_chla, cv_sst, cv_chla, source_effort) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm, species_commonname, mean_sst, mean_chla, cv_sst, cv_chla, source_effort, latitude, longitude, year) %>% 
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
      ntrees = trees 
      
      class_model_tune <- rand_forest(trees = ntrees, 
                                 mtry = tune(),  # tune the hyperparameters
                                 min_n = tune()) %>% 
        set_engine("ranger", 
                   importance = "impurity") %>% 
        set_mode("classification") %>% 
        translate()
      
      # Set Workflow
      class_wflow_tune <- 
        workflow() %>% 
        add_model(class_model_tune) %>% 
        add_recipe(class_recipe)
      
      # Cross validation - https://juliasilge.com/blog/sf-trees-random-tuning/
      ## EB stopped here
      set.seed(1234)
      class_folds <- vfold_cv(train_class)
      
      doParallel::registerDoParallel()
      set.seed(1234)
      tune_class <- tune_grid(class_wflow_tune, 
                              resamples = class_folds, 
                              grid = 25)
      
      
      
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
      reg_model <- rand_forest(trees = ntrees, 
                               mtry = tune(), # tune the hyperparameters 
                               min_n = tune()) %>% 
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
      
      # Predict Final Testing
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
        dplyr::select(-.estimator) 
      
      # Save outputs
      output_fit <- list("class_fit", class_fit, 
                         "reg_fit" = reg_fit, 
                         "class_predict" = final_class, 
                         "reg_predict" = final_reg, 
                         "final_predict" = final_predict, 
                         "metrics" = final_metrics)
      
      writeRDS(output_fit, paste0(save_loc, rfmo, "_trained_model.rds"))
    }
  }
} 
