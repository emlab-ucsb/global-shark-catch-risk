## Run All ML models - one for each RFMO 

#' Function that runs a wide variety of model scenarios for the WCPFC region
#' @param data dataframe of catch values
#' @param save_loc output file loc to save r2
#' @param trees a numeric value for number of trees to use in the modeling framework
#' 

train_all_models <- function(data, save_loc){
  
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
      
      set.seed(1234)
      class_folds <- vfold_cv(train_class)
      
      class_model_tune <- rand_forest(trees = 1000, 
                                 mtry = tune(),  # tune the hyperparameters
                                 min_n = tune()) %>% 
        set_engine("ranger", num.threads = 8) %>% 
        set_mode("classification")
      
      class_grid  <- grid_regular(
        min_n(range = c(10, 30)),
        mtry(range = c(2, 8)),
        levels = 5
      )
      
      # Set Workflow
      class_wflow_tune <- 
        workflow() %>% 
        add_model(class_model_tune) %>% 
        add_recipe(class_recipe)
      
      # Run tuning
      cluster <- parallel::makePSOCKcluster(8)
      doParallel::registerDoParallel(cluster)
      
      set.seed(1234)
      tune_res_rf_class <-  class_wflow_tune %>% finetune::tune_race_anova(resamples = class_folds, grid = class_grid)
      
      best_tune_class <- tune_res_rf_class %>% 
        tune::select_best(metric = "roc_auc")
      
      class_tuned <- class_wflow_tune %>% 
        finalize_workflow(best_tune_class) %>% 
        fit(data = train_class)
      
      ###
      # Regression Model
      ###
      set.seed(1234)
      reg_folds <- vfold_cv(train_reg)
      
      reg_model_tune <- rand_forest(trees = 1000, 
                               mtry = tune(), # tune the hyperparameters 
                               min_n = tune()) %>% 
        set_engine("ranger", num.threads = 8) %>% 
        set_mode("regression") %>% 
        translate()
      
      reg_grid <- grid_regular(
        min_n(range = c(10, 30)),
        mtry(range = c(2, 8)),
        levels = 5
      )
      
      # Set Workflow
      reg_wflow_tune <- 
        workflow() %>% 
        add_model(reg_model_tune) %>% 
        add_recipe(reg_recipe)
      
      cluster <- parallel::makePSOCKcluster(8)
      doParallel::registerDoParallel(cluster)
      
      set.seed(1235)
      tune_res_rf_reg <- reg_wflow_tune %>% finetune::tune_race_anova(resamples = reg_folds, grid = reg_grid)
      
      best_tune_reg <- tune_res_rf_reg %>% 
        tune::select_best(metric = "rsq")
      
      reg_tuned <- reg_wflow_tune %>% 
        finalize_workflow(best_tune_reg) %>% 
        fit(data = train_reg)
      
      class_tuned_pred <- class_tuned %>% 
        predict(final_test)
      
      # Predict Final Testing
      final_class <- predict(class_tuned, final_test) %>% 
        mutate(.pred_class = ifelse(.pred_class == "absent", 0, 1)) %>% 
        bind_cols(final_test)
      
      final_reg <- predict(reg_tuned, final_test) %>% 
        bind_cols(final_test)
      
      final_predict <- final_class %>% 
        left_join(final_reg) %>% 
        mutate(.final_pred = .pred_class*.pred) 
      
      final_metrics <- metrics(truth = catch, estimate = .final_pred, data = final_predict)
      
      final_metrics <- final_metrics %>% 
        dplyr::select(-.estimator) 
      
      # Save outputs
      output_fit <- list("class_tuned", class_tuned, 
                         "reg_tuned" = reg_tuned, 
                         "class_predict" = final_class, 
                         "reg_predict" = final_reg, 
                         "final_predict" = final_predict, 
                         "metrics" = final_metrics)
      
      writeRDS(output_fit, paste0(save_loc, rfmo, "_trained_model.rds"))
    }
  }
} 
