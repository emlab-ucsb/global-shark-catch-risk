# prep data
prep_rf_data <- function(data, geartype, source) { 
  if(source == "gfw") {
    # Grab data and clean up a bit
    data <- data %>% 
      filter(gear == geartype) %>% 
      mutate_if(is_character, as.factor) %>% 
      arrange(year, latitude, longitude) %>% 
      fill(mean_sst:cv_chla) %>% 
      mutate(pres_abs = factor(ifelse(catch > 0, "present", "absent")), 
             tot_kwh = day_total_fishing_kwh + night_total_fishing_kwh, 
             zone = paste(latitude, longitude, sep = "|"))  
    
    data <- data %>% 
      # remove columns of all zeros (any gt groups, countries, etc)
      dplyr::select(-colnames(data %>% 
                         dplyr::select_if(grepl("gt", colnames(.)) ) %>% 
                         dplyr::select_if(colSums(., na.rm = T) == 0))) %>% 
      filter(!is.na(tot_kwh))
    
    # Give consistent names
    data <- data %>% 
      mutate(source_effort = tot_kwh)
  }
  if(source == "observer") {
    # Grab data and clean up a bit
    data <- data %>% 
      filter(gear == geartype) %>% 
      dplyr::select(year:observer_coverage, median_price, median_tuna_price, median_tunalike_price) %>%
      mutate_if(is_character, as.factor) %>% 
      arrange(year, latitude, longitude) %>% 
      fill(mean_sst:cv_chla) %>% 
      mutate(pres_abs = factor(ifelse(catch > 0, "present", "absent")), 
             zone = paste(latitude, longitude, sep = "|")) %>% 
      filter(!is.na(effort)) 
    
    # Give consistent names
    data <- data %>% 
      mutate(source_effort = effort)
    
  }
  if(source %in% c("flagquarter", "flagyear", "month")) {
    # Grab data and clean up a bit
    data <- data %>% 
      filter(gear == geartype) %>% 
      mutate_if(is_character, as.factor) %>% 
      arrange(year, latitude, longitude) %>% 
      fill(mean_sst:cv_chla) %>% 
      mutate(pres_abs = factor(ifelse(catch > 0, "present", "absent")), 
             zone = paste(latitude, longitude, sep = "|")) 
    
    data <- data %>% 
      # remove columns of all zeros
      dplyr::select(-colnames(data %>% 
                         dplyr::select_if(grepl("target_effort_|target_catch_", colnames(.)) ) %>% 
                         dplyr::select_if(colSums(.) == 0))) %>% 
      filter(!is.na(target_effort))
    
    # Give consistent names
    data <- data %>% 
      mutate(source_effort = target_effort)
  }
  
    # Add spatial groups
    unique_pts <- data %>% 
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
    
    # Use clustering to group locations into 6 groups
    xy$clust <- cutree(hc, k = 6)
    
    unique_pts$location_cluster <- xy$clust
    
    # Add groups back to data
    data <- data %>% 
      left_join(unique_pts)
    
    # Clean column names 
    data <- janitor::clean_names(data)
    
    # If GFW data, remove NAs
    #if(source %in% c("gfw", "observer")) { 
     # data <- na.omit(data) 
    #} 
    
    # Create different datasets for marine mammals, sharks and rays, and turtles
    # Marine Mammals
    mammals <- data %>% 
      filter(species_group == "marine mammals") %>%
      group_by(species_commonname) %>% 
      mutate(scaled_catch = round(scale(catch) %>% as.vector()),
             scaled_source_effort = scale(source_effort) %>% as.vector(), 
             scaled_sdm = scale(sdm_probability) %>% as.vector()) %>%
      filter(!is.na(scaled_catch)) %>% 
      mutate(scaled_catch = ifelse(scaled_catch == -1, 0, scaled_catch), 
             species_commonname = droplevels(species_commonname)) %>% 
      ungroup()
    
    # Sharks and Rays
    sharks <- data %>% 
      filter(species_group == "sharks and rays") %>%
      group_by(species_commonname) %>% 
      mutate(scaled_catch = round(scale(catch) %>% as.vector()),
             scaled_source_effort = scale(source_effort) %>% as.vector(), 
             scaled_sdm = scale(sdm_probability) %>% as.vector()) %>% 
      filter(!is.na(scaled_catch)) %>% 
      mutate(scaled_catch = ifelse(scaled_catch == -1, 0, scaled_catch), 
             species_commonname = droplevels(species_commonname)) %>% 
      ungroup()
    
    # Turtles
    turtles <- data %>% 
      filter(species_group == "turtles") %>%
      group_by(species_commonname) %>% 
      mutate(scaled_catch = round(scale(catch) %>% as.vector()),
             scaled_source_effort = scale(source_effort) %>% as.vector(), 
             scaled_sdm = scale(sdm_probability) %>% as.vector()) %>% 
      filter(!is.na(scaled_catch)) %>% 
      mutate(scaled_catch = ifelse(scaled_catch == -1, 0, scaled_catch), 
             species_commonname = droplevels(species_commonname)) %>% 
      ungroup()
    
    # Tunas
    tunas <- data %>% 
      filter(species_group == "tunas") %>%
      group_by(species_commonname) %>% 
      mutate(scaled_catch = round(scale(catch) %>% as.vector()),
             scaled_source_effort = scale(source_effort) %>% as.vector(), 
             scaled_sdm = scale(sdm_probability) %>% as.vector()) %>% 
      filter(!is.na(scaled_catch)) %>% 
      mutate(scaled_catch = ifelse(scaled_catch == -1, 0, scaled_catch), 
             species_commonname = droplevels(species_commonname)) %>% 
      ungroup()
    
    # Tuna-Like-Species
    tunalike <- data %>% 
      filter(species_group == "tuna-like species") %>%
      group_by(species_commonname) %>% 
      mutate(scaled_catch = round(scale(catch) %>% as.vector()),
             scaled_source_effort = scale(source_effort) %>% as.vector(), 
             scaled_sdm = scale(sdm_probability) %>% as.vector()) %>% 
      filter(!is.na(scaled_catch)) %>% 
      mutate(scaled_catch = ifelse(scaled_catch == -1, 0, scaled_catch), 
             species_commonname = droplevels(species_commonname)) %>% 
      ungroup()
    
    return(list(mammals, sharks, turtles, tunas, tunalike))
  
}

prep_rf_data_final <- function(data, geartype, source) { 
  if(source == "gfw") {
    # Grab data and clean up a bit
    data <- data %>% 
      filter(gear == geartype) %>% 
      mutate_if(is_character, as.factor) %>% 
      arrange(year, latitude, longitude) %>% 
      fill(mean_sst:cv_chla) %>% 
      mutate(pres_abs = factor(ifelse(catch > 0, "present", "absent")), 
             tot_kwh = day_total_fishing_kwh + night_total_fishing_kwh)  
    
    data <- data %>% 
      # remove columns of all zeros (any gt groups, countries, etc)
      dplyr::select(-colnames(data %>% 
                                dplyr::select_if(grepl("gt", colnames(.)) ) %>% 
                                dplyr::select_if(colSums(., na.rm = T) == 0))) 
    
    # Give consistent names
    data <- data %>% 
      mutate(source_effort = tot_kwh)
  }
  if(source == "observer") {
    # Grab data and clean up a bit
    data <- data %>% 
      filter(gear == geartype) %>% 
      dplyr::select(year:observer_coverage, median_price, median_tuna_price, median_tunalike_price) %>%
      mutate_if(is_character, as.factor) %>% 
      arrange(year, latitude, longitude) %>% 
      fill(mean_sst:cv_chla) %>% 
      mutate(pres_abs = factor(ifelse(catch > 0, "present", "absent")))
    
    # Give consistent names
    data <- data %>% 
      mutate(source_effort = effort)
    
  }
  if(source %in% c("flagquarter", "flagyear", "month")) {
    # Grab data and clean up a bit
    data <- data %>% 
      filter(gear == geartype) %>% 
      mutate_if(is_character, as.factor) %>% 
      arrange(year, latitude, longitude) %>% 
      fill(mean_sst:cv_chla) %>% 
      mutate(pres_abs = factor(ifelse(catch > 0, "present", "absent"))) 
    
    data <- data %>% 
      # remove columns of all zeros
      dplyr::select(-colnames(data %>% 
                                dplyr::select_if(grepl("target_effort_|target_catch_", colnames(.)) ) %>% 
                                dplyr::select_if(colSums(., na.rm = T) == 0))) 
    
    # Give consistent names
    data <- data %>% 
      mutate(source_effort = target_effort)
  }
  
  # Clean column names 
  data <- janitor::clean_names(data)
  
  # Create different datasets for marine mammals, sharks and rays, and turtles
  # Marine Mammals
  mammals <- data %>% 
    filter(species_group == "marine mammals")
  
  # Sharks and Rays
  sharks <- data %>% 
    filter(species_group == "sharks and rays") 
  
  # Turtles
  turtles <- data %>% 
    filter(species_group == "turtles") 
  
  # Tunas
  tunas <- data %>% 
    filter(species_group == "tunas") 
  
  # Tuna-Like-Species
  tunalike <- data %>% 
    filter(species_group == "tuna-like species") 
  
  return(list(mammals, sharks, turtles, tunas, tunalike))
  
}
 

# Run models

run_rf_models <- function(data, model_run, source, include_price = FALSE, ntrees = 100, 
                          save_model = model_run) { 
  
  # Assign training and testing data
  
  train_data_class_orig <- data$training
  
  test_data_class_orig <- data$testing
  
  train_data_reg_orig <- data$training %>% 
    filter(pres_abs == "present")
  
  test_data_reg_orig <- data$testing %>% 
    filter(pres_abs == "present")
  
  # SDM models should be the same for all sources... 
  if(model_run == "null-sdm") { 
    
    predictors <- "sdm"
    
    # Grab data
    train_class <- train_data_class_orig %>% 
      dplyr::select(pres_abs, sdm_probability, species_commonname) %>% 
      distinct_all()
    
    train_reg <- train_data_reg_orig %>% 
      dplyr::select(catch, sdm_probability, species_commonname) %>% 
      distinct_all()
    
    final_test <- test_data_class_orig %>% 
      dplyr::select(pres_abs, catch, sdm_probability, species_commonname, latitude, longitude, year) %>% 
      distinct_all()
    
    # Start Recipes
    class_recipe <- recipe(pres_abs ~ sdm_probability, data = train_class) %>%
      themis::step_upsample(pres_abs, over_ratio = 1) %>% 
      step_center(all_numeric(), -all_outcomes()) %>% 
      step_scale(all_numeric(), -all_outcomes()) %>%
      step_zv(all_predictors())
    
    reg_recipe <- recipe(catch ~ sdm_probability, data = train_reg) %>% 
      step_center(all_numeric(), -all_outcomes()) %>% 
      step_scale(all_numeric(), -all_outcomes()) %>%
      step_zv(all_predictors())
  }
  
  # Different parameters for the GFW data... 
  if(include_price == FALSE) { 
  if(source == "gfw") { 
    if(model_run == "null-effort") { 
        
      predictors <- "tot kwh"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, tot_kwh, species_commonname) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, tot_kwh, species_commonname) %>% 
        distinct_all()
        
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, tot_kwh, species_commonname, latitude, longitude, year) %>% 
        distinct_all()
      
      # Start Recipes
      class_recipe <- recipe(pres_abs ~ tot_kwh, data = train_class) %>%
        themis::step_upsample(pres_abs, over_ratio = 1) %>% 
        step_center(all_numeric(), -all_outcomes()) %>% 
        step_scale(all_numeric(), -all_outcomes()) %>%
        step_zv(all_predictors())
      
      reg_recipe <- recipe(catch ~ tot_kwh, data = train_reg) %>% 
        step_center(all_numeric(), -all_outcomes()) %>% 
        step_scale(all_numeric(), -all_outcomes()) %>%
        step_zv(all_predictors())
    }
    if(model_run == "rf1") { 
      
      predictors <- "sdm, species name, sst cv, chla cv, tot kwh"
      
      # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, tot_kwh) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, tot_kwh) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, tot_kwh, latitude, longitude, year) %>% 
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
    }
    if(model_run == "rf2") { 
      
      predictors <- "sdm, species name, sst cv, chla cv, tot kwh (by country type, day/night, and gross tonnage group)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, 
               colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
        dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, 
               colnames(train_data_reg_orig)[grepl("fishing_kwh", colnames(train_data_reg_orig))]) %>% 
        dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, latitude, longitude, year, 
               colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
        dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
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
    }
    if(model_run == "rf3") { 
      
      predictors <- "sdm, species name, sst mean, chla mean, tot kwh"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, tot_kwh) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, tot_kwh) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, tot_kwh, latitude, longitude, year) %>% 
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
    }
    if(model_run == "rf4") { 
      
      predictors <- "sdm, species name, sst median, chla median, tot kwh"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, tot_kwh) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, tot_kwh) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, tot_kwh, latitude, longitude, year) %>% 
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
    }
    if(model_run == "rf5") { 
      
      predictors <- "sdm, species name, sst mean, chla mean, tot kwh (by country type, day/night, and gross tonnage group)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, 
               colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
        dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, 
               colnames(train_data_reg_orig)[grepl("fishing_kwh", colnames(train_data_reg_orig))]) %>% 
        dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, latitude, longitude, year,
               colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
        dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
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
    }
    if(model_run == "rf6") { 
      
      predictors <- "sdm, species name, sst median, chla median, tot kwh (by country type, day/night, and gross tonnage group)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, 
               colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
        dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, 
               colnames(train_data_reg_orig)[grepl("fishing_kwh", colnames(train_data_reg_orig))]) %>% 
        dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, latitude, longitude, year,
               colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
        dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
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
    }
  }
  if(source == "observer") { 
    if(model_run == "null-effort") { 
      
      predictors <- "effort (longines = hooks, purse seine = sets)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, effort, species_commonname) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, effort, species_commonname) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, effort, species_commonname, latitude, longitude, year) %>% 
        distinct_all()
      
      # Start Recipes
      class_recipe <- recipe(pres_abs ~ effort, data = train_class) %>%
        themis::step_upsample(pres_abs, over_ratio = 1) %>% 
        step_center(all_numeric(), -all_outcomes()) %>% 
        step_scale(all_numeric(), -all_outcomes()) %>%
        step_zv(all_predictors())
      
      reg_recipe <- recipe(catch ~ effort, data = train_reg) %>% 
        step_center(all_numeric(), -all_outcomes()) %>% 
        step_scale(all_numeric(), -all_outcomes()) %>%
        step_zv(all_predictors())
    }
    
    if(model_run == "rf1") { 
      
      predictors <- "sdm, species name, sst cv, chla cv, effort (longines = hooks, purse seine = sets)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, effort) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, effort) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, effort, latitude, longitude, year) %>% 
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
    }
    if(model_run == "rf3") { 
      
      predictors <- "sdm, species name, sst mean, chla mean, effort (longines = hooks, purse seine = sets)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, effort) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, effort) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, effort, latitude, longitude, year) %>% 
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
    }
    if(model_run == "rf4") { 
      
      predictors <- "sdm, species name, sst median, chla median, effort (longines = hooks, purse seine = sets)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, effort) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, effort) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, effort, latitude, longitude, year) %>% 
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
    }
  }
  if(source %in% c("flagquarter", "flagyear")) { 
    if(model_run == "null-effort") { 
      
      predictors <- "target effort (longines = hooks, purse seine = sets)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, target_effort, species_commonname) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>%  
        dplyr::select(catch, target_effort, species_commonname) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, target_effort, species_commonname, latitude, longitude, year) %>% 
        distinct_all()
      
      # Start Recipes
      class_recipe <- recipe(pres_abs ~ target_effort, data = train_class) %>%
        themis::step_upsample(pres_abs, over_ratio = 1) %>% 
        step_center(all_numeric(), -all_outcomes()) %>% 
        step_scale(all_numeric(), -all_outcomes()) %>%
        step_zv(all_predictors())
      
      reg_recipe <- recipe(catch ~ target_effort, data = train_reg) %>% 
        step_center(all_numeric(), -all_outcomes()) %>% 
        step_scale(all_numeric(), -all_outcomes()) %>%
        step_zv(all_predictors())
    }    
    if(model_run == "null-targetcatch") { 
      
      predictors <- "target catch"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        mutate(total_target_catch = rowSums(train_data_class_orig[,colnames(train_data_class_orig)[grepl("target_catch_tonnes", colnames(train_data_class_orig))]])) %>% 
        dplyr::select(pres_abs, total_target_catch, species_commonname) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        mutate(total_target_catch = rowSums(train_data_reg_orig[,colnames(train_data_reg_orig)[grepl("target_catch_tonnes", colnames(train_data_reg_orig))]])) %>% 
        dplyr::select(catch, total_target_catch, species_commonname) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        mutate(total_target_catch = rowSums(test_data_class_orig[,colnames(test_data_class_orig)[grepl("target_catch_tonnes", colnames(test_data_class_orig))]])) %>% 
        dplyr::select(pres_abs, catch, total_target_catch, species_commonname, latitude, longitude, year) %>% 
        distinct_all()
      
      # Start Recipes
      class_recipe <- recipe(pres_abs ~ total_target_catch, data = train_class) %>%
        themis::step_upsample(pres_abs, over_ratio = 1) %>% 
        step_center(all_numeric(), -all_outcomes()) %>% 
        step_scale(all_numeric(), -all_outcomes()) %>%
        step_zv(all_predictors())
      
      reg_recipe <- recipe(catch ~ total_target_catch, data = train_reg) %>% 
        step_center(all_numeric(), -all_outcomes()) %>% 
        step_scale(all_numeric(), -all_outcomes()) %>%
        step_zv(all_predictors())
    }
    if(model_run == "rf1") { 
      
      predictors <- "sdm, species name, sst cv, chla cv, target effort (longines = hooks, purse seine = sets)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort)  %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort, latitude, longitude, year) %>% 
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
    }
    if(model_run == "rf2") { 
      
      predictors <- "sdm, species name, sst cv, chla cv, target effort (longines = hooks, purse seine = sets) by country type"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, 
               colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))]) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, 
               colnames(train_data_reg_orig)[grepl("target_effort_", colnames(train_data_reg_orig))]) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, latitude, longitude, year,
               colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))])  %>% 
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
    }
    if(model_run == "rf3") { 
      
      predictors <- "sdm, species name, sst mean, chla mean, target effort (longines = hooks, purse seine = sets)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort, latitude, longitude, year) %>% 
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
    }
    if(model_run == "rf4") { 
      
      predictors <- "sdm, species name, sst median, chla median, target_effort (longines = hooks, purse seine = sets)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, target_effort) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, target_effort) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, target_effort, latitude, longitude, year) %>% 
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
    }
    if(model_run == "rf5") { 
      
      predictors <- "sdm, species name, sst mean, chla mean, target effort (longines = hooks, purse seine = sets) by country type"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, 
               colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))]) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, 
               colnames(train_data_reg_orig)[grepl("target_effort_", colnames(train_data_reg_orig))]) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, latitude, longitude, year,
               colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))]) %>% 
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
    }
    if(model_run == "rf6") { 
      
      predictors <- "sdm, species name, sst median, chla median, target effort (longines = hooks, purse seine = sets) by country type"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, 
               colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))]) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, 
               colnames(train_data_reg_orig)[grepl("target_effort_", colnames(train_data_reg_orig))]) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, latitude, longitude, year,
               colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))])  %>% 
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
    }
    if(model_run == "rf7") { 
      
      predictors <- "sdm, species name, sst cv, chla cv, target effort (longines = hooks, purse seine = sets) by country type, target catch (metric tonnes)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, 
               colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
               colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                            grepl("tonnes", colnames(train_data_class_orig)) & 
                                            !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, 
               colnames(train_data_reg_orig)[grepl("target_effort_", colnames(train_data_reg_orig))], 
               colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                                 grepl("tonnes", colnames(train_data_reg_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_reg_orig))]) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, latitude, longitude, year,
               colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
               colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                 grepl("tonnes", colnames(train_data_class_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
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
    }
    if(model_run == "rf8") { 
      
      predictors <- "sdm, species name, sst mean, chla mean, target effort (longines = hooks, purse seine = sets) by country type, target catch (metric tonnes)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, 
               colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
               colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                 grepl("tonnes", colnames(train_data_class_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, 
               colnames(train_data_reg_orig)[grepl("target_effort_", colnames(train_data_reg_orig))], 
               colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                               grepl("tonnes", colnames(train_data_reg_orig)) & 
                                               !grepl("target_catch_tonnes", colnames(train_data_reg_orig))]) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, latitude, longitude, year,
               colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
               colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                 grepl("tonnes", colnames(train_data_class_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
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
    }
    if(model_run == "rf9") { 
      
      predictors <- "sdm, species name, sst median, chla median, target effort (longines = hooks, purse seine = sets) by country type, target catch (metric tonnes)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, 
               colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
               colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                 grepl("tonnes", colnames(train_data_class_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, 
               colnames(train_data_reg_orig)[grepl("target_effort_", colnames(train_data_reg_orig))], 
               colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                               grepl("tonnes", colnames(train_data_reg_orig)) & 
                                               !grepl("target_catch_tonnes", colnames(train_data_reg_orig))]) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, latitude, longitude, year,
               colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
               colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                 grepl("tonnes", colnames(train_data_class_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
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
    }
  }
  if(source == "month") { 
    if(model_run == "null-effort") { 
      
      predictors <- "target effort (longines = hooks, purse seine = sets)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, target_effort, species_commonname)  %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, target_effort, species_commonname) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, target_effort, species_commonname, latitude, longitude, year) %>% 
        distinct_all()
      
      # Start Recipes
      class_recipe <- recipe(pres_abs ~ target_effort, data = train_class) %>%
        themis::step_upsample(pres_abs, over_ratio = 1) %>% 
        step_center(all_numeric(), -all_outcomes()) %>% 
        step_scale(all_numeric(), -all_outcomes()) %>%
        step_zv(all_predictors())
      
      reg_recipe <- recipe(catch ~ target_effort, data = train_reg) %>% 
        step_center(all_numeric(), -all_outcomes()) %>% 
        step_scale(all_numeric(), -all_outcomes()) %>%
        step_zv(all_predictors())
    }    
    if(model_run == "null-targetcatch") { 
      
      predictors <- "target catch"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        mutate(total_target_catch = rowSums(train_data_class_orig[,colnames(train_data_class_orig)[grepl("target_catch_tonnes", colnames(train_data_class_orig))]])) %>% 
        dplyr::select(pres_abs, total_target_catch, species_commonname) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        mutate(total_target_catch = rowSums(train_data_reg_orig[,colnames(train_data_reg_orig)[grepl("target_catch_tonnes", colnames(train_data_reg_orig))]])) %>% 
        dplyr::select(catch, total_target_catch, species_commonname) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        mutate(total_target_catch = rowSums(test_data_class_orig[,colnames(test_data_class_orig)[grepl("target_catch_tonnes", colnames(test_data_class_orig))]])) %>% 
        dplyr::select(pres_abs, catch, total_target_catch, species_commonname, latitude, longitude, year) %>% 
        distinct_all()
      
      # Start Recipes
      class_recipe <- recipe(pres_abs ~ total_target_catch, data = train_class) %>%
        themis::step_upsample(pres_abs, over_ratio = 1) %>% 
        step_center(all_numeric(), -all_outcomes()) %>% 
        step_scale(all_numeric(), -all_outcomes()) %>%
        step_zv(all_predictors())
      
      reg_recipe <- recipe(catch ~ total_target_catch, data = train_reg) %>% 
        step_center(all_numeric(), -all_outcomes()) %>% 
        step_scale(all_numeric(), -all_outcomes()) %>%
        step_zv(all_predictors())
    }
    if(model_run == "rf1") { 
      
      predictors <- "sdm, species name, sst cv, chla cv, target effort (longines = hooks, purse seine = sets)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort, latitude, longitude, year) %>% 
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
    }
    if(model_run == "rf3") { 
      
      predictors <- "sdm, species name, sst mean, chla mean, target effort (longines = hooks, purse seine = sets)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort, latitude, longitude, year) %>% 
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
    }
    if(model_run == "rf4") { 
      
      predictors <- "sdm, species name, sst median, chla median, target_effort (longines = hooks, purse seine = sets)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, target_effort) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, target_effort) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, target_effort, latitude, longitude, year) %>% 
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
    }
    if(model_run == "rf7") { 
      
      predictors <- "sdm, species name, sst cv, chla cv, target effort (longines = hooks, purse seine = sets), target catch (metric tonnes)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort,
               colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                 grepl("tonnes", colnames(train_data_class_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort,
               colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                               grepl("tonnes", colnames(train_data_reg_orig)) & 
                                               !grepl("target_catch_tonnes", colnames(train_data_reg_orig))]) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort,
               latitude, longitude, year,
               colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                 grepl("tonnes", colnames(train_data_class_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
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
    }
    if(model_run == "rf8") { 
      
      predictors <- "sdm, species name, sst mean, chla mean, target effort (longines = hooks, purse seine = sets), target catch (metric tonnes)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort,
               colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                 grepl("tonnes", colnames(train_data_class_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort,
               colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                               grepl("tonnes", colnames(train_data_reg_orig)) & 
                                               !grepl("target_catch_tonnes", colnames(train_data_reg_orig))]) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort,
               latitude, longitude, year,
               colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                 grepl("tonnes", colnames(train_data_class_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
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
    }
    if(model_run == "rf9") { 
      
      predictors <- "sdm, species name, sst median, chla median, target effort (longines = hooks, purse seine = sets), target catch (metric tonnes)"
      
      # Grab data
      train_class <- train_data_class_orig %>% 
        dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, target_effort, 
               colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                 grepl("tonnes", colnames(train_data_class_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
        distinct_all()
      
      train_reg <- train_data_reg_orig %>% 
        dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, target_effort, 
               colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                               grepl("tonnes", colnames(train_data_reg_orig)) & 
                                               !grepl("target_catch_tonnes", colnames(train_data_reg_orig))]) %>% 
        distinct_all()
      
      final_test <- test_data_class_orig %>% 
        dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, target_effort, 
               latitude, longitude, year,
               colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                 grepl("tonnes", colnames(train_data_class_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
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
    }
  } } 
  
  # Different parameters for the GFW data... 
  if(include_price == TRUE) { 
    if(source == "gfw") { 
      if(model_run == "null-effort") { 
        
        predictors <- "tot kwh"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, tot_kwh, species_commonname) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, tot_kwh, species_commonname) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, tot_kwh, species_commonname, latitude, longitude, year) %>% 
          distinct_all()
        
        # Start Recipes
        class_recipe <- recipe(pres_abs ~ tot_kwh, data = train_class) %>%
          themis::step_upsample(pres_abs, over_ratio = 1) %>% 
          step_center(all_numeric(), -all_outcomes()) %>% 
          step_scale(all_numeric(), -all_outcomes()) %>%
          step_zv(all_predictors())
        
        reg_recipe <- recipe(catch ~ tot_kwh, data = train_reg) %>% 
          step_center(all_numeric(), -all_outcomes()) %>% 
          step_scale(all_numeric(), -all_outcomes()) %>%
          step_zv(all_predictors())
      }
      if(model_run == "rf1") { 
        
        predictors <- "sdm, species name, sst cv, chla cv, tot kwh, median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, tot_kwh, median_price) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, tot_kwh, median_price) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, median_price, 
                 tot_kwh, latitude, longitude, year) %>% 
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
      }
      if(model_run == "rf2") { 
        
        predictors <- "sdm, species name, sst cv, chla cv, tot kwh (by country type, day/night, and gross tonnage group), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, median_price,
                 colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, median_price,
                 colnames(train_data_reg_orig)[grepl("fishing_kwh", colnames(train_data_reg_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, latitude, longitude, year, median_price,
                 colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
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
      }
      if(model_run == "rf3") { 
        
        predictors <- "sdm, species name, sst mean, chla mean, tot kwh, median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, tot_kwh, median_price) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, tot_kwh, median_price) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, tot_kwh, median_price, latitude, longitude, year) %>% 
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
      }
      if(model_run == "rf4") { 
        
        predictors <- "sdm, species name, sst median, chla median, tot kwh, median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, tot_kwh, median_price) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, tot_kwh, median_price) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, tot_kwh, median_price, latitude, longitude, year) %>% 
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
      }
      if(model_run == "rf5") { 
        
        predictors <- "sdm, species name, sst mean, chla mean, tot kwh (by country type, day/night, and gross tonnage group), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, median_price, 
                 colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, median_price, 
                 colnames(train_data_reg_orig)[grepl("fishing_kwh", colnames(train_data_reg_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, median_price, latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
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
      }
      if(model_run == "rf6") { 
        
        predictors <- "sdm, species name, sst median, chla median, tot kwh (by country type, day/night, and gross tonnage group), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, median_price, 
                 colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, median_price, 
                 colnames(train_data_reg_orig)[grepl("fishing_kwh", colnames(train_data_reg_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, median_price, latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
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
      }
      if(model_run == "rf13") { 
        
        predictors <- "sdm, species name, sst median, chla median, tot kwh (by country type, day/night, and gross tonnage group), median bycatch price, median target price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, median_price, median_tuna_price, median_tunalike_price,
                 colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, median_price, median_tuna_price, median_tunalike_price,
                 colnames(train_data_reg_orig)[grepl("fishing_kwh", colnames(train_data_reg_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, median_price, median_tuna_price, median_tunalike_price, latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
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
      }
      if(model_run == "rf14") { 
        
        predictors <- "sdm, species name, sst mean, chla mean, tot kwh (by country type, day/night, and gross tonnage group), median bycatch price, median target price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, median_price, median_tuna_price, median_tunalike_price, 
                 colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, median_price, median_tuna_price, median_tunalike_price, 
                 colnames(train_data_reg_orig)[grepl("fishing_kwh", colnames(train_data_reg_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, median_price, median_tuna_price, median_tunalike_price, latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
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
      }
      if(model_run == "rf15") { 
        
        predictors <- "sdm, species name, sst cv, chla cv, tot kwh (by country type, day/night, and gross tonnage group), median bycatch price, median target price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, median_price, median_tuna_price, median_tunalike_price, 
                 colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, median_price, median_tuna_price, median_tunalike_price, 
                 colnames(train_data_reg_orig)[grepl("fishing_kwh", colnames(train_data_reg_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, median_price, median_tuna_price, median_tunalike_price, latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("fishing_kwh", colnames(train_data_class_orig))]) %>% 
          dplyr::select(-day_total_fishing_kwh, -night_total_fishing_kwh) %>% 
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
      }
    }
    if(source == "observer") { 
      if(model_run == "null-effort") { 
        
        predictors <- "effort (longines = hooks, purse seine = sets)"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, effort, species_commonname) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, effort, species_commonname) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>%  
          dplyr::select(pres_abs, catch, effort, species_commonname, latitude, longitude, year) %>% 
          distinct_all()
        
        # Start Recipes
        class_recipe <- recipe(pres_abs ~ effort, data = train_class) %>%
          themis::step_upsample(pres_abs, over_ratio = 1) %>% 
          step_center(all_numeric(), -all_outcomes()) %>% 
          step_scale(all_numeric(), -all_outcomes()) %>%
          step_zv(all_predictors())
        
        reg_recipe <- recipe(catch ~ effort, data = train_reg) %>% 
          step_center(all_numeric(), -all_outcomes()) %>% 
          step_scale(all_numeric(), -all_outcomes()) %>%
          step_zv(all_predictors())
      }
      
      if(model_run == "rf1") { 
        
        predictors <- "sdm, species name, sst cv, chla cv, effort (longines = hooks, purse seine = sets), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, effort, median_price) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, effort, median_price) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, effort, median_price, latitude, longitude, year) %>% 
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
      }
      if(model_run == "rf3") { 
        
        predictors <- "sdm, species name, sst mean, chla mean, effort (longines = hooks, purse seine = sets), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, effort, median_price) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, effort, median_price) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, effort, median_price, latitude, longitude, year) %>% 
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
      }
      if(model_run == "rf4") { 
        
        predictors <- "sdm, species name, sst median, chla median, effort (longines = hooks, purse seine = sets), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, effort, median_price) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, effort, median_price) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, effort, latitude, longitude, year, median_price) %>% 
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
      }
      if(model_run == "rf13") { 
        
        predictors <- "sdm, species name, sst median, chla median, effort (longines = hooks, purse seine = sets), median bycatch price, median target price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, effort, median_price, median_tuna_price, median_tunalike_price) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, effort, median_price, median_tuna_price, median_tunalike_price) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, effort, latitude, longitude, year, median_price, median_tuna_price, median_tunalike_price) %>% 
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
      }
      if(model_run == "rf14") { 
        
        predictors <- "sdm, species name, sst mean, chla mean, effort (longines = hooks, purse seine = sets), median bycatch price, median target price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, effort, median_price, median_tuna_price, median_tunalike_price) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, effort, median_price, median_tuna_price, median_tunalike_price) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, effort, latitude, longitude, year, median_price, median_tuna_price, median_tunalike_price) %>% 
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
      }
      if(model_run == "rf15") { 
        
        predictors <- "sdm, species name, sst cv, chla cv, effort (longines = hooks, purse seine = sets), median bycatch price, median target price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, effort, median_price, median_tuna_price, median_tunalike_price) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, effort, median_price, median_tuna_price, median_tunalike_price) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, effort, latitude, longitude, year, median_price, median_tuna_price, median_tunalike_price) %>% 
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
      }
    }
    if(source %in% c("flagquarter", "flagyear")) { 
      if(model_run == "null-effort") { 
        
        predictors <- "target effort (longines = hooks, purse seine = sets)"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, target_effort, species_commonname) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, target_effort, species_commonname) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, target_effort, species_commonname, latitude, longitude, year) %>% 
          distinct_all()
        
        # Start Recipes
        class_recipe <- recipe(pres_abs ~ target_effort, data = train_class) %>%
          themis::step_upsample(pres_abs, over_ratio = 1) %>% 
          step_center(all_numeric(), -all_outcomes()) %>% 
          step_scale(all_numeric(), -all_outcomes()) %>%
          step_zv(all_predictors())
        
        reg_recipe <- recipe(catch ~ target_effort, data = train_reg) %>% 
          step_center(all_numeric(), -all_outcomes()) %>% 
          step_scale(all_numeric(), -all_outcomes()) %>%
          step_zv(all_predictors())
      }    
      if(model_run == "null-targetcatch") { 
        
        predictors <- "target catch"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          mutate(total_target_catch = rowSums(train_data_class_orig[,colnames(train_data_class_orig)[grepl("target_catch_tonnes", colnames(train_data_class_orig))]])) %>% 
          dplyr::select(pres_abs, total_target_catch, species_commonname) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          mutate(total_target_catch = rowSums(train_data_reg_orig[,colnames(train_data_reg_orig)[grepl("target_catch_tonnes", colnames(train_data_reg_orig))]])) %>% 
          dplyr::select(catch, total_target_catch, species_commonname) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          mutate(total_target_catch = rowSums(test_data_class_orig[,colnames(test_data_class_orig)[grepl("target_catch_tonnes", colnames(test_data_class_orig))]])) %>% 
          dplyr::select(pres_abs, catch, total_target_catch, species_commonname, latitude, longitude, year) %>% 
          distinct_all()
        
        # Start Recipes
        class_recipe <- recipe(pres_abs ~ total_target_catch, data = train_class) %>%
          themis::step_upsample(pres_abs, over_ratio = 1) %>% 
          step_center(all_numeric(), -all_outcomes()) %>% 
          step_scale(all_numeric(), -all_outcomes()) %>%
          step_zv(all_predictors())
        
        reg_recipe <- recipe(catch ~ total_target_catch, data = train_reg) %>% 
          step_center(all_numeric(), -all_outcomes()) %>% 
          step_scale(all_numeric(), -all_outcomes()) %>%
          step_zv(all_predictors())
      }
      if(model_run == "rf1") { 
        
        predictors <- "sdm, species name, sst cv, chla cv, target effort (longines = hooks, purse seine = sets), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort, median_price) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort, median_price) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort, median_price, latitude, longitude, year) %>% 
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
      }
      if(model_run == "rf2") { 
        
        predictors <- "sdm, species name, sst cv, chla cv, target effort (longines = hooks, purse seine = sets) by country type, median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, median_price, 
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))]) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, median_price, 
                 colnames(train_data_reg_orig)[grepl("target_effort_", colnames(train_data_reg_orig))]) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, median_price, latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))])  %>% 
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
      }
      if(model_run == "rf3") { 
        
        predictors <- "sdm, species name, sst mean, chla mean, target effort (longines = hooks, purse seine = sets), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort, median_price) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort, median_price) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort, median_price, latitude, longitude, year) %>% 
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
      }
      if(model_run == "rf4") { 
        
        predictors <- "sdm, species name, sst median, chla median, target_effort (longines = hooks, purse seine = sets), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, target_effort, median_price) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, target_effort, median_price) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, target_effort, median_price, latitude, longitude, year) %>% 
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
      }
      if(model_run == "rf5") { 
        
        predictors <- "sdm, species name, sst mean, chla mean, target effort (longines = hooks, purse seine = sets) by country type, median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, median_price, 
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))]) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, median_price, 
                 colnames(train_data_reg_orig)[grepl("target_effort_", colnames(train_data_reg_orig))]) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, median_price, latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))]) %>% 
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
      }
      if(model_run == "rf6") { 
        
        predictors <- "sdm, species name, sst median, chla median, target effort (longines = hooks, purse seine = sets) by country type, median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, median_price, 
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))]) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, median_price, 
                 colnames(train_data_reg_orig)[grepl("target_effort_", colnames(train_data_reg_orig))]) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, median_price, latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))])  %>% 
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
      }
      if(model_run == "rf7") { 
        
        predictors <- "sdm, species name, sst cv, chla cv, target effort (longines = hooks, purse seine = sets) by country type, target catch (metric tonnes), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, median_price, 
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, median_price, 
                 colnames(train_data_reg_orig)[grepl("target_effort_", colnames(train_data_reg_orig))], 
                 colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                                 grepl("tonnes", colnames(train_data_reg_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_reg_orig))]) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, median_price, latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
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
      }
      if(model_run == "rf8") { 
        
        predictors <- "sdm, species name, sst mean, chla mean, target effort (longines = hooks, purse seine = sets) by country type, target catch (metric tonnes), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, median_price, 
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, median_price, 
                 colnames(train_data_reg_orig)[grepl("target_effort_", colnames(train_data_reg_orig))], 
                 colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                                 grepl("tonnes", colnames(train_data_reg_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_reg_orig))]) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, median_price, latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
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
      }
      if(model_run == "rf9") { 
        
        predictors <- "sdm, species name, sst median, chla median, target effort (longines = hooks, purse seine = sets) by country type, target catch (metric tonnes), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, median_price, 
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, median_price, 
                 colnames(train_data_reg_orig)[grepl("target_effort_", colnames(train_data_reg_orig))], 
                 colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                                 grepl("tonnes", colnames(train_data_reg_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_reg_orig))]) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, median_price, latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
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
      }
      if(model_run == "rf10") { 
        
        predictors <- "sdm, species name, sst median, chla median, target effort (longines = hooks, purse seine = sets) by country type, target catch (metric tonnes), median bycatch price, median target price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, median_price, median_tuna_price, median_tunalike_price, 
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, median_price, median_tuna_price, median_tunalike_price,
                 colnames(train_data_reg_orig)[grepl("target_effort_", colnames(train_data_reg_orig))], 
                 colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                                 grepl("tonnes", colnames(train_data_reg_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_reg_orig))]) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, median_price, median_tuna_price, median_tunalike_price, latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
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
      }
      if(model_run == "rf11") { 
        
        predictors <- "sdm, species name, sst mean, chla mean, target effort (longines = hooks, purse seine = sets) by country type, target catch (metric tonnes), median bycatch price, median target price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, median_price,  median_tuna_price, median_tunalike_price,
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, median_price,  median_tuna_price, median_tunalike_price,
                 colnames(train_data_reg_orig)[grepl("target_effort_", colnames(train_data_reg_orig))], 
                 colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                                 grepl("tonnes", colnames(train_data_reg_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_reg_orig))]) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, median_price, median_tuna_price, median_tunalike_price, latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
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
      }
      if(model_run == "rf12") { 
        
        predictors <- "sdm, species name, sst cv, chla cv, target effort (longines = hooks, purse seine = sets) by country type, target catch (metric tonnes), median bycatch price, median target price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, median_price, median_tuna_price, median_tunalike_price, 
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, median_price,  median_tuna_price, median_tunalike_price,
                 colnames(train_data_reg_orig)[grepl("target_effort_", colnames(train_data_reg_orig))], 
                 colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                                 grepl("tonnes", colnames(train_data_reg_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_reg_orig))]) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, median_price, median_tuna_price, median_tunalike_price, latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("target_effort_", colnames(train_data_class_orig))], 
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
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
      }
    }
    if(source == "month") { 
      if(model_run == "null-effort") { 
        
        predictors <- "target effort (longines = hooks, purse seine = sets)"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, target_effort, species_commonname) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, target_effort, species_commonname) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, target_effort, species_commonname, latitude, longitude, year) %>% 
          distinct_all()
        
        # Start Recipes
        class_recipe <- recipe(pres_abs ~ target_effort, data = train_class) %>%
          themis::step_upsample(pres_abs, over_ratio = 1) %>% 
          step_center(all_numeric(), -all_outcomes()) %>% 
          step_scale(all_numeric(), -all_outcomes()) %>%
          step_zv(all_predictors())
        
        reg_recipe <- recipe(catch ~ target_effort, data = train_reg) %>% 
          step_center(all_numeric(), -all_outcomes()) %>% 
          step_scale(all_numeric(), -all_outcomes()) %>%
          step_zv(all_predictors())
      }    
      if(model_run == "null-targetcatch") { 
        
        predictors <- "target catch"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          mutate(total_target_catch = rowSums(train_data_class_orig[,colnames(train_data_class_orig)[grepl("target_catch_tonnes", colnames(train_data_class_orig))]])) %>% 
          dplyr::select(pres_abs, total_target_catch, species_commonname)  %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          mutate(total_target_catch = rowSums(train_data_reg_orig[,colnames(train_data_reg_orig)[grepl("target_catch_tonnes", colnames(train_data_reg_orig))]])) %>% 
          dplyr::select(catch, total_target_catch, species_commonname)  %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          mutate(total_target_catch = rowSums(test_data_class_orig[,colnames(test_data_class_orig)[grepl("target_catch_tonnes", colnames(test_data_class_orig))]])) %>% 
          dplyr::select(pres_abs, catch, total_target_catch, species_commonname, latitude, longitude, year) %>% 
          distinct_all()
        
        # Start Recipes
        class_recipe <- recipe(pres_abs ~ total_target_catch, data = train_class) %>%
          themis::step_upsample(pres_abs, over_ratio = 1) %>% 
          step_center(all_numeric(), -all_outcomes()) %>% 
          step_scale(all_numeric(), -all_outcomes()) %>%
          step_zv(all_predictors())
        
        reg_recipe <- recipe(catch ~ total_target_catch, data = train_reg) %>% 
          step_center(all_numeric(), -all_outcomes()) %>% 
          step_scale(all_numeric(), -all_outcomes()) %>%
          step_zv(all_predictors())
      }
      if(model_run == "rf1") { 
        
        predictors <- "sdm, species name, sst cv, chla cv, target effort (longines = hooks, purse seine = sets), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort, median_price)  %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort, median_price) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort, median_price, latitude, longitude, year) %>% 
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
      }
      if(model_run == "rf3") { 
        
        predictors <- "sdm, species name, sst mean, chla mean, target effort (longines = hooks, purse seine = sets), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort, median_price) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort, median_price) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort, median_price, latitude, longitude, year) %>% 
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
      }
      if(model_run == "rf4") { 
        
        predictors <- "sdm, species name, sst median, chla median, target_effort (longines = hooks, purse seine = sets), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, target_effort, median_price) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, target_effort, median_price) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, target_effort, median_price, latitude, longitude, year) %>% 
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
      }
      if(model_run == "rf7") { 
        
        predictors <- "sdm, species name, sst cv, chla cv, target effort (longines = hooks, purse seine = sets), target catch (metric tonnes), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort, median_price,
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort, median_price,
                 colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                                 grepl("tonnes", colnames(train_data_reg_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_reg_orig))]) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort, median_price,
                 latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
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
      }
      if(model_run == "rf8") { 
        
        predictors <- "sdm, species name, sst mean, chla mean, target effort (longines = hooks, purse seine = sets), target catch (metric tonnes), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort, median_price,
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort, median_price,
                 colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                                 grepl("tonnes", colnames(train_data_reg_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_reg_orig))]) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort, median_price,
                 latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
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
      }
      if(model_run == "rf9") { 
        
        predictors <- "sdm, species name, sst median, chla median, target effort (longines = hooks, purse seine = sets), target catch (metric tonnes), median bycatch price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, target_effort, median_price, 
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
          distinct_all()
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, target_effort, median_price, 
                 colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                                 grepl("tonnes", colnames(train_data_reg_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_reg_orig))]) %>% 
          distinct_all()
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, target_effort, median_price, 
                 latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))]) %>% 
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
      }
      if(model_run == "rf10") { 
        
        predictors <- "sdm, species name, sst median, chla median, target effort (longines = hooks, purse seine = sets), target catch (metric tonnes), median bycatch price, median target price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, median_sst, median_chla, target_effort, median_price, median_tuna_price, median_tunalike_price,
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))])
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, median_sst, median_chla, target_effort, median_price, median_tuna_price, median_tunalike_price,
                 colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                                 grepl("tonnes", colnames(train_data_reg_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_reg_orig))])
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, median_sst, median_chla, target_effort, median_price, median_tuna_price, median_tunalike_price,
                 latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))])
        
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
      }
      if(model_run == "rf11") { 
        
        predictors <- "sdm, species name, sst mean, chla mean, target effort (longines = hooks, purse seine = sets), target catch (metric tonnes), median bycatch price, median target price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort, median_price, median_tuna_price, median_tunalike_price,
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))])
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort, median_price, median_tuna_price, median_tunalike_price,
                 colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                                 grepl("tonnes", colnames(train_data_reg_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_reg_orig))])
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, mean_sst, mean_chla, target_effort, median_price, median_tuna_price, median_tunalike_price,
                 latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))])
        
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
      }
      if(model_run == "rf12") { 
        
        predictors <- "sdm, species name, sst cv, chla cv, target effort (longines = hooks, purse seine = sets), target catch (metric tonnes), median bycatch price, median target price"
        
        # Grab data
        train_class <- train_data_class_orig %>% 
          dplyr::select(pres_abs, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort, median_price, median_tuna_price, median_tunalike_price,
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))])
        
        train_reg <- train_data_reg_orig %>% 
          dplyr::select(catch, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort, median_price, median_tuna_price, median_tunalike_price,
                 colnames(train_data_reg_orig)[grepl("target_catch", colnames(train_data_reg_orig)) & 
                                                 grepl("tonnes", colnames(train_data_reg_orig)) & 
                                                 !grepl("target_catch_tonnes", colnames(train_data_reg_orig))])
        
        final_test <- test_data_class_orig %>% 
          dplyr::select(pres_abs, catch, sdm_probability, species_commonname, cv_sst, cv_chla, target_effort, median_price, median_tuna_price, median_tunalike_price,
                 latitude, longitude, year,
                 colnames(train_data_class_orig)[grepl("target_catch", colnames(train_data_class_orig)) & 
                                                   grepl("tonnes", colnames(train_data_class_orig)) & 
                                                   !grepl("target_catch_tonnes", colnames(train_data_class_orig))])
        
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
      }
    } } 
  
  # Classification
  # Set Model
  class_model <- rand_forest(trees = ntrees) %>% 
    set_engine("ranger") %>% 
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
  
  # Save fitted model
  write_rds(class_fit, paste0(save_model, "class.rds"))  
  
  # Regression
  # Set Model
  reg_model <- rand_forest(trees = ntrees) %>% 
    set_engine("ranger") %>% 
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
  
  # Save fitted model
  write_rds(reg_fit, paste0(save_model, "reg.rds")) 
  
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
    dplyr::select(-.estimator) %>% 
    mutate(model = model_run, 
           predictor = predictors) %>% 
    pivot_wider(names_from = .metric, values_from = .estimate)
  
  return(list(final_metrics, final_predict))
  
}
  
# function to run the multivariate rf (please ensure effort predictor to be used is called "effort")
run_multi_rf <- function(data, seed = 1234, split_prop = 0.25, n_tree = 100, m_feature = 1, min_leaf = 1) { 
  
  # Organize data
  multi_data <- data %>% 
    dplyr::select(gear, latitude, longitude, effort, year, median_tuna_price, median_tunalike_price,
                  marine_mammal_sdm, sharks_sdm, tunalike_sdm, tunas_sdm, turtles_sdm, cv_chla, cv_sst, target_catch_count_tunas, 
                  target_catch_count_tunalike) %>% 
    distinct_all() %>% 
    left_join(data %>% 
                dplyr::select(latitude, longitude, year, species_group, median_price) %>%
                distinct_all() %>% 
                pivot_wider(names_from = species_group, values_from = median_price) %>% 
                rename(median_marinemammal_price = `marine mammals`, 
                       median_sharks_price = `sharks and rays`, 
                       median_turtles_price = `turtles`), 
              by = c("latitude", "longitude", "year")) %>% 
    left_join(data %>% 
                dplyr::select(species_group, latitude, longitude, year, catch, median_price) %>% 
                distinct_all() %>% 
                group_by(species_group, latitude, longitude, year) %>% 
                summarise(total_bycatch = sum(catch)) %>% 
                ungroup() %>% 
                pivot_wider(names_from = species_group, values_from = total_bycatch) %>% 
                rename(marine_mammal_catch = `marine mammals`, 
                       sharks_catch = `sharks and rays`, 
                       turtles_catch = turtles), 
              by = c("year", "longitude", "latitude")) %>% 
    # below may change depending on what we're calculating
    mutate(scaled_sharks_catch = sharks_catch / (sharks_catch + target_catch_count_tunas + target_catch_count_tunalike), 
           scaled_tuna_catch = target_catch_count_tunas / (sharks_catch + target_catch_count_tunas + target_catch_count_tunalike), 
           scaled_tunalike_catch = target_catch_count_tunalike / (sharks_catch + target_catch_count_tunas +
                                                                    target_catch_count_tunalike), 
           zone = paste(latitude, longitude, sep = "|")) %>% 
    filter(year < 2018) # price data only goes up to 2017
  
  # Add spatial groups
  unique_pts <- multi_data %>% 
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
  
  # Use clustering to group locations into 6 groups
  xy$clust <- cutree(hc, k = 6)
  
  unique_pts$location_cluster <- xy$clust
  
  # Add groups back to data
  multi_data <- multi_data %>% 
    left_join(unique_pts)
  
  # Split data for model
  
  for(geartype in c("longline", "purse-seine")) { 
    
    split_gear <- multi_data %>% filter(gear == geartype)
    
    if(nrow(split_gear) > 1 &
       length(unique(split_gear$target_catch_count_tunas)) > 1) { 
    
    split_data <- split_it(data = split_gear, 
                           loc_group = "location_cluster", seed = seed, prop = split_prop)
    
    # Build training and testing datasets
    train_x <- split_data$training %>% 
      dplyr::select(effort, median_sharks_price, median_tuna_price, median_tunalike_price, sharks_sdm, tunas_sdm, tunalike_sdm,
                    cv_chla, cv_sst)
    
    train_y <- split_data$training %>% 
      dplyr::select(scaled_sharks_catch, scaled_tuna_catch)
    
    test_x <- split_data$testing %>% 
      dplyr::select(effort, median_sharks_price, median_tuna_price, median_tunalike_price, sharks_sdm, tunas_sdm, tunalike_sdm,
                    cv_chla, cv_sst)
    
    test_y <- split_data$testing %>% 
      dplyr::select(scaled_sharks_catch, scaled_tuna_catch)
    
    # Run multivariate random forest prediction
    prediction <- build_forest_predict(trainX = as.matrix(train_x), trainY = as.matrix(train_y), n_tree = n_tree, 
                                       m_feature = m_feature, min_leaf = min_leaf, testX = as.matrix(test_x))
    
    prediction <- as.data.frame(prediction)
    colnames(prediction)  <- c("estimate_sharks", "estimate_tunas") 
    
    # Add prediction back into the testing/training dataset
    assign(paste0("testing_output_", gsub("-", "", geartype)), 
           split_data$testing %>% 
             mutate(predicted_sharks = prediction$estimate_sharks, 
                    predicted_tunas = prediction$estimate_tunas, 
                    predicted_tunalike = 1 - (prediction$estimate_sharks + prediction$estimate_tunas)))
    
    # Evaluate how well the model did
    eval_data <- data.frame("truth_sharks" = test_y$scaled_sharks_catch, 
                            "truth_tunas" = test_y$scaled_tuna_catch, 
                            "truth_tunalike" = (split_data$training)$scaled_tunalike_catch, 
                            "estimate_sharks" = prediction$estimate_sharks, 
                            "estimate_tunas" = prediction$estimate_tunas, 
                            "estimate_tunalike" = 1 - (prediction$estimate_sharks + prediction$estimate_tunas))
    
    # I think this works by adding all three? Then we can see how the model works as a whole?
    assign(paste0("evaluation_", gsub("-", "", geartype)),
           data.frame(
             "rmse" = (rmse(data = eval_data, truth = c(truth_sharks, truth_tunas, truth_tunalike), 
                            estimate = c(estimate_sharks, estimate_tunas, estimate_tunalike)))$.estimate,
             "rsq" = (rsq(data = eval_data, truth = c(truth_sharks, truth_tunas, truth_tunalike), 
                          estimate = c(estimate_sharks, estimate_tunas, estimate_tunalike)))$.estimate ))
    
    } else { 
      assign(paste0("testing_output_", gsub("-", "", geartype)), NULL)
      assign(paste0("evaluation_", gsub("-", "", geartype)),
             NULL)
      }
  } 
    
  # Provide the output and the rmse
  return_list <- list(testing_output_longline, evaluation_longline, testing_output_purseseine, evaluation_purseseine)
  names(return_list) <- c("testing_output_longline", "evaluation_longline", 
                          "testing_output_purseseine", "evaluation_purseseine")
  return(return_list)
  
  }
  
  