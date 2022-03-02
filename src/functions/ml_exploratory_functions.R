# Functions for ML models

# Presence/Absence T-tests for explanatory variables across gear types and species groups
# data = dataset name
# variable = predictor variable to test
# title = label for predictor variable that is being tested; can be changed if column specified is not clear enough
pres_abs_test <- function(data, variable, gear = "both", title = variable) { 
  
  if(gear == "both") { 
  # Longlines (general)
  l1 <- t.test(data %>% 
                 filter(gear == "longline" & pres_abs == "present") %>% 
                 select(all_of(variable)), 
               data %>% 
                 filter(gear == "longline" & pres_abs == "absent") %>% 
                 select(all_of(variable)))
  
  # Purse seines (general)
  p1 <- t.test(data %>% 
                 filter(gear == "purse-seine" & pres_abs == "present") %>% 
                 select(all_of(variable)), 
               data %>% 
                 filter(gear == "purse-seine" & pres_abs == "absent") %>% 
                 select(all_of(variable)))
  
  # Save results
  result <- data.frame("gear_type" = c("longline", "purse-seine"),
                       "presence_estimate" = c(l1$estimate[1], p1$estimate[1]), 
                       "absence_estimate" = c(l1$estimate[2], p1$estimate[2]), 
                       "p" = c(l1$p.value, p1$p.value)) %>% 
    mutate(p = round(p/0.001)*0.001)
  
  } else { 
    t1 <- t.test(data %>% 
                   filter(gear == gear & pres_abs == "present") %>% 
                   select(all_of(variable)), 
                 data %>% 
                   filter(gear == gear & pres_abs == "absent") %>% 
                   select(all_of(variable)))
    result <- data.frame("gear_type" = gear, 
                         "presence_estimate" = t1$estimate[1], 
                         "absence_estimate" = t1$estimate[2], 
                         "p" = t1$p.value) %>% 
      mutate(p = round(p/0.001)*0.001)
    
    }
  
  # Return results table
  return(result)
  
} 

nonzero_lm <- function(data, variable, gear = "both", title = variable) { 
  # Rename the dataframe var
  data <- data %>% 
    rename(test = all_of(variable))
  
  if( gear == "both") { 
  # Longline (general)
  l1 <- tidy(lm(catch ~ test, data = data %>% filter(gear == "longline" & pres_abs == "present")))
  
  # All Purse Seine
  p1 <- tidy(lm(catch ~ test, data = data %>% filter(gear == "purse-seine" & pres_abs == "present")))
  
  # Save results
  result <- data.frame("gear_type" = c("longline", "purse-seine"),
                       # "predictor_variable" = title,
                       "intercept_estimate" = c(l1$estimate[1], p1$estimate[1]), 
                       "predictor_estimate" = c(l1$estimate[2], p1$estimate[2]), 
                       "predictor_p" = c(l1$p.value[2], p1$p.value[2])) %>% 
    mutate(predictor_p = round(predictor_p/0.001)*0.001)
  
  } else { 
    t1 <- tidy(lm(catch ~ test, data = data %>% filter(gear == gear & pres_abs == "present")))
    
    # Save results
    result <- data.frame("gear_type" = gear,
                         "intercept_estimate" = t1$estimate[1],
                         "predictor_estimate" = t1$estimate[2],
                         "predictor_p" = t1$p.value[2]) %>% 
      mutate(predictor_p = round(predictor_p/0.001)*0.001)
    
    }
  
  # Return results table
  return(result)
} 

# Split data function (for now)
split_it <- function(data, loc_group, seed, prop = 0.25) { 
  
  # Make sure location group is properly defined
  data <- data %>% 
    rename(location_cluster = loc_group)
  
  # Figure out how many rows we need to get the split proportion we desire
  nrow_break <- round(nrow(data)*prop)
  
  # Set seed and sample location group
  set.seed(seed)
  loc_sample <- sample(unique(data$location_cluster), 100, replace = TRUE)
  
  # Set seed and sample year group
  set.seed(seed+1)
  year_sample <- sample(unique(data$year), 100, replace = TRUE)
  
  # Create dataframe with unique combinations (so we don't get duplicate data)
  samples <- data.frame("location" = loc_sample, "year" = year_sample) %>% distinct_all()
  
  # Start adding to test dataset
  test_data <- NULL
  for(i in 1:nrow(samples)) { # For each sample row...
    # Filter the dataset
    test_data_temp <- data %>% filter(location_cluster == samples$location[i] & year == samples$year[i])
    
    # Add to the testing data list and remove duplicates
    test_data <- test_data %>% 
      bind_rows(test_data_temp) %>% 
      distinct_all()
    
    # If we exceed our breakpoint, put some back... 
    if(nrow(test_data) > nrow_break) { 
      test_data <- test_data[c(0:nrow_break),]
      break() # and break out of the loop
    } # end for breakpoint
  } # end for each sample row
  
  # The training data is whatever is not in the test dataset
  train_data <- anti_join(data, test_data, by = colnames(data))
  
  # And we return both datasets
  return(list("training" = train_data, "testing" = test_data))
}
