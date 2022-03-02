## Train/test split for ML

#' @param data dataframe to split
#' @param loc_group column name of location cluster
#' @param seed set seed so that we can reproduce results
#' @param prop set the proportion of data that you'd like to keep for the testing dataset

# Split data function (for now)
train_test_split <- function(data, loc_group, seed = 1234, prop = 0.25) { 
  
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
