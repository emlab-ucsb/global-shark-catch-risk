# Housing for functions for converting rfmo and gfw data into molleweide projection and resampling

##############################
#### rfmo_raster_blaster() ### 
##############################

# DESCRIPTION: Takes rfmo dataset with point data, converts to polygons, transforms to mollweide 
# projection, rasterizes, resamples, and produces raster object

# VARIABLES:
# data : variable name that contains dataset
# raster_layer : variable name that contains raster layer for masking
# species : character value - all capitals of species name
# geartype : character value - all lowercase of gear_group
# unit: character value - all lowercase of catch units
# time_res : character value - all lowercase of time_res ("yearly", "quarterly", "monthly")
# mm in format of c(month value, ...) in case you wanted to filter for multiple
# yy in format of c(year value, ...) in case you wanted to filter for multiple
# flag_code in format of c(flag code, ...) in case you wanted to filter for multiple
# include_generated : TRUE/FALSE - whether to include generated data
# spatial_res : interger value - spatial resolution of grid cell (only works if x and y have identical resolution)
#   ex: "center of 5x5 cell" in spatial_notes means spatial_res == 5
# summary_var : character value - which column to use for raster creation
# summary_func : character value - which function to use for raster collating if cells are overlapping

rfmo_raster_blaster <- function(data, raster_layer, species, geartype, unit, time_res, mm = NA, yy = NA,
                           flag_code = NA, include_generated = FALSE, spatial_res, summary_var = "catch", 
                           summary_func = "sum") { 
  
  # Pre-filter for all things apart from spatial resolution
  data <- data %>% 
    filter(species_sciname == species & 
             gear_group == geartype & 
             time_period == time_res)
  
  # Check
  if(nrow(data) == 0) { 
    stop("Specifications too constricted, cannot rasterize - try different species, geartype, time_res or remove some specifications")
  }
  
  # Pre filter for units
  if(summary_var == "effort" ) { 
    data <- data %>% filter(effort_units == unit)
  }
  if(summary_var == "catch" | summary_var == "rel_catch") { 
    data <- data %>% filter(grepl(unit, catch_units))
  }
  
  # Continue pre-filter but only pre-filter if values are not NAs
  # Month
  if(!is.na(mm)) { 
    data <- data %>% filter(month %in% mm)
    if(nrow(data) == 0) { 
      stop("Specifications too constricted, cannot rasterize - try different month or remove month specification")
    }
  }
  
  # Year
  if(!is.na(yy)) { 
    data <- data %>% filter(year %in% yy)
    if(nrow(data) == 0) { 
      stop("Specifications too constricted, cannot rasterize - try different year or remove year specification")
    }
  }
  
  # Flag
  if(!is.na(flag_code)) { 
    data <- data %>% filter(country %in% flag_code)
    if(nrow(data) == 0) { 
      stop("Specifications too constricted, cannot rasterize - try different flag_code or remove flag specification")
    }
  }
  
  # Generated Data
  if(include_generated == FALSE) { 
    data <- data %>% filter(was_generated == "no")
    if(nrow(data) == 0) { 
      stop("Specifications too constricted, cannot rasterize - try allowing for generated data")
    }
  }
  
  # Spatial Check
  if(nrow(data[grepl(spatial_res, data$spatial_notes),]) == 0) { 
    stop("Specifications too constricted, cannot rasterize - try different spatial resolution or changing other parameters")
  }
  
  # Convert to polygons for rasterize/fasterize
  # First step - figuring out spatial resolution of the points
  data_polygons <- data %>%
    # remove redundancy in spatial notes
    mutate(spatial_notes = gsub("center of | cell|~", "", spatial_notes)) %>%
    # convert to the range for the x and y axes
    mutate(x_range = reshape::colsplit(spatial_notes, 
                                       "x", c("x_range","y_range"))[,1], 
           y_range = reshape::colsplit(spatial_notes, 
                                       "x", c("x_range","y_range"))[,2]) %>% 
    # convert to numeric
    mutate(x_range = as.numeric(as.character(x_range)), 
           y_range = as.numeric(as.character(y_range))) %>% 
    # only keep data for particular spatial scale
    filter( if(spatial_res == "other") { 
      x_range != 1 & x_range != 5 & y_range != 1 & y_range != 5 
    } else {x_range == spatial_res & y_range == spatial_res}) %>%
    # create bounding box for each point
    mutate(x_min = longitude - (x_range/2), 
           x_max = longitude + (x_range/2), 
           y_min = latitude - (y_range/2), 
           y_max = latitude + (y_range/2)) 
  
  # Next - create list for polygons (each row = new polygon)
  # (Had to do in separate chunk because R keeps crashing)
  data_polygons <- data_polygons %>% 
    rowwise() %>%
    mutate(polygon_list = list(st_polygon(list(cbind(c(x_min, x_min, x_max, x_max, x_min), 
                                                     c(y_min, y_max, y_max, y_min, y_min)))))) %>% 
    st_as_sf(crs = 4326)
  
  # Reproject to mollweide
  data_moll <- st_transform(data_polygons, crs = proj4string(raster_layer))
  
  # Rasterize
  data_rasterize <- data_moll %>% 
    fasterize(raster_layer, field = summary_var, fun = summary_func)
  
  # Remove points on land
  output <- mask(data_rasterize, raster_layer > 0.99, maskvalue = 1)
  
  return(output)
}


###############################
#### model_raster_blaster() ### 
###############################

# DESCRIPTION: Takes rfmo dataset with point data, converts to polygons, transforms to mollweide 
# projection, rasterizes, resamples, and produces raster object

# VARIABLES:
# data : variable name that contains dataset
# raster_layer : variable name that contains raster layer for masking
# species : character value - all capitals of species name
# geartype : character value - all lowercase of gear_group
# unit: character value - all lowercase of catch units
# time_res : character value - all lowercase of time_res ("yearly", "quarterly", "monthly")
# mm in format of c(month value, ...) in case you wanted to filter for multiple
# yy in format of c(year value, ...) in case you wanted to filter for multiple
# flag_code in format of c(flag code, ...) in case you wanted to filter for multiple
# include_generated : TRUE/FALSE - whether to include generated data
# spatial_res : interger value - spatial resolution of grid cell (only works if x and y have identical resolution)
#   ex: "center of 5x5 cell" in spatial_notes means spatial_res == 5
# summary_var : character value - which column to use for raster creation
# summary_func : character value - which function to use for raster collating if cells are overlapping

model_raster_blaster <- function(data, raster_layer, species, geartype, unit, yy = NA,
                                include_generated = FALSE, spatial_res, summary_var = "catch", 
                                summary_func = "sum") { 
  
  # Pre-filter for all things apart from spatial resolution
  data <- data %>% 
    filter(species_sciname == species & 
             gear == geartype)
  
  # Check
  if(nrow(data) == 0) { 
    stop("Specifications too constricted, cannot rasterize - try different species, geartype, time_res or remove some specifications")
  }
  
  # Pre filter for units
  if(summary_var == "effort" ) { 
    data <- data %>% filter(effort_units == unit)
  }
  if(summary_var == "catch" | summary_var == "rel_catch") { 
    data <- data %>% filter(grepl(unit, catch_units))
  }
  
  # Continue pre-filter but only pre-filter if values are not NAs
  # Year
  if(!is.na(yy)) { 
    data <- data %>% filter(year %in% yy)
    if(nrow(data) == 0) { 
      stop("Specifications too constricted, cannot rasterize - try different year or remove year specification")
    }
  }
  
  # Spatial Check
  if(nrow(data[grepl(spatial_res, data$spatial_notes),]) == 0) { 
    stop("Specifications too constricted, cannot rasterize - try different spatial resolution or changing other parameters")
  }
  
  # Convert to polygons for rasterize/fasterize
  # First step - figuring out spatial resolution of the points
  data_polygons <- data %>%
    # remove redundancy in spatial notes
    mutate(spatial_notes = gsub("center of | cell|~", "", spatial_notes)) %>%
    # convert to the range for the x and y axes
    mutate(x_range = reshape::colsplit(spatial_notes, 
                                       "x", c("x_range","y_range"))[,1], 
           y_range = reshape::colsplit(spatial_notes, 
                                       "x", c("x_range","y_range"))[,2]) %>% 
    # convert to numeric
    mutate(x_range = as.numeric(as.character(x_range)), 
           y_range = as.numeric(as.character(y_range))) %>% 
    # only keep data for particular spatial scale
    filter( if(spatial_res == "other") { 
      x_range != 1 & x_range != 5 & y_range != 1 & y_range != 5 
    } else {x_range == spatial_res & y_range == spatial_res}) %>%
    # create bounding box for each point
    mutate(x_min = longitude - (x_range/2), 
           x_max = longitude + (x_range/2), 
           y_min = latitude - (y_range/2), 
           y_max = latitude + (y_range/2)) 
  
  # Next - create list for polygons (each row = new polygon)
  # (Had to do in separate chunk because R keeps crashing)
  data_polygons <- data_polygons %>% 
    rowwise() %>%
    mutate(polygon_list = list(st_polygon(list(cbind(c(x_min, x_min, x_max, x_max, x_min), 
                                                     c(y_min, y_max, y_max, y_min, y_min)))))) %>% 
    st_as_sf(crs = 4326)
  
  # Reproject to mollweide
  data_moll <- st_transform(data_polygons, crs = proj4string(raster_layer))
  
  # Rasterize
  data_rasterize <- data_moll %>% 
    fasterize(raster_layer, field = summary_var, fun = summary_func)
  
  # Remove points on land
  output <- mask(data_rasterize, raster_layer > 0.99, maskvalue = 1)
  
  return(output)
}

############################
### gfw_raster_blaster() ###
############################

# DESCRIPTION: Takes gfw queried dataset with point data, converts to polygons (1x1 degree), transforms to mollweide 
# projection, rasterizes, resamples, and produces raster object

# VARIABLES:
# data : variable name that contains dataset
# raster_layer : variable name that contains raster layer for masking
# tonnage_layer : variable name that contains the min/max values for each gear group
# geartype : character value - all lowercase of gear
# nighttime: TRUE/FALSE - should layer be for night or day hours
# mm in format of c(month value, ...) in case you wanted to filter for multiple
# yy in format of c(year value, ...) in case you wanted to filter for multiple
# rfmo_code in format of c(rfmo, ...) in case you wanted to filter for different rfmos
# flag_code in format of c(flag code, ...) in case you wanted to filter for multiple
# tonnage_group : integer value between 0 and 4 where 0 is the smallest binned value and 4 is the largest
# summary_func : character value - which function to use for raster collating if cells are overlapping

gfw_raster_blaster <- function(data, raster_layer, tonnage_layer  = NA, geartype, nighttime = FALSE, mm = NA, yy = NA, 
                               rfmo_code = c("ICCAT", "CCSBT", "IOTC", "IATTC", "WCPFC", "NAFO"), spatial_scale = 1,
                               flag_code = NA, tonnage_group = NA, summary_func = "sum") { 

  
  # Pre-filter for all things
  data <- data %>% 
    filter(gear == geartype & 
           night == nighttime &
           rfmo %in% rfmo_code)
  
  # Check
  if(nrow(data) == 0) { 
    stop("Specifications too constricted, cannot rasterize - try different geartype, nighttime, or rfmo_code or remove some specifications")
  }
  
  # Continue pre-filter but only pre-filter if values are not NAs
  # Month
  if(!is.na(mm)) { 
    data <- data %>% filter(month %in% mm)
    if(nrow(data) == 0) { 
      stop("Specifications too constricted, cannot rasterize - try different month or remove month specification")
    }
  }
  
  # Year
  if(!is.na(yy)) { 
    data <- data %>% filter(year %in% yy)
    if(nrow(data) == 0) { 
      stop("Specifications too constricted, cannot rasterize - try different year or remove year specification")
    }
  }
  
  # Flag
  if(!is.na(flag_code)) { 
    data <- data %>% filter(flag %in% flag_code)
    if(nrow(data) == 0) { 
      stop("Specifications too constricted, cannot rasterize - try different flag_code or remove flag specification")
    }
  }
  
  # Tonnage
  if(!is.na(tonnage_group)) { 
  tonnage_layer <- tonnage_layer %>% 
    filter(gear == geartype) %>%
    mutate(binning_value = (ceiling(max_tonnage)-floor(min_tonnage))/4)
    
    ton_value <- (tonnage_layer$binning_value)*tonnage_group
    data <- data %>% filter(binned_gt == ton_value)
    if(nrow(data) == 0) { 
      stop("Specifications too constricted, cannot rasterize - try different tonnage or remove tonnage specification")
    }
  }
  
  # Convert to polygons for rasterize/fasterize
  # First step - figuring out spatial resolution of the points
  data_polygons <- data %>%
    # create bounding box for each point
    mutate(x_min = lon_bin - (0.5), 
           x_max = lon_bin + (0.5), 
           y_min = lat_bin - (0.5), 
           y_max = lat_bin + (0.5)) 
  
  if(spatial_scale != 1) { 
    # Convert to polygons for rasterize/fasterize
    # First step - figuring out spatial resolution of the points
    data_polygons <- data %>%
      # create bounding box for each point
      mutate(lon_bin = round(lon_bin/spatial_scale)*spatial_scale, 
             lat_bin = round(lat_bin/spatial_scale)*spatial_scale) %>% 
      mutate(x_min = ifelse(lon_bin != -180, lon_bin - (spatial_scale/2), lon_bin), 
             x_max = ifelse(lon_bin != 180, lon_bin + (spatial_scale/2), lon_bin),
             y_min = lat_bin - (spatial_scale/2), 
             y_max = lat_bin + (spatial_scale/2))
  }
  
  # Next - create list for polygons (each row = new polygon)
  # (Had to do in separate chunk because R keeps crashing)
  data_polygons <- data_polygons %>% 
    rowwise() %>%
    mutate(polygon_list = list(st_polygon(list(cbind(c(x_min, x_min, x_max, x_max, x_min), 
                                                     c(y_min, y_max, y_max, y_min, y_min)))))) %>% 
    st_as_sf(crs = 4326)
  
  # Reproject to mollweide
  data_moll <- st_transform(data_polygons, crs = proj4string(raster_layer))
  
  # Rasterize
  data_rasterize <- data_moll %>% 
    fasterize(raster_layer, field = "total_fishing_kwh", fun = summary_func)
  
  # Remove points on land
  output <- mask(data_rasterize, raster_layer > 0.99, maskvalue = 1)
  
  return(output)
}

# for sdm data
sdm_raster_blaster <- function(data, raster_layer, species_group, scaled = TRUE, pacific = TRUE) { 
  sdm <- data 
  
  sdm <- sdm %>% filter(group == species_group) %>% 
    group_by(CenterLat, CenterLong) %>% 
    summarise(probability = sum(probability)) %>% 
    ungroup()
  
  sdm <- sdm %>% 
    rowwise() %>%
    mutate(polygon_list = list(st_polygon(list(cbind(c(CenterLong-0.25, CenterLong-0.25, 
                                                       CenterLong+0.25, CenterLong+0.25, CenterLong-0.25), 
                                                     c(CenterLat-0.25, CenterLat+0.25, CenterLat+0.25,
                                                       CenterLat-0.25, CenterLat-0.25)))))) %>% 
    st_as_sf(crs = 4326)
  
  # Sample up to 5x5 degree cells...
  sdm_polygons <- sdm %>% 
    mutate(CenterLong5x5 = round(CenterLong/5)*5,
           CenterLat5x5 = round(CenterLat/5)*5) %>% 
    group_by(CenterLong5x5, CenterLat5x5) %>% 
    mutate(sum_probability = sum(probability), 
           sum_probability = ifelse(CenterLong5x5 == 180 | CenterLong5x5 == -180, 
                                    sum_probability*2, sum_probability)) %>% 
    select(CenterLat5x5, CenterLong5x5, sum_probability) %>%
    distinct_all() %>% 
    rowwise() %>% 
    mutate(polygon_list = ifelse(CenterLong5x5 != -180 & CenterLong5x5 != 180,
                                 list(st_polygon(list(cbind(c(CenterLong5x5-2.5, CenterLong5x5-2.5, 
                                                              CenterLong5x5+2.5, CenterLong5x5+2.5, 
                                                              CenterLong5x5-2.5), 
                                                            c(CenterLat5x5-2.5, CenterLat5x5+2.5, 
                                                              CenterLat5x5+2.5, CenterLat5x5-2.5,
                                                              CenterLat5x5-2.5))))),
                                 ifelse(CenterLong5x5 == -180, 
                                        list(st_polygon(list(cbind(c(CenterLong5x5, CenterLong5x5, 
                                                                     CenterLong5x5+2.5, CenterLong5x5+2.5, 
                                                                     CenterLong5x5), 
                                                                   c(CenterLat5x5-2.5, CenterLat5x5+2.5, 
                                                                     CenterLat5x5+2.5, CenterLat5x5-2.5,
                                                                     CenterLat5x5-2.5))))),
                                        list(st_polygon(list(cbind(c(CenterLong5x5-2.5, CenterLong5x5-2.5, 
                                                                     CenterLong5x5, CenterLong5x5, 
                                                                     CenterLong5x5-2.5), 
                                                                   c(CenterLat5x5-2.5, CenterLat5x5+2.5, 
                                                                     CenterLat5x5+2.5, CenterLat5x5-2.5,
                                                                     CenterLat5x5-2.5)))))))) %>%
    st_as_sf(crs = 4326)
  
  # Reproject to mollweide
  sdm_moll <- st_transform(sdm_polygons, crs = proj4string(raster_layer))
  
  # Rasterize
  sdm_rasterize <- sdm_moll %>% 
    fasterize(raster_layer > 0.99, field = "sum_probability", fun = "sum")
  
  # Remove points on land
  sdm_masked <- mask(sdm_rasterize, raster_layer > 0.99, maskvalue = 1) 
  
  # Rescale
  if(pacific == TRUE) { 
    sdm_masked <- sdm_masked %>% 
    projectRaster(crs = CRS("+proj=moll +lon_0=-18035095 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")) %>% 
    extend(raster_layer)
  } 
  
  if(scaled == TRUE) { 
    sdm_masked = sdm_masked/sdm_masked@data@max
    }
  
  return(sdm_masked)
} 

# Simple raster
simple_raster_blaster <- function(data, raster_layer, summary_var, scale = 5, summary_func = "sum", pacific = TRUE) { 
  
  data_polygons <- data %>% 
    rowwise() %>%
    mutate(polygon_list = list(st_polygon(list(cbind(c(longitude-scale/2, longitude-scale/2, longitude+scale/2, longitude+scale/2, longitude-scale/2), 
                                                     c(latitude-scale/2, latitude+scale/2, latitude+scale/2, latitude-scale/2, latitude-scale/2)))))) %>% 
    st_as_sf(crs = 4326)
  
  # Reproject to mollweide
  data_moll <- st_transform(data_polygons, crs = proj4string(raster_layer))
  
  # Rasterize
  data_rasterize <- data_moll %>% 
    fasterize(raster_layer, field = summary_var, fun = summary_func)
  
  # Remove points on land
  output <- mask(data_rasterize, raster_layer > 0.99, maskvalue = 1)
  
  if(pacific == TRUE) { 
    output <- output %>% 
      projectRaster(crs = CRS("+proj=moll +lon_0=-18035095 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")) %>% 
      extend(raster_layer)
    }
  
  return(output)
  
  
  }
