# Plotting defaults for figures

# Basemap
basemap_raster <- raster::raster(file.path(here::here(), "data/mapping-templates/land_low_res_moll.tif"))
basemap_df <- raster::as.data.frame(basemap_raster, xy = T) %>% 
  mutate(land_low_res_moll = ifelse(land_low_res_moll >= 0.9, 1, NA))

# Four options for rasterizing RFMO data (depends on RFMO location reporting)
whole_numbers <- raster::raster(xmn = -179, xmx = 179, ymn = -89, ymx = 89, crs = 4326, resolution = 1)
whole_numbers_lon <- raster::raster(xmn = -179, xmx = 179, ymn = -89.5, ymx = 89.5, crs = 4326, resolution = 1)
whole_numbers_lat <- raster::raster(xmn = -179.5, xmx = 179.5, ymn = -89, ymx = 89, crs = 4326, resolution = 1)
fraction_numbers <- raster::raster(xmn = -179.5, xmx = 179.5, ymn = -89.5, ymx = 89.5, crs = 4326, resolution = 1)

# map theme
custom_theme <- theme_void() 
