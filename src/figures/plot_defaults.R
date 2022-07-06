# Plotting defaults for figures

# Basemap
basemap_raster <- raster::raster(file.path(here::here(), "data/mapping-templates/land_low_res_moll.tif")) %>% 
  projectRaster(., crs = 4326)
basemap_df <- raster::as.data.frame(basemap_raster, xy = T) %>% 
  mutate(land_low_res_moll = ifelse(land_low_res_moll >= 0.9, 1, NA))

# Four options for rasterizing RFMO data (depends on RFMO location reporting)
whole_numbers <- raster::raster(xmn = -179, xmx = 179, ymn = -89, ymx = 89, crs = 4326, resolution = 1) #%>% 
  #projectRaster(., crs = basemap_raster)
whole_numbers_lon <- raster::raster(xmn = -179, xmx = 179, ymn = -89.5, ymx = 89.5, crs = 4326, resolution = 1) #%>% 
  #projectRaster(., crs = basemap_raster)
whole_numbers_lat <- raster::raster(xmn = -179.5, xmx = 179.5, ymn = -89, ymx = 89, crs = 4326, resolution = 1) #%>% 
  #projectRaster(., crs = basemap_raster)
fraction_numbers <- raster::raster(xmn = -179.5, xmx = 179.5, ymn = -89.5, ymx = 89.5, crs = 4326, resolution = 1) #%>% 
  #projectRaster(., crs = basemap_raster)

# RFMO boundaries
wcpfc_boundary <- sf::st_read(file.path(here::here(), "data/rfmo-boundaries/RFB_WCPFC/RFB_WCPFC.shp"))
iotc_boundary <- sf::st_read(file.path(here::here(), "data/rfmo-boundaries/RFB_IOTC/RFB_IOTC.shp"))
iccat_boundary <- sf::st_read(file.path(here::here(), "data/rfmo-boundaries/RFB_ICCAT/RFB_ICCAT.shp"))
iattc_boundary <- sf::st_read(file.path(here::here(), "data/rfmo-boundaries/RFB_IATTC/RFB_IATTC.shp"))

# map theme
custom_theme <- theme_void() 
