# Load libraries
library(sf)
library(terra)
library(gstat)
library(dplyr)
library(spatstat)
library(sp)
library(stars)

#Load the data

data.directory <- 'Z://public//CTD//'  # may differ for OSX users
df <- read.csv(paste0(data.directory, 'mgl2506_castmetadata.csv'))

spatial.plot <- function(param = 'z_1026_.m.', name = '1026 isopycnal depth'){

  # Rename columns for convenience (optional)
  
  df <- df[,c('lon', 'lat', param)]
  colnames(df) <- c('lon', 'lat', 'z')
  df <- na.omit(df)

  # Convert to sf object

  points_sf <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)

  # Project to a suitable planar CRS for kriging (e.g., UTM zone)

  utm_crs <- 32610  # Change based on your region
  points_sf_proj <- st_transform(points_sf, crs = utm_crs)

  # Create prediction grid using terra

  bbox <- st_bbox(points_sf_proj)

  r <- rast(xmin = bbox["xmin"] - 2e5, xmax = bbox["xmax"],
            ymin = bbox["ymin"], ymax = bbox["ymax"],
            resolution = 1000)  # adjust resolution (in meters)

  grid_points <- as.points(r)

  # Convert to sf for easier manipulation

  grid_sf <- st_as_sf(grid_points)

  # Rotate grid using affine transformation

  angle_deg <- 32
  angle_rad <- angle_deg * pi / 180

  # Rotation matrix

  rot_mat <- matrix(c(
    cos(angle_rad), -sin(angle_rad),
    sin(angle_rad),  cos(angle_rad)
  ), ncol = 2, byrow = TRUE)

  # Extract coords and rotate
  
  coords <- st_coordinates(grid_sf)
  
  center_x <- mean(coords[, 1])
  center_y <- mean(coords[, 2])  # rotate around center of grid
  
  # Translate coordinates to origin
  
  translated <- coords
  translated[, 1] <- coords[, 1] - center_x
  translated[, 2] <- coords[, 2] - center_y
  
  rotated <- t(rot_mat %*% t(translated))
  
  rotated[, 1] <- rotated[, 1] + center_x
  rotated[, 2] <- rotated[, 2] + center_y
  
  grid_sf_rotated <- st_as_sf(
    data.frame(x = rotated[, 1], y = rotated[, 2]),
    coords = c("x", "y"),
    crs = st_crs(grid_sf)
  )
  
  # Convert to Spatial for kriging
  
  grid_sp_rotated <- as(grid_sf_rotated, "Spatial")
  
  # Set projection
  
  proj4string(grid_sp_rotated) <- CRS(SRS_string = "EPSG:32610")
  
  points_sp <- as(points_sf_proj, "Spatial")
  
  # Fit variogram and perform ordinary kriging
  
  vgm_model <- variogram(z ~ 1, data = points_sp)
  fitted_vgm <- fit.variogram(vgm_model, model = vgm(c("Exp", "Sph", "Gau", "Mat")))
  
  kriging_result <- krige(z ~ 1, locations = points_sp, newdata = grid_sp_rotated, model = fitted_vgm)
  
  # Convert kriging result to sf with coordinates and prediction values
  
  krig_sf <- st_as_sf(kriging_result)
  
  # Rasterize star object and then convert to df
  
  krig_rast <- st_rasterize(krig_sf)
  
  krig_df <- as.data.frame(krig_rast)[,1:3]
  colnames(krig_df) <- c('x', 'y', 'z')
  
  # Change utm_crs to match your data location (e.g., UTM zones 10–11 for the eastern Pacific).
  # Adjust resolution = 1000 for finer or coarser grids.
  # Use z ~ x + y if you want to try universal kriging with a trend.
  
  library(sf)
  library(terra)
  library(ggplot2)
  library(rnaturalearth)
  library(ggspatial)
  library(viridis)
  
  utm_crs <- 32610
  
  # Transform coastline and kriging raster to same CRS
  
  coastline <- ne_coastline(scale = "medium", returnclass = "sf")
  coastline <- ne_states()
  coastline_proj <- st_transform(coastline, crs = utm_crs)
  
  # Plot with proper projected CRS
  
  print(
  
    ggplot() +
      geom_raster(data = krig_df, aes(x = x, y = y, fill = z)) +
      scale_fill_viridis(option = "plasma", name = "(m)", direction = -1, guide = guide_colorbar(reverse = TRUE), na.value = 'transparent') +
      geom_sf(data = coastline_proj, color = "black", fill = 'gray') +
      geom_sf(data = points_sf_proj, color = "black", size = 1.5) +
      coord_sf(xlim = c(300000, 1122252), ylim = c(3275791, 4297791), crs = st_crs(utm_crs)) +
      theme_minimal() +
      labs(title = name, x = 'lon', y = 'lat')
  )
  
}

param.names <- list()

param.names["ml_depth_.m."] <- 'Mixed layer depth (m)'
param.names["maxc_.mg.m.3."] <- 'Max chlorophyll (mg m^-3)'
param.names["comp_depth_.m."] <- 'Compensation depth (m)'
param.names["z_1026_.m."] <- '1026 Isopycnal depth (m)'

pdf(paste0(data.directory, 'spatial_plots.pdf'),
    width = 8,
    height = 8)

for(param in names(param.names)){
  name <- param.names[param]
  spatial.plot(param, name)
}

dev.off()