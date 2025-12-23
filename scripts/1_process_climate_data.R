# Process environmental and geographic data

library(dplyr)
library(terra)
library(sf)
library(ncdf4)

# import raster template
my_ras <- readRDS("data/processed/my_ras.rds")

# breeding months
breed_N <- 3:6
breed_S <- 9:12
breed_T <- 1:12

# ------------------------------------------------------------------
# Compile CHELSA monthly data (1981-2010) into breeding season data
# ------------------------------------------------------------------
# Data source: CHELSA v2.1 monthly data (https://www.chelsa-climate.org/datasets/chelsa_monthly)
# Original data are not included in this code package due to their large file size.

compile_breeding_climate <- function(source_folder, filename) {
  stack_monthly <- c(rast(list.files(path = source_folder, pattern = ".tif", full.names = TRUE)))
  dt_n <- stack_monthly[[breed_N]] %>% crop(ext(-180, 180, 23.5, 90)) %>% mean # Nothern hemisphere
  dt_s <- stack_monthly[[breed_S]] %>% crop(ext(-180, 180, -90, -23.5)) %>% mean # Southern hemisphere 
  dt_t <- stack_monthly[[breed_T]] %>% crop(ext(-180, 180, -23.5, 23.5)) %>% mean # Tropical region
  dt_breed <- merge(dt_n, dt_s, dt_t) %>% project(my_ras, method = "bilinear")
  plot(dt_breed)
  writeRaster(dt_breed, filename, overwrite =  T)
}

path <- "data/processed/breeding_clim"

compile_breeding_climate("data/input/climate/chelsa_hurs", file.path(path, "hurs.tif"))
compile_breeding_climate("data/input/climate/chelsa_rsds", file.path(path, "rsds.tif"))
compile_breeding_climate("data/input/climate/chelsa_sfcWind", file.path(path, "sfcWind.tif"))
compile_breeding_climate("data/input/climate/chelsa_vpd", file.path(path, "vpd.tif"))
compile_breeding_climate("data/input/climate/chelsa_pr", file.path(path, "pr.tif"))
compile_breeding_climate("data/input/climate/chelsa_tas", file.path(path, "tmean.tif"))
compile_breeding_climate("data/input/climate/chelsa_tasmax", file.path(path, "tmax.tif"))
compile_breeding_climate("data/input/climate/chelsa_tasmin", file.path(path, "tmin.tif"))

# compute diurnal temperature range (DTR)
dtr <- rast("data/breeding_clim/tmax.tif") - rast("data/breeding_clim/tmin.tif") # tmax - tmin
writeRaster(dtr, "data/breeding_clim/dtr.tif")

# ------------------------------------------------------------------
# Visualization
# ------------------------------------------------------------------
library(ggplot2)

land <- rast("data/input/land.tif")
idx.land <- which(values(land) == 1)

plot_clim <- function(file, var_name) {
  data <- rast(file)
  values(data)[-idx.land] <- NA
  names(data) <- var_name
  ggplot() + 
    geom_raster(data = as.data.frame(data, xy = T), aes(x = x, y = y, fill = !!sym(var_name))) +
    scale_fill_viridis_c(option = "inferno", name = var_name) + 
    coord_fixed() + 
    theme_void()
}

plot_clim(file.path(path, "hurs.tif"), "HURS")
plot_clim(file.path(path, "rsds.tif"), "RSDS")
plot_clim(file.path(path, "sfcWind.tif"), "sfcWind")
plot_clim(file.path(path, "vpd.tif"), "VPD")
plot_clim(file.path(path, "pr.tif"), "Precp")
plot_clim(file.path(path, "tmean.tif"), "TMean")
plot_clim(file.path(path, "dtr.tif"), "DTR")