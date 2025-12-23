# Extract environmental variables for each species

library(dplyr)
library(terra)
library(Matrix)

# Versions
North <- FALSE # Global
North <- TRUE # Northern Hemisphere only

if (TRUE) {
  # import environmental datasets
  pred <- rast("data/processed/predator/predatorSR.tif")
  tmean <- rast("data/processed/breeding_clim/tmean.tif")
  dtr <- rast("data/processed/breeding_clim/dtr.tif")
  hurs <- rast("data/processed/breeding_clim/hurs.tif")
  rsds <- rast("data/processed/breeding_clim/rsds.tif")
  vpd <- rast("data/processed/breeding_clim/vpd.tif")
  precp <- rast("data/processed/breeding_clim/pr.tif")
  wind <- rast("data/processed/breeding_clim/sfcWind.tif")
  if (North) {
    sp_occur <- readRDS("data/processed/sp_north.rds") # sparse matrix
  } else {
    sp_occur <- readRDS("data/processed/sp_occur.rds") # sparse matrix
  }
  
  ## subset environmental and distribution datasets to keep shared non-NA grid cells
  # extract terrestrial and non-NA cell indices across all datasets
  land <- rast("data/input/land.tif")
  idx_land <- which(values(land) == 1)
  if (North) idx_land <- idx_land[which(idx_land <= 1001 * ncol(land))] #### for norther hemisphere analysis
  clim_mtx <- cbind(values(pred), values(tmean), values(dtr), values(hurs), values(rsds), values(vpd), values(precp), values(wind))[idx_land, ]
  idx_keep <- which(rowSums(is.na(clim_mtx)) == 0)
  clim_mtx <- clim_mtx[idx_keep, ]
  
  # compute environmental values for each species (mean value across the species' distribution range)
  range_size <- rowSums(sp_occur)
  sp_clim <- sp_occur[, idx_keep] %*% clim_mtx %>% # matrix multiplication of [species x cell] * [cell x envi]
    sweep(1, range_size, "/") %>% # divide each species value by the number of cells of its distribution range
    as.matrix
  colnames(sp_clim) <- c("pred", "tmean", "dtr", "hurs", "rsds", "vpd", "precp", "wind")
  
  if (North) {
    dir.create("data/processed/trait_north", showWarnings = FALSE, recursive = TRUE)
    write.csv(sp_clim, "data/processed/trait_north/sp_envi.csv")
  } else {
    dir.create("data/processed/trait", showWarnings = FALSE, recursive = TRUE)
    write.csv(sp_clim, "data/processed/trait/sp_envi.csv")
  }
}

