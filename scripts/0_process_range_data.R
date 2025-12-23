# Process distribution data and export breeding range maps for each species

# Data source: Area of Habitat Maps (https://datadryad.org/dataset/doi:10.5061/dryad.02v6wwq48)
# Original data are not included in this code package due to their large file size.
# Preparation: data downloaded and unzipped, migratory and non-migratory files separated into two folders

library(dplyr)
library(terra)
library(raster)
library(parallel) # for multiprocessing, only works on Unix OS

# create raster template
extent <- c(-180, 180, -90, 90) * 2 * pi * 6371000 / 360 # degree converting to meter (earth radius 6371 km)
my_ras <- rast(extent = extent, resolution = 10000, crs = "+proj=moll +datum=WGS84") %>%
  raster() # convert to object of the raster package for faster reprojection/resampling
saveRDS(my_ras, "data/processed/my_ras.rds")

# -------------------------------------------------------
# Migratory species (merge breeding and resident ranges)
# -------------------------------------------------------
merge_resident_breeding_range <- function (sppath) {
  # extract breeding and resident map files
  files <- paste0(path, sppath, c("_B.tif", "_R.tif"))
  
  # reproject them to raster template
  rasters <- lapply(files, function (x) if (file.exists(x)) raster(x) %>% 
    projectRaster(to = my_ras, method = "ngb"))
  rasters <- lapply(files, FUN = reproject_map)
  rasters <- rasters[!sapply(rasters, is.null)] # exclude null results in the list
  
  # merge maps (if both breeding and resident maps are present, merge; if only one present, use that one)
  if (length(rasters) == 2) {
    range_merge <- do.call("merge", rasters)
  } else if (length(rasters) == 1) {
    range_merge <- rasters[[1]]
  }
  
  # mark raster cells within the distribution as 1 and others as 0
  values(range_merge)[is.na(values(range_merge))] <- 0
  
  # save map
  spname <- gsub(".*/", "", sppath)
  writeRaster(range_merge, paste0("data/processed/breeding_map/AOH_mig/", spname, ".tif"))
}

# extract species list to process
path <- "data/input/AOH_bird/migratory/"
filelist <- list.files(path = path, pattern = "\\.tif$", recursive = TRUE)
sp.paths <- sapply(strsplit(filelist, "_"), function(x) paste(x[1:4], collapse = "_")) %>% unique

mclapply(sp.paths, FUN = merge_resident_breeding_range, mc.cores = 7) # ~1 hr with 7 cores on PC

# -------------------------------------------------------
# Non-migratory species
# -------------------------------------------------------
reproject_map <- function (sppath) {
  rs <- raster(paste0(path, sppath)) %>% 
    projectRaster(to = my_ras, method = "ngb")
  values(rs)[is.na(values(rs))] <- 0

  filename <- gsub(".*/", "", sppath)
  writeRaster(rs, paste0("data/processed/breeding_map/AOH_nonmig/", filename))
}
path <- "data/AOH_bird/nonmigratory/"
files <- list.files(path = path, pattern = "\\.tif$", recursive = TRUE)

mclapply(files, FUN = reproject_map, mc.cores = 7) # ~2 hr with 7 cores on PC
